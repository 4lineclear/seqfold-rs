//! Predict nucleic acid secondary structure

use std::{
    fmt::Display,
    ops::{Index, IndexMut},
};

use crate::{
    Energies, TmCache,
    util::{ByteStr, ToIsize, round1, round2},
};

#[cfg(test)]
mod test;

// TODO: explore branchless code techniques

/// A single structure with a free energy, description, and inward children.
#[derive(Debug, Clone)]
pub struct Value {
    pub e: f64,
    pub desc: Desc,
    pub ij: Indices,
}

// Vec<[T; self.len]>
// struct Storage {
//     len: usize,
//     values: Vec<T>,
// }

impl PartialEq for Value {
    fn eq(&self, other: &Self) -> bool {
        self.e == other.e
    }
}

pub type Indices = Vec<(usize, usize)>;

#[derive(Clone)]
pub enum Desc {
    Empty,
    Bifurcation(usize, usize),
    InteriorLoop(usize, usize),
    Hairpin([u8; 5]),
    Stack(Vec<u8>),
    StackDe([u8; 5]),
    Bulge(usize),
}

impl std::fmt::Debug for Desc {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        Display::fmt(self, f)
    }
}

impl Display for Desc {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Desc::Empty => write!(f, ""),
            Desc::Bifurcation(unpaired, branches_count) => {
                write!(f, "BIFURCATION:{unpaired}n/{branches_count}h")
            }
            Desc::InteriorLoop(i, j) => write!(f, "INTERIOR_LOOP:{i}/{j}"),
            Desc::Hairpin(pair) => write!(f, "HAIRPIN:{}", ByteStr(pair)),
            Desc::Stack(pair) => write!(f, "STACK:{}", ByteStr(pair)),
            Desc::StackDe(pair) => write!(f, "STACK_DE:{}", ByteStr(pair)),
            Desc::Bulge(i) => write!(f, "BULGE:{i}"),
        }
    }
}

impl From<f64> for Value {
    fn from(value: f64) -> Self {
        Self::new(value, Desc::Empty, Indices::new())
    }
}

impl Value {
    pub const NULL: Self = Self::empty(f64::INFINITY, Desc::Empty);
    pub const DEFAULT: Self = Self::empty(f64::NEG_INFINITY, Desc::Empty);

    pub const fn empty(e: f64, desc: Desc) -> Value {
        Self::new(e, desc, Indices::new())
    }

    pub const fn new(e: f64, desc: Desc, ij: Indices) -> Self {
        Self { e, desc, ij }
    }

    pub fn with_ij(mut self, ij: Indices) -> Self {
        self.ij = ij;
        self
    }

    pub fn valid(&self) -> bool {
        self.e != f64::INFINITY && self.e != f64::NEG_INFINITY
    }
}

/// A map from i, j tuple to a min free energy Struct.
pub struct Context {
    /// - v_cache: Values of energies where V(i,j) bond
    v: Matrix,
    /// - w_cache: Values of min energy of substructures between W(i,j)
    w: Matrix,
}

impl Context {
    pub fn new(n: usize) -> Self {
        Context {
            v: Matrix::new(n),
            w: Matrix::new(n),
        }
    }
}

pub struct Matrix {
    values: Box<[Value]>,
    len: usize,
}

impl Matrix {
    pub fn new(len: usize) -> Self {
        let values = vec![Value::DEFAULT; len * len].into_boxed_slice();
        Self { values, len }
    }
}

impl Index<usize> for Matrix {
    type Output = [Value];

    fn index(&self, index: usize) -> &Self::Output {
        &self.values[index * self.len..(index + 1) * self.len]
    }
}

impl IndexMut<usize> for Matrix {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.values[index * self.len..(index + 1) * self.len]
    }
}

pub type SeqResult<T> = Result<T, SeqError>;

/// An invalid sequence was inputted.
#[derive(Debug)]
pub enum SeqError {
    /// Both T and U in sequence. Provide one or the other for DNA OR RNA.
    DnaOrRna,
    /// Unknown bp encountered.
    UnknownBp(String),
}

impl Display for SeqError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            SeqError::DnaOrRna => write!(
                f,
                "Both T and U in sequence. Provide one or the other for DNA OR RNA."
            ),
            SeqError::UnknownBp(s) => write!(f, "Unknown bp encountered: {s}"),
        }
    }
}

impl std::error::Error for SeqError {}

/// Fold the DNA sequence and return the lowest free energy score.
///
/// Based on the approach described in:
/// Zuker and Stiegler, 1981
/// <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC326673/pdf/nar00394-0137.pdf>
///
/// If the sequence is 50 or more bp long, "isolated" matching bp
/// are ignored in V(i,j). This is based on an approach described in:
/// Mathews, Sabina, Zuker and Turner, 1999
/// <https://www.ncbi.nlm.nih.gov/pubmed/10329189>
///
/// # Args
///
/// - seq: The sequence to fold
/// - temp: The temperature the fold takes place in, in Celcius
///
/// # Returns
///
/// - [`Vec<Value>`]: A list of structures. Stacks, bulges, hairpins, etc.
pub fn fold(seq: &[u8], temp: Option<f64>) -> Vec<Value> {
    try_fold(seq, temp).expect("Invalid Sequence Inputted")
}

/// A non-panicing version of [`fold`]
pub fn try_fold(seq: &[u8], temp: Option<f64>) -> SeqResult<Vec<Value>> {
    let ctx = cache(seq, temp)?;
    let n = seq.len();

    // get the minimum free energy structure out of the cache
    Ok(traceback(0, n - 1, &ctx))
}

/// Fold the sequence and return just the delta G of the structure
///
/// # Args
///
/// - seq: The sequence to fold
/// - temp: The temperature to fold at
///
/// # Returns
///
/// - [`f64`]: The minimum free energy of the folded sequence
pub fn dg(seq: &[u8], temp: Option<f64>) -> f64 {
    try_dg(seq, temp).expect("Invalid Sequence Inputted")
}

/// A non-panicing version of [`dg`]
pub fn try_dg(seq: &[u8], temp: Option<f64>) -> SeqResult<f64> {
    let values = try_fold(seq, temp)?;
    let dg_sum = values.iter().map(|v| v.e).sum();
    Ok(round2(dg_sum))
}

/// Fold a nucleic acid sequence and return the estimated dg of each (i,j) pairing.
///
/// # Args
///
/// - seq: The nucleic acid sequence to fold
/// - temp: The temperature to fold at
///
/// # Returns
///
/// - [`Cache`]: A 2D matrix where each (i, j) pairing corresponds to the
///   minimum free energy between i and j
pub fn dg_cache(seq: &[u8], temp: Option<f64>) -> TmCache {
    try_dg_cache(seq, temp).expect("Invalid Sequence Inputted")
}

/// A non-panicing version of [`dg_cache`]
pub fn try_dg_cache(seq: &[u8], temp: Option<f64>) -> SeqResult<TmCache> {
    let cache = cache(seq, temp)?;
    Ok((0..cache.w.len)
        .map(|i| cache.w[i].iter().map(|v| v.e).collect())
        .collect())
}

/// Get the dot bracket notation for a secondary structure.
///
/// # Args
///
/// - structs: A list of structs, usually from the fold function
///
/// # Returns
///
/// - [`Vec<u8>`]: The dot bracket notation of the secondary structure
pub fn dot_bracket(seq: &[u8], values: &[Value]) -> Vec<u8> {
    let mut result = vec![b'.'; seq.len()];
    for v in values {
        if v.ij.len() == 1 {
            let (i, j) = v.ij[0];
            result[i] = b'(';
            result[j] = b')';
        }
    }
    result
}

/// Create caches for the w_cache and v_cache
///
/// The Values are useful for gathering many possible energies
/// between a series of (i,j) combinations.
///
/// # Args
///
/// - seq: The sequence to fold
/// - temp: The temperature to fold at
///
/// # Returns
///
/// - [`(Values, Values)`]: The w_cache and the v_cache for traversal later
pub fn cache(seq: &[u8], temp: Option<f64>) -> SeqResult<Context> {
    let temp = temp.unwrap_or(37.0) + 273.15; // kelvin

    let (seq, emap) = parse_sequence(seq)?;

    let n = seq.len();
    let mut cache = Context::new(n);
    // fill the cache
    w(&seq, 0, n - 1, temp, &mut cache, emap);
    Ok(cache)
}

fn parse_sequence(seq: &[u8]) -> SeqResult<(Vec<u8>, &'static Energies)> {
    let seq = seq.to_ascii_uppercase();
    // figure out whether it's DNA or RNA, choose energy map
    let mut dna = true;

    let bps: rustc_hash::FxHashSet<_> = seq.iter().copied().collect();

    if bps.contains(&b'U') && bps.contains(&b'T') {
        return Err(SeqError::DnaOrRna);
    }

    if bps.iter().all(|b| b"AUCG".contains(b)) {
        dna = false;
    } else if bps.iter().any(|b| !b"ATGC".contains(b)) {
        return Err(SeqError::UnknownBp(
            bps.iter()
                .filter(|b| !b"ATGC".contains(b))
                .map(|b| *b as char)
                .collect::<String>(),
        ));
    }

    Ok((seq, if dna { crate::dna() } else { crate::rna() }))
}

/// Find and return the lowest free energy structure in Sij subsequence
///
/// Figure 2B in Zuker and Stiegler, 1981
///
/// # Args
///
/// - seq: The sequence being folded
/// - i: The start index
/// - j: The end index (inclusive)
/// - temp: The temperature in Kelvin
/// - ctx: Combined free energy cache, etc
/// - emap: Energy map for DNA/RNA
///
/// # Returns
///
/// - [`f64`]: The free energy for the subsequence from i to j
pub fn w(
    seq: &[u8],
    i: usize,
    j: usize,
    temp: f64,
    ctx: &mut Context,
    emap: &Energies,
) -> (usize, usize) {
    if ctx.w[i][j] != Value::DEFAULT {
        return (i, j);
    }

    if j - i < 4 {
        ctx.w[i][j] = Value::NULL;
        return (i, j);
    }

    let w1 = w(seq, i + 1, j, temp, ctx, emap);
    let w2 = w(seq, i, j - 1, temp, ctx, emap);
    let w3 = v(seq, i, j, temp, ctx, emap);

    let mut w4 = Value::NULL;
    for k in i + 1..j - 1 {
        let w4_test = multi_branch(seq, i, k, j, temp, ctx, emap, false);

        if w4_test.valid() && w4_test.e < w4.e {
            w4 = w4_test;
        }
    }

    ctx.w[i][j] = min_value([
        &ctx.w[w1.0][w1.1],
        &ctx.w[w2.0][w2.1],
        &ctx.v[w3.0][w3.1],
        &w4,
    ]);

    (i, j)
}

/// Find, store and return the minimum free energy of the structure between i and j
///
/// If i and j don't bp, store and return INF.
/// See: Figure 2B of Zuker, 1981
///
/// # Args
///
/// - seq: The sequence being folded
/// - i: The start index
/// - j: The end index (inclusive)
/// - temp: The temperature in Kelvin
/// - ctx: Combined free energy cache, etc
/// - emap: Energy map for DNA/RNA
///
/// # Returns
///
/// - [`f64`]: The minimum energy folding structure possible between i and j on seq
pub fn v(
    seq: &[u8],
    i: usize,
    j: usize,
    temp: f64,
    ctx: &mut Context,
    emap: &Energies,
) -> (usize, usize) {
    if ctx.v[i][j] != Value::DEFAULT {
        return (i, j);
    }

    // the ends must basepair for V(i,j)
    if emap.complement[seq[i]] != seq[j] {
        ctx.v[i][j] = Value::NULL;
        return (i, j);
    }

    // if the basepair is isolated, and the seq large, penalize at 1,600 kcal/mol
    // heuristic for speeding this up
    // from https://www.ncbi.nlm.nih.gov/pubmed/10329189
    let isolated_outer = i == 0 || j >= seq.len() - 1 || emap.complement[seq[i - 1]] != seq[j + 1];
    let isolated_inner = emap.complement[seq[i + 1]] != seq[j - 1];

    if isolated_outer && isolated_inner {
        ctx.v[i][j] = Value::from(1600.0);
        return (i, j);
    }

    // E1 = FH(i, j); hairpin
    let pair = calc_pair(seq, i, i + 1, j, j - 1);
    let e1 = Value::empty(hairpin(seq, i, j, temp, emap), Desc::Hairpin(pair));
    if j - i == 4 {
        // small hairpin; 4bp
        ctx.v[i][j] = e1.clone();
        ctx.w[i][j] = e1;
        return (i, j);
    }

    // E2 = min{FL(i, j, i', j') + V(i', j')}, i<i'<j'<j
    // stacking region or bulge or interior loop; Figure 2A(2)
    // j-i=d>4; various pairs i',j' for j'-i'<d
    let n = seq.len();
    let mut e2 = Value::from(f64::INFINITY);
    for i1 in i + 1..j - 4 {
        for j1 in i1 + 4..j {
            // i1 and j1 must match
            if emap.complement[seq[i1]] != seq[j1] {
                continue;
            }

            let pair = calc_pair(seq, i, i1, j, j1);
            let pair_left = calc_pair(seq, i, i + 1, j, j - 1);
            let pair_right = calc_pair(seq, i1 - 1, i1, j1 + 1, j1);
            let pair_inner = emap.nn.contains_key(pair_left) || emap.nn.contains_key(pair_right);

            let stack = i1 == i + 1 && j1 == j - 1;
            let bulge_left = i1 > i + 1;
            let bulge_right = j1 < j - 1;

            let (mut e2_test, mut e2_test_type);
            if stack {
                // it's a neighboring/stacking pair in a helix
                e2_test = calc_stack(seq, i, i1, j, j1, temp, emap);
                e2_test_type = Desc::Stack(pair.to_vec());

                if i > 0 && j == n - 1 || i == 0 && j < n - 1 {
                    // there's a dangling end
                    e2_test_type = Desc::StackDe(pair);
                }
            } else if bulge_left && bulge_right && !pair_inner {
                // it's an interior loop
                e2_test = internal_loop(seq, i, i1, j, j1, temp, emap);
                e2_test_type = Desc::InteriorLoop(i1 - i, j - j1);

                if i1 - i == 2 && j - j1 == 2 {
                    let mut stack = Vec::new();
                    stack.extend(&seq[i..i1 + 1]); // loop_left
                    stack.push(b'/');
                    stack.extend(seq[j1..j + 1].iter().rev().copied()); // loop_right

                    // technically an interior loop of 1. really 1bp mismatch
                    e2_test_type = Desc::Stack(stack);
                }
            } else if bulge_left && !bulge_right {
                // it's a bulge on the left side
                e2_test = bulge(seq, i, i1, j, j1, temp, emap);
                e2_test_type = Desc::Bulge(i1 - i);
            } else if !bulge_left && bulge_right {
                // it's a bulge on the right side
                e2_test = bulge(seq, i, i1, j, j1, temp, emap);
                e2_test_type = Desc::Bulge(j - j1);
            } else {
                // it's basically a hairpin, only outside bp match
                continue;
            }

            let (i, j) = v(seq, i1, j1, temp, ctx, emap);
            e2_test += ctx.v[i][j].e;
            if e2_test != f64::NEG_INFINITY && e2_test < e2.e {
                e2 = Value::new(e2_test, e2_test_type, vec![(i1, j1)]);
            }
        }
    }

    // E3 = min{W(i+1,i') + W(i'+1,j-1)}, i+1<i'<j-2
    let mut e3 = Value::NULL;
    if isolated_outer || i != 0 || j == seq.len() - 1 {
        for k in i + 1..j - 1 {
            let e3_test = multi_branch(seq, i, k, j, temp, ctx, emap, true);

            if e3_test.valid() && e3_test.e < e3.e {
                e3 = e3_test
            }
        }
    }

    ctx.v[i][j] = min_value([&e1, &e2, &e3]);
    (i, j)
}

/// Return a stack representation, a key for the NN maps
///
/// # Args
///
/// - s: Sequence being folded
/// - i: leftmost index
/// - i1: index to right of i
/// - j: rightmost index
/// - j1: index to left of j
///
/// # Returns
///
/// - [`[u8; 5]`]: string representation of the pair
pub fn calc_pair(
    s: &[u8],
    i: impl ToIsize,
    i1: impl ToIsize,
    j: impl ToIsize,
    j1: impl ToIsize,
) -> [u8; 5] {
    let (i, i1, j, j1) = (i.to_isize(), i1.to_isize(), j.to_isize(), j1.to_isize());
    [
        if i >= 0 { s[i as usize] } else { b'.' },
        if i1 >= 0 { s[i1 as usize] } else { b'.' },
        b'/',
        if j >= 0 { s[j as usize] } else { b'.' },
        if j1 >= 0 { s[j1 as usize] } else { b'.' },
    ]
}

/// Return the struct with the lowest free energy that isn't -inf (undef)
///
/// # Args
///
/// - values: Values being compared
///
/// # Returns
///
/// - [`Value`]: The min free energy structure
pub fn min_value<const N: usize>(values: [&Value; N]) -> Value {
    let mut vt = &Value::NULL;
    for v in values {
        if v.e != f64::NEG_INFINITY && v.e < vt.e {
            vt = v;
        }
    }
    vt.clone()
}

/// Find the free energy given delta h, s and temp
///
/// # Args
///
/// - d_h: The enthalpy increment in kcal / mol
/// - d_s: The entropy increment in cal / mol
/// - temp: The temperature in Kelvin
///
/// # Returns
///
/// - The free energy increment in kcal / (mol x K)
pub fn calc_d_g(d_h: f64, d_s: f64, temp: f64) -> f64 {
    d_h - temp * (d_s / 1000.0)
}

/// Estimate the free energy of length query_len based on one of length known_len.
///
/// The Jacobson-Stockmayer entry extrapolation formula is used
/// for bulges, hairpins, etc that fall outside the 30nt upper limit
/// for pre-calculated free-energies. See SantaLucia and Hicks (2004).
///
/// # Args
///
/// - query_len: Length of element without known free energy value
/// - known_len: Length of element with known free energy value (d_g_x)
/// - d_g_x: The free energy of the element known_len
/// - temp: Temperature in Kelvin
///
/// # Returns
///
/// - [`f64`]: The free energy for a structure of length query_len
pub fn calc_j_s(query_len: usize, known_len: usize, d_g_x: f64, temp: f64) -> f64 {
    let gas_constant = 1.9872e-3;
    d_g_x + 2.44 * gas_constant * temp * (query_len as f64 / known_len as f64).ln()
}

/// Get the free energy for a stack.
///
/// Using the indexes i and j, check whether it's at the end of
/// the sequence or internal. Then check whether it's a match
/// or mismatch, and return.
///
/// Two edge-cases are terminal mismatches and dangling ends.
/// The energy of a dangling end is added to the energy of a pair
/// where i XOR j is at the sequence's end.
///
/// # Args
///
/// - seq: The full folding sequence
/// - i: The start index on left side of the pair/stack
/// - i1: The index to the right of i
/// - j: The end index on right side of the pair/stack
/// - j1: The index to the left of j
/// - temp: Temperature in Kelvin
/// - emap: Map of energies
///
/// # Returns
///
/// - [`f64`]: The free energy of the NN pairing
pub fn calc_stack(
    seq: &[u8],
    i: impl ToIsize,
    i1: impl ToIsize,
    j: impl ToIsize,
    j1: impl ToIsize,
    temp: f64,
    emap: &Energies,
) -> f64 {
    let (i, i1, j, j1) = (i.to_isize(), i1.to_isize(), j.to_isize(), j1.to_isize());
    if [i, i1, j, j1].iter().any(|&x| x >= seq.len() as isize) {
        return 0.0;
    }

    let pair = calc_pair(seq, i, i1, j, j1);

    if [i, i1, j, j1].contains(&-1) {
        // it's a dangling end
        let (d_h, d_s) = emap.de[pair];
        return calc_d_g(d_h, d_s, temp);
    }
    let (i, j) = (i as usize, j as usize);

    if i > 0 && j < seq.len() - 1 {
        // it's internal
        let (d_h, d_s) = emap.nn.get(pair).unwrap_or_else(|| emap.internal_mm[pair]);
        // d_h, d_s = emap.NN[pair] if pair in emap.NN else emap.INTERNAL_MM[pair]
        return calc_d_g(d_h, d_s, temp);
    }

    let (d_h, d_s) = emap.nn.get(pair).unwrap_or_else(|| emap.terminal_mm[pair]);
    let mut d_g = calc_d_g(d_h, d_s, temp);

    if i == 0 && j == seq.len() - 1 {
        // it's terminal
        d_g
    } else if i > 0 && j == seq.len() - 1 {
        // it's dangling on left
        let pair_de = [seq[i - 1], seq[i], b'/', b'.', seq[j]];
        if let Some((d_h, d_s)) = emap.de.get(pair_de) {
            d_g += calc_d_g(d_h, d_s, temp);
        }
        d_g
    } else if i == 0 && j < seq.len() - 1 {
        // it's dangling on right
        let pair_de = [b'.', seq[i], b'/', seq[j + 1], seq[j]];
        if let Some((d_h, d_s)) = emap.de.get(pair_de) {
            d_g += calc_d_g(d_h, d_s, temp);
        }
        d_g
    } else {
        0.0
    }
}

/// Calculate the free energy of a hairpin.
///
/// # Args
///
/// - seq: The sequence we're folding
/// - i: The index of start of hairpin
/// - j: The index of end of hairpin
/// - temp: Temperature in Kelvin
/// - emap: Map of energies
///
/// # Returns
///
/// - [`f64`]: The free energy increment from the hairpin structure
pub fn hairpin(seq: &[u8], i: usize, j: usize, temp: f64, emap: &Energies) -> f64 {
    if j - i < 4 {
        return f64::INFINITY;
    }

    let hairpin = &seq[i..=j];
    let hairpin_len = hairpin.len() - 2;
    let pair = calc_pair(seq, i, i + 1, j, j - 1);

    assert_eq!(emap.complement[hairpin[0]], hairpin[hairpin.len() - 1]);

    let mut d_g = emap
        .tri_tetra_loops
        .as_ref()
        .and_then(|ttl| ttl.get(<[u8; 6]>::try_from(hairpin).ok()?))
        .map_or(0.0, |(d_h, d_s)| calc_d_g(d_h, d_s, temp));

    // add penalty based on size
    d_g += emap.hairpin_loops.get(hairpin_len, temp);

    // add penalty for a terminal mismatch
    if hairpin_len > 3
        && let Some((d_h, d_s)) = emap.terminal_mm.get(pair)
    {
        d_g += calc_d_g(d_h, d_s, temp);
    }

    // add penalty if length 3 and AT closing, formula 8 from SantaLucia, 2004
    if hairpin_len == 3 && (hairpin[0] == b'A' || hairpin[hairpin.len() - 1] == b'A') {
        d_g += 0.5; //  convert to entropy
    }

    d_g
}

/// Calculate the free energy associated with a bulge.
///
/// # Args
///
/// - seq: The full folding DNA sequence
/// - i: The start index of the bulge
/// - i1: The index to the right of i
/// - j: The end index of the bulge
/// - j1: The index to the left of j
/// - loop: The sequence of the bulge
/// - temp: Temperature in Kelvin
/// - emap: Map to DNA/RNA energies
///
/// # Returns
///
/// - [`f64`]: The increment in free energy from the bulge
pub fn bulge(
    seq: &[u8],
    i: usize,
    i1: usize,
    j: usize,
    j1: usize,
    temp: f64,
    emap: &Energies,
) -> f64 {
    let loop_len = (i1 - i - 1).max(j - j1 - 1);
    assert!(loop_len > 0, "Loop len can't be <= 0: {loop_len}");

    // add penalty based on size
    let mut d_g = emap.bulge_loops.get(loop_len, temp);

    if loop_len == 1 {
        // if len 1, include the delta G of intervening NN (SantaLucia 2004)
        let pair = calc_pair(seq, i, i1, j, j1);
        assert!(emap.nn.contains_key(pair));
        d_g += calc_stack(seq, i, i1, j, j1, temp, emap);
    }

    // penalize AT terminal bonds
    if [i, i1, j, j1].iter().any(|&k| seq[k] == b'A') {
        d_g += 0.5;
    }

    d_g
}

/// Calculate the free energy of an internal loop.
///
/// The first and last bp of both left and right sequences
/// are not themselves parts of the loop, but are the terminal
/// bp on either side of it. They are needed for when there's
/// a single internal looping bp (where just the mismatching
/// free energies are used)
///
/// Note that both left and right sequences are in 5' to 3' direction
///
/// This is adapted from the "Internal Loops" section of SantaLucia/Hicks, 2004
///
/// # Args
///
/// - seq: The sequence we're folding
/// - i: The index of the start of structure on left side
/// - i1: The index to the right of i
/// - j: The index of the end of structure on right side
/// - j1: The index to the left of j
/// - temp: Temperature in Kelvin
/// - emap: Dictionary mapping to energies for DNA/RNA
///
/// # Returns
///
/// - [`f64`]: The free energy associated with the internal loop
pub fn internal_loop(
    seq: &[u8],
    i: usize,
    i1: usize,
    j: usize,
    j1: usize,
    temp: f64,
    emap: &Energies,
) -> f64 {
    let loop_left = i1 - i - 1;
    let loop_right = j - j1 - 1;
    let loop_len = loop_left + loop_right;

    assert!(loop_left >= 1 && loop_right >= 1);

    // single bp mismatch, sum up the two single mismatch pairs
    if loop_left == 1 && loop_right == 1 {
        let mm_left = calc_stack(seq, i, i1, j, j1, temp, emap);
        let mm_right = calc_stack(seq, i1 - 1, i1, j1 + 1, j1, temp, emap);
        return mm_left + mm_right;
    }

    // apply a penalty based on loop size
    let mut d_g = emap.internal_loops.get(loop_len, temp);

    // apply an asymmetry penalty
    let loop_asymmetry = loop_left.abs_diff(loop_right);
    d_g += 0.3 * loop_asymmetry as f64;

    // apply penalty based on the mismatching pairs on either side of the loop
    let pair_left_mm = calc_pair(seq, i, i + 1, j, j - 1);
    let (d_h, d_s) = emap.terminal_mm[pair_left_mm];
    d_g += calc_d_g(d_h, d_s, temp);

    let pair_right_mm = calc_pair(seq, i1 - 1, i1, j1 + 1, j1);
    let (d_h, d_s) = emap.terminal_mm[pair_right_mm];
    d_g += calc_d_g(d_h, d_s, temp);

    d_g
}

pub fn add_branch(
    (i, j): (usize, usize),
    seq: &[u8],
    temp: f64,
    ctx: &mut Context,
    emap: &Energies,
    branches: &mut Indices,
) {
    let this = &ctx.w[i][j];

    if !this.valid() || this.ij.is_empty() {
        return;
    }

    if this.ij.len() == 1 {
        branches.push(this.ij[0]);
        return;
    }

    for (i1, j1) in this.ij.clone() {
        let index = w(seq, i1, j1, temp, ctx, emap);
        add_branch(index, seq, temp, ctx, emap, branches);
    }
}

/// Calculate a multi-branch energy penalty using a linear formula.
///
/// From Jaeger, Turner, and Zuker, 1989.
/// Found to be better than logarithmic in Ward, et al. 2017
///
/// # Args
///
/// - seq: The sequence being folded
/// - i: The left starting index
/// - k: The mid-point in the search
/// - j: The right ending index
/// - temp: Folding temp
/// - ctx: Combined free energy cache, etc
/// - emap: Map to DNA/RNA energies
/// - helix: Whether this multibranch is enclosed by a helix
///
/// # Returns
///
/// - [`Value`]: A multi-branch structure
pub fn multi_branch(
    seq: &[u8],
    i: usize,
    k: usize,
    j: usize,
    temp: f64,
    ctx: &mut Context,
    emap: &Energies,
    helix: bool,
) -> Value {
    let (li, rj) = if helix { (i + 1, j - 1) } else { (i, j) };
    let (li, lj) = w(seq, li, k, temp, ctx, emap);
    let (ri, rj) = w(seq, k + 1, rj, temp, ctx, emap);

    if !ctx.w[li][lj].valid() || !ctx.w[ri][rj].valid() {
        return Value::NULL;
    }

    // gather all branches of this multi-branch structure
    let mut branches = Indices::new();

    add_branch((li, lj), seq, temp, ctx, emap, &mut branches);
    add_branch((ri, rj), seq, temp, ctx, emap, &mut branches);

    // this isn't multi-branched
    if branches.len() < 2 {
        return Value::NULL;
    }

    // if there's a helix, i,j counts as well
    if helix {
        branches.push((i, j));
    }

    // count up unpaired bp and asymmetry
    let branches_count = branches.len();
    let mut unpaired = 0;
    let mut e_sum = 0.0;
    for (index, &(i2, j2)) in branches.iter().enumerate() {
        let (_, j1) = branches[(index as isize - 1).rem_euclid(branches.len() as isize) as usize];
        let (i3, j3) = branches[(index as isize + 1).rem_euclid(branches.len() as isize) as usize];

        // add energy from unpaired bp to the right
        // of the helix as though it was a dangling end
        // if there's only one bp, it goes to whichever
        // helix (this or the next) has the more favorable energy
        let unpaired_left;
        let mut unpaired_right = 0;
        let mut de = 0.0;
        let [i2, i3, j1, j2, j3] = [i2, i3, j1, j2, j3].map(|n| n as isize);

        if index == branches.len() - 1 && !helix {
            // do nothing
        } else if (i3, j3) == (i as isize, j as isize) {
            unpaired_left = i2 - j1 - 1;
            unpaired_right = j3 - j2 - 1;

            if unpaired_left != 0 && unpaired_right != 0 {
                de = calc_stack(seq, i2 - 1, i2, j2 + 1, j2, temp, emap);
            } else if unpaired_right != 0 {
                de = calc_stack(seq, -1, i2, j2 + 1, j2, temp, emap);
                if unpaired_right == 1 {
                    de = calc_stack(seq, i3, -1, j3, j3 - 1, temp, emap).min(de);
                }
            }
        } else if (i2, j2) == (i as isize, j as isize) {
            unpaired_left = j2 - j1 - 1;
            unpaired_right = i3 - i2 - 1;

            if unpaired_left != 0 && unpaired_right != 0 {
                de = calc_stack(seq, i2 - 1, i2, j2 + 1, j2, temp, emap);
            } else if unpaired_right != 0 {
                de = calc_stack(seq, i2, i2 + 1, j2, -1, temp, emap);
                if unpaired_right == 1 {
                    de = calc_stack(seq, i3 - 1, i3, -1, j3, temp, emap).min(de);
                }
            }
        } else {
            unpaired_left = i2 - j1 - 1;
            unpaired_right = i3 - j2 - 1;

            if unpaired_left != 0 && unpaired_right != 0 {
                de = calc_stack(seq, i2 - 1, i2, j2 + 1, j2, temp, emap);
            } else if unpaired_right != 0 {
                de = calc_stack(seq, -1, i2, j2 + 1, j2, temp, emap);
                if unpaired_right == 1 {
                    de = calc_stack(seq, i2 - 1, i2, j2 + 1, j2, temp, emap).min(de);
                }
            }
        }

        e_sum += de;
        unpaired += unpaired_right;
        assert!(unpaired_right >= 0);

        if (i2, j2) != (i as isize, j as isize) {
            // add energy
            let (i, j) = w(seq, i2 as usize, j2 as usize, temp, ctx, emap);
            e_sum += ctx.w[i][j].e;
        }
    }

    assert!(unpaired >= 0);

    // penalty for unmatched bp and multi-branch
    let (a, b, c, d) = emap.multibranch;
    let e_multibranch = if unpaired == 0 {
        a + d
    } else {
        a + b * branches.len() as f64 + c * unpaired as f64
    };

    // energy of min-energy neighbors
    let e = e_multibranch + e_sum;

    // pointer to next structures
    if helix {
        branches.pop();
    }

    Value::new(
        e,
        Desc::Bifurcation(unpaired as usize, branches_count),
        branches,
    )
}

/// Traceback thru the V(i,j) and W(i,j) caches to find the structure
///
/// For each step, get to the lowest energy W(i,j) within that block
/// Store the structure in W(i,j)
/// Inc i and j
/// If the next structure is viable according to V(i,j), store as well
/// Repeat
///
/// # Args
///
/// - i: The leftmost index to start searching in
/// - j: The rightmost index to start searching in
/// - ctx: Combined free energy cache, etc
///
/// # Returns
///
/// - A list of Values in the final secondary structure
pub fn traceback(mut i: usize, mut j: usize, ctx: &Context) -> Vec<Value> {
    // move i,j down-left to start coordinates
    let s_w = &ctx.w[i][j];
    if !matches!(s_w.desc, Desc::Hairpin(_)) {
        while &ctx.w[i + 1][j] == s_w {
            i += 1;
        }
        while &ctx.w[i][j - 1] == s_w {
            j -= 1;
        }
    }

    let mut values = Vec::new();
    loop {
        // multibrach structures are only in the w_cache.
        let s_v = if s_w.ij.len() > 1 { s_w } else { &ctx.v[i][j] };

        values.push(s_v.clone().with_ij(vec![(i, j)]));

        // it's a multibranch
        if s_v.ij.len() > 1 {
            let mut e_sum = 0.0;
            let mut values = trackback_energy(&values);
            let last_index = values.len() - 1;

            let mut branches = Vec::new();
            for &(i1, j1) in &s_v.ij {
                let tb = traceback(i1, j1, ctx);

                if let Some(&(i2, j2)) = tb.first().and_then(|v| v.ij.first()) {
                    e_sum += ctx.w[i2][j2].e;
                    branches.extend(tb);
                }
            }

            values[last_index].e = round1(values[last_index].e - e_sum);
            values.extend(branches);

            break values;
        }

        // it's a stack, bulge, etc
        // there's another single structure beyond this
        if s_v.ij.len() == 1 {
            (i, j) = s_v.ij[0];
        } else {
            // it's a hairpin, end of structure
            // set the energy of everything relative to the hairpin
            break trackback_energy(&values);
        }
    }
}

/// Add energy to each structure, based on how it's W(i,j) differs from the one after
///
/// # Args
///
/// - vals: The structures for whom energy is being calculated
///
/// # Returns
///
/// - [`Vec<Value>`]: Values in the folded DNA with energy
pub fn trackback_energy(vals: &[Value]) -> Vec<Value> {
    vals.iter()
        .enumerate()
        .map(|(index, v)| {
            let l = vals.len() - 1;
            let e_next = if index == l { 0.0 } else { vals[index + 1].e };
            let e_corrected = round1(v.e - e_next);
            Value::new(e_corrected, v.desc.clone(), v.ij.clone())
        })
        .collect()
}
