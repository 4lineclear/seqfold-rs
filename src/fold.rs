//! Predict nucleic acid secondary structure

use std::{collections::HashSet, fmt::Display};

use crate::{
    Cache, Energies,
    util::{ByteStr, Matrix, ToIsize, round1, round2},
};

#[cfg(test)]
mod test;

/// A single structure with a free energy, description, and inward children.
#[derive(Debug, Clone)]
pub struct Value {
    pub e: f64,
    pub desc: Desc,
    pub ij: Vec<(usize, usize)>,
}

#[derive(Clone, PartialEq, Eq)]
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

/// A map from i, j tuple to a min free energy Struct.
pub type Values = Vec<Vec<Value>>;

impl From<f64> for Value {
    fn from(value: f64) -> Self {
        Self::new(value, Desc::Empty, Vec::new())
    }
}

impl Value {
    pub const NULL: Self = Self::empty(f64::INFINITY, Desc::Empty);
    pub const DEFAULT: Self = Self::empty(f64::NEG_INFINITY, Desc::Empty);

    pub const fn empty(e: f64, desc: Desc) -> Value {
        Self::new(e, desc, Vec::new())
    }

    pub const fn new(e: f64, desc: Desc, ij: Vec<(usize, usize)>) -> Self {
        Self { e, desc, ij }
    }

    pub fn with_ij(mut self, ij: Vec<(usize, usize)>) -> Self {
        self.ij = ij;
        self
    }

    pub fn valid(&self) -> bool {
        self.e != f64::INFINITY && self.e != f64::NEG_INFINITY
    }
}

impl PartialEq for Value {
    fn eq(&self, other: &Self) -> bool {
        self.e == other.e && self.ij == other.ij
    }
}

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
    let (v_cache, w_cache) = cache(seq, temp);
    let n = seq.len();

    // get the minimum free energy structure out of the cache
    traceback(0, n - 1, &v_cache, &w_cache)
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
    let values: Vec<Value> = fold(seq, temp);
    let dg_sum = values.iter().map(|v| v.e).sum();
    round2(dg_sum)
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
///              minimum free energy between i and j
pub fn dg_cache(seq: &[u8], temp: Option<f64>) -> Cache {
    cache(seq, temp)
        .1
        .into_iter()
        .map(|row| row.into_iter().map(|v| v.e).collect())
        .collect()
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
pub fn cache(seq: &[u8], temp: Option<f64>) -> (Values, Values) {
    let temp = temp.unwrap_or(37.0);

    let seq = seq.to_ascii_uppercase();
    let temp = temp + 273.15; // kelvin

    // figure out whether it's DNA or RNA, choose energy map
    let mut dna = true;
    let bps = seq.iter().copied().collect::<HashSet<_>>();
    if bps.contains(&b'U') && bps.contains(&b'T') {
        panic!("Both T and U in sequence. Provide one or the other for DNA OR RNA.");
    }
    if bps.iter().all(|b| b"AUCG".contains(b)) {
        dna = false;
    } else if bps.iter().any(|b| !b"ATGC".contains(b)) {
        panic!(
            "Unknown bp: {}. Only DNA/RNA foldable",
            bps.iter()
                .filter(|b| !b"ATGC".contains(b))
                .map(|b| *b as char)
                .collect::<String>()
        );
    }
    let emap = if dna { crate::dna() } else { crate::rna() };

    let n = seq.len();
    let mut v_cache = Values::new();
    let mut w_cache = Values::new();
    for _ in 0..n {
        v_cache.push(vec![Value::DEFAULT; n]);
        w_cache.push(vec![Value::DEFAULT; n]);
    }

    // fill the cache
    w(&seq, 0, n - 1, temp, &mut v_cache, &mut w_cache, emap);

    (v_cache, w_cache)
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
/// - v_cache: Free energy cache for if i and j bp
/// - w_cache: Free energy cache for lowest energy structure from i to j. 0 otherwise
///
/// # Returns
///
/// - [`f64`]: The free energy for the subsequence from i to j
pub fn w(
    seq: &[u8],
    i: usize,
    j: usize,
    temp: f64,
    v_cache: &mut Values,
    w_cache: &mut Values,
    emap: &Energies,
) -> (usize, usize) {
    if w_cache[i][j] != Value::DEFAULT {
        return (i, j);
    }

    if j - i < 4 {
        w_cache[i][j] = Value::NULL;
        return (i, j);
    }

    let w1 = w(seq, i + 1, j, temp, v_cache, w_cache, emap);
    let w2 = w(seq, i, j - 1, temp, v_cache, w_cache, emap);
    let w3 = v(seq, i, j, temp, v_cache, w_cache, emap);

    let mut w4 = Value::NULL;
    for k in i + 1..j - 1 {
        let w4_test = multi_branch(seq, i, k, j, temp, v_cache, w_cache, emap, Some(false));

        if w4_test.valid() && w4_test.e < w4.e {
            w4 = w4_test;
        }
    }

    let w = Matrix(&*w_cache);
    let v = Matrix(&*v_cache);
    let min_value = min_value([&w[w1], &w[w2], &v[w3], &w4]);
    w_cache[i][j] = min_value;
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
/// - v_cache: Free energy cache for if i and j bp. INF otherwise
/// - w_cache: Free energy cache for lowest energy structure from i to j. 0 otherwise
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
    v_cache: &mut Values,
    w_cache: &mut Values,
    emap: &Energies,
) -> (usize, usize) {
    if v_cache[i][j] != Value::DEFAULT {
        return (i, j);
    }

    // the ends must basepair for V(i,j)
    if emap.complement[&seq[i]] != seq[j] {
        v_cache[i][j] = Value::NULL;
        return (i, j);
    }

    // if the basepair is isolated, and the seq large, penalize at 1,600 kcal/mol
    // heuristic for speeding this up
    // from https://www.ncbi.nlm.nih.gov/pubmed/10329189
    let mut isolated_outer = true;
    if i != 0 && j < seq.len() - 1 {
        isolated_outer = emap.complement[&seq[i - 1]] != seq[j + 1];
    }
    let isolated_inner = emap.complement[&seq[i + 1]] != seq[j - 1];

    if isolated_outer && isolated_inner {
        v_cache[i][j] = Value::from(1600.0);
        return (i, j);
    }

    // E1 = FH(i, j); hairpin
    let pair = calc_pair(seq, i, i + 1, j, j - 1);
    let e1 = Value::empty(hairpin(seq, i, j, temp, emap), Desc::Hairpin(pair));
    if j - i == 4 {
        // small hairpin; 4bp
        v_cache[i][j] = e1.clone();
        w_cache[i][j] = e1;
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
            if emap.complement[&seq[i1]] != seq[j1] {
                continue;
            }

            let pair = calc_pair(seq, i, i1, j, j1);
            let pair_left = calc_pair(seq, i, i + 1, j, j - 1);
            let pair_right = calc_pair(seq, i1 - 1, i1, j1 + 1, j1);
            let pair_left = pair_left.as_slice();
            let pair_right = pair_right.as_slice();
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
                    stack.extend((&seq[j1..j + 1]).iter().rev().copied()); // loop_right

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

            let (i, j) = v(seq, i1, j1, temp, v_cache, w_cache, emap);
            e2_test += v_cache[i][j].e;
            if e2_test != f64::NEG_INFINITY && e2_test < e2.e {
                e2 = Value::new(e2_test, e2_test_type, vec![(i1, j1)]);
            }
        }
    }

    // E3 = min{W(i+1,i') + W(i'+1,j-1)}, i+1<i'<j-2
    let mut e3 = Value::NULL;
    if isolated_outer || i != 0 || j == seq.len() - 1 {
        for k in i + 1..j - 1 {
            let e3_test = multi_branch(seq, i, k, j, temp, v_cache, w_cache, emap, Some(true));

            if e3_test.valid() && e3_test.e < e3.e {
                e3 = e3_test
            }
        }
    }

    v_cache[i][j] = min_value([&e1, &e2, &e3]);
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
pub fn min_value<'a>(values: impl IntoIterator<Item = &'a Value>) -> Value {
    let mut value = Value::NULL;
    for v in values {
        if v.e != f64::NEG_INFINITY && v.e < value.e {
            value = v.clone();
        }
    }
    value
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
fn calc_j_s(query_len: usize, known_len: usize, d_g_x: f64, temp: f64) -> f64 {
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
///
/// # Returns
///
/// - [`f64`]: The free energy of the NN pairing
fn calc_stack(
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

    if [i, i1, j, j1].iter().any(|&x| x as isize == -1) {
        // it's a dangling end
        let (d_h, d_s) = emap.de[pair.as_slice()];
        return calc_d_g(d_h, d_s, temp);
    }
    let (i, j) = (i as usize, j as usize);

    if i > 0 && j < seq.len() - 1 {
        // it's internal
        let &(d_h, d_s) = emap
            .nn
            .get(pair.as_slice())
            .unwrap_or_else(|| &emap.internal_mm[pair.as_slice()]);
        // d_h, d_s = emap.NN[pair] if pair in emap.NN else emap.INTERNAL_MM[pair]
        return calc_d_g(d_h, d_s, temp);
    }

    if i == 0 && j == seq.len() - 1 {
        // it's terminal
        let &(d_h, d_s) = emap
            .nn
            .get(pair.as_slice())
            .unwrap_or_else(|| &emap.terminal_mm[pair.as_slice()]);
        return calc_d_g(d_h, d_s, temp);
    }

    if i > 0 && j == seq.len() - 1 {
        // it's dangling on left
        let &(d_h, d_s) = emap
            .nn
            .get(pair.as_slice())
            .unwrap_or_else(|| &emap.terminal_mm[pair.as_slice()]);
        let mut d_g = calc_d_g(d_h, d_s, temp);

        let pair_de = [seq[i - 1], seq[i], b'/', b'.', seq[j]];
        if emap.de.contains_key(pair_de.as_slice()) {
            let (d_h, d_s) = emap.de[pair_de.as_slice()];
            d_g += calc_d_g(d_h, d_s, temp);
        }
        return d_g;
    }

    if i == 0 && j < seq.len() - 1 {
        // it's dangling on right
        let &(d_h, d_s) = emap
            .nn
            .get(pair.as_slice())
            .unwrap_or_else(|| &emap.terminal_mm[pair.as_slice()]);
        let mut d_g = calc_d_g(d_h, d_s, temp);

        let pair_de = [b'.', seq[i], b'/', seq[j + 1], seq[j]];
        if emap.de.contains_key(pair_de.as_slice()) {
            let (d_h, d_s) = emap.de[pair_de.as_slice()];
            d_g += calc_d_g(d_h, d_s, temp);
        }
        return d_g;
    }

    0.0
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

    let hairpin = &seq[i..j + 1];
    let hairpin_len = hairpin.len() - 2;
    let pair = calc_pair(seq, i, i + 1, j, j - 1);

    if emap.complement[&hairpin[0]] != hairpin[hairpin.len() - 1] {
        // not known terminal pair, nothing to close "hairpin"
        panic!()
    }

    let mut d_g = emap
        .tri_tetra_loops
        .as_ref()
        .and_then(|ttl| ttl.get(hairpin))
        .map_or(0.0, |&(d_h, d_s)| calc_d_g(d_h, d_s, temp));

    // add penalty based on size
    if emap.hairpin_loops.contains_key(&hairpin_len) {
        let (d_h, d_s) = emap.hairpin_loops[&hairpin_len];
        d_g += calc_d_g(d_h, d_s, temp);
    } else {
        // it's too large, extrapolate
        let (d_h, d_s) = emap.hairpin_loops[&30];
        let d_g_inc = calc_d_g(d_h, d_s, temp);
        d_g += calc_j_s(hairpin_len, 30, d_g_inc, temp);
    }

    // add penalty for a terminal mismatch
    if hairpin_len > 3 && emap.terminal_mm.contains_key(pair.as_slice()) {
        if emap.terminal_mm.contains_key(pair.as_slice()) {
            let (d_h, d_s) = emap.terminal_mm[pair.as_slice()];
            d_g += calc_d_g(d_h, d_s, temp);
        }
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
    i: impl ToIsize,
    i1: impl ToIsize,
    j: impl ToIsize,
    j1: impl ToIsize,
    temp: f64,
    emap: &Energies,
) -> f64 {
    let (i, i1, j, j1) = (i.to_isize(), i1.to_isize(), j.to_isize(), j1.to_isize());
    let loop_len = (i1 - i - 1).max(j - j1 - 1);
    assert!(loop_len > 0, "Loop len can't be <= 0: {loop_len}");
    let loop_len = loop_len as usize;

    // add penalty based on size
    let mut d_g = if emap.bulge_loops.contains_key(&loop_len) {
        let (d_h, d_s) = emap.bulge_loops[&loop_len];
        calc_d_g(d_h, d_s, temp)
    } else {
        // it's too large for pre-calculated list, extrapolate
        let (d_h, d_s) = emap.bulge_loops[&30];
        let d_g = calc_d_g(d_h, d_s, temp);
        calc_j_s(loop_len, 30, d_g, temp)
    };

    if loop_len == 1 {
        // if len 1, include the delta G of intervening NN (SantaLucia 2004)
        let pair = calc_pair(seq, i, i1, j, j1);
        assert!(emap.nn.contains_key(pair.as_slice()));
        d_g += calc_stack(seq, i, i1, j, j1, temp, emap);
    }

    // penalize AT terminal bonds
    if [i, i1, j, j1].iter().any(|&k| seq[k as usize] == b'A') {
        d_g += 0.5
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
fn internal_loop(
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

    if loop_left < 1 || loop_right < 1 {
        panic!();
    }

    // single bp mismatch, sum up the two single mismatch pairs
    if loop_left == 1 && loop_right == 1 {
        let mm_left = calc_stack(seq, i, i1, j, j1, temp, emap);
        let mm_right = calc_stack(seq, i1 - 1, i1, j1 + 1, j1, temp, emap);
        return mm_left + mm_right;
    }

    // apply a penalty based on loop size
    let mut d_g = if emap.internal_loops.contains_key(&loop_len) {
        let (d_h, d_s) = emap.internal_loops[&loop_len];
        calc_d_g(d_h, d_s, temp)
    } else {
        // it's too large an internal loop, extrapolate
        let (d_h, d_s) = emap.internal_loops[&30];
        let d_g = calc_d_g(d_h, d_s, temp);
        calc_j_s(loop_len, 30, d_g, temp)
    };

    // apply an asymmetry penalty
    let loop_asymmetry = loop_left.abs_diff(loop_right);
    d_g += 0.3 * loop_asymmetry as f64;

    // apply penalty based on the mismatching pairs on either side of the loop
    let pair_left_mm = calc_pair(seq, i, i + 1, j, j - 1);
    let (d_h, d_s) = emap.terminal_mm[pair_left_mm.as_slice()];
    d_g += calc_d_g(d_h, d_s, temp);

    let pair_right_mm = calc_pair(seq, i1 - 1, i1, j1 + 1, j1);
    let (d_h, d_s) = emap.terminal_mm[pair_right_mm.as_slice()];
    d_g += calc_d_g(d_h, d_s, temp);

    d_g
}

pub fn add_branch(
    (i, j): (usize, usize),
    seq: &[u8],
    temp: f64,
    v_cache: &mut Values,
    w_cache: &mut Values,
    emap: &Energies,
    branches: &mut Vec<(usize, usize)>,
) {
    let this = &w_cache[i][j];

    if !this.valid() || this.ij.is_empty() {
        return;
    }

    if this.ij.len() == 1 {
        branches.push(this.ij[0]);
        return;
    }

    for (i1, j1) in this.ij.clone() {
        let index = w(seq, i1, j1, temp, v_cache, w_cache, emap);
        add_branch(index, seq, temp, v_cache, w_cache, emap, branches);
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
/// - v_cache: Values of energies where V(i,j) bond
/// - w_cache: Values of min energy of substructures between W(i,j)
/// - helix: Whether this multibranch is enclosed by a helix
/// - emap: Map to DNA/RNA energies
/// - helix: Whether V(i, j) bond with one another in a helix
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
    v_cache: &mut Values,
    w_cache: &mut Values,
    emap: &Energies,
    helix: Option<bool>,
) -> Value {
    let helix = helix.unwrap_or(false);
    let (li, lj, ri, rj) = if helix {
        (i + 1, k, k + 1, j - 1)
    } else {
        (i, k, k + 1, j)
    };
    let (li, lj) = w(seq, li, lj, temp, v_cache, w_cache, emap);
    let (ri, rj) = w(seq, ri, rj, temp, v_cache, w_cache, emap);

    let left = &w_cache[li][lj];
    let right = &w_cache[ri][rj];

    if !left.valid() || !right.valid() {
        return Value::NULL;
    }

    // gather all branches of this multi-branch structure
    let mut branches = Vec::new();

    add_branch((li, lj), seq, temp, v_cache, w_cache, emap, &mut branches);
    add_branch((ri, rj), seq, temp, v_cache, w_cache, emap, &mut branches);

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
            (); // do nothing
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
            let (i, j) = w(seq, i2 as usize, j2 as usize, temp, v_cache, w_cache, emap);
            e_sum += w_cache[i][j].e;
        }
    }

    assert!(unpaired >= 0);

    // penalty for unmatched bp and multi-branch
    let (a, b, c, d) = emap.multibranch;
    let mut e_multibranch = a + b * branches.len() as f64 + c * unpaired as f64;

    if unpaired == 0 {
        e_multibranch = a + d;
    }

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
/// - v_cache: Energies where i and j bond
/// - w_cache: Energies/sub-structures between or with i and j
///
/// # Returns
///
/// - A list of Values in the final secondary structure
pub fn traceback(mut i: usize, mut j: usize, v_cache: &Values, w_cache: &Values) -> Vec<Value> {
    // move i,j down-left to start coordinates
    let s_w = &w_cache[i][j];
    if !matches!(s_w.desc, Desc::Hairpin(_)) {
        while &w_cache[i + 1][j] == s_w {
            i += 1
        }
        while &w_cache[i][j - 1] == s_w {
            j -= 1
        }
    }

    let mut values = Vec::new();
    loop {
        // multibrach structures are only in the w_cache.
        let s_v = if s_w.ij.len() > 1 {
            s_w
        } else {
            &v_cache[i][j]
        };

        values.push(s_v.clone().with_ij(vec![(i, j)]));

        // it's a multibranch
        if s_v.ij.len() > 1 {
            let mut e_sum = 0.0;
            let mut values = trackback_energy(&values);
            let last_index = values.len() - 1;

            let mut branches = Vec::new();
            for &(i1, j1) in s_v.ij.iter() {
                let tb = traceback(i1, j1, &v_cache, &w_cache);

                if let Some(&(i2, j2)) = tb.get(0).and_then(|v| v.ij.first()) {
                    e_sum += w_cache[i2][j2].e;
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
            continue;
        }

        // it's a hairpin, end of structure
        // set the energy of everything relative to the hairpin
        break trackback_energy(&values);
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
