//! Calculate the tm of a DNA sequence

use crate::{dna::dna, types::Cache, util::round1};

/// Calculate the annealing temperature between seq1 and seq2.
///
/// If seq2 is not provided, its exact complement is used.
/// In otherwords, it's assumed to be an exact match. This tm
/// calculate does not account for pseudoknots or anything other
/// than an exact, unpadded alignment between seq1 and seq2.
///
/// This is largley influenced by Bio.SeqUtils.MeltingTemp with
/// some different defaults. Here, the reaction mixture is assumed to
/// be PCR and concentrations for Mg, Tris, K, and dNTPs are included
/// that match a typical PCR reaction according to Thermo and NEB. Additionally,
/// the salt correction formula from IDT's Owczarzy et al. (2008) is used.
///
/// NEB: <https://www.neb.com/tools-and-resources/usage-guidelines/guidelines-for-pcr-optimization-with-taq-dna-polymerase>
/// ThermoFisher: <https://www.thermofisher.com/order/catalog/product/18067017?SID=srch-srp-18067017>
///
/// NOTE: Sequences are assumed not to be symmetrical. Oligo not binding to self.
///
/// # Args
///
/// - seq1: The seq whose tm is calculated
/// - seq2: The seq that seq1 anneals to in 3' -> 5' direction
/// - pcr: Whether tm is being calculated for the oligo is in a
///     - PCR reaction mixture. If so, ion and Tris concentrations
///     - that match a typical NEB/ThermoFisher PCR mixture are used
///
/// # Returns
///
/// - [`f64`]: The estimated tm as a float
pub fn tm(seq1: &[u8], seq2: Option<&[u8]>, pcr: Option<bool>) -> f64 {
    let (seq1, seq2) = parse_input(seq1, seq2);

    // sum enthalpy and entropy. Enthaply is first value of each tuple and
    // entropy is the second value of each tuple in:
    // SantaLucia & Hicks (2004), Annu. Rev. Biophys. Biomol. Struct 33: 415-440

    // start with initiation enthalpy and entropy
    let (mut dh, mut ds) = dna().nn[b"init".as_slice()];

    // add in initial A/T and initial G/Cs
    let init = [seq1[0], seq1[seq1.len() - 1]];
    let [a, t, g, c] = counts(&init, [b'A', b'T', b'G', b'C']);
    let init_at = (a + t) as f64;
    let init_gc = (g + c) as f64;
    let (init_at_h, init_at_s) = dna().nn[b"init_A/T".as_slice()];
    let (init_gc_h, init_gc_s) = dna().nn[b"init_G/C".as_slice()];
    dh += init_at * init_at_h + init_gc * init_gc_h;
    ds += init_at * init_at_s + init_gc * init_gc_s;

    // work through each nearest neighbor pair
    for i in 0..seq1.len() - 1 {
        let pair = [seq1[i], seq1[i + 1], b'/', seq2[i], seq2[i + 1]];
        let pair = pair.as_slice();

        // assuming internal neighbor pair
        let (mut pair_dh, mut pair_ds) = (0.0, 0.0);
        if dna().nn.contains_key(pair) {
            (pair_dh, pair_ds) = dna().nn[pair];
        } else if dna().internal_mm.contains_key(pair) {
            (pair_dh, pair_ds) = dna().internal_mm[pair];
        }

        // overwrite if it's a terminal pair
        if [0, seq1.len() - 2].contains(&i) && dna().terminal_mm.contains_key(pair) {
            (pair_dh, pair_ds) = dna().terminal_mm[pair];
        }

        dh += pair_dh;
        ds += pair_ds;
    }

    let gc = gc(&seq1);
    calc_tm(dh, ds, pcr.unwrap_or(true), gc, seq1.len())
}

/// Return a TmCache where each (i, j) returns the Tm for that subspan.
///
/// 1. Build up the 2D matrixes for the tm calculation:
///     - dh
///     - ds
///     - gc count
/// 2. Fill in each cell, (i, j), with the estimated tm for that range from i to j
///
/// # Args
///
/// - seq1: The seq whose tm is calculated
/// - seq2: The seq that seq1 anneals to in 3' -> 5' direction
/// - pcr:  Whether tm is being calculated for the oligo is in a
///         PCR reaction mixture. If so, ion and Tris concentrations
///         that match a typical NEB/ThermoFisher PCR mixture are used
///
/// # Returns
///
/// - [`Cache`]: Where querying the cache with (i, j) returns the tm of the
///              subsequence starting with i and ending with j, inclusive
pub fn tm_cache(seq1: &[u8], seq2: Option<&[u8]>, pcr: Option<bool>) -> Cache {
    let pcr = pcr.unwrap_or(true);
    let (seq1, seq2) = parse_input(seq1, seq2);
    let n = seq1.len(); // using nearest neighbors, -1

    let arr_gc = gc_cache(&seq1);
    let mut arr_dh = Vec::new();
    let mut arr_ds = Vec::new();
    let mut arr_tm = Vec::new();

    for _ in 0..n {
        arr_dh.push(vec![0.0; n]);
        arr_ds.push(vec![0.0; n]);
        arr_tm.push(vec![f64::INFINITY; n]);
    }

    // fill in the diagonal
    for i in 0..n {
        if i == n - 1 {
            // hackish
            arr_dh[i][i] = arr_dh[i - 1][i - 1];
            arr_ds[i][i] = arr_ds[i - 1][i - 1];
            continue;
        }

        let pair = [seq1[i], seq1[i + 1], b'/', seq2[i], seq2[i + 1]];
        let pair = pair.as_slice();
        let &(dh, ds) = dna()
            .nn
            .get(pair)
            .unwrap_or_else(|| &dna().internal_mm[pair]);

        arr_dh[i][i] = dh;
        arr_ds[i][i] = ds;

        // fill in the tm array
        for i in 0..n {
            for j in i + 1..n {
                arr_dh[i][j] = arr_dh[i][j - 1] + arr_dh[j][j];
                arr_ds[i][j] = arr_ds[i][j - 1] + arr_ds[j][j];
                arr_tm[i][j] = calc_tm(arr_dh[i][j], arr_ds[i][j], pcr, arr_gc[i][j], j - i + 1);
            }
        }
    }

    arr_tm
}

/// Return the GC ratio of each range, between i and j, in the sequence
///
/// # Args
///
/// - seq: The sequence whose tm we're querying
///
/// # Returns
///
/// - [`Cache`]: A cache for GC ratio lookup
pub fn gc_cache(seq: &[u8]) -> Cache {
    let n = seq.len();
    let mut arr_gc = vec![vec![f64::INFINITY; n]; n];

    // fill in the diagonal
    for i in 0..n {
        if i == n - 1 {
            // hackish
            arr_gc[i][i] = arr_gc[i - 1][i - 1];
            continue;
        }

        arr_gc[i][i] = b"GC".contains(&seq[i]) as u64 as f64;

        if i == n - 2 && arr_gc[i][i] == 0.0 {
            // don't ignore last pair
            arr_gc[i][i] = b"GC".contains(&seq[i + 1]) as u64 as f64;
        }
    }

    // fill in the upper right of the array
    for i in 0..n {
        for j in i + 1..n {
            arr_gc[i][j] = arr_gc[i][j - 1] + arr_gc[j][j];
        }
    }

    // convert to ratios
    for i in 0..n {
        for j in i..n {
            arr_gc[i][j] = round1(arr_gc[i][j] / (j - i + 1) as f64);
        }
    }

    arr_gc
}

// TODO: turn this into a result

/// Parse and prepare the input sequences. Throw if there's an issue.
///
/// # Args
///
/// - seq1: The main sequence whose tm is being calculated
/// - seq2: The second sequence that seq2 is annealing to
///
/// # Returns
///
/// - `(Vec<u8>, Vec<u8>)`: The sequences to use for tm calculation
pub fn parse_input(seq1: &[u8], seq2: Option<&[u8]>) -> (Vec<u8>, Vec<u8>) {
    let seq1 = seq1.to_ascii_uppercase();
    let seq2 = seq2.map_or_else(
        || seq1.iter().map(|b| dna().complement[b]).collect(),
        Vec::from,
    );

    assert_eq!(
        seq1.len(),
        seq2.len(),
        "Length mismatch between seq1 and seq2"
    );
    assert!(
        seq1.len() >= 2,
        "Sequence, {}bp, is too short for tm calculation",
        seq1.len()
    );

    return (seq1, seq2);
}

/// Apply the correction formula to estimate Tm
///
/// # Args
///
/// - dh: Accumulated entropy
/// - ds: Accumulated enthalpy
/// - pcr: Whether this is for PCR or not
/// - gc: The GC% of the sequence
/// - seq_len: The length of the sequence
///
/// # Returns
///
/// - [`f64`]: The estimated tm
#[expect(non_snake_case)]
pub fn calc_tm(dh: f64, ds: f64, pcr: bool, gc: f64, seq_len: usize) -> f64 {
    // adjust salt based on mode
    let (seq1_conc, seq2_conc, Na, K, Tris, Mg, dNTPs) = if pcr {
        // see Thermo for tris
        // see NEB for Mg & dNTPS
        (250.0, 0.0, 0.0, 50.0, 2.0, 1.5, 0.2)
    } else {
        (25.0, 25.0, 50.0, 0.0, 0.0, 0.0, 0.0)
    };

    // salt correction for deltaS
    // copied-pasted from Bio.SeqUtils' use of a decision tree by:
    // Owczarzy et al. (2008), Biochemi[u8]y 4 7: 5336-5353
    let Mon = Na + K + Tris / 2.0; // monovalent ions;
    let mut mg = Mg * 1e-3; // Lowercase ions (mg, mon, dntps) are molar;
    let mon = Mon * 1e-3;

    // coefficients to a multi-variate from the paper
    let (mut a, b, c, mut d, e, f, mut g) = (3.92, -0.911, 6.26, 1.42, -48.2, 52.5, 8.31);

    if dNTPs > 0.0 {
        let dntps = dNTPs * 1e-3;
        let ka = 3e4; // Dissociation constant for Mg:dNTP
        // Free Mg2+ calculation:
        mg = (-(ka * dntps - ka * mg + 1.0)
            + ((ka * dntps - ka * mg + 1.0f64).powi(2) + 4.0 * ka * mg).sqrt())
            / (2.0 * ka);
    }
    if Mon > 0.0 {
        let R = mg.sqrt() / mon;
        if R < 0.22 {
            return (4.29 * gc / 100.0 - 3.95) * 1e-5 * mon.ln() + 9.4e-6 * mon.ln().powi(2);
        } else if R < 6.0 {
            a = 3.92 * (0.843 - 0.352 * mon.sqrt() * mon.ln());
            d = 1.42 * (1.279 - 4.03e-3 * mon.ln() - 8.03e-3 * mon.ln().powi(2));
            g = 8.31 * (0.486 - 0.258 * mon.ln() + 5.25e-3 * mon.ln().powi(3));
        }
    }
    let corr = (a
        + b * mg.ln()
        + (gc / 100.0) * (c + d * mg.ln())
        + (1.0 / (2.0 * (seq_len as f64 - 1.0))) * (e + f * mg.ln() + g * mg.ln().powi(2)))
        * 1e-5;

    // tm with concentration consideration
    let k = (seq1_conc - (seq2_conc / 2.0f64)) * 1e-9;
    let R = 1.9872;
    let est = (dh * 1000.0) / (ds + R * k.ln()) - 273.15;

    // add in salt correction
    let est = 1.0 / (1.0 / (est + 273.15) + corr) - 273.1;

    round1(est)
}

/// Return the GC ratio of a sequence.
pub fn gc(seq: &[u8]) -> f64 {
    let [g, c] = counts(seq, [b'G', b'C']);
    ((g + c) as f64) / (seq.len() as f64)
}

fn counts<const N: usize>(seq: &[u8], symbols: [u8; N]) -> [usize; N] {
    let mut indices = [0; 256];
    let mut counts = [0; N];

    for (i, s) in symbols.into_iter().enumerate() {
        indices[s as usize] = i + 1;
    }

    for &b in seq {
        if indices[b as usize] != 0 {
            counts[indices[b as usize] - 1] += 1;
        }
    }

    counts
}

#[cfg(test)]
mod test {
    //! Test Tm calculation

    use approx::assert_relative_eq;

    #[test]
    fn counts() {
        let counts = super::counts(b"AcBdBBrBeAlBaThakaZAAzaaBBAA", [b'a', b'A', b'B']);
        assert_eq!(counts, [5, 6, 7]);
    }

    /// Oligo tm calculation.
    #[test]
    fn test_calc_tm() {
        // values are from Table 1 of IDT's:
        // Owczarzy et al. (2008), Biochemistry 4 7: 5336-5353
        // with a 1.5mM Mg concentration which looks typical according to NEB
        let experimental_tms = [
            ("GGGACCGCCT", 51.9),
            ("CCATTGCTACC", 42.7),
            ("GCAGTGGATGTGAGA", 55.1),
            ("CTGGTCTGGATCTGAGAACTTCAGG", 67.7),
            ("CTTAAGATATGAGAACTTCAACTAATGTGT", 59.7),
            ("AGTCTGGTCTGGATCTGAGAACTTCAGGCT", 71.6),
        ];

        for (seq, actual) in experimental_tms {
            let calc = super::tm(seq.as_bytes(), None, None);
            assert_relative_eq!(calc, actual, epsilon = 7.0);
        }
    }

    /// Create a cache for tms over ranges in the sequence.
    #[test]
    fn test_tm_cache() {
        let seq = "AGTCTGGTCTGGATCTGAGAACTTCAGGCT";
        let n = seq.len();

        let cache = super::tm_cache(seq.as_bytes(), None, None);

        assert_relative_eq!(cache[0][n - 1], 71.6, epsilon = 3.0);
        assert!(cache[3][n - 3] < 71.6);
        assert_eq!(f64::INFINITY, cache[5][5]);
        assert_eq!(f64::INFINITY, cache[5][1]);
    }

    /// Create a cache of GC ratios from i to j.
    #[test]
    fn test_gc_cache() {
        let seq = "GGATTACCCAGATAGATAGAT";
        let ranges = [(0, seq.len() - 1), (5, 9), (3, 15)];

        let cache = super::gc_cache(seq.as_bytes());

        for (s, e) in ranges {
            let est = cache[s][e];
            let ss = &seq[s..e + 1];
            let [g, c] = super::counts(ss.as_bytes(), [b'G', b'C']);
            let gc_count = g + c;
            let gc_actual = (gc_count as f64) / ((e - s + 1) as f64);

            assert_relative_eq!(gc_actual, est, epsilon = 0.02);
        }
    }
}
