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
