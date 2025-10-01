//! Test DNA/RNA folding.

use approx::assert_relative_eq;

use super::{Value, Values};
use crate::{dna, fold::Desc, rna};

/// Fold function.
#[test]
#[should_panic]
fn test_fold_p1() {
    // should throw if a nonsense sequence is provided
    // with self.assertRaises(RuntimeError):
    super::dg(b"EASFEASFAST", Some(37.0));
}

#[test]
#[should_panic]
fn test_fold_p2() {
    // U and T, mix of RNA and DNA
    // with self.assertRaises(RuntimeError):
    super::dg(b"ATGCATGACGATUU", Some(37.0));
}

#[test]
fn test_fold_s1() {
    // not throw
    super::dg(b"ATGGATTTAGATAGAT", None);
}

/// Gather a cache of the folded structure.
#[test]
fn test_fold_cache() {
    let seq = b"ATGGATTTAGATAGAT";
    let cache = super::dg_cache(seq, None);
    let seq_dg = super::dg(seq.as_slice(), None);

    assert_relative_eq!(seq_dg, cache[0][seq.len() - 1], epsilon = 1.0);
}

/// DNA folding to find min energy secondary structure.
#[test]
fn test_fold_dna() {
    //'s estimates for free energy estimates of DNA oligos
    let unafold_dgs = [
        ("GGGAGGTCGTTACATCTGGGTAACACCGGTACTGATCCGGTGACCTCCC", -10.94), // branched structure
        (
            "GGGAGGTCGCTCCAGCTGGGAGGAGCGTTGGGGGTATATACCCCCAACACCGGTACTGATCCGGTGACCTCCC",
            -23.4,
        ), // branched structure
        ("CGCAGGGAUACCCGCG", -3.8),
        ("TAGCTCAGCTGGGAGAGCGCCTGCTTTGCACGCAGGAGGT", -6.85),
        (
            "GGGGGCATAGCTCAGCTGGGAGAGCGCCTGCTTTGCACGCAGGAGGTCTGCGGTTCGATCCCGCGCGCTCCCACCA",
            -15.50,
        ),
        ("TGAGACGGAAGGGGATGATTGTCCCCTTCCGTCTCA", -18.10),
        ("ACCCCCTCCTTCCTTGGATCAAGGGGCTCAA", -3.65),
    ];

    for (seq, ufold) in unafold_dgs {
        let d = super::dg(seq.as_bytes(), Some(37.0));

        // a 60% difference
        let delta = (0.6 * d.min(ufold)).abs();
        assert_relative_eq!(d, ufold, epsilon = delta);
    }
}

/// RNA folding to find min energy secondary structure.
#[test]
fn test_fold_rna() {
    //'s estimates for free energy estimates of RNA oligos
    // tests available at https://github.com/jaswindersingh2/SPOT-RNA/blob/master/sample_inputs/batch_seq.fasta
    let unafold_dgs = [
        ("ACCCCCUCCUUCCUUGGAUCAAGGGGCUCAA", -9.5),
        ("AAGGGGUUGGUCGCCUCGACUAAGCGGCUUGGAAUUCC", -10.1),
        ("UUGGAGUACACAACCUGUACACUCUUUC", -4.3),
        ("AGGGAAAAUCCC", -3.3),
        ("GCUUACGAGCAAGUUAAGCAAC", -4.6),
        (
            "UGGGAGGUCGUCUAACGGUAGGACGGCGGACUCUGGAUCCGCUGGUGGAGGUUCGAGUCCUCCCCUCCCAGCCA",
            -32.8,
        ),
        (
            "GGGCGAUGAGGCCCGCCCAAACUGCCCUGAAAAGGGCUGAUGGCCUCUACUG",
            -20.7,
        ),
        (
            "GGGGGCAUAGCUCAGCUGGGAGAGCGCCUGCUUUGCACGCAGGAGGUCUGCGGUUCGAUCCCGCGCGCUCCCACCA",
            -31.4,
        ),
        (
            "CAGCGCGGCGGGCGGGAGUCCGGCGCGCCCUCCAUCCCCGGCGGCGUCGGCAAGGAGUAG",
            -18.26,
        ),
    ];

    for (seq, ufold) in unafold_dgs {
        let d = super::dg(seq.as_bytes(), Some(37.0));

        // a 30% difference
        let delta = (0.3 * d.min(ufold)).abs();
        assert_relative_eq!(d, ufold, epsilon = delta);
    }
}

/// Get the dot bracket notation for a folded structure.
#[test]
fn test_dot_bracket() {
    let seq = b"GGGAGGTCGTTACATCTGGGTAACACCGGTACTGATCCGGTGACCTCCC";
    let structs = super::fold(seq, None);

    assert_eq!(
        b"((((((((.((((......))))..((((.......)))).))))))))".as_slice(),
        super::dot_bracket(seq, &structs),
    );

    let seq = b"ACGCTCACCGTGCCCAGTGAGCGA";
    let structs = super::fold(seq, None);
    assert_eq!(seq.len(), super::dot_bracket(seq, &structs).len());
}

/// Fold a multibranch structure.
#[test]
fn test_multibranch() {
    let seq = b"GGGAGGTCGTTACATCTGGGTAACACCGGTACTGATCCGGTGACCTCCC"; // branch

    let structs = super::fold(seq, None);
    assert!(
        structs
            .iter()
            .any(|v| matches!(v.desc, Desc::Bifurcation(_, _)) && v.ij.contains(&(7, 41)))
    );

    let seq = b"CAGCGCGGCGGGCGGGAGUCCGGCGCGCCCUCCAUCCCCGGCGGCGUCGGCAAGGAGUAG";

    let structs = super::fold(seq, None);
    assert!(
        structs
            .iter()
            .any(|v| matches!(v.desc, Desc::Bifurcation(_, _)) && v.ij.contains(&(2, 56)))
    );
}

/// Create a pair for stack checking.
#[test]
fn test_pair() {
    let seq = b"ATGGAATAGTG";
    assert_eq!(super::calc_pair(seq, 0, 1, 9, 10).as_slice(), b"AT/TG");
}

/// Calc delta G of a stack.
#[test]
fn test_stack() {
    let seq = b"GCUCAGCUGGGAGAGC";
    let temp = 310.15;

    assert_relative_eq!(
        super::calc_stack(seq, 1, 2, 14, 13, temp, rna()),
        -2.1,
        epsilon = 0.1
    );
}

/// Calc delta G calc of a bulge.
#[test]
fn test_bulge() {
    // bulge of CAT on one side and AG on other
    // pg 429 of SantaLucia, 2004
    let seq = b"ACCCCCATCCTTCCTTGAGTCAAGGGGCTCAA";

    let pair_dg = super::bulge(seq, 5, 7, 18, 17, 310.15, dna());
    assert_relative_eq!(3.22, pair_dg, epsilon = 0.4);
}

/// Calc delta G of a hairpin structure.
#[test]
fn test_hairpin() {
    // = b"CCTTGG"
    let seq = b"ACCCCCTCCTTCCTTGGATCAAGGGGCTCAA";
    let i = 11;
    let j = 16;
    let temp = 310.15;
    let hairpin_dg = super::hairpin(seq, i, j, temp, dna());
    // differs from Unafold
    assert_relative_eq!(hairpin_dg, 4.3, epsilon = 1.0);

    // page 428 of SantaLucia, 2004
    // = b"CGCAAG"
    let seq = b"ACCCGCAAGCCCTCCTTCCTTGGATCAAGGGGCTCAA";
    let i = 3;
    let j = 8;
    let hairpin_dg = super::hairpin(seq, i, j, temp, dna());
    assert_relative_eq!(0.67, hairpin_dg, epsilon = 0.1);

    let seq = b"CUUUGCACG";
    let i = 0;
    let j = 8;
    let hairpin_dg = super::hairpin(seq, i, j, temp, rna());
    assert_relative_eq!(4.5, hairpin_dg, epsilon = 0.2);
}

/// Calc dg of an internal loop.
#[test]
fn test_internal_loop() {
    let seq = b"ACCCCCTCCTTCCTTGGATCAAGGGGCTCAA";
    let i = 6;
    let j = 21;
    let temp = 310.15;
    let dg = super::internal_loop(seq, i, i + 4, j, j - 4, temp, dna());
    assert_relative_eq!(dg, 3.5, epsilon = 0.1);
}

/// Calculate _w over some range.
#[test]
fn test_w() {
    let seq = b"GCUCAGCUGGGAGAGC";
    let i = 0;
    let j = 15;
    let temp = 310.15;
    let mut v_cache = Values::new();
    let mut w_cache = Values::new();
    for _ in 0..seq.len() {
        v_cache.push(vec![Value::DEFAULT; seq.len()]);
        w_cache.push(vec![Value::DEFAULT; seq.len()]);
    }
    let (i, j) = super::w(seq, i, j, temp, &mut v_cache, &mut w_cache, rna());
    let value = &w_cache[i][j];
    assert_relative_eq!(value.e, -3.8, epsilon = 0.2);

    let seq = b"CCUGCUUUGCACGCAGG";
    let i = 0;
    let j = 16;
    let temp = 310.15;
    let mut v_cache = Values::new();
    let mut w_cache = Values::new();
    for _ in 0..seq.len() {
        v_cache.push(vec![Value::DEFAULT; seq.len()]);
        w_cache.push(vec![Value::DEFAULT; seq.len()]);
    }
    let (i, j) = super::w(seq, i, j, temp, &mut v_cache, &mut w_cache, rna());
    let value = &w_cache[i][j];
    assert_relative_eq!(value.e, -6.4, epsilon = 0.2);

    let seq = b"GCGGUUCGAUCCCGC";
    let i = 0;
    let j = 14;
    let mut v_cache = Values::new();
    let mut w_cache = Values::new();
    for _ in 0..seq.len() {
        v_cache.push(vec![Value::DEFAULT; seq.len()]);
        w_cache.push(vec![Value::DEFAULT; seq.len()]);
    }
    let (i, j) = super::w(seq, i, j, temp, &mut v_cache, &mut w_cache, rna());
    let value = &w_cache[i][j];
    assert_relative_eq!(value.e, -4.2, epsilon = 0.2);
}
