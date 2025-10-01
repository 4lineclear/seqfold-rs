//! DNA enthalpy and entropy change parameters.

use crate::{BpEnergy, Comp, Energies, LoopEnergy, MultiBranch};

use std::sync::LazyLock;

// TODO: consider bpenergy computation in a LazyLock

/// a, b, c, d in a linear multi-branch energy change function.
///
/// Inferred from:
/// Supplemental Material: Annu.Rev.Biophs.Biomol.Struct.33:415-40
/// doi: 10.1146/annurev.biophys.32.110601.141800
/// The Thermodynamics of DNA Structural Motifs
/// SantaLucia and Hicks, 2004
pub const fn multibranch() -> MultiBranch {
    (2.6, 0.2, 0.2, 2.0)
}

pub fn complement() -> Comp {
    pub static RAW_COMPLEMENT: [(u8, u8); 5] = [
        (b'A', b'T'),
        (b'T', b'A'),
        (b'G', b'C'),
        (b'C', b'G'),
        (b'N', b'N'),
    ];

    Comp::from_iter(RAW_COMPLEMENT)
}

/// The Thermodynamics of DNA Structural Motifs
/// SantaLucia and Hicks, 2004
pub fn nn() -> BpEnergy {
    pub static RAW_NN: [(&[u8], (f64, f64)); 14] = [
        (b"init", (0.2, -5.7)),
        (b"init_G/C", (0.0, 0.0)),
        (b"init_A/T", (2.2, 6.9)),
        (b"sym", (0.0, -1.4)),
        (b"AA/TT", (-7.6, -21.3)),
        (b"AT/TA", (-7.2, -20.4)),
        (b"TA/AT", (-7.2, -21.3)),
        (b"CA/GT", (-8.5, -22.7)),
        (b"GT/CA", (-8.4, -22.4)),
        (b"CT/GA", (-7.8, -21.0)),
        (b"GA/CT", (-8.2, -22.2)),
        (b"CG/GC", (-10.6, -27.2)),
        (b"GC/CG", (-9.8, -24.4)),
        (b"GG/CC", (-8.0, -19.9)),
    ];

    let mut nn = BpEnergy::default();
    for (s, v) in RAW_NN {
        nn.insert(s.to_owned(), v);
        nn.insert(s.iter().copied().rev().collect(), v);
    }
    nn
}

/// Internal mismatch table (DNA)
/// Allawi & SantaLucia (1997), Biochemistry 36: 10581-10594
/// Allawi & SantaLucia (1998), Biochemistry 37: 9435-9444
/// Allawi & SantaLucia (1998), Biochemistry 37: 2170-2179 *
/// Allawi & SantaLucia (1998), Nucl Acids Res 26: 2694-2701 *
/// Peyret et al. (1999), Biochemistry 38: 3468-3477 *
pub fn internal_mm() -> BpEnergy {
    pub static RAW_INTERNAL_MM: [(&[u8], (f64, f64)); 51] = [
        (b"AG/TT", (1.0, 0.9)),
        (b"AT/TG", (-2.5, -8.3)),
        (b"CG/GT", (-4.1, -11.7)),
        (b"CT/GG", (-2.8, -8.0)),
        (b"GG/CT", (3.3, 10.4)),
        (b"GG/TT", (5.8, 16.3)),
        (b"GT/CG", (-4.4, -12.3)),
        (b"GT/TG", (4.1, 9.5)),
        (b"TG/AT", (-0.1, -1.7)),
        (b"TG/GT", (-1.4, -6.2)),
        (b"TT/AG", (-1.3, -5.3)),
        (b"AA/TG", (-0.6, -2.3)),
        (b"AG/TA", (-0.7, -2.3)),
        (b"CA/GG", (-0.7, -2.3)),
        (b"CG/GA", (-4.0, -13.2)),
        (b"GA/CG", (-0.6, -1.0)),
        (b"GG/CA", (0.5, 3.2)),
        (b"TA/AG", (0.7, 0.7)),
        (b"TG/AA", (3.0, 7.4)),
        (b"AC/TT", (0.7, 0.2)),
        (b"AT/TC", (-1.2, -6.2)),
        (b"CC/GT", (-0.8, -4.5)),
        (b"CT/GC", (-1.5, -6.1)),
        (b"GC/CT", (2.3, 5.4)),
        (b"GT/CC", (5.2, 13.5)),
        (b"TC/AT", (1.2, 0.7)),
        (b"TT/AC", (1.0, 0.7)),
        (b"AA/TC", (2.3, 4.6)),
        (b"AC/TA", (5.3, 14.6)),
        (b"CA/GC", (1.9, 3.7)),
        (b"CC/GA", (0.6, -0.6)),
        (b"GA/CC", (5.2, 14.2)),
        (b"GC/CA", (-0.7, -3.8)),
        (b"TA/AC", (3.4, 8.0)),
        (b"TC/AA", (7.6, 20.2)),
        (b"AA/TA", (1.2, 1.7)),
        (b"CA/GA", (-0.9, -4.2)),
        (b"GA/CA", (-2.9, -9.8)),
        (b"TA/AA", (4.7, 12.9)),
        (b"AC/TC", (0.0, -4.4)),
        (b"CC/GC", (-1.5, -7.2)),
        (b"GC/CC", (3.6, 8.9)),
        (b"TC/AC", (6.1, 16.4)),
        (b"AG/TG", (-3.1, -9.5)),
        (b"CG/GG", (-4.9, -15.3)),
        (b"GG/CG", (-6.0, -15.8)),
        (b"TG/AG", (1.6, 3.6)),
        (b"AT/TT", (-2.7, -10.8)),
        (b"CT/GT", (-5.0, -15.8)),
        (b"GT/CT", (-2.2, -8.4)),
        (b"TT/AT", (0.2, -1.5)),
    ];

    let mut internal_mm = BpEnergy::default();
    for (s, v) in RAW_INTERNAL_MM {
        internal_mm.insert(s.to_owned(), v);
        let srev = s.iter().copied().rev().collect();
        if !internal_mm.contains_key(&srev) {
            internal_mm.insert(srev, v);
        }
    }
    internal_mm
}

/// Terminal mismatch table (DNA)
/// SantaLucia & Peyret (2001) Patent Application WO 01/94611
pub fn terminal_mm() -> BpEnergy {
    pub static RAW_TERMINAL_MM: [(&[u8], (f64, f64)); 48] = [
        (b"AA/TA", (-3.1, -7.8)),
        (b"TA/AA", (-2.5, -6.3)),
        (b"CA/GA", (-4.3, -10.7)),
        (b"GA/CA", (-8.0, -22.5)),
        (b"AC/TC", (-0.1, 0.5)),
        (b"TC/AC", (-0.7, -1.3)),
        (b"CC/GC", (-2.1, -5.1)),
        (b"GC/CC", (-3.9, -10.6)),
        (b"AG/TG", (-1.1, -2.1)),
        (b"TG/AG", (-1.1, -2.7)),
        (b"CG/GG", (-3.8, -9.5)),
        (b"GG/CG", (-0.7, -19.2)),
        (b"AT/TT", (-2.4, -6.5)),
        (b"TT/AT", (-3.2, -8.9)),
        (b"CT/GT", (-6.1, -16.9)),
        (b"GT/CT", (-7.4, -21.2)),
        (b"AA/TC", (-1.6, -4.0)),
        (b"AC/TA", (-1.8, -3.8)),
        (b"CA/GC", (-2.6, -5.9)),
        (b"CC/GA", (-2.7, -6.0)),
        (b"GA/CC", (-5.0, -13.8)),
        (b"GC/CA", (-3.2, -7.1)),
        (b"TA/AC", (-2.3, -5.9)),
        (b"TC/AA", (-2.7, -7.0)),
        (b"AC/TT", (-0.9, -1.7)),
        (b"AT/TC", (-2.3, -6.3)),
        (b"CC/GT", (-3.2, -8.0)),
        (b"CT/GC", (-3.9, -10.6)),
        (b"GC/CT", (-4.9, -13.5)),
        (b"GT/CC", (-3.0, -7.8)),
        (b"TC/AT", (-2.5, -6.3)),
        (b"TT/AC", (-0.7, -1.2)),
        (b"AA/TG", (-1.9, -4.4)),
        (b"AG/TA", (-2.5, -5.9)),
        (b"CA/GG", (-3.9, -9.6)),
        (b"CG/GA", (-6.0, -15.5)),
        (b"GA/CG", (-4.3, -11.1)),
        (b"GG/CA", (-4.6, -11.4)),
        (b"TA/AG", (-2.0, -4.7)),
        (b"TG/AA", (-2.4, -5.8)),
        (b"AG/TT", (-3.2, -8.7)),
        (b"AT/TG", (-3.5, -9.4)),
        (b"CG/GT", (-3.8, -9.0)),
        (b"CT/GG", (-6.6, -18.7)),
        (b"GG/CT", (-5.7, -15.9)),
        (b"GT/CG", (-5.9, -16.1)),
        (b"TG/AT", (-3.9, -10.5)),
        (b"TT/AG", (-3.6, -9.8)),
    ];

    let mut terminal_mm = BpEnergy::default();
    for (s, v) in RAW_TERMINAL_MM {
        terminal_mm.insert(s.to_owned(), v);
        let srev = s.iter().copied().rev().collect();
        if !terminal_mm.contains_key(&srev) {
            terminal_mm.insert(srev, v);
        }
    }
    terminal_mm
}

/// DNA dangling ends
///
/// Bommarito et al. (2000), Nucl Acids Res 28: 1929-1934
pub fn de() -> BpEnergy {
    pub static RAW_DE: [(&[u8], (f64, f64)); 32] = [
        (b"AA/.T", (0.2, 2.3)),
        (b"AC/.G", (-6.3, -17.1)),
        (b"AG/.C", (-3.7, -10.0)),
        (b"AT/.A", (-2.9, -7.6)),
        (b"CA/.T", (0.6, 3.3)),
        (b"CC/.G", (-4.4, -12.6)),
        (b"CG/.C", (-4.0, -11.9)),
        (b"CT/.A", (-4.1, -13.0)),
        (b"GA/.T", (-1.1, -1.6)),
        (b"GC/.G", (-5.1, -14.0)),
        (b"GG/.C", (-3.9, -10.9)),
        (b"GT/.A", (-4.2, -15.0)),
        (b"TA/.T", (-6.9, -20.0)),
        (b"TC/.G", (-4.0, -10.9)),
        (b"TG/.C", (-4.9, -13.8)),
        (b"TT/.A", (-0.2, -0.5)),
        (b".A/AT", (-0.7, -0.8)),
        (b".C/AG", (-2.1, -3.9)),
        (b".G/AC", (-5.9, -16.5)),
        (b".T/AA", (-0.5, -1.1)),
        (b".A/CT", (4.4, 14.9)),
        (b".C/CG", (-0.2, -0.1)),
        (b".G/CC", (-2.6, -7.4)),
        (b".T/CA", (4.7, 14.2)),
        (b".A/GT", (-1.6, -3.6)),
        (b".C/GG", (-3.9, -11.2)),
        (b".G/GC", (-3.2, -10.4)),
        (b".T/GA", (-4.1, -13.1)),
        (b".A/TT", (2.9, 10.4)),
        (b".C/TG", (-4.4, -13.1)),
        (b".G/TC", (-5.2, -15.0)),
        (b".T/TA", (-3.8, -12.6)),
    ];

    let mut de = BpEnergy::default();
    for (s, v) in RAW_DE {
        de.insert(s.to_owned(), v);
        let srev = s.iter().copied().rev().collect();
        if !de.contains_key(&srev) {
            de.insert(srev, v);
        }
    }
    de
}

/// Experimental delta H and delta S for tri/tetra loops
///
/// Supplemental Material: Annu.Rev.Biophs.Biomol.Struct.33:415-40
/// doi: 10.1146/annurev.biophys.32.110601.141800
/// The Thermodynamics of DNA Structural Motifs
/// SantaLucia and Hicks, 2004
///
/// delta S was computed using delta G and delta H and is in cal / (K x mol)
/// (versus delta H in kcal / mol)
pub fn tri_tetra_loops() -> BpEnergy {
    pub static RAW_TRI_TETRA_LOOPS: [(&[u8], (f64, f64)); 131] = [
        (b"AGAAT", (-1.5, 0.0)),
        (b"AGCAT", (-1.5, 0.0)),
        (b"AGGAT", (-1.5, 0.0)),
        (b"AGTAT", (-1.5, 0.0)),
        (b"CGAAG", (-2.0, 0.0)),
        (b"CGCAG", (-2.0, 0.0)),
        (b"CGGAG", (-2.0, 0.0)),
        (b"CGTAG", (-2.0, 0.0)),
        (b"GGAAC", (-2.0, 0.0)),
        (b"GGCAC", (-2.0, 0.0)),
        (b"GGGAC", (-2.0, 0.0)),
        (b"GGTAC", (-2.0, 0.0)),
        (b"TGAAA", (-1.5, 0.0)),
        (b"TGCAA", (-1.5, 0.0)),
        (b"TGGAA", (-1.5, 0.0)),
        (b"TGTAA", (-1.5, 0.0)),
        (b"AAAAAT", (0.5, 0.6)),
        (b"AAAACT", (0.7, -1.6)),
        (b"AAACAT", (1.0, -1.6)),
        (b"ACTTGT", (0.0, -4.2)),
        (b"AGAAAT", (-1.1, -1.6)),
        (b"AGAGAT", (-1.1, -1.6)),
        (b"AGATAT", (-1.5, -1.6)),
        (b"AGCAAT", (-1.6, -1.6)),
        (b"AGCGAT", (-1.1, -1.6)),
        (b"AGCTTT", (0.2, -1.6)),
        (b"AGGAAT", (-1.1, -1.6)),
        (b"AGGGAT", (-1.1, -1.6)),
        (b"AGGGGT", (0.5, -0.6)),
        (b"AGTAAT", (-1.6, -1.6)),
        (b"AGTGAT", (-1.1, -1.6)),
        (b"AGTTCT", (0.8, -1.6)),
        (b"ATTCGT", (-0.2, -1.6)),
        (b"ATTTGT", (0.0, -1.6)),
        (b"ATTTTT", (-0.5, -1.6)),
        (b"CAAAAG", (0.5, 1.3)),
        (b"CAAACG", (0.7, 0.0)),
        (b"CAACAG", (1.0, 0.0)),
        (b"CAACCG", (0.0, 0.0)),
        (b"CCTTGG", (0.0, -2.6)),
        (b"CGAAAG", (-1.1, 0.0)),
        (b"CGAGAG", (-1.1, 0.0)),
        (b"CGATAG", (-1.5, 0.0)),
        (b"CGCAAG", (-1.6, 0.0)),
        (b"CGCGAG", (-1.1, 0.0)),
        (b"CGCTTG", (0.2, 0.0)),
        (b"CGGAAG", (-1.1, 0.0)),
        (b"CGGGAG", (-1.0, 0.0)),
        (b"CGGGGG", (0.5, 1.0)),
        (b"CGTAAG", (-1.6, 0.0)),
        (b"CGTGAG", (-1.1, 0.0)),
        (b"CGTTCG", (0.8, 0.0)),
        (b"CTTCGG", (-0.2, 0.0)),
        (b"CTTTGG", (0.0, 0.0)),
        (b"CTTTTG", (-0.5, 0.0)),
        (b"GAAAAC", (0.5, 3.2)),
        (b"GAAACC", (0.7, 0.0)),
        (b"GAACAC", (1.0, 0.0)),
        (b"GCTTGC", (0.0, -2.6)),
        (b"GGAAAC", (-1.1, 0.0)),
        (b"GGAGAC", (-1.1, 0.0)),
        (b"GGATAC", (-1.6, 0.0)),
        (b"GGCAAC", (-1.6, 0.0)),
        (b"GGCGAC", (-1.1, 0.0)),
        (b"GGCTTC", (0.2, 0.0)),
        (b"GGGAAC", (-1.1, 0.0)),
        (b"GGGGAC", (-1.1, 0.0)),
        (b"GGGGGC", (0.5, 1.0)),
        (b"GGTAAC", (-1.6, 0.0)),
        (b"GGTGAC", (-1.1, 0.0)),
        (b"GGTTCC", (0.8, 0.0)),
        (b"GTTCGC", (-0.2, 0.0)),
        (b"GTTTGC", (0.0, 0.0)),
        (b"GTTTTC", (-0.5, 0.0)),
        (b"GAAAAT", (0.5, 3.2)),
        (b"GAAACT", (1.0, 0.0)),
        (b"GAACAT", (1.0, 0.0)),
        (b"GCTTGT", (0.0, -1.6)),
        (b"GGAAAT", (-1.1, 0.0)),
        (b"GGAGAT", (-1.1, 0.0)),
        (b"GGATAT", (-1.6, 0.0)),
        (b"GGCAAT", (-1.6, 0.0)),
        (b"GGCGAT", (-1.1, 0.0)),
        (b"GGCTTT", (-0.1, 0.0)),
        (b"GGGAAT", (-1.1, 0.0)),
        (b"GGGGAT", (-1.1, 0.0)),
        (b"GGGGGT", (0.5, 1.0)),
        (b"GGTAAT", (-1.6, 0.0)),
        (b"GGTGAT", (-1.1, 0.0)),
        (b"GTATAT", (-0.5, 0.0)),
        (b"GTTCGT", (-0.4, 0.0)),
        (b"GTTTGT", (-0.4, 0.0)),
        (b"GTTTTT", (-0.5, 0.0)),
        (b"TAAAAA", (0.5, -0.3)),
        (b"TAAACA", (0.7, -1.6)),
        (b"TAACAA", (1.0, -1.6)),
        (b"TCTTGA", (0.0, -4.2)),
        (b"TGAAAA", (-1.1, -1.6)),
        (b"TGAGAA", (-1.1, -1.6)),
        (b"TGATAA", (-1.6, -1.6)),
        (b"TGCAAA", (-1.6, -1.6)),
        (b"TGCGAA", (-1.1, -1.6)),
        (b"TGCTTA", (0.2, -1.6)),
        (b"TGGAAA", (-1.1, -1.6)),
        (b"TGGGAA", (-1.1, -1.6)),
        (b"TGGGGA", (0.5, -0.6)),
        (b"TGTAAA", (-1.6, -1.6)),
        (b"TGTGAA", (-1.1, -1.6)),
        (b"TGTTCA", (0.8, -1.6)),
        (b"TTTCGA", (-0.2, -1.6)),
        (b"TTTTGA", (0.0, -1.6)),
        (b"TTTTTA", (-0.5, -1.6)),
        (b"TAAAAG", (0.5, 1.6)),
        (b"TAAACG", (1.0, -1.6)),
        (b"TAACAG", (1.0, -1.6)),
        (b"TCTTGG", (0.0, -3.2)),
        (b"TGAAAG", (-1.0, -1.6)),
        (b"TGAGAG", (-1.0, -1.6)),
        (b"TGATAG", (-1.5, -1.6)),
        (b"TGCAAG", (-1.5, -1.6)),
        (b"TGCGAG", (-1.0, -1.6)),
        (b"TGCTTG", (-0.1, -1.6)),
        (b"TGGAAG", (-1.0, -1.6)),
        (b"TGGGAG", (-1.0, -1.6)),
        (b"TGGGGG", (0.5, -0.6)),
        (b"TGTAAG", (-1.5, -1.6)),
        (b"TGTGAG", (-1.0, -1.6)),
        (b"TTTCGG", (-0.4, -1.6)),
        (b"TTTTAG", (-1.0, -1.6)),
        (b"TTTTGG", (-0.4, -1.6)),
        (b"TTTTTG", (-0.5, -1.6)),
    ];

    BpEnergy::from_iter(
        RAW_TRI_TETRA_LOOPS
            .into_iter()
            .map(|(s, v)| (s.to_owned(), v)),
    )
}

/// Enthalpy and entropy increments for length dependence of internal loops
///
/// Were calculated from delta G Table 4 of SantaLucia, 2004:
///
/// Annu.Rev.Biophs.Biomol.Struct.33:415-40
/// doi: 10.1146/annurev.biophys.32.110601.141800
/// The Thermodynamics of DNA Structural Motifs
/// SantaLucia and Hicks, 2004
///
/// Additional loop sizes are accounted for with the Jacobson-Stockmayer
/// entry extrapolation formula in paper:
/// delta G (loop-n) = delta G (loop-x) + 2.44 x R x 310.15 x ln(n / x)
///
/// Additional correction is applied for asymmetric loops in paper:
/// delta G (asymmetry) = |length A - length B| x 0.3 (kcal / mol)
/// where A and B are lengths of both sides of loop
pub fn internal_loops() -> LoopEnergy {
    pub static RAW_INTERNAL_LOOPS: [(usize, (f64, f64)); 30] = [
        (1, (0.0, 0.0)),
        (2, (0.0, 0.0)),
        (3, (0.0, -10.3)),
        (4, (0.0, -11.6)),
        (5, (0.0, -12.9)),
        (6, (0.0, -14.2)),
        (7, (0.0, -14.8)),
        (8, (0.0, -15.5)),
        (9, (0.0, -15.8)),
        (10, (0.0, -15.8)),
        (11, (0.0, -16.1)),
        (12, (0.0, -16.8)),
        (13, (0.0, -16.4)),
        (14, (0.0, -17.4)),
        (15, (0.0, -17.7)),
        (16, (0.0, -18.1)),
        (17, (0.0, -18.4)),
        (18, (0.0, -18.7)),
        (19, (0.0, -18.7)),
        (20, (0.0, -19.0)),
        (21, (0.0, -19.0)),
        (22, (0.0, -19.3)),
        (23, (0.0, -19.7)),
        (24, (0.0, -20.0)),
        (25, (0.0, -20.3)),
        (26, (0.0, -20.3)),
        (27, (0.0, -20.6)),
        (28, (0.0, -21.0)),
        (29, (0.0, -21.0)),
        (30, (0.0, -21.3)),
    ];

    LoopEnergy::from_iter(RAW_INTERNAL_LOOPS)
}

/// Enthalpy and entropy increments for length depedence of bulge loops
///
/// Were calculated from delta G Table 4 of SantaLucia, 2004:
///
/// Annu.Rev.Biophs.Biomol.Struct.33:415-40
/// doi: 10.1146/annurev.biophys.32.110601.141800
/// The Thermodynamics of DNA Structural Motifs
/// SantaLucia and Hicks, 2004
///
/// For bulge loops of size 1, the intervening NN energy is used.
/// Closing AT penalty is applied on both sides
pub fn bulge_loops() -> LoopEnergy {
    pub static RAW_BULGE_LOOPS: [(usize, (f64, f64)); 30] = [
        (1, (0.0, -12.9)),
        (2, (0.0, -9.4)),
        (3, (0.0, -10.0)),
        (4, (0.0, -10.3)),
        (5, (0.0, -10.6)),
        (6, (0.0, -11.3)),
        (7, (0.0, -11.9)),
        (8, (0.0, -12.6)),
        (9, (0.0, -13.2)),
        (10, (0.0, -13.9)),
        (11, (0.0, -14.2)),
        (12, (0.0, -14.5)),
        (13, (0.0, -14.8)),
        (14, (0.0, -15.5)),
        (15, (0.0, -15.8)),
        (16, (0.0, -16.1)),
        (17, (0.0, -16.4)),
        (18, (0.0, -16.8)),
        (19, (0.0, -16.8)),
        (20, (0.0, -17.1)),
        (21, (0.0, -17.4)),
        (22, (0.0, -17.4)),
        (23, (0.0, -17.7)),
        (24, (0.0, -17.7)),
        (25, (0.0, -18.1)),
        (26, (0.0, -18.1)),
        (27, (0.0, -18.4)),
        (28, (0.0, -18.7)),
        (29, (0.0, -18.7)),
        (30, (0.0, -19.0)),
    ];

    LoopEnergy::from_iter(RAW_BULGE_LOOPS)
}

/// Enthalpy and entropy increments for length depedence of hairpin loops
///
/// Were calculated from delta G Table 4 of SantaLucia, 2004:
///
/// Annu.Rev.Biophs.Biomol.Struct.33:415-40
/// doi: 10.1146/annurev.biophys.32.110601.141800
/// The Thermodynamics of DNA Structural Motifs
/// SantaLucia and Hicks, 2004
///
/// For hairpins of length 3 and 4, the entropy values are looked up
/// in the DNA_TRI_TETRA_LOOPS Dict
///
/// From formula 8-9 of the paper:
/// An additional 1.6 delta entropy penalty if the hairpin is closed by AT
fn raw_hairpin_loops() -> LoopEnergy {
    pub static RAW_HAIRPIN_LOOPS: [(usize, (f64, f64)); 30] = [
        (1, (0.0, 0.0)),
        (2, (0.0, 0.0)),
        (3, (0.0, -11.3)),
        (4, (0.0, -11.3)),
        (5, (0.0, -10.6)),
        (6, (0.0, -12.9)),
        (7, (0.0, -13.5)),
        (8, (0.0, -13.9)),
        (9, (0.0, -14.5)),
        (10, (0.0, -14.8)),
        (11, (0.0, -15.5)),
        (12, (0.0, -16.1)),
        (13, (0.0, -16.1)),
        (14, (0.0, -16.4)),
        (15, (0.0, -16.8)),
        (16, (0.0, -17.1)),
        (17, (0.0, -17.4)),
        (18, (0.0, -17.7)),
        (19, (0.0, -18.1)),
        (20, (0.0, -18.4)),
        (21, (0.0, -18.7)),
        (22, (0.0, -18.7)),
        (23, (0.0, -19.0)),
        (24, (0.0, -19.3)),
        (25, (0.0, -19.7)),
        (26, (0.0, -19.7)),
        (27, (0.0, -19.7)),
        (28, (0.0, -20.0)),
        (29, (0.0, -20.0)),
        (30, (0.0, -20.3)),
    ];

    LoopEnergy::from_iter(RAW_HAIRPIN_LOOPS)
}

pub fn calc_dna() -> Energies {
    Energies {
        bulge_loops: bulge_loops(),
        complement: complement(),
        de: de(),
        hairpin_loops: raw_hairpin_loops(),
        multibranch: multibranch(),
        internal_loops: internal_loops(),
        internal_mm: internal_mm(),
        nn: nn(),
        terminal_mm: terminal_mm(),
        tri_tetra_loops: Some(tri_tetra_loops()),
    }
}

pub fn dna() -> &'static Energies {
    pub static DNA: LazyLock<Energies> = LazyLock::new(|| calc_dna());
    &DNA
}
