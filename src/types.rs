//! Types shared between dna.py and rna.py

use std::collections::HashMap;

// TODO: see if we can move away from statics

pub type Cache = Vec<Vec<f64>>;
pub type MultiBranch = (f64, f64, f64, f64);

pub type Comp = HashMap<u8, u8>;
pub type BpEnergy = HashMap<String, (f64, f64)>;
pub type LoopEnergy = HashMap<u64, (f64, f64)>;

#[derive(Debug)]
pub struct Engergies {
    pub bulge_loops: &'static LoopEnergy,
    pub complement: &'static Comp,
    pub de: &'static BpEnergy,
    pub hairpin_loops: &'static LoopEnergy,
    pub multibranch: MultiBranch,
    pub internal_loops: &'static LoopEnergy,
    pub internal_mm: &'static BpEnergy,
    pub nn: &'static BpEnergy,
    pub terminal_mm: &'static BpEnergy,
    pub tri_tetra_loops: Option<&'static BpEnergy>,
}
