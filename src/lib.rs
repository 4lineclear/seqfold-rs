#![doc = include_str!("../README.md")]

mod util;

pub mod dna;
pub mod fold;
pub mod rna;
pub mod tm;

pub use dna::dna;
pub use rna::rna;

pub use fold::fold;

use std::collections::HashMap;

// TODO: see if we can move away from statics

pub type Cache = Vec<Vec<f64>>;
pub type MultiBranch = (f64, f64, f64, f64);

pub type Comp = HashMap<u8, u8>;
pub type BpEnergy = HashMap<Vec<u8>, (f64, f64)>;
pub type LoopEnergy = HashMap<usize, (f64, f64)>;

#[derive(Debug)]
pub struct Energies {
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
