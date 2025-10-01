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
    pub bulge_loops: LoopEnergy,
    pub complement: Comp,
    pub de: BpEnergy,
    pub hairpin_loops: LoopEnergy,
    pub multibranch: MultiBranch,
    pub internal_loops: LoopEnergy,
    pub internal_mm: BpEnergy,
    pub nn: BpEnergy,
    pub terminal_mm: BpEnergy,
    pub tri_tetra_loops: Option<BpEnergy>,
}
