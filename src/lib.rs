#![doc = include_str!("../README.md")]
#![allow(unsafe_code)]
// #![deny(
//     clippy::all,
//     clippy::pedantic,
//     // clippy::cargo,
//     clippy::nursery,
//     // missing_docs,
//     // rustdoc::all,
//     future_incompatible
// )]
// #![allow(clippy::must_use_candidate)]
// #![allow(clippy::missing_panics_doc)]
// #![allow(clippy::cast_precision_loss)]

mod util;

pub mod dna;
pub mod fold;
pub mod rna;
pub mod tm;

use std::ops::Index;

pub use dna::dna;
pub use rna::rna;

pub use fold::fold;

pub type Cache = Vec<Vec<f64>>;
pub type MultiBranch = (f64, f64, f64, f64);

use rustc_hash::FxHashMap as HashMap;

// TODO: consider moving from hashmap to phf

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

#[derive(Debug)]
pub struct Comp([Option<u8>; 256]);

impl FromIterator<(u8, u8)> for Comp {
    fn from_iter<T: IntoIterator<Item = (u8, u8)>>(iter: T) -> Self {
        let mut this = Self([None; 256]);
        for (i, b) in iter {
            this.0[i as usize] = Some(b);
        }
        this
    }
}

impl Index<u8> for Comp {
    type Output = u8;

    fn index(&self, index: u8) -> &Self::Output {
        self.0[index as usize].as_ref().unwrap()
    }
}

#[derive(Debug)]
pub struct LoopEnergy([(f64, f64); 30]);

impl LoopEnergy {
    fn get(&self, index: usize, temp: f64) -> f64 {
        assert!(index != 0);

        if let Some((d_h, d_s)) = self.0.get(index - 1).copied() {
            fold::calc_d_g(d_h, d_s, temp)
        } else {
            // it's too large, extrapolate
            let (d_h, d_s) = self.0[self.0.len() - 1];
            let d_g_inc = fold::calc_d_g(d_h, d_s, temp);
            fold::calc_j_s(index, self.0.len(), d_g_inc, temp)
        }
    }
}

impl From<[(f64, f64); 30]> for LoopEnergy {
    fn from(value: [(f64, f64); 30]) -> Self {
        Self(value)
    }
}

#[derive(Debug, Default)]
pub struct BpEnergy {
    values: HashMap<u64, (f64, f64)>,
}

// TODO: create a byte-slice trait

pub fn interpret_bytes<'a>(mut b: impl Iterator<Item = &'a u8>) -> u64 {
    u64::from_be_bytes(std::array::from_fn(|_| b.next().copied().unwrap_or(0)))
}

impl BpEnergy {
    pub fn new<'a, I>(values: I) -> Self
    where
        I: IntoIterator<Item = (&'a [u8], (f64, f64))>,
    {
        let values = values
            .into_iter()
            .map(|(b, v)| {
                assert!(b.len() <= 8);
                (interpret_bytes(b.iter()), v)
            })
            .collect();
        Self { values }
    }

    pub fn build<'a, I>(replace: bool, iter: I) -> Self
    where
        I: IntoIterator<Item = (&'a [u8], (f64, f64))>,
    {
        let mut values = HashMap::default();

        for (b, v) in iter {
            assert!(b.len() <= 8);
            let forward = interpret_bytes(b.iter());
            let backward = interpret_bytes(b.iter().rev());

            values.insert(forward, v);
            if replace || !values.contains_key(&backward) {
                values.insert(backward, v);
            }
        }

        Self { values }
    }

    pub fn get<const N: usize>(&self, i: [u8; N]) -> Option<(f64, f64)> {
        self.get_ref(i).copied()
    }

    pub fn get_ref<const N: usize>(&self, b: [u8; N]) -> Option<&(f64, f64)> {
        self.values.get(&interpret_bytes(b.iter()))
    }

    pub fn contains_key<const N: usize>(&self, i: [u8; N]) -> bool {
        self.get(i).is_some()
    }
}

impl<const N: usize> Index<[u8; N]> for BpEnergy {
    type Output = (f64, f64);

    fn index(&self, i: [u8; N]) -> &Self::Output {
        self.get_ref(i).unwrap()
    }
}
