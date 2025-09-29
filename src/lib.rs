#![doc = include_str!("../README.md")]

pub mod dna;
pub mod fold;
pub mod rna;
pub mod tm;
pub mod types;

pub use dna::dna;
pub use rna::rna;

mod util {
    /// Round to one decimal
    pub fn round1(n: f64) -> f64 {
        (n * 10.0).round() / 10.0
    }

    /// Round to two decimal
    pub fn round2(n: f64) -> f64 {
        (n * 100.0).round() / 100.0
    }
}
