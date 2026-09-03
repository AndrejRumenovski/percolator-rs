//! Reusable scientific pipeline and benchmark support.

pub mod benchmark_comparison;
pub mod benchmark_manifest;
pub mod benchmark_result;
pub mod competition;
pub mod mlp;
pub mod output;
pub mod peptide;
pub mod percolator;
pub mod pin;
pub mod pipeline;
mod preprocessing;
#[cfg(feature = "profiling")]
pub mod profile;
pub mod protein;
pub mod protein_bayes;
pub mod rt;
pub mod simd;
pub mod stats;
pub mod svm;
pub mod tiebreak;
