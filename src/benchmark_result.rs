//! Stable, machine-readable records emitted by benchmark tooling.

use serde::{Deserialize, Serialize};

/// Increment only for incompatible changes to [`BenchmarkResult`].
pub const RESULT_SCHEMA_VERSION: u32 = 1;

/// Complete result for one implementation on one dataset.
#[derive(Debug, Serialize, Deserialize, PartialEq)]
pub struct BenchmarkResult {
    pub schema_version: u32,
    pub benchmark_timestamp_unix_seconds: u64,
    pub dataset_id: String,
    pub dataset_accession: Option<String>,
    pub implementation: String,
    pub percolator_rs_git_commit: Option<String>,
    pub rust_compiler_version: Option<String>,
    pub cpp_percolator_version: Option<String>,
    pub os: Option<String>,
    pub cpu: Option<String>,
    pub available_threads: Option<usize>,
    pub configured_concurrency: usize,
    pub random_seed: u64,
    /// Exact Percolator command lines, one for each attempted input file.
    pub command_line_arguments: Vec<Vec<String>>,
    pub wall_seconds: Option<f64>,
    pub peak_rss_kb: Option<u64>,
    pub files_attempted: usize,
    pub files_successful: usize,
    pub failed_files: Vec<FailedFile>,
    pub psms_q_lt_0_01: Option<u64>,
    pub peptides_q_lt_0_01: Option<u64>,
    pub proteins_q_lt_0_01: Option<u64>,
    pub per_file_results: Vec<PerFileResult>,
}

/// A failed process or unreadable expected result file.
#[derive(Debug, Serialize, Deserialize, PartialEq)]
pub struct FailedFile {
    pub input: String,
    pub exit_status: Option<i32>,
    pub termination: Option<String>,
    pub stderr_log: String,
    pub reason: String,
}

/// Metrics and provenance for one attempted input file.
#[derive(Debug, Serialize, Deserialize, PartialEq)]
pub struct PerFileResult {
    pub input: String,
    pub output_dir: String,
    pub command_line_arguments: Vec<String>,
    pub exit_status: Option<i32>,
    pub termination: Option<String>,
    pub wall_seconds: Option<f64>,
    pub peak_rss_kb: Option<u64>,
    pub psms_q_lt_0_01: Option<u64>,
    pub peptides_q_lt_0_01: Option<u64>,
    pub proteins_q_lt_0_01: Option<u64>,
    pub failure: Option<String>,
}
