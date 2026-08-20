//! Comparison of compatible Rust and C++ benchmark result artifacts.

use crate::benchmark_result::{BenchmarkResult, PerFileResult, RESULT_SCHEMA_VERSION};
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};
use std::fmt;

/// Increment only for incompatible changes to [`BenchmarkComparison`].
pub const COMPARISON_SCHEMA_VERSION: u32 = 1;

/// A versioned comparison suitable for machine ingestion.
#[derive(Debug, Serialize, Deserialize, PartialEq)]
pub struct BenchmarkComparison {
    pub schema_version: u32,
    pub comparison_timestamp_unix_seconds: u64,
    pub dataset_id: String,
    pub dataset_accession: Option<String>,
    pub rust_run: ComparisonRunMetadata,
    pub cpp_run: ComparisonRunMetadata,
    pub metrics: ComparisonMetrics,
    /// Yield language is intentional: count differences do not establish FDR calibration or accuracy.
    pub interpretation: String,
}

/// Provenance retained so consumers can audit whether a comparison was fair.
#[derive(Debug, Serialize, Deserialize, PartialEq)]
pub struct ComparisonRunMetadata {
    pub result_schema_version: u32,
    pub implementation: String,
    pub percolator_rs_git_commit: Option<String>,
    pub rust_compiler_version: Option<String>,
    pub cpp_percolator_version: Option<String>,
    pub os: Option<String>,
    pub cpu: Option<String>,
    pub available_threads: Option<usize>,
    pub configured_concurrency: usize,
    pub random_seed: u64,
    pub command_line_arguments: Vec<Vec<String>>,
}

/// Aggregate and per-file identification/yield differences.
#[derive(Debug, Serialize, Deserialize, PartialEq)]
pub struct ComparisonMetrics {
    /// C++ wall time divided by Rust wall time; `None` when either measurement is unusable.
    pub runtime_speedup: Option<f64>,
    /// C++ peak RSS divided by Rust peak RSS; `None` when either measurement is unusable.
    pub peak_memory_ratio: Option<f64>,
    pub psms: CountDifference,
    pub peptides: CountDifference,
    pub proteins: Option<CountDifference>,
    pub files_rust_identifies_more: usize,
    pub files_cpp_identifies_more: usize,
    pub files_tied: usize,
    pub median_per_file_psm_difference: Option<f64>,
    pub median_per_file_peptide_difference: Option<f64>,
    /// Most negative Rust-minus-C++ per-file PSM difference, if Rust is lower on any file.
    pub worst_rust_psm_deficit: Option<i64>,
    /// Largest positive Rust-minus-C++ per-file PSM difference, if Rust is higher on any file.
    pub largest_rust_psm_gain: Option<i64>,
}

/// Rust-minus-C++ count difference. Percentage uses C++ as the denominator.
#[derive(Debug, Serialize, Deserialize, PartialEq)]
pub struct CountDifference {
    pub rust_count: Option<u64>,
    pub cpp_count: Option<u64>,
    pub absolute_difference: Option<i64>,
    pub percentage_difference: Option<f64>,
}

/// A comparison cannot proceed safely.
#[derive(Debug, PartialEq, Eq)]
pub struct ComparisonError {
    reasons: Vec<String>,
}

impl fmt::Display for ComparisonError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "cannot compare incompatible benchmark results:")?;
        for reason in &self.reasons {
            writeln!(f, "  - {reason}")?;
        }
        Ok(())
    }
}

impl std::error::Error for ComparisonError {}

/// Compare two complete, compatible results. Rust-minus-C++ differences are used throughout.
pub fn compare(
    rust: &BenchmarkResult,
    cpp: &BenchmarkResult,
    comparison_timestamp_unix_seconds: u64,
) -> Result<BenchmarkComparison, ComparisonError> {
    let reasons = compatibility_errors(rust, cpp);
    if !reasons.is_empty() {
        return Err(ComparisonError { reasons });
    }

    let rust_files = file_map(rust).expect("validated Rust result");
    let cpp_files = file_map(cpp).expect("validated C++ result");
    let mut psm_differences = Vec::new();
    let mut peptide_differences = Vec::new();
    let mut rust_more = 0;
    let mut cpp_more = 0;
    let mut ties = 0;
    for (input, rust_file) in &rust_files {
        let cpp_file = cpp_files.get(input).expect("validated matching input set");
        let psm_difference = signed_difference(
            rust_file.psms_q_lt_0_01.expect("validated PSM result"),
            cpp_file.psms_q_lt_0_01.expect("validated PSM result"),
        );
        let peptide_difference = signed_difference(
            rust_file
                .peptides_q_lt_0_01
                .expect("validated peptide result"),
            cpp_file
                .peptides_q_lt_0_01
                .expect("validated peptide result"),
        );
        match psm_difference.cmp(&0) {
            std::cmp::Ordering::Greater => rust_more += 1,
            std::cmp::Ordering::Less => cpp_more += 1,
            std::cmp::Ordering::Equal => ties += 1,
        }
        psm_differences.push(psm_difference);
        peptide_differences.push(peptide_difference);
    }

    let proteins = match (rust.proteins_q_lt_0_01, cpp.proteins_q_lt_0_01) {
        (Some(rust_count), Some(cpp_count)) => Some(count_difference(rust_count, cpp_count)),
        _ => None,
    };
    let metrics = ComparisonMetrics {
        runtime_speedup: ratio(cpp.wall_seconds, rust.wall_seconds),
        peak_memory_ratio: ratio_u64(cpp.peak_rss_kb, rust.peak_rss_kb),
        psms: optional_count_difference(rust.psms_q_lt_0_01, cpp.psms_q_lt_0_01),
        peptides: optional_count_difference(rust.peptides_q_lt_0_01, cpp.peptides_q_lt_0_01),
        proteins,
        files_rust_identifies_more: rust_more,
        files_cpp_identifies_more: cpp_more,
        files_tied: ties,
        median_per_file_psm_difference: median(&mut psm_differences),
        median_per_file_peptide_difference: median(&mut peptide_differences),
        worst_rust_psm_deficit: psm_differences
            .iter()
            .copied()
            .filter(|difference| *difference < 0)
            .min(),
        largest_rust_psm_gain: psm_differences
            .iter()
            .copied()
            .filter(|difference| *difference > 0)
            .max(),
    };
    Ok(BenchmarkComparison {
        schema_version: COMPARISON_SCHEMA_VERSION,
        comparison_timestamp_unix_seconds,
        dataset_id: rust.dataset_id.clone(),
        dataset_accession: rust.dataset_accession.clone(),
        rust_run: metadata(rust),
        cpp_run: metadata(cpp),
        metrics,
        interpretation: "Identification/yield differences at reported q < 0.01; they are not evidence of higher accuracy or calibrated FDR without separate calibration evidence.".to_owned(),
    })
}

fn compatibility_errors(rust: &BenchmarkResult, cpp: &BenchmarkResult) -> Vec<String> {
    let mut reasons = Vec::new();
    if rust.implementation != "rust" {
        reasons.push(format!(
            "first result implementation must be 'rust', got {:?}",
            rust.implementation
        ));
    }
    if cpp.implementation != "cpp" {
        reasons.push(format!(
            "second result implementation must be 'cpp', got {:?}",
            cpp.implementation
        ));
    }
    compare_required(
        &mut reasons,
        "dataset ID",
        &rust.dataset_id,
        &cpp.dataset_id,
    );
    compare_optional(
        &mut reasons,
        "dataset accession",
        &rust.dataset_accession,
        &cpp.dataset_accession,
    );
    compare_required(
        &mut reasons,
        "random seed",
        &rust.random_seed,
        &cpp.random_seed,
    );
    compare_required(
        &mut reasons,
        "configured concurrency",
        &rust.configured_concurrency,
        &cpp.configured_concurrency,
    );
    compare_optional(
        &mut reasons,
        "percolator-rs git commit",
        &rust.percolator_rs_git_commit,
        &cpp.percolator_rs_git_commit,
    );
    compare_optional(
        &mut reasons,
        "Rust compiler version",
        &rust.rust_compiler_version,
        &cpp.rust_compiler_version,
    );
    compare_optional(
        &mut reasons,
        "C++ Percolator version",
        &rust.cpp_percolator_version,
        &cpp.cpp_percolator_version,
    );
    compare_optional(&mut reasons, "OS", &rust.os, &cpp.os);
    compare_optional(&mut reasons, "CPU", &rust.cpu, &cpp.cpu);
    compare_optional(
        &mut reasons,
        "available threads",
        &rust.available_threads,
        &cpp.available_threads,
    );
    validate_completed(&mut reasons, "Rust", rust);
    validate_completed(&mut reasons, "C++", cpp);
    match (file_map(rust), file_map(cpp)) {
        (Ok(rust_files), Ok(cpp_files)) => {
            let rust_inputs: BTreeSet<_> = rust_files.keys().collect();
            let cpp_inputs: BTreeSet<_> = cpp_files.keys().collect();
            let only_rust: Vec<_> = rust_inputs
                .difference(&cpp_inputs)
                .take(3)
                .map(|input| **input)
                .collect();
            let only_cpp: Vec<_> = cpp_inputs
                .difference(&rust_inputs)
                .take(3)
                .map(|input| **input)
                .collect();
            if !only_rust.is_empty() || !only_cpp.is_empty() {
                reasons.push(format!(
                    "input-file sets differ (only Rust: {}; only C++: {})",
                    only_rust.join(", "),
                    only_cpp.join(", ")
                ));
            }
            for (input, rust_file) in rust_files {
                if let Some(cpp_file) = cpp_files.get(input) {
                    require_metric(
                        &mut reasons,
                        input,
                        "PSM",
                        rust_file.psms_q_lt_0_01,
                        cpp_file.psms_q_lt_0_01,
                    );
                    require_metric(
                        &mut reasons,
                        input,
                        "peptide",
                        rust_file.peptides_q_lt_0_01,
                        cpp_file.peptides_q_lt_0_01,
                    );
                }
            }
        }
        (Err(reason), _) | (_, Err(reason)) => reasons.push(reason),
    }
    reasons
}

fn validate_completed(reasons: &mut Vec<String>, label: &str, result: &BenchmarkResult) {
    if result.schema_version != RESULT_SCHEMA_VERSION {
        reasons.push(format!(
            "{label} result schema version {} is unsupported (expected {RESULT_SCHEMA_VERSION})",
            result.schema_version
        ));
    }
    if result.files_attempted != result.per_file_results.len() {
        reasons.push(format!(
            "{label} files_attempted does not match per_file_results length"
        ));
    }
    if result.files_successful != result.files_attempted || !result.failed_files.is_empty() {
        reasons.push(format!(
            "{label} result contains failed or incomplete files"
        ));
    }
    if result
        .per_file_results
        .iter()
        .any(|file| file.failure.is_some() || file.exit_status != Some(0))
    {
        reasons.push(format!("{label} per-file results contain a failed process"));
    }
    if result.psms_q_lt_0_01.is_none() || result.peptides_q_lt_0_01.is_none() {
        reasons.push(format!(
            "{label} result is missing aggregate PSM or peptide counts"
        ));
    }
}

fn file_map(result: &BenchmarkResult) -> Result<BTreeMap<&str, &PerFileResult>, String> {
    let mut files = BTreeMap::new();
    for file in &result.per_file_results {
        if files.insert(file.input.as_str(), file).is_some() {
            return Err(format!(
                "{} result contains duplicate input file {:?}",
                result.implementation, file.input
            ));
        }
    }
    Ok(files)
}

fn require_metric(
    reasons: &mut Vec<String>,
    input: &str,
    name: &str,
    rust: Option<u64>,
    cpp: Option<u64>,
) {
    if rust.is_none() || cpp.is_none() {
        reasons.push(format!("{name} count is missing for input {input:?}"));
    }
}

fn compare_required<T: PartialEq + fmt::Debug>(
    reasons: &mut Vec<String>,
    name: &str,
    rust: &T,
    cpp: &T,
) {
    if rust != cpp {
        reasons.push(format!("{name} differs (Rust {rust:?}; C++ {cpp:?})"));
    }
}

fn compare_optional<T: PartialEq + fmt::Debug>(
    reasons: &mut Vec<String>,
    name: &str,
    rust: &Option<T>,
    cpp: &Option<T>,
) {
    if let (Some(rust), Some(cpp)) = (rust, cpp) {
        compare_required(reasons, name, rust, cpp);
    }
}

fn metadata(result: &BenchmarkResult) -> ComparisonRunMetadata {
    ComparisonRunMetadata {
        result_schema_version: result.schema_version,
        implementation: result.implementation.clone(),
        percolator_rs_git_commit: result.percolator_rs_git_commit.clone(),
        rust_compiler_version: result.rust_compiler_version.clone(),
        cpp_percolator_version: result.cpp_percolator_version.clone(),
        os: result.os.clone(),
        cpu: result.cpu.clone(),
        available_threads: result.available_threads,
        configured_concurrency: result.configured_concurrency,
        random_seed: result.random_seed,
        command_line_arguments: result.command_line_arguments.clone(),
    }
}

fn optional_count_difference(rust: Option<u64>, cpp: Option<u64>) -> CountDifference {
    match (rust, cpp) {
        (Some(rust), Some(cpp)) => count_difference(rust, cpp),
        _ => CountDifference {
            rust_count: rust,
            cpp_count: cpp,
            absolute_difference: None,
            percentage_difference: None,
        },
    }
}

fn count_difference(rust: u64, cpp: u64) -> CountDifference {
    let absolute_difference = signed_difference(rust, cpp);
    CountDifference {
        rust_count: Some(rust),
        cpp_count: Some(cpp),
        absolute_difference: Some(absolute_difference),
        percentage_difference: (cpp != 0).then(|| absolute_difference as f64 / cpp as f64 * 100.0),
    }
}

fn signed_difference(rust: u64, cpp: u64) -> i64 {
    i64::try_from(rust).unwrap_or(i64::MAX) - i64::try_from(cpp).unwrap_or(i64::MAX)
}

fn ratio(numerator: Option<f64>, denominator: Option<f64>) -> Option<f64> {
    match (numerator, denominator) {
        (Some(numerator), Some(denominator)) if denominator != 0.0 => Some(numerator / denominator),
        _ => None,
    }
}

fn ratio_u64(numerator: Option<u64>, denominator: Option<u64>) -> Option<f64> {
    match (numerator, denominator) {
        (Some(numerator), Some(denominator)) if denominator != 0 => {
            Some(numerator as f64 / denominator as f64)
        }
        _ => None,
    }
}

fn median(values: &mut [i64]) -> Option<f64> {
    if values.is_empty() {
        return None;
    }
    values.sort_unstable();
    let middle = values.len() / 2;
    if values.len() % 2 == 1 {
        Some(values[middle] as f64)
    } else {
        Some((values[middle - 1] as f64 + values[middle] as f64) / 2.0)
    }
}
