use percolator_rs::benchmark_comparison::{compare, BenchmarkComparison};
use percolator_rs::benchmark_result::{BenchmarkResult, PerFileResult, RESULT_SCHEMA_VERSION};

fn file(input: &str, psms: u64, peptides: u64) -> PerFileResult {
    PerFileResult {
        input: input.to_owned(),
        output_dir: format!("out/{input}"),
        command_line_arguments: vec!["percolator".to_owned(), input.to_owned()],
        exit_status: Some(0),
        termination: None,
        wall_seconds: Some(1.0),
        peak_rss_kb: Some(10),
        psms_q_lt_0_01: Some(psms),
        peptides_q_lt_0_01: Some(peptides),
        proteins_q_lt_0_01: None,
        failure: None,
    }
}

fn result(implementation: &str, files: Vec<PerFileResult>) -> BenchmarkResult {
    let attempted = files.len();
    BenchmarkResult {
        schema_version: RESULT_SCHEMA_VERSION,
        benchmark_timestamp_unix_seconds: 1,
        dataset_id: "fixture".to_owned(),
        dataset_accession: Some("PXDfixture".to_owned()),
        implementation: implementation.to_owned(),
        percolator_rs_git_commit: Some("commit".to_owned()),
        rust_compiler_version: Some("rustc".to_owned()),
        cpp_percolator_version: Some("percolator".to_owned()),
        os: Some("linux".to_owned()),
        cpu: Some("cpu".to_owned()),
        available_threads: Some(4),
        configured_concurrency: 1,
        random_seed: 1,
        command_line_arguments: files
            .iter()
            .map(|file| file.command_line_arguments.clone())
            .collect(),
        wall_seconds: Some(if implementation == "rust" { 10.0 } else { 20.0 }),
        peak_rss_kb: Some(if implementation == "rust" { 100 } else { 200 }),
        files_attempted: attempted,
        files_successful: attempted,
        failed_files: Vec::new(),
        psms_q_lt_0_01: Some(files.iter().map(|file| file.psms_q_lt_0_01.unwrap()).sum()),
        peptides_q_lt_0_01: Some(
            files
                .iter()
                .map(|file| file.peptides_q_lt_0_01.unwrap())
                .sum(),
        ),
        proteins_q_lt_0_01: None,
        per_file_results: files,
    }
}

#[test]
fn compares_complete_matching_results_and_round_trips_json() {
    let mut rust = result(
        "rust",
        vec![
            file("a.pin", 110, 10),
            file("b.pin", 90, 20),
            file("c.pin", 100, 30),
        ],
    );
    let mut cpp = result(
        "cpp",
        vec![
            file("a.pin", 100, 9),
            file("b.pin", 100, 21),
            file("c.pin", 100, 30),
        ],
    );
    rust.proteins_q_lt_0_01 = Some(3);
    cpp.proteins_q_lt_0_01 = Some(2);
    let comparison = compare(&rust, &cpp, 42).unwrap();

    assert_eq!(comparison.dataset_id, "fixture");
    assert_eq!(comparison.metrics.runtime_speedup, Some(2.0));
    assert_eq!(comparison.metrics.peak_memory_ratio, Some(2.0));
    assert_eq!(comparison.metrics.psms.absolute_difference, Some(0));
    assert_eq!(comparison.metrics.psms.percentage_difference, Some(0.0));
    assert_eq!(comparison.metrics.peptides.absolute_difference, Some(0));
    assert_eq!(
        comparison
            .metrics
            .proteins
            .as_ref()
            .unwrap()
            .absolute_difference,
        Some(1)
    );
    assert_eq!(comparison.metrics.files_rust_identifies_more, 1);
    assert_eq!(comparison.metrics.files_cpp_identifies_more, 1);
    assert_eq!(comparison.metrics.files_tied, 1);
    assert_eq!(comparison.metrics.median_per_file_psm_difference, Some(0.0));
    assert_eq!(
        comparison.metrics.median_per_file_peptide_difference,
        Some(0.0)
    );
    assert_eq!(comparison.metrics.worst_rust_psm_deficit, Some(-10));
    assert_eq!(comparison.metrics.largest_rust_psm_gain, Some(10));
    assert!(comparison
        .interpretation
        .contains("not evidence of higher accuracy"));

    let json = serde_json::to_string(&comparison).unwrap();
    let decoded: BenchmarkComparison = serde_json::from_str(&json).unwrap();
    assert_eq!(decoded, comparison);
}

#[test]
fn handles_zero_denominators_and_missing_proteins_explicitly() {
    let rust = result("rust", vec![file("a.pin", 2, 1)]);
    let mut cpp = result("cpp", vec![file("a.pin", 0, 1)]);
    cpp.psms_q_lt_0_01 = Some(0);
    cpp.per_file_results[0].psms_q_lt_0_01 = Some(0);
    let comparison = compare(&rust, &cpp, 42).unwrap();
    assert_eq!(comparison.metrics.psms.absolute_difference, Some(2));
    assert_eq!(comparison.metrics.psms.percentage_difference, None);
    assert_eq!(comparison.metrics.proteins, None);
}

#[test]
fn rejects_dataset_seed_input_and_failure_mismatches() {
    let rust = result("rust", vec![file("a.pin", 1, 1)]);
    let mut cpp = result("cpp", vec![file("b.pin", 1, 1)]);
    cpp.dataset_id = "other".to_owned();
    cpp.random_seed = 2;
    cpp.files_successful = 0;
    cpp.per_file_results[0].exit_status = Some(7);
    cpp.per_file_results[0].failure = Some("failed".to_owned());
    let error = compare(&rust, &cpp, 42).expect_err("incompatible results must fail");
    let message = error.to_string();
    assert!(message.contains("dataset ID differs"));
    assert!(message.contains("random seed differs"));
    assert!(message.contains("failed or incomplete files"));
    assert!(message.contains("input-file sets differ"));
}
