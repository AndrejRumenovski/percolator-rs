use percolator_rs::benchmark_result::{
    BenchmarkResult, FailedFile, PerFileResult, RESULT_SCHEMA_VERSION,
};

#[test]
fn result_json_round_trips_without_losing_missing_values() {
    let result = BenchmarkResult {
        schema_version: RESULT_SCHEMA_VERSION,
        benchmark_timestamp_unix_seconds: 1_700_000_000,
        dataset_id: "fixture".to_owned(),
        dataset_accession: None,
        implementation: "rust".to_owned(),
        percolator_rs_git_commit: Some("abc123".to_owned()),
        rust_compiler_version: Some("rustc fixture".to_owned()),
        cpp_percolator_version: None,
        os: Some("fixture-os".to_owned()),
        cpu: None,
        available_threads: Some(4),
        configured_concurrency: 1,
        random_seed: 1,
        command_line_arguments: vec![vec!["percolator-rs".to_owned(), "input.pin".to_owned()]],
        wall_seconds: Some(1.5),
        peak_rss_kb: None,
        files_attempted: 1,
        files_successful: 0,
        failed_files: vec![FailedFile {
            input: "input.pin".to_owned(),
            exit_status: Some(7),
            termination: None,
            stderr_log: "stderr.log".to_owned(),
            reason: "fixture failure".to_owned(),
        }],
        psms_q_lt_0_01: None,
        peptides_q_lt_0_01: None,
        proteins_q_lt_0_01: None,
        per_file_results: vec![PerFileResult {
            input: "input.pin".to_owned(),
            output_dir: "output".to_owned(),
            command_line_arguments: vec!["percolator-rs".to_owned(), "input.pin".to_owned()],
            exit_status: Some(7),
            termination: None,
            wall_seconds: Some(1.5),
            peak_rss_kb: None,
            psms_q_lt_0_01: None,
            peptides_q_lt_0_01: None,
            proteins_q_lt_0_01: None,
            failure: Some("fixture failure".to_owned()),
        }],
    };

    let json = serde_json::to_string(&result).unwrap();
    assert!(json.contains("\"peak_rss_kb\":null"));
    let decoded: BenchmarkResult = serde_json::from_str(&json).unwrap();
    assert_eq!(decoded, result);
}
