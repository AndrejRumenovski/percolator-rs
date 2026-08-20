#![cfg(unix)]

use std::fs;
use std::os::unix::fs::PermissionsExt;
use std::path::{Path, PathBuf};
use std::process::Command;

use percolator_rs::benchmark_result::BenchmarkResult;

struct TempDir(PathBuf);

impl TempDir {
    fn new(name: &str) -> Self {
        let path =
            std::env::temp_dir().join(format!("percolator-rs-{name}-{}", std::process::id()));
        let _ = fs::remove_dir_all(&path);
        fs::create_dir_all(&path).unwrap();
        Self(path)
    }

    fn path(&self) -> &Path {
        &self.0
    }
}

impl Drop for TempDir {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

fn write_executable(path: &Path) {
    fs::write(
        path,
        r#"#!/bin/sh
while [ "$#" -gt 0 ]; do
  case "$1" in
    --results-psms) psm="$2"; shift 2 ;;
    --results-peptides) pep="$2"; shift 2 ;;
    --results-proteins) protein="$2"; shift 2 ;;
    --seed|--num-threads|--search-input|--decoy-results-psms|--decoy-results-peptides|--decoy-results-proteins) shift 2 ;;
    --canonical) shift ;;
    *) input="$1"; shift ;;
  esac
done
if [ "$(basename "$0")" = fake-cpp ] && [ "$(basename "$input")" = fail.pin ]; then
  echo "intentional mock failure" >&2
  exit 7
fi
printf 'PSMId\tq-value\nfirst\t0.009\nedge\t0.01\n' > "$psm"
printf 'peptide\tq-value\nfirst\t0.009\nedge\t0.01\n' > "$pep"
printf 'ProteinId\tq-value\nfirst\t0.009\nedge\t0.01\n' > "$protein"
"#,
    )
    .unwrap();
    fs::set_permissions(path, fs::Permissions::from_mode(0o755)).unwrap();
}

fn write_manifest(path: &Path) {
    fs::write(
        path,
        r#"version = 1
[[datasets]]
id = "mock"
source = "fixture"
organism = "fixture"
experiment_type = "fixture"
search_engine = "fixture"
pin_path = "${BENCHMARK_RUNNER_TEST_DATA}/*.pin"
protein_level_evaluation = true
notes = "fixture"
reference_search_input = "concatenated"
"#,
    )
    .unwrap();
}

fn runner() -> String {
    std::env::var("CARGO_BIN_EXE_benchmark-dataset").expect("Cargo supplies runner path")
}

#[test]
fn previews_commands_and_records_mock_failure_without_discarding_it() {
    let temp = TempDir::new("benchmark-runner");
    let inputs = temp.path().join("inputs");
    fs::create_dir(&inputs).unwrap();
    fs::write(inputs.join("ok.pin"), "fixture").unwrap();
    fs::write(inputs.join("fail.pin"), "fixture").unwrap();
    let manifest = temp.path().join("datasets.toml");
    write_manifest(&manifest);
    let rust = temp.path().join("fake-rust");
    let cpp = temp.path().join("fake-cpp");
    write_executable(&rust);
    write_executable(&cpp);

    let preview_root = temp.path().join("preview");
    let preview = Command::new(runner())
        .args([
            "--manifest",
            manifest.to_str().unwrap(),
            "--dataset",
            "mock",
            "--output",
            preview_root.to_str().unwrap(),
            "--rust",
            rust.to_str().unwrap(),
            "--percolator",
            cpp.to_str().unwrap(),
            "--dry-run",
        ])
        .env("BENCHMARK_RUNNER_TEST_DATA", &inputs)
        .output()
        .unwrap();
    assert!(preview.status.success());
    let stdout = String::from_utf8(preview.stdout).unwrap();
    assert_eq!(stdout.lines().count(), 4);
    assert!(stdout.contains("--search-input concatenated"));
    assert!(stdout.contains("--canonical"));
    assert!(!preview_root.exists());

    let output_root = temp.path().join("results");
    let run = Command::new(runner())
        .args([
            "--manifest",
            manifest.to_str().unwrap(),
            "--dataset",
            "mock",
            "--output",
            output_root.to_str().unwrap(),
            "--rust",
            rust.to_str().unwrap(),
            "--percolator",
            cpp.to_str().unwrap(),
        ])
        .env("BENCHMARK_RUNNER_TEST_DATA", &inputs)
        .output()
        .unwrap();
    assert_eq!(run.status.code(), Some(2));
    assert!(String::from_utf8(run.stderr)
        .unwrap()
        .contains("1 file/implementation run(s) failed"));

    let result = output_root.join("mock");
    let rust_summary = fs::read_to_string(result.join("rust-summary.tsv")).unwrap();
    assert!(rust_summary.contains("\t2\t2\t0\t2\t2\t2\n"));
    let cpp_summary = fs::read_to_string(result.join("cpp-summary.tsv")).unwrap();
    assert!(cpp_summary.contains("\tnonzero\t2\t1\t1\tNA\tNA\tNA\n"));
    let per_file = fs::read_to_string(result.join("per-file.tsv")).unwrap();
    assert_eq!(per_file.lines().count(), 5);
    assert!(per_file.contains("cpp\t"));
    let failures = fs::read_to_string(result.join("failures.tsv")).unwrap();
    assert!(failures.contains("cpp"));
    assert!(failures.contains("\t7\t"));
    assert!(result.join("rust/0001/target.psms.tsv").exists());
    assert!(result.join("cpp/0002/target.psms.tsv").exists());

    let rust_result: BenchmarkResult =
        serde_json::from_str(&fs::read_to_string(result.join("rust-result.json")).unwrap())
            .unwrap();
    assert_eq!(rust_result.schema_version, 1);
    assert_eq!(rust_result.dataset_id, "mock");
    assert_eq!(rust_result.implementation, "rust");
    assert!(rust_result.benchmark_timestamp_unix_seconds > 0);
    assert_eq!(rust_result.files_attempted, 2);
    assert_eq!(rust_result.files_successful, 2);
    assert_eq!(rust_result.psms_q_lt_0_01, Some(2));
    assert_eq!(rust_result.peptides_q_lt_0_01, Some(2));
    assert_eq!(rust_result.proteins_q_lt_0_01, Some(2));
    assert_eq!(rust_result.command_line_arguments.len(), 2);
    assert_eq!(rust_result.per_file_results[0].psms_q_lt_0_01, Some(1));

    let cpp_result: BenchmarkResult =
        serde_json::from_str(&fs::read_to_string(result.join("cpp-result.json")).unwrap()).unwrap();
    assert_eq!(cpp_result.files_successful, 1);
    assert_eq!(cpp_result.failed_files.len(), 1);
    assert_eq!(cpp_result.failed_files[0].exit_status, Some(7));
    assert_eq!(cpp_result.psms_q_lt_0_01, None);
}
