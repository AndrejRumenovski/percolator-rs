//! Write a strict comparison artifact from Rust and C++ benchmark result JSON.

use percolator_rs::benchmark_comparison::compare;
use percolator_rs::benchmark_result::BenchmarkResult;
use std::path::PathBuf;
use std::time::{SystemTime, UNIX_EPOCH};

fn usage() -> ! {
    eprintln!("usage: compare-benchmark-results --rust RUST_RESULT.json --cpp CPP_RESULT.json --output COMPARISON.json");
    std::process::exit(2);
}

fn main() {
    let mut rust = None;
    let mut cpp = None;
    let mut output = None;
    let mut arguments = std::env::args().skip(1);
    while let Some(argument) = arguments.next() {
        let mut value = || PathBuf::from(arguments.next().unwrap_or_else(|| usage()));
        match argument.as_str() {
            "--rust" => rust = Some(value()),
            "--cpp" => cpp = Some(value()),
            "--output" => output = Some(value()),
            "--help" | "-h" => usage(),
            _ => usage(),
        }
    }
    let rust = load(rust.unwrap_or_else(|| usage()));
    let cpp = load(cpp.unwrap_or_else(|| usage()));
    let timestamp = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap_or_else(|_| usage())
        .as_secs();
    let comparison = compare(&rust, &cpp, timestamp).unwrap_or_else(|error| {
        eprintln!("compare-benchmark-results: {error}");
        std::process::exit(2);
    });
    let output = output.unwrap_or_else(|| usage());
    if output.exists() {
        eprintln!(
            "compare-benchmark-results: output {} already exists",
            output.display()
        );
        std::process::exit(2);
    }
    let file = std::fs::File::create(&output).unwrap_or_else(|error| {
        eprintln!(
            "compare-benchmark-results: could not create {}: {error}",
            output.display()
        );
        std::process::exit(2);
    });
    serde_json::to_writer_pretty(file, &comparison).unwrap_or_else(|error| {
        eprintln!(
            "compare-benchmark-results: could not write {}: {error}",
            output.display()
        );
        std::process::exit(2);
    });
}

fn load(path: PathBuf) -> BenchmarkResult {
    let file = std::fs::File::open(&path).unwrap_or_else(|error| {
        eprintln!(
            "compare-benchmark-results: could not read {}: {error}",
            path.display()
        );
        std::process::exit(2);
    });
    serde_json::from_reader(file).unwrap_or_else(|error| {
        eprintln!(
            "compare-benchmark-results: invalid result {}: {error}",
            path.display()
        );
        std::process::exit(2);
    })
}
