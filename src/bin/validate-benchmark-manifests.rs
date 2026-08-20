//! Validate benchmark dataset metadata without running a benchmark.

use percolator_rs::benchmark_manifest::DatasetRegistry;

fn main() {
    let path = std::env::args()
        .nth(1)
        .unwrap_or_else(|| "bench/datasets.toml".to_owned());
    match DatasetRegistry::load(&path) {
        Ok(registry) => println!(
            "valid dataset manifest: {} dataset(s) in {path}",
            registry.datasets.len()
        ),
        Err(error) => {
            eprintln!("{error}");
            std::process::exit(2);
        }
    }
}
