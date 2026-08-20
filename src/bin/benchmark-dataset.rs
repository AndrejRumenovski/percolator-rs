//! Run both Percolator implementations for one configured benchmark dataset.

use glob::glob;
use percolator_rs::benchmark_manifest::{expand_environment_templates, Dataset, DatasetRegistry};
use std::fs::{self, File};
use std::io::{self, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::process::{Command, ExitStatus, Stdio};
use std::time::Instant;

const SEED: &str = "1";

struct Args {
    manifest: PathBuf,
    dataset: String,
    output: PathBuf,
    rust: PathBuf,
    reference: PathBuf,
    dry_run: bool,
}

struct Implementation<'a> {
    name: &'a str,
    binary: &'a Path,
    reference: bool,
}

struct FileResult {
    implementation: String,
    input: PathBuf,
    output_dir: PathBuf,
    exit_status: String,
    wall_seconds: Option<f64>,
    peak_rss_kb: Option<u64>,
    psms: Option<u64>,
    peptides: Option<u64>,
    proteins: Option<u64>,
    failure: Option<String>,
}

fn usage() -> ! {
    eprintln!(
        "usage: benchmark-dataset --dataset ID --output DIR --percolator PATH [options]\n\
         \noptions:\n\
           --manifest PATH       Registry path (default: bench/datasets.toml)\n\
           --rust PATH           percolator-rs binary (default: target/release/percolator-rs)\n\
           --percolator PATH     Reference C++ Percolator binary\n\
           --dry-run             Print exact commands without creating outputs or executing\n\
           --help                Show this help"
    );
    std::process::exit(2);
}

fn parse_args() -> Args {
    let mut manifest = PathBuf::from("bench/datasets.toml");
    let mut dataset = None;
    let mut output = None;
    let mut rust = PathBuf::from("target/release/percolator-rs");
    let mut reference = None;
    let mut dry_run = false;
    let mut args = std::env::args().skip(1);
    while let Some(arg) = args.next() {
        match arg.as_str() {
            "--manifest" => manifest = PathBuf::from(args.next().unwrap_or_else(|| usage())),
            "--dataset" => dataset = args.next().unwrap_or_else(|| usage()).into(),
            "--output" => output = PathBuf::from(args.next().unwrap_or_else(|| usage())).into(),
            "--rust" => rust = PathBuf::from(args.next().unwrap_or_else(|| usage())),
            "--percolator" => {
                reference = PathBuf::from(args.next().unwrap_or_else(|| usage())).into()
            }
            "--dry-run" => dry_run = true,
            "--help" | "-h" => usage(),
            _ => usage(),
        }
    }
    Args {
        manifest,
        dataset: dataset.unwrap_or_else(|| usage()),
        output: output.unwrap_or_else(|| usage()),
        rust,
        reference: reference.unwrap_or_else(|| usage()),
        dry_run,
    }
}

fn main() {
    let args = parse_args();
    if let Err(error) = run(args) {
        eprintln!("benchmark-dataset: {error}");
        std::process::exit(2);
    }
}

fn run(args: Args) -> Result<(), String> {
    let registry = DatasetRegistry::load(&args.manifest).map_err(|error| error.to_string())?;
    let dataset = registry
        .datasets
        .iter()
        .find(|dataset| dataset.id == args.dataset)
        .ok_or_else(|| {
            format!(
                "dataset {:?} is not in {}",
                args.dataset,
                args.manifest.display()
            )
        })?;
    let inputs = discover_inputs(dataset)?;
    if let Some(expected) = dataset.file_count {
        if inputs.len() != expected {
            return Err(format!(
                "dataset {} expects {expected} PIN files, but {} matched {}; refusing to run a partial benchmark",
                dataset.id,
                inputs.len(),
                dataset.pin_path
            ));
        }
    }

    let implementations = [
        Implementation {
            name: "rust",
            binary: &args.rust,
            reference: false,
        },
        Implementation {
            name: "cpp",
            binary: &args.reference,
            reference: true,
        },
    ];

    if args.dry_run {
        for implementation in &implementations {
            for (index, input) in inputs.iter().enumerate() {
                let command = command_for(
                    implementation,
                    dataset,
                    input,
                    &args
                        .output
                        .join(&dataset.id)
                        .join(implementation.name)
                        .join(format!("{:04}", index + 1)),
                );
                println!(
                    "{}",
                    display_command(
                        &command,
                        &args
                            .output
                            .join(&dataset.id)
                            .join(implementation.name)
                            .join(format!("{:04}", index + 1))
                            .join("time.tsv")
                    )
                );
            }
        }
        return Ok(());
    }

    let root = args.output.join(&dataset.id);
    if root.exists() {
        return Err(format!(
            "output directory {} already exists; choose a new --output directory to avoid mixing runs",
            root.display()
        ));
    }
    fs::create_dir_all(&root).map_err(io_error)?;
    let mut all_results = Vec::new();
    for implementation in &implementations {
        let start = Instant::now();
        for (index, input) in inputs.iter().enumerate() {
            let output_dir = root
                .join(implementation.name)
                .join(format!("{:04}", index + 1));
            let result = execute_one(implementation, dataset, input, output_dir)?;
            all_results.push(result);
        }
        let wall = start.elapsed().as_secs_f64();
        write_summary(&root, implementation.name, wall, &all_results)?;
    }
    write_per_file(&root, &all_results)?;
    write_failures(&root, &all_results)?;

    let failed = all_results
        .iter()
        .filter(|result| result.failure.is_some())
        .count();
    if failed > 0 {
        return Err(format!(
            "{failed} file/implementation run(s) failed; inspect {}/failures.tsv",
            root.display()
        ));
    }
    Ok(())
}

fn discover_inputs(dataset: &Dataset) -> Result<Vec<PathBuf>, String> {
    let pattern = expand_environment_templates(&dataset.pin_path)?;
    let mut inputs = Vec::new();
    for path in glob(&pattern).map_err(|error| format!("invalid PIN glob {pattern:?}: {error}"))? {
        let path =
            path.map_err(|error| format!("could not expand PIN glob {pattern:?}: {error}"))?;
        if path.is_file() {
            inputs.push(path);
        }
    }
    inputs.sort();
    if inputs.is_empty() {
        return Err(format!("PIN glob {pattern:?} matched no files"));
    }
    Ok(inputs)
}

fn command_for(
    implementation: &Implementation<'_>,
    dataset: &Dataset,
    input: &Path,
    output: &Path,
) -> Vec<String> {
    let mut command = vec![
        implementation.binary.display().to_string(),
        "--seed".to_owned(),
        SEED.to_owned(),
    ];
    if implementation.reference {
        command.extend(["--num-threads".to_owned(), "1".to_owned()]);
        if let Some(search_input) = &dataset.reference_search_input {
            command.extend([
                "--search-input".to_owned(),
                search_input.as_str().to_owned(),
            ]);
        }
    } else {
        command.extend([
            "--canonical".to_owned(),
            "--num-threads".to_owned(),
            "1".to_owned(),
        ]);
    }
    command.extend([
        "--results-psms".to_owned(),
        output.join("target.psms.tsv").display().to_string(),
        "--decoy-results-psms".to_owned(),
        output.join("decoy.psms.tsv").display().to_string(),
        "--results-peptides".to_owned(),
        output.join("target.peptides.tsv").display().to_string(),
        "--decoy-results-peptides".to_owned(),
        output.join("decoy.peptides.tsv").display().to_string(),
    ]);
    if dataset.protein_level_evaluation {
        command.extend([
            "--results-proteins".to_owned(),
            output.join("target.proteins.tsv").display().to_string(),
            "--decoy-results-proteins".to_owned(),
            output.join("decoy.proteins.tsv").display().to_string(),
        ]);
    }
    command.push(input.display().to_string());
    command
}

fn execute_one(
    implementation: &Implementation<'_>,
    dataset: &Dataset,
    input: &Path,
    output_dir: PathBuf,
) -> Result<FileResult, String> {
    fs::create_dir_all(&output_dir).map_err(io_error)?;
    let command = command_for(implementation, dataset, input, &output_dir);
    let stdout = File::create(output_dir.join("stdout.log")).map_err(io_error)?;
    let stderr = File::create(output_dir.join("stderr.log")).map_err(io_error)?;
    let time_path = output_dir.join("time.tsv");

    let status = Command::new("/usr/bin/time")
        .args(["-f", "%e\\t%M", "-o"])
        .arg(&time_path)
        .arg(&command[0])
        .args(&command[1..])
        .stdout(Stdio::from(stdout))
        .stderr(Stdio::from(stderr))
        .status()
        .map_err(|error| {
            format!(
                "could not start {} for {}: {error}",
                implementation.name,
                input.display()
            )
        })?;

    let mut result = FileResult {
        implementation: implementation.name.to_owned(),
        input: input.to_path_buf(),
        output_dir: output_dir.clone(),
        exit_status: status_label(&status),
        wall_seconds: None,
        peak_rss_kb: None,
        psms: None,
        peptides: None,
        proteins: None,
        failure: None,
    };
    match read_time(&time_path) {
        Ok((wall, rss)) => {
            result.wall_seconds = Some(wall);
            result.peak_rss_kb = Some(rss);
        }
        Err(error) => result.failure = Some(error),
    }
    if !status.success() {
        result
            .failure
            .get_or_insert_with(|| format!("process exited with {}", result.exit_status));
        return Ok(result);
    }
    match count_q_lt_001(&output_dir.join("target.psms.tsv")) {
        Ok(count) => result.psms = Some(count),
        Err(error) => result.failure = Some(format!("{error} ({})", input.display())),
    }
    match count_q_lt_001(&output_dir.join("target.peptides.tsv")) {
        Ok(count) => result.peptides = Some(count),
        Err(error) => {
            result.failure.get_or_insert(error);
        }
    }
    if result.psms.is_none() || result.peptides.is_none() {
        result.failure.get_or_insert_with(|| {
            "successful process did not produce readable PSM and peptide outputs".to_owned()
        });
    }
    if dataset.protein_level_evaluation {
        match count_q_lt_001(&output_dir.join("target.proteins.tsv")) {
            Ok(count) => result.proteins = Some(count),
            Err(error) => {
                result.failure.get_or_insert(error);
            }
        }
    }
    Ok(result)
}

fn read_time(path: &Path) -> Result<(f64, u64), String> {
    let text = fs::read_to_string(path).map_err(io_error)?;
    let mut fields = text.trim().split('\t');
    let wall = fields
        .next()
        .ok_or_else(|| format!("{} has no wall-time value", path.display()))?
        .parse()
        .map_err(|_| format!("{} has invalid wall-time value", path.display()))?;
    let rss = fields
        .next()
        .ok_or_else(|| format!("{} has no peak-RSS value", path.display()))?
        .parse()
        .map_err(|_| format!("{} has invalid peak-RSS value", path.display()))?;
    Ok((wall, rss))
}

fn count_q_lt_001(path: &Path) -> Result<u64, String> {
    let text = fs::read_to_string(path).map_err(io_error)?;
    let mut rows = text.lines();
    let header = rows
        .next()
        .ok_or_else(|| format!("{} is empty", path.display()))?;
    let q_column = header
        .split('\t')
        .position(|name| name == "q-value")
        .ok_or_else(|| format!("{} has no q-value column", path.display()))?;
    let mut count = 0;
    for (line_number, row) in rows.enumerate() {
        let q = row
            .split('\t')
            .nth(q_column)
            .ok_or_else(|| format!("{}:{} has no q-value", path.display(), line_number + 2))?;
        let q: f64 = q.parse().map_err(|_| {
            format!(
                "{}:{} has invalid q-value {q:?}",
                path.display(),
                line_number + 2
            )
        })?;
        if q < 0.01 {
            count += 1;
        }
    }
    Ok(count)
}

fn write_summary(
    root: &Path,
    implementation: &str,
    wall_seconds: f64,
    results: &[FileResult],
) -> Result<(), String> {
    let matching: Vec<_> = results
        .iter()
        .filter(|result| result.implementation == implementation)
        .collect();
    let successful: Vec<_> = matching
        .iter()
        .copied()
        .filter(|result| result.failure.is_none())
        .collect();
    let failed = matching.len() - successful.len();
    let peak = matching
        .iter()
        .filter_map(|result| result.peak_rss_kb)
        .max()
        .map(|value| value.to_string())
        .unwrap_or_else(|| "NA".to_owned());
    let total = |counts: Vec<Option<u64>>| {
        if successful.len() == matching.len() && counts.iter().any(Option::is_some) {
            counts.into_iter().flatten().sum::<u64>().to_string()
        } else {
            "NA".to_owned()
        }
    };
    let path = root.join(format!("{implementation}-summary.tsv"));
    let mut out = BufWriter::new(File::create(path).map_err(io_error)?);
    writeln!(out, "implementation\twall_seconds\tpeak_rss_kb\toverall_exit_status\tfiles_attempted\tfiles_succeeded\tfiles_failed\tpsms_q_lt_0.01\tpeptides_q_lt_0.01\tproteins_q_lt_0.01") .map_err(io_error)?;
    writeln!(
        out,
        "{implementation}\t{wall_seconds:.6}\t{peak}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        if failed == 0 { "0" } else { "nonzero" },
        matching.len(),
        successful.len(),
        failed,
        total(successful.iter().map(|result| result.psms).collect()),
        total(successful.iter().map(|result| result.peptides).collect()),
        total(successful.iter().map(|result| result.proteins).collect())
    )
    .map_err(io_error)
}

fn write_per_file(root: &Path, results: &[FileResult]) -> Result<(), String> {
    let mut out = BufWriter::new(File::create(root.join("per-file.tsv")).map_err(io_error)?);
    writeln!(out, "implementation\tinput\toutput_dir\texit_status\twall_seconds\tpeak_rss_kb\tpsms_q_lt_0.01\tpeptides_q_lt_0.01\tproteins_q_lt_0.01\tfailure").map_err(io_error)?;
    for result in results {
        writeln!(
            out,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            result.implementation,
            result.input.display(),
            result.output_dir.display(),
            result.exit_status,
            optional_display(result.wall_seconds),
            optional_display(result.peak_rss_kb),
            optional_display(result.psms),
            optional_display(result.peptides),
            optional_display(result.proteins),
            result.failure.as_deref().unwrap_or("")
        )
        .map_err(io_error)?;
    }
    Ok(())
}

fn write_failures(root: &Path, results: &[FileResult]) -> Result<(), String> {
    let mut out = BufWriter::new(File::create(root.join("failures.tsv")).map_err(io_error)?);
    writeln!(
        out,
        "implementation\tinput\texit_status\tstderr_log\tfailure"
    )
    .map_err(io_error)?;
    for result in results.iter().filter(|result| result.failure.is_some()) {
        writeln!(
            out,
            "{}\t{}\t{}\t{}\t{}",
            result.implementation,
            result.input.display(),
            result.exit_status,
            result.output_dir.join("stderr.log").display(),
            result.failure.as_deref().unwrap()
        )
        .map_err(io_error)?;
    }
    Ok(())
}

fn optional_display<T: std::fmt::Display>(value: Option<T>) -> String {
    value
        .map(|value| value.to_string())
        .unwrap_or_else(|| "NA".to_owned())
}

fn status_label(status: &ExitStatus) -> String {
    status
        .code()
        .map(|code| code.to_string())
        .unwrap_or_else(|| "signal".to_owned())
}

fn display_command(command: &[String], time_path: &Path) -> String {
    let mut shown = vec![
        "/usr/bin/time".to_owned(),
        "-f".to_owned(),
        shell_quote("%e\\t%M"),
        "-o".to_owned(),
        shell_quote(&time_path.display().to_string()),
    ];
    shown.extend(command.iter().map(|argument| shell_quote(argument)));
    shown.join(" ")
}

fn shell_quote(value: &str) -> String {
    if value
        .bytes()
        .all(|byte| byte.is_ascii_alphanumeric() || b"_+-./=".contains(&byte))
    {
        value.to_owned()
    } else {
        format!("'{}'", value.replace('\'', "'\\\"'\\\"'"))
    }
}

fn io_error(error: io::Error) -> String {
    error.to_string()
}
