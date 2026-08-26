//! Adversarial held-out-label attack on every advertised cross-validation mode.
//!
//! # The property under test
//!
//! A row's out-of-fold score is produced by a model that was trained without
//! that row's fold. Its own label therefore cannot reach it: flip the label of
//! one row, leave every feature byte-identical, and that row's reported score
//! must not move. Everything else in the run may move — the other folds see the
//! flipped row in their training data, and the q-values see a different label
//! composition — so only the attacked row's own score is compared.
//!
//! This formulation needs no knowledge of the fold assignment, which is what
//! makes it an attack rather than a restatement of the implementation.
//!
//! Each mode is tested independently. Passing on one mode says nothing about
//! another: `--select-c` used to fail this while the default passed, because its
//! hyperparameter search sat outside the folds.

use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::process::Command;

const ROWS: usize = 600;
/// Rows attacked one at a time. Fixed in advance, spread across the file so the
/// deterministic fold assignment cannot put them all in one fold.
const ATTACKED: [usize; 6] = [3, 97, 198, 301, 444, 577];

struct Row {
    label: i8,
    features: [f64; 4],
}

/// Deterministic separable-with-noise fixture; no dependence on the crate.
fn fixture() -> Vec<Row> {
    let mut state = 0x2545_F491_4F6C_DD1Du64;
    let mut next = || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        (state >> 11) as f64 / (1u64 << 53) as f64
    };
    (0..ROWS)
        .map(|_| {
            let label = if next() < 0.65 { 1i8 } else { -1 };
            let signal = f64::from(label);
            Row {
                label,
                features: [
                    signal * 1.4 + (next() - 0.5) * 3.0,
                    signal * 0.6 + (next() - 0.5) * 2.0,
                    (next() - 0.5) * 2.0,
                    (next() - 0.5) * 2.0,
                ],
            }
        })
        .collect()
}

/// `engine` selects which half of the rows this file reports, so an ensemble run
/// gets two files that overlap on scan numbers.
fn write_pin(path: &Path, rows: &[Row], flipped: Option<usize>, half: Option<usize>) {
    let mut text =
        String::from("SpecId\tLabel\tScanNr\tExpMass\tf0\tf1\tf2\tf3\tPeptide\tProteins\n");
    for (index, row) in rows.iter().enumerate() {
        if let Some(half) = half {
            // Engine 0 reports the first two thirds, engine 1 the last two
            // thirds, so the middle third is reported by both.
            let keep = if half == 0 {
                index < ROWS * 2 / 3
            } else {
                index >= ROWS / 3
            };
            if !keep {
                continue;
            }
        }
        let label = if Some(index) == flipped {
            -row.label
        } else {
            row.label
        };
        text.push_str(&format!(
            "p{index}\t{label}\t{index}\t500.0\t{}\t{}\t{}\t{}\tK.PEP{index}.R\tPROT{index}\n",
            row.features[0], row.features[1], row.features[2], row.features[3]
        ));
    }
    std::fs::write(path, text).expect("fixture should write");
}

fn read_scores(paths: &[PathBuf]) -> HashMap<String, String> {
    let mut scores = HashMap::new();
    for path in paths {
        let text = std::fs::read_to_string(path).expect("result file should exist");
        for line in text.lines().skip(1) {
            let mut fields = line.split('\t');
            let id = fields.next().unwrap_or_default().to_string();
            let score = fields.next().unwrap_or_default().to_string();
            scores.insert(id, score);
        }
    }
    scores
}

fn run(mode: &str, directory: &Path, tag: &str, inputs: &[PathBuf]) -> HashMap<String, String> {
    let targets = directory.join(format!("{tag}.targets.tsv"));
    let decoys = directory.join(format!("{tag}.decoys.tsv"));
    let mut command = Command::new(env!("CARGO_BIN_EXE_percolator-rs"));
    command
        .arg("--canonical")
        .args(["--seed", "1"])
        .args(["--num-threads", "1"])
        .arg("--no-psm-competition")
        .args(["--results-psms", targets.to_str().unwrap()])
        .args(["--decoy-results-psms", decoys.to_str().unwrap()]);
    match mode {
        "fixed-c" => {
            command.arg("--no-select-c");
            command.arg(inputs[0].to_str().unwrap());
        }
        "select-c" => {
            command.arg("--select-c");
            command.arg(inputs[0].to_str().unwrap());
        }
        "ensemble" => {
            command.arg("--no-select-c").arg("--ensemble");
            command.arg(format!("engineA={}", inputs[0].display()));
            command.arg(format!("engineB={}", inputs[1].display()));
        }
        other => panic!("unknown mode {other}"),
    }
    let output = command.output().expect("percolator-rs should run");
    assert!(
        output.status.success(),
        "{mode}/{tag} failed:\n{}",
        String::from_utf8_lossy(&output.stderr)
    );
    read_scores(&[targets, decoys])
}

fn attack(mode: &str, ensemble: bool) {
    let directory = std::env::temp_dir().join(format!(
        "percolator-rs-leakage-{mode}-{}",
        std::process::id()
    ));
    std::fs::create_dir_all(&directory).expect("temp directory");
    let rows = fixture();

    let inputs: Vec<PathBuf> = if ensemble {
        vec![directory.join("clean.a.pin"), directory.join("clean.b.pin")]
    } else {
        vec![directory.join("clean.pin")]
    };
    for (half, path) in inputs.iter().enumerate() {
        write_pin(path, &rows, None, ensemble.then_some(half));
    }
    let clean = run(mode, &directory, "clean", &inputs);

    for &attacked in &ATTACKED {
        let dirty: Vec<PathBuf> = if ensemble {
            vec![
                directory.join(format!("dirty{attacked}.a.pin")),
                directory.join(format!("dirty{attacked}.b.pin")),
            ]
        } else {
            vec![directory.join(format!("dirty{attacked}.pin"))]
        };
        for (half, path) in dirty.iter().enumerate() {
            write_pin(path, &rows, Some(attacked), ensemble.then_some(half));
        }
        let scores = run(mode, &directory, &format!("dirty{attacked}"), &dirty);

        let mut compared = 0;
        for key in clean.keys() {
            // The attacked row's own identifier, under either engine prefix.
            let own = key == &format!("p{attacked}")
                || key.ends_with(&format!(":p{attacked}"))
                    && key.rsplit(':').next() == Some(&format!("p{attacked}"));
            if !own {
                continue;
            }
            compared += 1;
            assert_eq!(
                scores.get(key),
                clean.get(key),
                "{mode}: flipping the label of row {attacked} changed its own held-out score \
                 ({key}); the model that scored it saw its label"
            );
        }
        assert!(
            compared > 0,
            "{mode}: attacked row {attacked} was not present in the output"
        );
    }
    std::fs::remove_dir_all(&directory).ok();
}

#[test]
fn fixed_c_cross_validation_is_leakage_free() {
    attack("fixed-c", false);
}

/// **Frozen failure case (independent audit M5).** `--select-c` chose its class
/// weights on the same out-of-fold predictions it then reported.
#[test]
fn selected_c_cross_validation_is_leakage_free() {
    attack("select-c", false);
}

/// **Frozen failure case (independent audit M5).** `--ensemble` built a
/// whole-dataset agreement feature keyed on `(ScanNr, Label, Peptide)` before
/// the folds existed.
#[test]
fn ensemble_cross_validation_is_leakage_free() {
    attack("ensemble", true);
}
