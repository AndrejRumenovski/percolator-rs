//! End-to-end regression attacks for joined-input permutation invariance.

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicU64, Ordering};

static NEXT_DIRECTORY: AtomicU64 = AtomicU64::new(0);

#[derive(Clone)]
struct PinRow {
    id: String,
    label: i8,
    scan: i64,
    mass: f64,
    features: [f64; 3],
    peptide: String,
    proteins: String,
}

fn row(id: &str, label: i8, scan: i64, f0: f64) -> PinRow {
    PinRow {
        id: id.to_string(),
        label,
        scan,
        mass: 500.0,
        features: [f0, (scan % 7) as f64, (scan % 11) as f64],
        peptide: format!("K.{id}.R"),
        proteins: format!("{}P_{id}", if label < 0 { "DECOY_" } else { "" }),
    }
}

fn temp_directory(test: &str, arm: &str) -> PathBuf {
    let serial = NEXT_DIRECTORY.fetch_add(1, Ordering::Relaxed);
    let path = std::env::temp_dir().join(format!(
        "percolator-rs-join-{test}-{arm}-{}-{serial}",
        std::process::id()
    ));
    std::fs::create_dir_all(&path).expect("temporary test directory");
    path
}

fn write_pin(path: &Path, rows: &[PinRow]) {
    let mut text = String::from("SpecId\tLabel\tScanNr\tExpMass\tf0\tf1\tf2\tPeptide\tProteins\n");
    for row in rows {
        text.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            row.id,
            row.label,
            row.scan,
            row.mass,
            row.features[0],
            row.features[1],
            row.features[2],
            row.peptide,
            row.proteins,
        ));
    }
    std::fs::write(path, text).expect("PIN fixture should write");
}

#[derive(Debug, PartialEq, Eq)]
struct Run {
    rows: BTreeMap<String, (i8, String, String, String)>, // label, score, q, PEP
    discoveries_q01: usize,
}

fn read_results(path: &Path, label: i8, output: &mut Run) {
    let text = std::fs::read_to_string(path).expect("result should exist");
    for line in text.lines().skip(1) {
        let fields: Vec<&str> = line.split('\t').collect();
        assert!(fields.len() >= 4, "short result row: {line}");
        if label > 0 && fields[2].parse::<f64>().unwrap() < 0.01 {
            output.discoveries_q01 += 1;
        }
        output.rows.insert(
            fields[0].to_string(),
            (
                label,
                fields[1].to_string(),
                fields[2].to_string(),
                fields[3].to_string(),
            ),
        );
    }
}

fn execute(
    test: &str,
    arm: &str,
    files: &[(&str, Vec<PinRow>)],
    argument_order: &[usize],
    maxiter: usize,
) -> Run {
    let directory = temp_directory(test, arm);
    let paths: Vec<PathBuf> = files
        .iter()
        .map(|(name, rows)| {
            let path = directory.join(name);
            write_pin(&path, rows);
            path
        })
        .collect();
    let targets = directory.join("targets.tsv");
    let decoys = directory.join("decoys.tsv");
    let mut command = Command::new(env!("CARGO_BIN_EXE_percolator-rs"));
    command
        .arg("--canonical")
        .arg("--no-select-c")
        .args(["--maxiter", &maxiter.to_string()])
        .arg("--join")
        .args(["--seed", "17"])
        .args(["--num-threads", "1"])
        .args(["--results-psms", targets.to_str().unwrap()])
        .args(["--decoy-results-psms", decoys.to_str().unwrap()]);
    for &index in argument_order {
        command.arg(&paths[index]);
    }
    let result = command.output().expect("percolator-rs should execute");
    assert!(
        result.status.success(),
        "{test}/{arm} failed:\n{}",
        String::from_utf8_lossy(&result.stderr)
    );
    let mut output = Run {
        rows: BTreeMap::new(),
        discoveries_q01: 0,
    };
    read_results(&targets, 1, &mut output);
    read_results(&decoys, -1, &mut output);
    std::fs::remove_dir_all(directory).ok();
    output
}

fn shuffled(mut rows: Vec<PinRow>, mut state: u64) -> Vec<PinRow> {
    for index in (1..rows.len()).rev() {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        rows.swap(index, (state % (index as u64 + 1)) as usize);
    }
    rows
}

/// The three-row minimum: two exactly tied candidates in one source and one
/// lower candidate in a second source. Current source-position keying flips T/D
/// when the two file arguments are reversed.
#[test]
fn minimal_joined_exact_tie_is_file_order_invariant() {
    let files = vec![
        (
            "alpha.pin",
            vec![row("T", 1, 10, 0.0), row("D", -1, 10, 0.0)],
        ),
        ("beta.pin", vec![row("LOW", -1, 20, -10.0)]),
    ];
    let original = execute("minimal", "alpha-beta", &files, &[0, 1], 0);
    let reversed = execute("minimal", "beta-alpha", &files, &[1, 0], 0);
    assert_eq!(
        reversed, original,
        "reversing joined files flipped an exact-tie winner"
    );
}

/// Independent-audit cutoff fixture. The test deliberately does not prescribe
/// whether the seeded exact tie is won by target or decoy; it prescribes that a
/// permutation cannot change that draw or the resulting rejection boundary.
#[test]
fn joined_cutoff_is_invariant_to_file_and_row_layout() {
    let mut alpha: Vec<PinRow> = (1..=100)
        .map(|scan| row(&format!("HIGH_T_{scan}"), 1, scan, 10.0))
        .collect();
    alpha.push(row("BOUNDARY_T", 1, 10_002, 0.0));
    alpha.push(row("BOUNDARY_D", -1, 10_002, 0.0));
    let beta: Vec<PinRow> = (20_001..=20_100)
        .map(|scan| row(&format!("LOW_D_{scan}"), -1, scan, -10.0))
        .collect();
    let files = vec![("alpha.pin", alpha.clone()), ("beta.pin", beta.clone())];
    let reference = execute("cutoff", "original", &files, &[0, 1], 0);

    let arms = [
        ("reversed-files", files.clone(), vec![1, 0]),
        (
            "reversed-rows",
            vec![
                ("alpha.pin", alpha.iter().cloned().rev().collect()),
                ("beta.pin", beta.iter().cloned().rev().collect()),
            ],
            vec![0, 1],
        ),
        (
            "target-first",
            vec![
                (
                    "alpha.pin",
                    alpha
                        .iter()
                        .filter(|row| row.label > 0)
                        .cloned()
                        .chain(alpha.iter().filter(|row| row.label < 0).cloned())
                        .collect(),
                ),
                ("beta.pin", beta.clone()),
            ],
            vec![0, 1],
        ),
        (
            "decoy-first",
            vec![
                (
                    "alpha.pin",
                    alpha
                        .iter()
                        .filter(|row| row.label < 0)
                        .cloned()
                        .chain(alpha.iter().filter(|row| row.label > 0).cloned())
                        .collect(),
                ),
                ("beta.pin", beta.clone()),
            ],
            vec![0, 1],
        ),
        (
            "shuffle-991",
            vec![
                ("alpha.pin", shuffled(alpha.clone(), 991)),
                ("beta.pin", shuffled(beta.clone(), 992)),
            ],
            vec![1, 0],
        ),
        (
            "shuffle-1777",
            vec![
                ("alpha.pin", shuffled(alpha, 1777)),
                ("beta.pin", shuffled(beta, 1778)),
            ],
            vec![0, 1],
        ),
    ];
    for (name, arm_files, order) in arms {
        let actual = execute("cutoff", name, &arm_files, &order, 0);
        assert_eq!(
            actual, reference,
            "{name} changed scores, q-values, PEPs, winners, or q<0.01 discoveries"
        );
    }
}

/// Three joined sources exercise iterative training plus an exact score tie and
/// a one-representable-step near tie. Near ties stay strict; no tolerance or
/// epsilon is permitted to turn them into exact ties.
#[test]
fn trained_three_file_join_preserves_exact_and_near_ties_under_permutation() {
    let mut alpha = Vec::new();
    let mut beta = Vec::new();
    let mut gamma = Vec::new();
    for scan in 1..=90i64 {
        let label = if scan % 3 == 0 { -1 } else { 1 };
        let value = label as f64 * 1.5 + (scan % 13) as f64 / 10.0;
        let destination = match scan % 3 {
            0 => &mut alpha,
            1 => &mut beta,
            _ => &mut gamma,
        };
        destination.push(row(&format!("BG_{scan}"), label, scan, value));
    }
    alpha.extend([
        row("EXACT_T", 1, 1_001, 3.0),
        row("EXACT_D", -1, 1_001, 3.0),
    ]);
    let near = 3.0f64;
    beta.extend([
        row("NEAR_T", 1, 1_002, near),
        row("NEAR_D", -1, 1_002, f64::from_bits(near.to_bits() + 1)),
    ]);
    let files = vec![
        ("alpha.pin", alpha.clone()),
        ("beta.pin", beta.clone()),
        ("gamma.pin", gamma.clone()),
    ];
    let reference = execute("trained", "original", &files, &[0, 1, 2], 2);
    let reversed = execute(
        "trained",
        "reversed",
        &[
            ("alpha.pin", alpha.into_iter().rev().collect()),
            ("beta.pin", beta.into_iter().rev().collect()),
            ("gamma.pin", gamma.into_iter().rev().collect()),
        ],
        &[2, 1, 0],
        2,
    );
    assert_eq!(
        reversed, reference,
        "three-file permutation changed the trained/scored statistical result"
    );
    assert_eq!(
        usize::from(reference.rows.contains_key("NEAR_T"))
            + usize::from(reference.rows.contains_key("NEAR_D")),
        1,
        "near-tie precursor did not produce exactly one competed winner"
    );
}
