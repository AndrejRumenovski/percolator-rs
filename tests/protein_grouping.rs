//! End-to-end protection for peptide-to-protein association and grouping.

use std::collections::BTreeSet;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicU64, Ordering};

static NEXT_DIRECTORY: AtomicU64 = AtomicU64::new(0);

#[derive(Clone)]
struct PinRow {
    id: &'static str,
    peptide: &'static str,
    proteins: &'static str,
    scan: i64,
    score: f64,
}

fn temp_directory(arm: &str) -> PathBuf {
    let serial = NEXT_DIRECTORY.fetch_add(1, Ordering::Relaxed);
    let path = std::env::temp_dir().join(format!(
        "percolator-rs-protein-grouping-{arm}-{}-{serial}",
        std::process::id()
    ));
    std::fs::create_dir_all(&path).expect("temporary directory");
    path
}

fn write_pin(path: &Path, rows: &[PinRow]) {
    let mut text =
        String::from("SpecId\tLabel\tScanNr\tExpMass\tscore\tf1\tf2\tPeptide\tProteins\n");
    for row in rows {
        text.push_str(&format!(
            "{}\t1\t{}\t500\t{}\t0\t1\t{}\t{}\n",
            row.id, row.scan, row.score, row.peptide, row.proteins
        ));
    }
    std::fs::write(path, text).expect("fixture should write");
}

fn groups_for(arm: &str, rows: &[PinRow]) -> BTreeSet<(String, String, String)> {
    let directory = temp_directory(arm);
    let pin = directory.join("fixture.pin");
    let targets = directory.join("targets.tsv");
    let decoys = directory.join("decoys.tsv");
    let proteins = directory.join("proteins.tsv");
    let decoy_proteins = directory.join("decoy-proteins.tsv");
    write_pin(&pin, rows);

    let output = Command::new(env!("CARGO_BIN_EXE_percolator-rs"))
        .arg("--canonical")
        .arg("--no-select-c")
        .args(["--maxiter", "0"])
        .arg("--no-psm-competition")
        .args(["--seed", "1"])
        .args(["--num-threads", "1"])
        .args(["--results-psms", targets.to_str().unwrap()])
        .args(["--decoy-results-psms", decoys.to_str().unwrap()])
        .args(["--results-proteins", proteins.to_str().unwrap()])
        .args(["--decoy-results-proteins", decoy_proteins.to_str().unwrap()])
        .arg(&pin)
        .output()
        .expect("percolator-rs should execute");
    assert!(
        output.status.success(),
        "{arm} failed:\n{}",
        String::from_utf8_lossy(&output.stderr)
    );

    let mut groups = BTreeSet::new();
    for path in [&proteins, &decoy_proteins] {
        let text = std::fs::read_to_string(path).expect("protein output should exist");
        for line in text.lines().skip(1) {
            let fields: Vec<&str> = line.split('\t').collect();
            groups.insert((
                fields[0].to_string(),
                fields[2].to_string(),
                fields[5].to_string(),
            ));
        }
    }
    std::fs::remove_dir_all(directory).ok();
    groups
}

fn shuffle(rows: &mut [PinRow], mut state: u64) {
    for index in (1..rows.len()).rev() {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        rows.swap(index, (state % (index as u64 + 1)) as usize);
    }
}

/// Four-row structural counterexample from the independent audit. One peptide
/// is observed twice with equal score and complementary protein mappings. Its
/// association is the union {PROT_A, PROT_B}; choosing the first representative
/// loses one edge and changes the distinguishable group set.
#[test]
fn repeated_peptide_mappings_are_unioned_in_every_input_order() {
    let rows = vec![
        PinRow {
            id: "AMB_A",
            peptide: "K.AMBIGUOUS.R",
            proteins: "PROT_A",
            scan: 1,
            score: 5.0,
        },
        PinRow {
            id: "AMB_B",
            peptide: "K.AMBIGUOUS.R",
            proteins: "PROT_B",
            scan: 1,
            score: 5.0,
        },
        PinRow {
            id: "UNIQUE_A",
            peptide: "K.UNIQUEA.R",
            proteins: "PROT_A",
            scan: 2,
            score: 4.0,
        },
        PinRow {
            id: "UNIQUE_C",
            peptide: "K.UNIQUEC.R",
            proteins: "PROT_C",
            scan: 3,
            score: 3.0,
        },
    ];
    let expected: BTreeSet<(String, String, String)> = ["PROT_A", "PROT_B", "PROT_C"]
        .into_iter()
        .map(|protein| (protein.to_string(), "NA".to_string(), protein.to_string()))
        .collect();

    let mut reversed = rows.clone();
    reversed.reverse();
    let mut shuffled_17 = rows.clone();
    shuffle(&mut shuffled_17, 17);
    let mut shuffled_91 = rows.clone();
    shuffle(&mut shuffled_91, 91);
    for (arm, arranged) in [
        ("original", rows),
        ("reversed", reversed),
        ("shuffle-17", shuffled_17),
        ("shuffle-91", shuffled_91),
    ] {
        assert_eq!(
            groups_for(arm, &arranged),
            expected,
            "{arm}: repeated peptide mapping changed the protein graph"
        );
    }
}
