//! Independent hand-graph probe of the repaired picked-protein implementation.
//!
//! This covers graph shapes different from the repair's unit fixtures and also
//! tests a target/decoy mixed-evidence case that the public PIN contract does not
//! reject.

#[path = "../src/protein.rs"]
mod protein;
#[path = "../src/stats.rs"]
mod stats;
#[path = "../src/tiebreak.rs"]
mod tiebreak;

use protein::{infer, ProtGroup};

fn entry(score: f64, pep: f64, proteins: &str) -> (f64, f64, String) {
    (score, pep, proteins.to_string())
}

fn groups(entries: &[(f64, f64, String)], seed: u64) -> Vec<Vec<String>> {
    let mut output: Vec<Vec<String>> = infer(entries, seed)
        .into_iter()
        .map(|group| group.proteins)
        .collect();
    output.sort();
    output
}

fn picked_signature(entries: &[(f64, f64, String)], seed: u64) -> Vec<(Vec<String>, bool, String)> {
    let mut output: Vec<_> = infer(entries, seed)
        .into_iter()
        .filter(|group| group.picked)
        .map(|group| {
            (
                group.proteins,
                group.is_decoy,
                format!("{:.12}", group.qval),
            )
        })
        .collect();
    output.sort();
    output
}

fn find<'a>(groups: &'a [ProtGroup], members: &[&str]) -> &'a ProtGroup {
    groups
        .iter()
        .find(|group| {
            group
                .proteins
                .iter()
                .map(String::as_str)
                .eq(members.iter().copied())
        })
        .unwrap_or_else(|| panic!("missing group {members:?}"))
}

fn main() {
    // Two proteins supported by the same two peptides are truly indistinguishable.
    let indistinguishable = vec![
        entry(12.0, 0.01, "IOTA KAPPA"),
        entry(9.0, 0.03, "KAPPA IOTA"),
    ];
    assert_eq!(
        groups(&indistinguishable, 31),
        vec![vec![String::from("IOTA"), String::from("KAPPA")]]
    );

    // A six-node sharing cycle with alternating unique evidence must remain six
    // groups. Connected components would collapse all six.
    let distinguishable = vec![
        entry(20.0, 0.01, "A B"),
        entry(19.0, 0.02, "B C"),
        entry(18.0, 0.03, "C D"),
        entry(17.0, 0.04, "D E"),
        entry(16.0, 0.05, "E F"),
        entry(15.0, 0.06, "F A"),
        entry(14.0, 0.07, "A"),
        entry(13.0, 0.08, "C"),
        entry(12.0, 0.09, "E"),
    ];
    let expected: Vec<Vec<String>> = ["A", "B", "C", "D", "E", "F"]
        .into_iter()
        .map(|name| vec![name.to_string()])
        .collect();
    assert_eq!(groups(&distinguishable, 31), expected);

    // Nested evidence sets are distinguishable and both shared peptides support
    // every mapped group (no hidden razor-peptide assignment).
    let subset = vec![
        entry(10.0, 0.1, "SUPER SUB"),
        entry(8.0, 0.2, "SUPER"),
        entry(7.0, 0.3, "SUPER OTHER"),
        entry(6.0, 0.4, "INDEPENDENT"),
    ];
    let inferred = infer(&subset, 31);
    assert_eq!(groups(&subset, 31).len(), 4);
    assert_eq!(find(&inferred, &["SUPER"]).n_peptides, 3);
    assert_eq!(find(&inferred, &["SUB"]).n_peptides, 1);
    assert_eq!(find(&inferred, &["OTHER"]).n_peptides, 1);
    assert_eq!(find(&inferred, &["INDEPENDENT"]).n_peptides, 1);
    assert!(inferred.iter().all(|group| group.pep.is_none()));

    // 317 target/decoy protein ties: entry reversal must preserve every pick,
    // while a different run seed must reflip at least one pair.
    let mut tied = Vec::new();
    for index in 0..317 {
        tied.push(entry(5.0, 0.123, &format!("PAIR_{index}")));
        tied.push(entry(5.0, 0.987, &format!("DECOY_PAIR_{index}")));
    }
    let signature = picked_signature(&tied, 31);
    let mut reversed = tied.clone();
    reversed.reverse();
    assert_eq!(picked_signature(&reversed, 31), signature);
    assert_ne!(picked_signature(&tied, 32), signature);
    let tied_groups = infer(&tied, 31);
    let picked_targets = tied_groups
        .iter()
        .filter(|group| group.picked && !group.is_decoy)
        .count();
    assert!(
        (picked_targets as isize - 158).abs() < 50,
        "unfair tied split: {picked_targets}/317"
    );
    assert!(tied_groups.iter().all(|group| group.pep.is_none()));

    // A peptide mapping to both a target and a decoy protein currently collapses
    // them into one group and marks the whole group decoy (`any(decoy)`).  This
    // is an explicit reproduced defect: the target disappears from target output
    // and contaminates a decoy group, even though target/decoy provenance is the
    // quantity protein competition is supposed to preserve.
    let mixed = infer(&[entry(11.0, 0.01, "MIXED DECOY_MIXED")], 31);
    assert_eq!(mixed.len(), 1);
    assert_eq!(mixed[0].proteins, vec!["DECOY_MIXED", "MIXED"]);
    assert!(mixed[0].is_decoy && mixed[0].picked);
    println!(
        "MIXED_TARGET_DECOY_COLLAPSE members={} classified_decoy={}",
        mixed[0].proteins.join(","),
        mixed[0].is_decoy
    );

    println!(
        "PASS: ordinary grouping and fair protein ties; tied target winners={picked_targets}/317; picked PEP is unavailable"
    );
}
