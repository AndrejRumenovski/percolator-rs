//! Independent structural probes for picked-protein inference.
//!
//! This is an audit harness, not a Cargo target or production change.  Build it
//! standalone:
//!
//! ```text
//! rustc --edition=2021 validation/adversarial_protein.rs -o /tmp/adversarial-protein
//! ```

#[path = "../src/protein.rs"]
mod protein;
#[path = "../src/stats.rs"]
mod stats;
#[path = "../src/tiebreak.rs"]
mod tiebreak;

fn entry(score: f64, pep: f64, proteins: &str) -> (f64, f64, String) {
    (score, pep, proteins.to_string())
}

fn members(groups: &[protein::ProtGroup]) -> Vec<Vec<String>> {
    let mut sets: Vec<Vec<String>> = groups.iter().map(|g| g.proteins.clone()).collect();
    sets.sort();
    sets
}

fn main() {
    // A and B are distinguishable: A has a unique observed peptide, B does not.
    // Connected-component grouping collapses them into one group; evidence-set
    // grouping keeps them apart.
    let distinguishable = vec![entry(10.0, 0.01, "A B"), entry(9.0, 0.02, "A")];
    let grouped = protein::infer(&distinguishable, 1);
    println!(
        "DISTINGUISHABLE_GROUPING groups={} members={:?}",
        grouped.len(),
        members(&grouped)
    );

    // A subset protein is distinguishable from its superset.
    let subset = vec![
        entry(10.0, 0.01, "SUB SUPER"),
        entry(9.0, 0.02, "SUB SUPER"),
        entry(8.0, 0.03, "SUPER"),
    ];
    let grouped = protein::infer(&subset, 1);
    println!(
        "SUBSET_GROUPING groups={} members={:?}",
        grouped.len(),
        members(&grouped)
    );

    // A chain of shared peptides must not collapse a whole component.
    let chain = vec![
        entry(10.0, 0.01, "A B"),
        entry(9.0, 0.02, "B C"),
        entry(8.0, 0.03, "C D"),
    ];
    let grouped = protein::infer(&chain, 1);
    println!(
        "SHARING_CHAIN groups={} members={:?}",
        grouped.len(),
        members(&grouped)
    );

    // Genuinely identical evidence still groups.
    let identical = vec![entry(10.0, 0.01, "A B"), entry(9.0, 0.02, "A B")];
    let grouped = protein::infer(&identical, 1);
    println!(
        "IDENTICAL_EVIDENCE groups={} members={:?}",
        grouped.len(),
        members(&grouped)
    );

    // One exactly tied target/decoy protein pair per bucket, under several
    // input orders and two seeds.
    let build = |reversed: bool, count: usize| -> Vec<(f64, f64, String)> {
        let mut tied = Vec::new();
        for index in 0..count {
            let target = entry(1.0, 0.5, &format!("P{index:03}"));
            let decoy = entry(1.0, 0.5, &format!("DECOY_P{index:03}"));
            if reversed {
                tied.push(decoy);
                tied.push(target);
            } else {
                tied.push(target);
                tied.push(decoy);
            }
        }
        tied
    };
    let report = |label: &str, entries: &[(f64, f64, String)], seed: u64| {
        let inferred = protein::infer(entries, seed);
        let targets: Vec<_> = inferred
            .iter()
            .filter(|group| group.picked && !group.is_decoy)
            .collect();
        let decoys = inferred
            .iter()
            .filter(|group| group.picked && group.is_decoy)
            .count();
        let mut signature: Vec<String> = inferred
            .iter()
            .filter(|group| group.picked)
            .map(|group| group.proteins.join("|"))
            .collect();
        signature.sort();
        println!(
            "PICKED_TIE_ATTACK arm={label} seed={seed} target_winners={} decoy_winners={} targets_q_lt_0.01={} min_target_q={:?} signature_hash={:016x}",
            targets.len(),
            decoys,
            targets.iter().filter(|group| group.qval < 0.01).count(),
            targets.iter().map(|group| group.qval).reduce(f64::min),
            signature
                .iter()
                .fold(0xcbf2_9ce4_8422_2325u64, |acc, name| {
                    name.bytes()
                        .fold(acc, |a, b| (a ^ u64::from(b)).wrapping_mul(0x1000_0000_01b3))
                })
        );
    };
    report("target_first", &build(false, 200), 1);
    report("decoy_first", &build(true, 200), 1);
    report("target_first", &build(false, 200), 2);

    // Picked mode must not emit a protein-level posterior it does not estimate.
    let one = protein::infer(&[entry(8.0, 0.123456, "ONLY_PROTEIN")], 1);
    println!(
        "PICKED_PEP input_peptide_pep=0.123456 output_protein_pep={:?}",
        one[0].pep
    );
}
