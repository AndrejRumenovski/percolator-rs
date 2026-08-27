//! Fresh hand-graph and competition probe for the final post-repair audit.
//!
//! It uses graph shapes not present in the repair fixtures, reverses every
//! insertion order, and includes a delimiter-collision attack on the picked
//! competition key.

#[path = "../src/protein.rs"]
mod protein;
#[path = "../src/stats.rs"]
mod stats;
#[path = "../src/tiebreak.rs"]
mod tiebreak;

use protein::{infer, ProtGroup};

fn entry(score: f64, proteins: &str) -> (f64, f64, String) {
    (score, 0.123, proteins.to_string())
}

#[derive(Debug, Eq, Ord, PartialEq, PartialOrd)]
struct Signature {
    proteins: Vec<String>,
    decoy: bool,
    score: u64,
    qvalue: u64,
    pep: Option<u64>,
    peptides: usize,
    picked: bool,
}

fn signature(groups: &[ProtGroup]) -> Vec<Signature> {
    let mut output: Vec<_> = groups
        .iter()
        .map(|group| Signature {
            proteins: group.proteins.clone(),
            decoy: group.is_decoy,
            score: group.score.to_bits(),
            qvalue: group.qval.to_bits(),
            pep: group.pep.map(f64::to_bits),
            peptides: group.n_peptides,
            picked: group.picked,
        })
        .collect();
    output.sort();
    output
}

fn members(groups: &[ProtGroup]) -> Vec<Vec<String>> {
    let mut output: Vec<_> = groups.iter().map(|group| group.proteins.clone()).collect();
    output.sort();
    output
}

fn shuffled<T: Clone>(values: &[T], mut state: u64) -> Vec<T> {
    let mut output = values.to_vec();
    for index in (1..output.len()).rev() {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        output.swap(index, (state % (index as u64 + 1)) as usize);
    }
    output
}

fn main() {
    let graph = vec![
        // IX/IY are indistinguishable: exactly the same two peptides.
        entry(20.0, "IX IY"),
        entry(19.0, "IY IX"),
        // A/B/C share a connected component but have distinct evidence sets.
        entry(18.0, "A B"),
        entry(17.0, "B C"),
        entry(16.0, "A"),
        // SUB is a strict subset of SUPER.
        entry(15.0, "SUB SUPER"),
        entry(14.0, "SUPER"),
        // One observed peptide maps to paired target and decoy proteins.
        entry(13.0, "PAIR DECOY_PAIR"),
        // An independent exact-score target/decoy pair.
        entry(12.0, "TIE"),
        entry(12.0, "DECOY_TIE"),
        entry(11.0, "ALONE"),
    ];
    let expected = vec![
        vec!["A"],
        vec!["ALONE"],
        vec!["B"],
        vec!["C"],
        vec!["DECOY_PAIR"],
        vec!["DECOY_TIE"],
        vec!["IX", "IY"],
        vec!["PAIR"],
        vec!["SUB"],
        vec!["SUPER"],
        vec!["TIE"],
    ]
    .into_iter()
    .map(|group| group.into_iter().map(str::to_string).collect::<Vec<_>>())
    .collect::<Vec<_>>();
    let reference = infer(&graph, 37);
    assert_eq!(members(&reference), expected);
    assert!(reference.iter().all(|group| group.pep.is_none()));
    for arranged in [
        graph.iter().cloned().rev().collect::<Vec<_>>(),
        shuffled(&graph, 0x51A7),
        shuffled(&graph, 0xA991),
    ] {
        assert_eq!(members(&infer(&arranged, 37)), expected);
        assert_eq!(signature(&infer(&arranged, 37)), signature(&reference));
    }

    // A new 509-pair exact-tie population checks fair target/decoy picking and
    // reversed insertion order independently of the repair's populations.
    let mut ties = Vec::new();
    for pair in 0..509usize {
        ties.push(entry(7.0, &format!("P{pair:04}")));
        ties.push(entry(7.0, &format!("DECOY_P{pair:04}")));
    }
    let tied = infer(&ties, 37);
    let target_wins = tied
        .iter()
        .filter(|group| group.picked && !group.is_decoy)
        .count();
    assert_eq!(tied.iter().filter(|group| group.picked).count(), 509);
    assert!((target_wins as isize - 254).abs() < 70);
    ties.reverse();
    assert_eq!(signature(&infer(&ties, 37)), signature(&tied));
    assert_ne!(signature(&infer(&ties, 38)), signature(&tied));
    let target_wins_100_seeds: usize = (1..=100u64)
        .map(|seed| {
            let groups = infer(&ties, seed);
            assert_eq!(groups.iter().filter(|group| group.picked).count(), 509);
            groups
                .iter()
                .filter(|group| group.picked && !group.is_decoy)
                .count()
        })
        .sum();
    // Binomial(50,900, 0.5): sigma is about 112.8.  This six-sigma gate is
    // intentionally diagnostic rather than tuned to one observed draw.
    assert!((target_wins_100_seeds as isize - 25_450).abs() < 677);

    // IMPLEMENTATION DEFECT: the pairing key serializes a protein set with an
    // unescaped '|'.  The distinct target groups {LEFT, RIGHT} and
    // {LEFT|RIGHT} therefore occupy the same target slot; one is silently
    // discarded before target/decoy competition.
    let collision_entries = vec![
        entry(10.0, "LEFT RIGHT"),
        entry(9.0, "LEFT|RIGHT"),
        entry(8.0, "DECOY_LEFT DECOY_RIGHT"),
    ];
    let collision = infer(&collision_entries, 37);
    assert_eq!(collision.len(), 3);
    let target_groups = collision
        .iter()
        .filter(|group| !group.is_decoy)
        .collect::<Vec<_>>();
    assert_eq!(target_groups.len(), 2);
    assert_eq!(
        target_groups.iter().filter(|group| group.picked).count(),
        1,
        "the collision should reproduce one silently dropped target group"
    );
    let lost = target_groups
        .iter()
        .find(|group| !group.picked)
        .expect("one colliding target group must be lost");

    println!(
        "PASS graph_groups={} permutations=4 tied_target_wins={target_wins}/509 tied_target_wins_100_seeds={target_wins_100_seeds}/50900 picked_pep=NA",
        reference.len()
    );
    println!(
        "FAILURE PICKING_KEY_COLLISION lost_group={} colliding_key=LEFT|RIGHT",
        lost.proteins.join(",")
    );
}
