//! Standalone adversarial probe for `src/stats.rs`.
//!
//! Compile with:
//! `rustc --edition=2021 validation/independent_stats_probe.rs -o /tmp/independent-stats-probe`
//!
//! This is deliberately outside the Cargo targets.  Its oracle is written from
//! the mathematical definition rather than by calling another production helper.

#[path = "../src/stats.rs"]
mod production;

use production::{qvalues, qvalues_and_peps, Tdc};

fn oracle(scores: &[f64], labels: &[i8], null_target_win_prob: f64) -> Vec<f64> {
    assert_eq!(scores.len(), labels.len());
    assert!((0.0..1.0).contains(&null_target_win_prob));
    let factor = null_target_win_prob / (1.0 - null_target_win_prob);
    let mut order: Vec<usize> = (0..scores.len()).collect();
    order.sort_by(
        |&left, &right| match (scores[left].is_nan(), scores[right].is_nan()) {
            (true, true) => std::cmp::Ordering::Equal,
            (true, false) => std::cmp::Ordering::Greater,
            (false, true) => std::cmp::Ordering::Less,
            (false, false) => scores[right].total_cmp(&scores[left]),
        },
    );

    let mut fdp = vec![1.0; scores.len()];
    let mut targets = 0usize;
    let mut decoys = 1usize;
    let mut group_start = 0usize;
    for rank in 0..order.len() {
        if labels[order[rank]] > 0 {
            targets += 1;
        } else {
            decoys += 1;
        }
        let group_ends = rank + 1 == order.len() || scores[order[rank + 1]] != scores[order[rank]];
        if group_ends {
            let value = ((decoys as f64) * factor / targets.max(1) as f64).min(1.0);
            fdp[group_start..=rank].fill(value);
            group_start = rank + 1;
        }
    }

    let mut result = vec![1.0; scores.len()];
    let mut minimum = 1.0_f64;
    for rank in (0..order.len()).rev() {
        minimum = minimum.min(fdp[rank]);
        result[order[rank]] = minimum;
    }
    result
}

fn assert_close(name: &str, actual: &[f64], expected: &[f64]) {
    assert_eq!(actual.len(), expected.len(), "{name}: length");
    for (index, (&actual, &expected)) in actual.iter().zip(expected).enumerate() {
        assert!(
            (actual - expected).abs() <= 1e-12,
            "{name}[{index}]: actual={actual:?}, expected={expected:?}"
        );
    }
}

fn main() {
    let cases: Vec<(&str, Vec<f64>, Vec<i8>, Vec<f64>)> = vec![
        ("empty", vec![], vec![], vec![]),
        ("one_target", vec![1.0], vec![1], vec![1.0]),
        ("one_decoy", vec![1.0], vec![-1], vec![1.0]),
        (
            "target_then_decoy",
            vec![2.0, 1.0],
            vec![1, -1],
            vec![1.0, 1.0],
        ),
        (
            "decoy_then_target",
            vec![2.0, 1.0],
            vec![-1, 1],
            vec![1.0, 1.0],
        ),
        (
            "all_targets",
            vec![4.0, 3.0, 2.0, 1.0],
            vec![1, 1, 1, 1],
            vec![0.25; 4],
        ),
        (
            "all_decoys",
            vec![4.0, 3.0, 2.0, 1.0],
            vec![-1, -1, -1, -1],
            vec![1.0; 4],
        ),
        (
            "alternating",
            vec![6.0, 5.0, 4.0, 3.0, 2.0, 1.0],
            vec![1, 1, -1, 1, -1, 1],
            vec![0.5, 0.5, 2.0 / 3.0, 2.0 / 3.0, 0.75, 0.75],
        ),
        (
            "mixed_ties",
            vec![2.0, 2.0, 1.0, 1.0],
            vec![1, -1, 1, 1],
            vec![2.0 / 3.0; 4],
        ),
        (
            "all_identical",
            vec![1.0; 6],
            vec![1, -1, 1, -1, 1, -1],
            vec![1.0; 6],
        ),
        (
            "extreme_finite",
            vec![f64::MAX, f64::MIN, 0.0, -0.0],
            vec![1, -1, 1, -1],
            vec![1.0; 4],
        ),
    ];

    for (name, scores, labels, hand_expected) in cases {
        let actual = qvalues(&scores, &labels, Tdc::reported(0.5));
        assert_close(name, &actual, &hand_expected);
        assert_close(name, &actual, &oracle(&scores, &labels, 0.5));
        println!("PASS {name}: {actual:?}");
    }

    // Highly imbalanced opportunity ratios: two decoys per target imply p=1/3
    // and therefore half the balanced-search decoy factor.
    let scores = vec![5.0, 4.0, 3.0, 2.0, 1.0, 0.0];
    let labels = vec![1, 1, 1, 1, 1, -1];
    let unbalanced = qvalues(&scores, &labels, Tdc::reported(1.0 / 3.0));
    assert_close(
        "opportunity_ratio",
        &unbalanced,
        &oracle(&scores, &labels, 1.0 / 3.0),
    );
    assert_close(
        "opportunity_ratio_hand",
        &unbalanced,
        &[0.1, 0.1, 0.1, 0.1, 0.1, 0.2],
    );
    println!("PASS opportunity_ratio: {unbalanced:?}");

    // Production rejects non-finite PIN features before they can become scores.
    // The estimator boundary used to accept them anyway, and a trailing NaN
    // target then changed every earlier q-value through the reverse cumulative
    // minimum. It must now refuse the input instead.
    let finite = qvalues(&[3.0, 2.0, 1.0], &[1, -1, 1], Tdc::reported(0.5));
    let refused = std::panic::catch_unwind(|| {
        qvalues(
            &[3.0, 2.0, 1.0, f64::NAN],
            &[1, -1, 1, 1],
            Tdc::reported(0.5),
        )
    });
    println!(
        "INTERNAL_NAN finite={finite:?} non_finite_refused={}",
        refused.is_err()
    );

    // The PEP identity, checked against counts enumerated here rather than
    // against another production helper: the reported PEPs of the top k targets
    // sum to the estimated number of false discoveries among them.
    let (q, pep) = qvalues_and_peps(&scores, &labels, Tdc::reported(0.5));
    println!("PEP_DIAGNOSTIC scores={scores:?} labels={labels:?} q={q:?} pep={pep:?}");
    for (name, scores, labels, expected_total) in [
        ("five_targets_one_decoy", vec![5.0, 4.0, 3.0, 2.0, 1.0, 0.0], vec![1i8, 1, 1, 1, 1, -1], 1.0),
        (
            "ten_one_ten",
            (0..21).map(|k| 100.0 - k as f64).collect::<Vec<f64>>(),
            (0..21).map(|k| if k == 10 { -1i8 } else { 1 }).collect::<Vec<i8>>(),
            2.0,
        ),
    ] {
        let (_, pep) = qvalues_and_peps(&scores, &labels, Tdc::reported(0.5));
        let total: f64 = labels
            .iter()
            .zip(&pep)
            .filter(|(&label, _)| label > 0)
            .map(|(_, &value)| value)
            .sum();
        let zeros = pep.iter().filter(|&&value| value == 0.0).count();
        println!(
            "PEP_IDENTITY {name} summed_target_pep={total} expected={expected_total} match={} exact_zeros={zeros}",
            (total - expected_total).abs() < 1e-9
        );
    }
}
