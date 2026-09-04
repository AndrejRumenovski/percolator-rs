//! Independent, hand-calculable post-repair probe for q-values and PEPs.
//!
//! This deliberately does not call any test helper from `src/stats.rs`.  The
//! q-value oracle below is a separate implementation of the documented TDC+
//! formula.  PEP cases use closed-form expectations, plus checks of the public
//! estimator's conservation and boundary claims.

#[path = "../src/stats.rs"]
mod stats;

use stats::{qvalues, qvalues_and_peps, Tdc};

const EPS: f64 = 1e-11;

fn approx(left: f64, right: f64) -> bool {
    (left - right).abs() <= EPS * (1.0 + left.abs().max(right.abs()))
}

fn assert_vec(actual: &[f64], expected: &[f64], name: &str) {
    assert_eq!(actual.len(), expected.len(), "{name}: length");
    for (index, (&left, &right)) in actual.iter().zip(expected).enumerate() {
        assert!(
            approx(left, right),
            "{name}[{index}]: actual={left:.17} expected={right:.17}"
        );
    }
}

/// Independent direct implementation of tie-grouped TDC+ q-values.
fn oracle_q(scores: &[f64], labels: &[i8], p: f64) -> Vec<f64> {
    assert_eq!(scores.len(), labels.len());
    let lambda = p / (1.0 - p);
    let mut order: Vec<usize> = (0..scores.len()).collect();
    order.sort_by(|&left, &right| scores[right].total_cmp(&scores[left]));
    let mut raw = vec![1.0; scores.len()];
    let mut targets = 0usize;
    let mut decoys = 1usize;
    let mut start = 0usize;
    while start < order.len() {
        let score = scores[order[start]];
        let mut end = start + 1;
        while end < order.len() && scores[order[end]] == score {
            end += 1;
        }
        for &row in &order[start..end] {
            if labels[row] > 0 {
                targets += 1;
            } else {
                decoys += 1;
            }
        }
        let fdp = (lambda * decoys as f64 / targets.max(1) as f64).min(1.0);
        raw[start..end].fill(fdp);
        start = end;
    }
    let mut q = vec![1.0; scores.len()];
    let mut running = 1.0f64;
    for rank in (0..order.len()).rev() {
        running = running.min(raw[rank]);
        q[order[rank]] = running;
    }
    q
}

fn check_q_case(name: &str, scores: Vec<f64>, labels: Vec<i8>) {
    let expected = oracle_q(&scores, &labels, 0.5);
    let actual = qvalues(&scores, &labels, Tdc::reported(0.5));
    assert_vec(&actual, &expected, name);
    for &value in &actual {
        assert!(
            value.is_finite() && (0.0..=1.0).contains(&value),
            "{name}: bounds"
        );
    }
    let mut order: Vec<usize> = (0..scores.len()).collect();
    order.sort_by(|&left, &right| scores[right].total_cmp(&scores[left]));
    for pair in order.windows(2) {
        assert!(
            actual[pair[0]] <= actual[pair[1]] + EPS,
            "{name}: nonmonotone"
        );
        if scores[pair[0]] == scores[pair[1]] {
            assert!(
                approx(actual[pair[0]], actual[pair[1]]),
                "{name}: split tie"
            );
        }
    }
}

fn permutations(values: &mut [usize], start: usize, output: &mut Vec<Vec<usize>>) {
    if start == values.len() {
        output.push(values.to_vec());
        return;
    }
    for index in start..values.len() {
        values.swap(start, index);
        permutations(values, start + 1, output);
        values.swap(start, index);
    }
}

fn threshold_boundary(threshold: f64, targets: usize) {
    // One observed decoy plus the finite-sample safeguard gives numerator 2.
    let mut scores = Vec::with_capacity(targets + 1);
    let mut labels = Vec::with_capacity(targets + 1);
    // Put the observed decoy first so every target-containing threshold includes
    // it; otherwise a target-only prefix legitimately attains 1/T instead.
    scores.push((targets + 1) as f64);
    labels.push(-1);
    for index in 0..targets {
        scores.push((targets - index) as f64);
        labels.push(1);
    }
    let q = qvalues(&scores, &labels, Tdc::reported(0.5));
    let accepted = q
        .iter()
        .zip(&labels)
        .filter(|(q, label)| **label > 0 && **q < threshold)
        .count();
    assert_eq!(
        accepted, 0,
        "q exactly {threshold} must fail strict q<{threshold}"
    );
    assert!(approx(q[1], threshold), "boundary {threshold}: q={}", q[1]);

    scores.push(-1.0);
    labels.push(1);
    let q = qvalues(&scores, &labels, Tdc::reported(0.5));
    let accepted = q
        .iter()
        .zip(&labels)
        .filter(|(q, label)| **label > 0 && **q < threshold)
        .count();
    assert_eq!(
        accepted,
        targets + 1,
        "one extra target must cross q<{threshold}"
    );
}

fn main() {
    let cases = vec![
        ("empty", vec![], vec![]),
        ("single-target", vec![3.0], vec![1]),
        ("single-decoy", vec![3.0], vec![-1]),
        ("target-decoy", vec![2.0, 1.0], vec![1, -1]),
        ("decoy-target", vec![2.0, 1.0], vec![-1, 1]),
        ("all-target", vec![6.0, 5.0, 4.0, 3.0, 2.0, 1.0], vec![1; 6]),
        ("all-decoy", vec![6.0, 5.0, 4.0, 3.0, 2.0, 1.0], vec![-1; 6]),
        (
            "alternating",
            vec![9.0, 8.0, 7.0, 6.0, 5.0, 4.0, 3.0],
            vec![1, -1, 1, -1, 1, -1, 1],
        ),
        (
            "mixed-tie",
            vec![9.0, 8.0, 8.0, 8.0, 7.0, 6.0],
            vec![1, 1, -1, 1, -1, 1],
        ),
        (
            "all-tied",
            vec![4.0; 9],
            vec![1, -1, 1, 1, -1, -1, 1, -1, 1],
        ),
        (
            "extreme-values",
            vec![f64::MAX, 1e300, 0.0, -0.0, -1e300, -f64::MAX],
            vec![1, -1, 1, -1, 1, -1],
        ),
        (
            "target-heavy",
            (0..101).map(|i| 200.0 - i as f64).collect(),
            (0..101).map(|i| if i == 100 { -1 } else { 1 }).collect(),
        ),
        (
            "decoy-heavy",
            (0..101).map(|i| 200.0 - i as f64).collect(),
            (0..101).map(|i| if i == 100 { 1 } else { -1 }).collect(),
        ),
    ];
    for (name, scores, labels) in cases {
        check_q_case(name, scores, labels);
    }

    // Every permutation of a mixed four-row score tie has the same values by identity.
    let base_scores = [5.0, 5.0, 5.0, 5.0, 1.0];
    let base_labels = [1, 1, 1, -1, -1];
    let mut seed = vec![0usize, 1, 2, 3];
    let mut orders = Vec::new();
    permutations(&mut seed, 0, &mut orders);
    let reference = qvalues(&base_scores, &base_labels, Tdc::reported(0.5));
    for order in orders {
        let mut scores: Vec<f64> = order.iter().map(|&row| base_scores[row]).collect();
        let mut labels: Vec<i8> = order.iter().map(|&row| base_labels[row]).collect();
        scores.push(base_scores[4]);
        labels.push(base_labels[4]);
        let q = qvalues(&scores, &labels, Tdc::reported(0.5));
        let mut restored = vec![0.0; 5];
        for (position, &row) in order.iter().enumerate() {
            restored[row] = q[position];
        }
        restored[4] = q[4];
        assert_vec(&restored, &reference, "24 tie permutations");
    }

    // Exact strict boundaries and just-below cases at every requested threshold.
    for (threshold, targets) in [
        (0.001, 2000usize),
        (0.005, 400),
        (0.01, 200),
        (0.02, 100),
        (0.05, 40),
        (0.10, 20),
    ] {
        threshold_boundary(threshold, targets);
    }

    // Closed-form PEP cases, independent of the implementation's internal PAVA.
    let (q, pep) = qvalues_and_peps(
        &[5.0, 4.0, 3.0, 2.0, 1.0, 0.0],
        &[1, 1, 1, 1, 1, -1],
        Tdc::reported(0.5),
    );
    assert_vec(&q[..5], &[0.2; 5], "five-target q");
    assert_vec(&pep[..5], &[0.2; 5], "five-target PEP");

    let (_, pep) = qvalues_and_peps(
        &[10.0, 9.0, 8.0, 7.0, 6.0, 5.0],
        &[1, 1, -1, 1, 1, 1],
        Tdc::reported(0.5),
    );
    let target_sum: f64 = pep
        .iter()
        .zip([1, 1, -1, 1, 1, 1])
        .filter(|(_, label)| *label > 0)
        .map(|(pep, _)| pep)
        .sum();
    assert!(
        approx(target_sum, 2.0),
        "two estimated false discoveries, got {target_sum}"
    );
    assert!(pep
        .iter()
        .all(|value| value.is_finite() && (0.0..=1.0).contains(value)));

    let (_, pep) = qvalues_and_peps(&[3.0, 2.0, 1.0], &[-1, -1, -1], Tdc::reported(0.5));
    assert_vec(&pep, &[1.0, 1.0, 1.0], "decoy-only PEP placeholder");

    let (_, pep) = qvalues_and_peps(&[2.0, 1.0], &[1, -1], Tdc::reported(0.5));
    assert_vec(&pep, &[1.0, 1.0], "one-target PEP");

    // Publicly accepted extreme p exposes the PEP floor: conservation and exact
    // opportunity-ratio scaling no longer hold, contrary to the broad docs.
    let tiny_p = 1e-15;
    let (q, pep) = qvalues_and_peps(&[5.0, 4.0, 3.0, 2.0, 1.0], &[1; 5], Tdc::reported(tiny_p));
    let q_false_count = 5.0 * q[0];
    let pep_sum: f64 = pep.iter().sum();
    assert!(
        pep_sum > q_false_count * 1000.0,
        "floor boundary was expected to break conservation"
    );
    println!("EXTREME_P_FLOOR q_false_count={q_false_count:.17e} pep_sum={pep_sum:.17e}");

    println!("PASS: independent q-value oracle cases, 24 tie permutations, six strict boundaries, and closed-form PEP cases");
}
