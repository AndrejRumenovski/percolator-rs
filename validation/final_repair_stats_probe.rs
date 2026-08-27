//! Fresh independent arithmetic probe for the final post-repair audit.
//!
//! This is intentionally outside the Cargo test graph.  It enumerates score
//! partitions and label patterns against a separately written TDC+ oracle,
//! checks the optimized count/mask paths against materialized q-values, and
//! preserves the known extreme-opportunity PEP counterexample.

#[path = "../src/stats.rs"]
mod stats;

use stats::{
    qvalues, qvalues_and_peps, qvalues_into, target_count_at_fdr_into,
    target_mask_at_fdr_into, QValueWorkspace, Tdc,
};

fn close(left: f64, right: f64) -> bool {
    (left - right).abs() <= 2e-13 * (1.0 + left.abs().max(right.abs()))
}

fn oracle(scores: &[f64], labels: &[i8], p: f64) -> Vec<f64> {
    assert_eq!(scores.len(), labels.len());
    let lambda = p / (1.0 - p);
    let mut order: Vec<usize> = (0..scores.len()).collect();
    order.sort_unstable_by(|&left, &right| scores[right].total_cmp(&scores[left]));
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
    let mut output = vec![1.0; scores.len()];
    let mut running = 1.0f64;
    for rank in (0..order.len()).rev() {
        running = running.min(raw[rank]);
        output[order[rank]] = running;
    }
    output
}

fn score_patterns(n: usize) -> Vec<Vec<f64>> {
    let mut output = Vec::new();
    output.push(vec![7.0; n]);
    output.push((0..n).map(|row| (n - row) as f64).collect());
    output.push((0..n).map(|row| ((n - row) / 2) as f64).collect());
    output.push(
        (0..n)
            .map(|row| match row % 5 {
                0 | 1 => 3.0,
                2 => 1.0,
                _ => -2.0,
            })
            .collect(),
    );
    if n >= 2 {
        let base = 4.0f64;
        output.push(
            (0..n)
                .map(|row| {
                    if row == 0 {
                        f64::from_bits(base.to_bits() + 1)
                    } else if row == 1 {
                        base
                    } else {
                        -(row as f64)
                    }
                })
                .collect(),
        );
    }
    output
}

fn main() {
    let thresholds = [0.001, 0.005, 0.01, 0.02, 0.05, 0.10, 0.37, 1.0];
    let probabilities = [0.2, 1.0 / 3.0, 0.5, 0.8];
    let mut cases = 0usize;
    let mut fast_path_checks = 0usize;

    for n in 0..=10usize {
        let label_patterns = if n == 0 { 1 } else { 1usize << n };
        for bits in 0..label_patterns {
            let labels: Vec<i8> = (0..n)
                .map(|row| if bits & (1usize << row) == 0 { -1 } else { 1 })
                .collect();
            for scores in score_patterns(n) {
                for p in probabilities {
                    let expected = oracle(&scores, &labels, p);
                    let actual = qvalues(&scores, &labels, Tdc::reported(p));
                    assert_eq!(actual.len(), expected.len());
                    for (row, (&observed, &wanted)) in actual.iter().zip(&expected).enumerate() {
                        assert!(
                            close(observed, wanted),
                            "oracle mismatch n={n} bits={bits} p={p} row={row}: {observed:.17} != {wanted:.17}"
                        );
                    }

                    let mut workspace = QValueWorkspace::default();
                    let mut reused = vec![f64::NAN; n + 7];
                    qvalues_into(
                        &scores,
                        &labels,
                        Tdc::reported(p),
                        &mut workspace,
                        &mut reused,
                    );
                    assert_eq!(reused, actual, "reused q-value buffer mismatch");

                    for threshold in thresholds {
                        let materialized = actual
                            .iter()
                            .zip(&labels)
                            .filter(|(qvalue, label)| **label > 0 && **qvalue < threshold)
                            .count();
                        let mut order = vec![usize::MAX; n + 11];
                        let fast = target_count_at_fdr_into(
                            &scores,
                            &labels,
                            Tdc::reported(p),
                            threshold,
                            &mut order,
                        );
                        assert_eq!(fast, materialized, "fast count mismatch");

                        let mut accepted = vec![9u8; n + 13];
                        target_mask_at_fdr_into(
                            &scores,
                            &labels,
                            Tdc::reported(p),
                            threshold,
                            &mut workspace,
                            &mut accepted,
                        );
                        assert_eq!(accepted.len(), n);
                        for row in 0..n {
                            assert_eq!(
                                accepted[row] != 0,
                                labels[row] > 0 && actual[row] < threshold,
                                "fast mask mismatch"
                            );
                        }
                        fast_path_checks += 1;
                    }
                    cases += 1;
                }
            }
        }
    }

    // Exact threshold boundaries: equality must not pass a strict threshold.
    for (threshold, targets) in [
        (0.001, 2_000usize),
        (0.005, 400),
        (0.01, 200),
        (0.02, 100),
        (0.05, 40),
        (0.10, 20),
    ] {
        let mut scores = vec![(targets + 1) as f64];
        let mut labels = vec![-1i8];
        for rank in 0..targets {
            scores.push((targets - rank) as f64);
            labels.push(1);
        }
        let qvalue = qvalues(&scores, &labels, Tdc::reported(0.5));
        assert!(close(qvalue[1], threshold));
        assert_eq!(
            qvalue
                .iter()
                .zip(&labels)
                .filter(|(q, label)| **label > 0 && **q < threshold)
                .count(),
            0
        );
    }

    // Ordinary PEP conservation holds away from the explicit numerical floor.
    let scores = [12.0, 11.0, 10.0, 9.0, 9.0, 8.0, 7.0, 6.0];
    let labels = [1, 1, -1, 1, -1, 1, 1, -1];
    let (_, pep) = qvalues_and_peps(&scores, &labels, Tdc::reported(0.5));
    let target_sum: f64 = pep
        .iter()
        .zip(labels)
        .filter(|(_, label)| *label > 0)
        .map(|(value, _)| *value)
        .sum();
    // The final score group is decoy-only, so its increment has no lower target
    // to carry it; the estimate at the worst target score is three.
    assert!(close(target_sum, 3.0), "ordinary PEP mass={target_sum}");
    assert!(pep.iter().all(|value| value.is_finite() && (0.0..=1.0).contains(value)));

    // The accepted public p range still falsifies the exact PEP-conservation
    // and exact opportunity-scaling claims because the 1e-12 floor dominates.
    let tiny_p = 1e-15;
    let (qvalue, pep) = qvalues_and_peps(
        &[5.0, 4.0, 3.0, 2.0, 1.0],
        &[1, 1, 1, 1, 1],
        Tdc::reported(tiny_p),
    );
    let estimated_false = 5.0 * qvalue[0];
    let pep_sum: f64 = pep.iter().sum();
    assert!(pep_sum > estimated_false * 1_000.0);

    println!(
        "PASS q_oracle_cases={cases} fast_path_checks={fast_path_checks} strict_boundaries=6 ordinary_pep_mass={target_sum:.17}"
    );
    println!(
        "FAILURE EXTREME_P_FLOOR p={tiny_p:.17e} estimated_false={estimated_false:.17e} pep_sum={pep_sum:.17e} ratio={:.3}",
        pep_sum / estimated_false
    );
}
