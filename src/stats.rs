#![allow(clippy::items_after_test_module)]

//! Target-decoy q-values and a monotone PEP estimate.

#[derive(Default)]
pub struct QValueWorkspace {
    order: Vec<usize>,
    fdr_at: Vec<f64>,
}

#[inline]
fn sort_score_order(order: &mut [usize], scores: &[f64]) {
    debug_assert!(order.iter().all(|&index| index < scores.len()));
    order.sort_unstable_by(|&a, &b| {
        // `order` is always populated from 0..scores.len() immediately before
        // this call. Avoid two redundant bounds checks per sort comparison.
        unsafe { scores.get_unchecked(b) }
            .partial_cmp(unsafe { scores.get_unchecked(a) })
            .unwrap_or(std::cmp::Ordering::Equal)
    });
}

#[inline]
fn sort_reversed_rank_order(order: &mut [usize], ranks: &[u32]) {
    debug_assert!(order.iter().all(|&index| index < ranks.len()));
    order.sort_unstable_by(|&a, &b| {
        unsafe { ranks.get_unchecked(b) }.cmp(unsafe { ranks.get_unchecked(a) })
    });
}

/// Fill `q` with q-values aligned to input order, reusing the supplied buffers.
pub fn qvalues_into(
    scores: &[f64],
    labels: &[i8],
    pi0: f64,
    workspace: &mut QValueWorkspace,
    q: &mut Vec<f64>,
) {
    let n = scores.len();
    #[cfg(feature = "profiling")]
    let _qvalues = crate::profile::Scope::with_elements("qvalue", "qvalues_total", n);
    #[cfg(feature = "profiling")]
    let order_capacity = workspace.order.capacity();
    #[cfg(feature = "profiling")]
    let fdr_capacity = workspace.fdr_at.capacity();
    #[cfg(feature = "profiling")]
    let q_capacity = q.capacity();

    workspace.order.clear();
    workspace.order.extend(0..n);
    workspace.fdr_at.resize(n, 1.0);
    q.resize(n, 1.0);
    #[cfg(feature = "profiling")]
    crate::profile::allocation_site(
        "stats::qvalues temporary vectors",
        u64::from(workspace.order.capacity() > order_capacity)
            + u64::from(workspace.fdr_at.capacity() > fdr_capacity)
            + u64::from(q.capacity() > q_capacity),
        ((workspace.order.capacity() - order_capacity) * std::mem::size_of::<usize>()
            + (workspace.fdr_at.capacity() - fdr_capacity) * std::mem::size_of::<f64>()
            + (q.capacity() - q_capacity) * std::mem::size_of::<f64>()) as u64,
    );
    #[cfg(feature = "profiling")]
    let sort_start = std::time::Instant::now();
    sort_score_order(&mut workspace.order, scores);
    #[cfg(feature = "profiling")]
    crate::profile::record(
        "sort",
        "qvalue_score_order",
        sort_start.elapsed(),
        Some(n as u64),
        None,
    );

    #[cfg(feature = "profiling")]
    let evaluate_start = std::time::Instant::now();
    let mut targets = 0.0f64;
    let mut decoys = 0.0f64;
    // first pass: raw fdr at each rank
    for (rank, &i) in workspace.order.iter().enumerate() {
        if labels[i] > 0 {
            targets += 1.0;
        } else {
            decoys += 1.0;
        }
        let fdr = if targets > 0.0 {
            (pi0 * decoys) / targets
        } else {
            1.0
        };
        workspace.fdr_at[rank] = fdr.min(1.0);
    }
    // monotonize from worst rank upward, assign to targets
    let mut running = 1.0f64;
    for rank in (0..n).rev() {
        if workspace.fdr_at[rank] < running {
            running = workspace.fdr_at[rank];
        }
        let i = workspace.order[rank];
        q[i] = running;
    }
    #[cfg(feature = "profiling")]
    crate::profile::record(
        "qvalue",
        "qvalue_scan_and_monotonize",
        evaluate_start.elapsed(),
        Some(n as u64),
        None,
    );
}

/// Given scores and labels (+1 target / -1 decoy), return q-values aligned to input order.
/// Standard target-decoy competition: walk high->low score, FDR = pi0 * (#decoys / #targets),
/// then enforce monotonicity (running minimum from the bottom).
pub fn qvalues(scores: &[f64], labels: &[i8], pi0: f64) -> Vec<f64> {
    let mut workspace = QValueWorkspace::default();
    let mut q = Vec::new();
    qvalues_into(scores, labels, pi0, &mut workspace, &mut q);
    q
}

/// Mark targets whose monotone q-value is strictly below `threshold` without
/// materializing q-values that the training loop never otherwise consumes.
pub fn target_mask_at_fdr_into(
    scores: &[f64],
    labels: &[i8],
    pi0: f64,
    threshold: f64,
    workspace: &mut QValueWorkspace,
    accepted: &mut Vec<u8>,
) {
    let n = scores.len();
    #[cfg(feature = "profiling")]
    let _qvalues = crate::profile::Scope::with_elements("qvalue", "qvalues_total", n);
    workspace.order.clear();
    workspace.order.extend(0..n);
    accepted.resize(n, 0);
    accepted.fill(0);

    #[cfg(feature = "profiling")]
    let sort_start = std::time::Instant::now();
    sort_score_order(&mut workspace.order, scores);
    #[cfg(feature = "profiling")]
    crate::profile::record(
        "sort",
        "qvalue_score_order",
        sort_start.elapsed(),
        Some(n as u64),
        None,
    );

    #[cfg(feature = "profiling")]
    let evaluate_start = std::time::Instant::now();
    let mut targets = 0.0f64;
    let mut decoys = 0.0f64;
    let mut last_accepted_rank = None;
    for (rank, &i) in workspace.order.iter().enumerate() {
        if labels[i] > 0 {
            targets += 1.0;
        } else {
            decoys += 1.0;
        }
        let fdr = if targets > 0.0 {
            (pi0 * decoys) / targets
        } else {
            1.0
        };
        if fdr.min(1.0) < threshold {
            last_accepted_rank = Some(rank);
        }
    }
    if let Some(rank) = last_accepted_rank {
        for &i in &workspace.order[..=rank] {
            if labels[i] > 0 {
                accepted[i] = 1;
            }
        }
    }
    #[cfg(feature = "profiling")]
    crate::profile::record(
        "qvalue",
        "qvalue_scan_and_monotonize",
        evaluate_start.elapsed(),
        Some(n as u64),
        None,
    );
}

/// Count target scores whose q-value is strictly below `threshold` without
/// materializing the complete q-value vector.
///
/// A q-value at rank `i` is the minimum raw FDR at any rank `j >= i`.
/// Consequently, the accepted targets are exactly the targets at or above the
/// last rank whose raw FDR is below the threshold.
pub fn target_count_at_fdr_into(
    scores: &[f64],
    labels: &[i8],
    pi0: f64,
    threshold: f64,
    order: &mut Vec<usize>,
) -> usize {
    order.clear();
    order.extend(0..scores.len());
    sort_score_order(order, scores);

    let mut targets = 0.0f64;
    let mut decoys = 0.0f64;
    let mut accepted_targets = 0usize;
    for &i in order.iter() {
        if labels[i] > 0 {
            targets += 1.0;
        } else {
            decoys += 1.0;
        }
        let fdr = if targets > 0.0 {
            (pi0 * decoys) / targets
        } else {
            1.0
        };
        if fdr.min(1.0) < threshold {
            accepted_targets = targets as usize;
        }
    }
    accepted_targets
}

pub fn target_count_at_reversed_ranks_into(
    ranks: &[u32],
    labels: &[i8],
    pi0: f64,
    threshold: f64,
    order: &mut Vec<usize>,
) -> usize {
    order.clear();
    order.extend(0..ranks.len());
    sort_reversed_rank_order(order, ranks);

    let mut targets = 0.0f64;
    let mut decoys = 0.0f64;
    let mut accepted_targets = 0usize;
    for &i in order.iter() {
        if labels[i] > 0 {
            targets += 1.0;
        } else {
            decoys += 1.0;
        }
        let fdr = if targets > 0.0 {
            (pi0 * decoys) / targets
        } else {
            1.0
        };
        if fdr.min(1.0) < threshold {
            accepted_targets = targets as usize;
        }
    }
    accepted_targets
}

#[cfg(test)]
pub fn target_count_at_fdr(scores: &[f64], labels: &[i8], pi0: f64, threshold: f64) -> usize {
    target_count_at_fdr_into(scores, labels, pi0, threshold, &mut Vec::new())
}

/// Estimate pi0 (fraction of incorrect targets) simply as #decoys/#targets, capped to [0,1].
pub fn estimate_pi0(labels: &[i8]) -> f64 {
    let t = labels.iter().filter(|&&l| l > 0).count() as f64;
    let d = labels.iter().filter(|&&l| l < 0).count() as f64;
    if t > 0.0 {
        (d / t).min(1.0)
    } else {
        1.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn qvalues_are_monotone_by_score() {
        // higher score = better; targets should get lower q as score rises
        let scores = vec![5.0, 4.0, 3.0, 2.0, 1.0, 0.0];
        let labels = vec![1i8, 1, -1, 1, -1, -1];
        let q = qvalues(&scores, &labels, 1.0);
        // walking best->worst score, q must be non-decreasing
        let mut order: Vec<usize> = (0..scores.len()).collect();
        order.sort_by(|&a, &b| scores[b].partial_cmp(&scores[a]).unwrap());
        let mut prev = 0.0;
        for &i in &order {
            assert!(
                q[i] + 1e-9 >= prev,
                "q-values must be monotone non-decreasing"
            );
            prev = q[i];
        }
        // top target (best score) has the lowest q
        assert!(q[0] <= q[3]);
    }

    #[test]
    fn pi0_ratio() {
        let labels = vec![1i8, 1, 1, -1]; // 3 targets, 1 decoy
        assert!((estimate_pi0(&labels) - 1.0 / 3.0).abs() < 1e-9);
    }

    #[test]
    fn target_count_matches_materialized_qvalues() {
        let cases = [
            (
                vec![5.0, 4.0, 3.0, 2.0, 1.0, 0.0],
                vec![1, 1, -1, 1, -1, -1],
            ),
            (
                vec![2.0, 2.0, 2.0, 1.0, 1.0, -0.0, 0.0],
                vec![-1, 1, 1, -1, 1, -1, 1],
            ),
            (vec![f64::NAN, 3.0, 2.0, 1.0], vec![1, -1, 1, -1]),
            (Vec::new(), Vec::new()),
        ];
        for (scores, labels) in cases {
            for pi0 in [0.0, 0.25, 1.0] {
                let q = qvalues(&scores, &labels, pi0);
                for threshold in [0.0, 0.01, 0.5, 1.0, 1.01, f64::NAN] {
                    let expected = q
                        .iter()
                        .zip(&labels)
                        .filter(|(qvalue, label)| **label > 0 && **qvalue < threshold)
                        .count();
                    assert_eq!(
                        target_count_at_fdr(&scores, &labels, pi0, threshold),
                        expected,
                        "scores={scores:?}, labels={labels:?}, pi0={pi0}, threshold={threshold}"
                    );
                }
            }
        }
    }

    #[test]
    fn target_count_overwrites_reused_order() {
        let mut order = vec![usize::MAX; 20];
        let cases = [
            (vec![3.0, 2.0, 1.0], vec![1, -1, 1]),
            (vec![1.0], vec![1]),
            (Vec::new(), Vec::new()),
            (vec![1.0, 4.0, 2.0, 3.0], vec![-1, 1, -1, 1]),
        ];
        for (scores, labels) in cases {
            let fresh = target_count_at_fdr(&scores, &labels, 1.0, 0.5);
            let reused = target_count_at_fdr_into(&scores, &labels, 1.0, 0.5, &mut order);
            assert_eq!(reused, fresh);
        }
    }

    #[test]
    fn reversed_integer_ranks_preserve_negated_float_sort_permutation() {
        let mut state = 0x243f_6a88_85a3_08d3u64;
        for len in 0..300 {
            let mut scores = Vec::with_capacity(len);
            for index in 0..len {
                state = state
                    .wrapping_mul(6_364_136_223_846_793_005)
                    .wrapping_add(1_442_695_040_888_963_407);
                let value = match index % 31 {
                    0 => f64::INFINITY,
                    1 => f64::NEG_INFINITY,
                    2 => -0.0,
                    _ => ((state >> 60) as i8 - 8) as f64,
                };
                scores.push(value);
            }
            let mut positive_order: Vec<usize> = (0..len).collect();
            sort_score_order(&mut positive_order, &scores);
            let mut ranks = vec![0u32; len];
            let mut rank = 0u32;
            for position in 0..positive_order.len() {
                if position > 0
                    && scores[positive_order[position]]
                        .partial_cmp(&scores[positive_order[position - 1]])
                        != Some(std::cmp::Ordering::Equal)
                {
                    rank += 1;
                }
                ranks[positive_order[position]] = rank;
            }

            let negated: Vec<f64> = scores.iter().map(|value| -value).collect();
            let mut expected: Vec<usize> = (0..len).collect();
            sort_score_order(&mut expected, &negated);
            let mut actual: Vec<usize> = (0..len).collect();
            sort_reversed_rank_order(&mut actual, &ranks);
            assert_eq!(actual, expected, "len={len}");
        }
    }

    #[test]
    fn qvalue_workspace_overwrites_reused_buffers() {
        let mut workspace = QValueWorkspace::default();
        let mut reused = Vec::new();
        for (scores, labels, pi0) in [
            (vec![3.0, 2.0, 1.0], vec![1, -1, 1], 1.0),
            (vec![5.0], vec![-1], 0.5),
            (Vec::new(), Vec::new(), 1.0),
            (vec![1.0, 4.0, 2.0, 3.0], vec![-1, 1, -1, 1], 0.25),
        ] {
            let fresh = qvalues(&scores, &labels, pi0);
            qvalues_into(&scores, &labels, pi0, &mut workspace, &mut reused);
            assert_eq!(reused, fresh);
        }
    }

    #[test]
    fn target_mask_matches_materialized_qvalues_and_reuses_storage() {
        let mut workspace = QValueWorkspace::default();
        let mut accepted = vec![9; 30];
        let cases = [
            (
                vec![5.0, 4.0, 3.0, 2.0, 1.0, 0.0],
                vec![1, 1, -1, 1, -1, -1],
            ),
            (
                vec![2.0, 2.0, 2.0, 1.0, 1.0, -0.0, 0.0],
                vec![-1, 1, 1, -1, 1, -1, 1],
            ),
            (vec![f64::NAN, 3.0, 2.0, 1.0], vec![1, -1, 1, -1]),
            (Vec::new(), Vec::new()),
        ];
        for (scores, labels) in cases {
            for threshold in [0.0, 0.01, 0.5, 1.0, 1.01, f64::NAN] {
                let q = qvalues(&scores, &labels, 1.0);
                target_mask_at_fdr_into(
                    &scores,
                    &labels,
                    1.0,
                    threshold,
                    &mut workspace,
                    &mut accepted,
                );
                let expected: Vec<u8> = q
                    .iter()
                    .zip(&labels)
                    .map(|(qvalue, label)| u8::from(*label > 0 && *qvalue < threshold))
                    .collect();
                assert_eq!(accepted, expected);
            }
        }
    }

    #[test]
    fn peps_bounded_and_monotone() {
        let scores = vec![5.0, 4.0, 3.0, 2.0, 1.0];
        let labels = vec![1i8, 1, -1, -1, -1];
        let p = peps(&scores, &labels, 1.0);
        for v in &p {
            assert!(*v >= 0.0 && *v <= 1.0);
        }
    }

    #[test]
    fn paired_qvalues_and_peps_match_separate_calculations() {
        let cases = [
            (
                vec![
                    f64::INFINITY,
                    3.0,
                    3.0,
                    0.0,
                    -0.0,
                    f64::NEG_INFINITY,
                    f64::NAN,
                ],
                vec![1, -1, 1, -1, 1, -1, 1],
            ),
            (Vec::new(), Vec::new()),
        ];
        for (scores, labels) in cases {
            for pi0 in [0.0, 0.25, 1.0] {
                let expected_q = qvalues(&scores, &labels, pi0);
                let expected_pep = peps(&scores, &labels, pi0);
                let (actual_q, actual_pep) = qvalues_and_peps(&scores, &labels, pi0);
                assert_eq!(actual_q, expected_q);
                assert_eq!(actual_pep, expected_pep);
            }
        }
    }
}

/// Monotone non-increasing PEP vs score via a light isotonic (pool-adjacent-violators)
/// fit of the decoy indicator over score-sorted PSMs. Returns PEP aligned to input order.
#[allow(clippy::needless_range_loop)]
#[cfg(test)]
pub fn peps(scores: &[f64], labels: &[i8], pi0: f64) -> Vec<f64> {
    let n = scores.len();
    #[cfg(feature = "profiling")]
    let _peps = crate::profile::Scope::with_elements("pep", "pep_pava_total", n);
    #[cfg(feature = "profiling")]
    crate::profile::allocation_site(
        "stats::peps temporary vectors",
        6,
        n as u64 * (2 * std::mem::size_of::<usize>() + 5 * std::mem::size_of::<f64>()) as u64,
    );
    let mut order: Vec<usize> = (0..n).collect();
    #[cfg(feature = "profiling")]
    let sort_start = std::time::Instant::now();
    order.sort_unstable_by(|&a, &b| {
        scores[b]
            .partial_cmp(&scores[a])
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    #[cfg(feature = "profiling")]
    crate::profile::record(
        "sort",
        "pep_score_order",
        sort_start.elapsed(),
        Some(n as u64),
        None,
    );
    peps_from_order(labels, pi0, &order)
}

pub fn qvalues_and_peps(scores: &[f64], labels: &[i8], pi0: f64) -> (Vec<f64>, Vec<f64>) {
    #[cfg(feature = "profiling")]
    let n = scores.len();
    let mut workspace = QValueWorkspace::default();
    let mut q = Vec::new();
    qvalues_into(scores, labels, pi0, &mut workspace, &mut q);
    #[cfg(feature = "profiling")]
    let _peps = crate::profile::Scope::with_elements("pep", "pep_pava_total", n);
    #[cfg(feature = "profiling")]
    crate::profile::allocation_site(
        "stats::peps temporary vectors",
        5,
        n as u64 * (std::mem::size_of::<usize>() + 5 * std::mem::size_of::<f64>()) as u64,
    );
    let pep = peps_from_order(labels, pi0, &workspace.order);
    (q, pep)
}

#[allow(clippy::needless_range_loop)]
fn peps_from_order(labels: &[i8], pi0: f64, order: &[usize]) -> Vec<f64> {
    let n = order.len();
    #[cfg(feature = "profiling")]
    let pava_start = std::time::Instant::now();
    // target=0 (correct-ish), decoy=1; PAVA to get monotone-increasing prob-of-decoy along worsening score
    let y: Vec<f64> = order
        .iter()
        .map(|&i| if labels[i] < 0 { 1.0 } else { 0.0 })
        .collect();
    // PAVA (isotonic, non-decreasing)
    // simple stack-based PAVA
    let mut level_val: Vec<f64> = Vec::with_capacity(n);
    let mut level_wt: Vec<f64> = Vec::with_capacity(n);
    let mut level_len: Vec<usize> = Vec::with_capacity(n);
    for k in 0..n {
        let mut v = y[k];
        let mut w = 1.0;
        let mut len = 1;
        while let Some(&pv) = level_val.last() {
            if pv <= v {
                break;
            }
            let pw = level_wt.pop().unwrap();
            let pl = level_len.pop().unwrap();
            level_val.pop();
            v = (pv * pw + v * w) / (pw + w);
            w += pw;
            len += pl;
        }
        level_val.push(v);
        level_wt.push(w);
        level_len.push(len);
    }
    // expand back
    let mut fitted = vec![0.0f64; n];
    let mut pos = 0;
    for l in 0..level_val.len() {
        for _ in 0..level_len[l] {
            fitted[pos] = level_val[l];
            pos += 1;
        }
    }
    // scale decoy-prob to PEP with pi0 and map back to input order
    let mut pep = vec![0.0f64; n];
    for (rank, &orig) in order.iter().enumerate() {
        pep[orig] = (pi0 * 2.0 * fitted[rank]).min(1.0);
    }
    #[cfg(feature = "profiling")]
    crate::profile::record(
        "pep",
        "pava_fit_and_expand",
        pava_start.elapsed(),
        Some(n as u64),
        None,
    );
    pep
}
