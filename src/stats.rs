#![allow(clippy::items_after_test_module)]

//! Target-decoy competition (TDC) q-values and posterior error probabilities.
//!
//! # Supported input contract
//!
//! These estimators are valid for a **concatenated** target-decoy search in which
//! the supplied rows are the winners of a spectrum-level competition against a
//! decoy database of the same size as the target database.  Under that design an
//! incorrect target beats its paired decoy with probability
//! [`Tdc::null_target_win_prob`] (0.5 for a 1:1 competition), which converts an
//! observed decoy count into an expected count of incorrect targets through the
//! *opportunity ratio* `p / (1 - p)`.
//!
//! Separate target/decoy searches (mix-max post-processing) are **not** supported:
//! `pi0` is fixed at 1, which is the conservative choice for direct TDC.  Feeding
//! a separate-database search into these estimators produces an unvalidated number,
//! not a calibrated q-value.
//!
//! # Estimator
//!
//! Walking the score-sorted list from best to worst and letting `D` and `T` be the
//! decoy and target counts at or above the current *score*, the false discovery
//! proportion is estimated as
//!
//! ```text
//! FDP(s) = min(1, pi0 * (D(s) + 1) * lambda / max(1, T(s)))     lambda = p / (1 - p)
//! ```
//!
//! The `+1` is the finite-sample safeguard of Storey-style TDC estimation: with a
//! finite list one must account for the decoy that would have been observed just
//! below the threshold, otherwise a leading run of targets receives an
//! FDP estimate of exactly zero regardless of how little evidence supports it.
//! See Käll et al. (2008) and Emery et al. (2020).  q-values are the reverse
//! cumulative minimum of `FDP`, so a q-value is the smallest estimated FDP of any
//! accepted set containing the PSM.
//!
//! All PSMs sharing an exact score form one *tie group*: `FDP` is evaluated once,
//! after the whole group has been counted, and assigned to every member.  A
//! q-value is therefore a function of the score threshold alone and cannot depend
//! on the arbitrary order in which equal-scoring rows happen to appear.
//!
//! # Posterior error probabilities
//!
//! PEPs are derived from the q-values through the identity of Käll et al. (2008),
//! "Posterior error probabilities and false discovery rates: two sides of the same
//! coin": for targets ranked best-first, `q_k` is the mean PEP over the top `k`,
//! so `sum_{i<=k} PEP_i = k * q_k` and
//!
//! ```text
//! raw PEP_k = k * q_k - (k - 1) * q_{k-1}
//! ```
//!
//! A Bayesian pseudocount of half a false discovery spread uniformly over the
//! target list is added before an isotonic (PAVA) fit enforces monotonicity in
//! score.  The pseudocount is what keeps the leading tail strictly positive; it is
//! prior mass, not a clamp applied after the fact.

use std::cmp::Ordering;

/// Configuration of the target-decoy competition estimator.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Tdc {
    /// Proportion of incorrect targets.  Fixed at 1 for direct TDC, which is the
    /// conservative choice and the only design these estimators support.
    pub pi0: f64,
    /// Probability that an incorrect target outranks its paired decoy under the
    /// null.  0.5 for a 1:1 concatenated competition.
    pub null_target_win_prob: f64,
    /// Drop the finite-sample `+1` decoy safeguard.
    ///
    /// Only ever set for the semi-supervised training heuristic, which selects
    /// probable-correct targets to feed the SVM and makes no error-rate claim.
    /// Reported statistics must always keep the safeguard.
    pub skip_decoys_plus_one: bool,
}

impl Default for Tdc {
    fn default() -> Self {
        Tdc {
            pi0: 1.0,
            null_target_win_prob: 0.5,
            skip_decoys_plus_one: false,
        }
    }
}

impl Tdc {
    /// Estimator used for every reported q-value and PEP.
    pub fn reported(null_target_win_prob: f64) -> Self {
        Tdc {
            pi0: 1.0,
            null_target_win_prob,
            skip_decoys_plus_one: false,
        }
    }

    /// Estimator used to pick positive training examples inside a fold.  The
    /// safeguard is dropped so that a small training partition can still nominate
    /// positives; this is a selection heuristic, never a reported error rate.
    pub fn training(null_target_win_prob: f64) -> Self {
        Tdc {
            pi0: 1.0,
            null_target_win_prob,
            skip_decoys_plus_one: true,
        }
    }

    /// Opportunity ratio `p / (1 - p)` converting decoy counts into expected
    /// counts of incorrect targets.
    #[inline]
    fn decoy_factor(&self) -> f64 {
        let p = self.null_target_win_prob;
        if !(0.0..1.0).contains(&p) {
            return 1.0;
        }
        p / (1.0 - p)
    }

    /// Decoy count to start the scan from: the finite-sample safeguard.
    #[inline]
    fn initial_decoys(&self) -> f64 {
        if self.skip_decoys_plus_one {
            0.0
        } else {
            1.0
        }
    }

    /// `min(1, pi0 * decoys * lambda / max(1, targets))`.
    #[inline]
    fn raw_fdp(&self, decoys: f64, targets: f64) -> f64 {
        let value = self.pi0 * decoys * self.decoy_factor() / targets.max(1.0);
        if value.is_nan() {
            1.0
        } else {
            value.clamp(0.0, 1.0)
        }
    }
}

#[derive(Default)]
pub struct QValueWorkspace {
    order: Vec<usize>,
    fdr_at: Vec<f64>,
    target_ranks: Vec<usize>,
    target_pep: Vec<f64>,
    pava_value: Vec<f64>,
    pava_weight: Vec<f64>,
    pava_len: Vec<usize>,
}

/// Descending score order with a deterministic total order.
///
/// `f64::total_cmp` is a total order, so the permutation is reproducible even
/// with `-0.0`/`0.0` or repeated values.  NaN is ordered *last* rather than
/// treated as equal to everything, which would break the strict weak ordering the
/// sort requires and let a NaN acquire an arbitrary rank.
#[inline]
fn score_cmp_desc(a: f64, b: f64) -> Ordering {
    match (a.is_nan(), b.is_nan()) {
        (true, true) => Ordering::Equal,
        (true, false) => Ordering::Greater,
        (false, true) => Ordering::Less,
        (false, false) => b.total_cmp(&a),
    }
}

#[inline]
fn sort_score_order(order: &mut [usize], scores: &[f64]) {
    debug_assert!(order.iter().all(|&index| index < scores.len()));
    order.sort_unstable_by(|&a, &b| {
        // `order` is always populated from 0..scores.len() immediately before
        // this call. Avoid two redundant bounds checks per sort comparison.
        let (a, b) = unsafe { (*scores.get_unchecked(a), *scores.get_unchecked(b)) };
        score_cmp_desc(a, b)
    });
}

#[inline]
fn sort_reversed_rank_order(order: &mut [usize], ranks: &[u32]) {
    debug_assert!(order.iter().all(|&index| index < ranks.len()));
    order.sort_unstable_by(|&a, &b| {
        unsafe { ranks.get_unchecked(b) }.cmp(unsafe { ranks.get_unchecked(a) })
    });
}

/// True when rank `rank` ends a run of exactly equal scores.
///
/// Numeric equality (not the sort's total order) defines a tie group, so `-0.0`
/// and `0.0` share a group while each NaN forms its own.
#[inline]
fn ends_score_group(order: &[usize], scores: &[f64], rank: usize) -> bool {
    rank + 1 == order.len() || scores[order[rank + 1]] != scores[order[rank]]
}

#[inline]
fn ends_rank_group(order: &[usize], ranks: &[u32], rank: usize) -> bool {
    rank + 1 == order.len() || ranks[order[rank + 1]] != ranks[order[rank]]
}

/// Fill `q` with q-values aligned to input order, reusing the supplied buffers.
pub fn qvalues_into(
    scores: &[f64],
    labels: &[i8],
    tdc: Tdc,
    workspace: &mut QValueWorkspace,
    q: &mut Vec<f64>,
) {
    let n = scores.len();
    #[cfg(feature = "profiling")]
    let _qvalues = crate::profile::Scope::with_elements("qvalue", "qvalues_total", n);

    workspace.order.clear();
    workspace.order.extend(0..n);
    workspace.fdr_at.clear();
    workspace.fdr_at.resize(n, 1.0);
    q.clear();
    q.resize(n, 1.0);

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
    // Forward pass: one raw FDP per tie group, shared by every member.
    let mut targets = 0.0f64;
    let mut decoys = tdc.initial_decoys();
    let mut group_start = 0usize;
    for rank in 0..n {
        let i = workspace.order[rank];
        if labels[i] > 0 {
            targets += 1.0;
        } else {
            decoys += 1.0;
        }
        if ends_score_group(&workspace.order, scores, rank) {
            let fdp = tdc.raw_fdp(decoys, targets);
            workspace.fdr_at[group_start..=rank].fill(fdp);
            group_start = rank + 1;
        }
    }
    // Reverse cumulative minimum: q is the best achievable FDP of any set that
    // still contains this PSM.
    let mut running = 1.0f64;
    for rank in (0..n).rev() {
        if workspace.fdr_at[rank] < running {
            running = workspace.fdr_at[rank];
        }
        q[workspace.order[rank]] = running;
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
pub fn qvalues(scores: &[f64], labels: &[i8], tdc: Tdc) -> Vec<f64> {
    let mut workspace = QValueWorkspace::default();
    let mut q = Vec::new();
    qvalues_into(scores, labels, tdc, &mut workspace, &mut q);
    q
}

/// Mark targets whose monotone q-value is strictly below `threshold` without
/// materializing q-values that the training loop never otherwise consumes.
pub fn target_mask_at_fdr_into(
    scores: &[f64],
    labels: &[i8],
    tdc: Tdc,
    threshold: f64,
    workspace: &mut QValueWorkspace,
    accepted: &mut Vec<u8>,
) {
    let n = scores.len();
    #[cfg(feature = "profiling")]
    let _qvalues = crate::profile::Scope::with_elements("qvalue", "qvalues_total", n);
    workspace.order.clear();
    workspace.order.extend(0..n);
    accepted.clear();
    accepted.resize(n, 0);

    sort_score_order(&mut workspace.order, scores);

    let mut targets = 0.0f64;
    let mut decoys = tdc.initial_decoys();
    let mut last_accepted_rank = None;
    for rank in 0..n {
        let i = workspace.order[rank];
        if labels[i] > 0 {
            targets += 1.0;
        } else {
            decoys += 1.0;
        }
        if ends_score_group(&workspace.order, scores, rank)
            && tdc.raw_fdp(decoys, targets) < threshold
        {
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
}

/// Count target scores whose q-value is strictly below `threshold` without
/// materializing the complete q-value vector.
///
/// A q-value at rank `i` is the minimum raw FDP at any rank `j >= i`, so the
/// accepted targets are exactly the targets at or above the last tie group whose
/// raw FDP is below the threshold.
pub fn target_count_at_fdr_into(
    scores: &[f64],
    labels: &[i8],
    tdc: Tdc,
    threshold: f64,
    order: &mut Vec<usize>,
) -> usize {
    order.clear();
    order.extend(0..scores.len());
    sort_score_order(order, scores);

    let mut targets = 0.0f64;
    let mut decoys = tdc.initial_decoys();
    let mut accepted_targets = 0usize;
    for rank in 0..order.len() {
        let i = order[rank];
        if labels[i] > 0 {
            targets += 1.0;
        } else {
            decoys += 1.0;
        }
        if ends_score_group(order, scores, rank) && tdc.raw_fdp(decoys, targets) < threshold {
            accepted_targets = targets as usize;
        }
    }
    accepted_targets
}

/// As [`target_count_at_fdr_into`] but over dense integer ranks that already
/// collapse ties (rank increments only when the underlying score changes).
pub fn target_count_at_reversed_ranks_into(
    ranks: &[u32],
    labels: &[i8],
    tdc: Tdc,
    threshold: f64,
    order: &mut Vec<usize>,
) -> usize {
    order.clear();
    order.extend(0..ranks.len());
    sort_reversed_rank_order(order, ranks);

    let mut targets = 0.0f64;
    let mut decoys = tdc.initial_decoys();
    let mut accepted_targets = 0usize;
    for rank in 0..order.len() {
        let i = order[rank];
        if labels[i] > 0 {
            targets += 1.0;
        } else {
            decoys += 1.0;
        }
        if ends_rank_group(order, ranks, rank) && tdc.raw_fdp(decoys, targets) < threshold {
            accepted_targets = targets as usize;
        }
    }
    accepted_targets
}

#[cfg(test)]
pub fn target_count_at_fdr(scores: &[f64], labels: &[i8], tdc: Tdc, threshold: f64) -> usize {
    target_count_at_fdr_into(scores, labels, tdc, threshold, &mut Vec::new())
}

/// Q-values and PEPs from one score sort.
pub fn qvalues_and_peps(scores: &[f64], labels: &[i8], tdc: Tdc) -> (Vec<f64>, Vec<f64>) {
    let mut workspace = QValueWorkspace::default();
    let mut q = Vec::new();
    let mut pep = Vec::new();
    qvalues_into(scores, labels, tdc, &mut workspace, &mut q);
    peps_from_qvalues_into(labels, &q, &mut workspace, &mut pep);
    (q, pep)
}

#[cfg(test)]
pub fn peps(scores: &[f64], labels: &[i8], tdc: Tdc) -> Vec<f64> {
    qvalues_and_peps(scores, labels, tdc).1
}

/// Isotonic (non-decreasing) PAVA fit in place over `values`, unit weights.
fn pava_non_decreasing(
    values: &mut [f64],
    level_value: &mut Vec<f64>,
    level_weight: &mut Vec<f64>,
    level_len: &mut Vec<usize>,
) {
    level_value.clear();
    level_weight.clear();
    level_len.clear();
    for &current in values.iter() {
        let mut v = current;
        let mut w = 1.0f64;
        let mut len = 1usize;
        while let Some(&previous) = level_value.last() {
            if previous <= v {
                break;
            }
            let pw = level_weight.pop().unwrap();
            let pl = level_len.pop().unwrap();
            level_value.pop();
            v = (previous * pw + v * w) / (pw + w);
            w += pw;
            len += pl;
        }
        level_value.push(v);
        level_weight.push(w);
        level_len.push(len);
    }
    let mut position = 0usize;
    for (level, &value) in level_value.iter().enumerate() {
        for _ in 0..level_len[level] {
            values[position] = value;
            position += 1;
        }
    }
}

/// Smallest PEP the estimator will report.  Guards only against floating-point
/// underflow in the pseudocount; it is never the value that removes an otherwise
/// exact zero.
const PEP_FLOOR: f64 = 1e-12;

/// PEPs from q-values, aligned to input order.
///
/// `workspace.order` must hold the descending score order that produced `q`.
fn peps_from_qvalues_into(
    labels: &[i8],
    q: &[f64],
    workspace: &mut QValueWorkspace,
    pep: &mut Vec<f64>,
) {
    let n = workspace.order.len();
    #[cfg(feature = "profiling")]
    let _peps = crate::profile::Scope::with_elements("pep", "pep_from_qvalues_total", n);
    #[cfg(feature = "profiling")]
    let pep_start = std::time::Instant::now();
    pep.clear();
    pep.resize(n, 1.0);

    // Targets, best-first.  q along this sequence is non-decreasing because it is
    // a reverse cumulative minimum over the full list.
    workspace.target_ranks.clear();
    workspace.target_pep.clear();
    for rank in 0..n {
        if labels[workspace.order[rank]] > 0 {
            workspace.target_ranks.push(rank);
            workspace.target_pep.push(0.0);
        }
    }
    let target_count = workspace.target_ranks.len();
    if target_count == 0 {
        #[cfg(feature = "profiling")]
        crate::profile::record(
            "pep",
            "pep_from_qvalues",
            pep_start.elapsed(),
            Some(n as u64),
            None,
        );
        return;
    }

    // raw PEP_k = k * q_k - (k - 1) * q_{k-1}, plus half a false discovery of
    // prior mass spread over the list.
    let pseudocount = 0.5 / target_count as f64;
    let mut previous_scaled = 0.0f64;
    for k in 0..target_count {
        let qk = q[workspace.order[workspace.target_ranks[k]]];
        let scaled = qk * (k + 1) as f64;
        let raw = (scaled - previous_scaled).max(0.0);
        workspace.target_pep[k] = raw + pseudocount;
        previous_scaled = scaled;
    }

    let mut values = std::mem::take(&mut workspace.target_pep);
    pava_non_decreasing(
        &mut values,
        &mut workspace.pava_value,
        &mut workspace.pava_weight,
        &mut workspace.pava_len,
    );
    for value in values.iter_mut() {
        *value = value.clamp(PEP_FLOOR, 1.0);
    }
    workspace.target_pep = values;

    for k in 0..target_count {
        pep[workspace.order[workspace.target_ranks[k]]] = workspace.target_pep[k];
    }

    // Decoys carry no error-rate claim, but the reported column should stay a
    // monotone function of score.  Interpolate each decoy in q between the
    // targets that bracket it, holding the end values outside that range.
    let mut next_target = 0usize;
    for rank in 0..n {
        let row = workspace.order[rank];
        if labels[row] > 0 {
            next_target += 1;
            continue;
        }
        let value = if next_target == 0 {
            workspace.target_pep[0]
        } else if next_target == target_count {
            workspace.target_pep[target_count - 1]
        } else {
            let low_q = q[workspace.order[workspace.target_ranks[next_target - 1]]];
            let high_q = q[workspace.order[workspace.target_ranks[next_target]]];
            let low_pep = workspace.target_pep[next_target - 1];
            let high_pep = workspace.target_pep[next_target];
            if high_q > low_q {
                let position = ((q[row] - low_q) / (high_q - low_q)).clamp(0.0, 1.0);
                low_pep + position * (high_pep - low_pep)
            } else {
                high_pep
            }
        };
        pep[row] = value.clamp(PEP_FLOOR, 1.0);
    }
    #[cfg(feature = "profiling")]
    crate::profile::record(
        "pep",
        "pep_from_qvalues",
        pep_start.elapsed(),
        Some(n as u64),
        None,
    );
}

#[cfg(test)]
mod tests {
    use super::*;

    fn reported() -> Tdc {
        Tdc::reported(0.5)
    }

    /// Descending score order used by the tests to check monotonicity.
    fn descending(scores: &[f64]) -> Vec<usize> {
        let mut order: Vec<usize> = (0..scores.len()).collect();
        sort_score_order(&mut order, scores);
        order
    }

    // ----- q-value invariants -------------------------------------------------

    #[test]
    fn qvalues_are_bounded_in_the_unit_interval() {
        let scores = vec![5.0, 4.0, 3.0, 2.0, 1.0, 0.0];
        let labels = vec![1i8, 1, -1, 1, -1, -1];
        for q in qvalues(&scores, &labels, reported()) {
            assert!((0.0..=1.0).contains(&q), "q out of range: {q}");
        }
    }

    #[test]
    fn qvalues_are_monotone_by_score() {
        let scores = vec![5.0, 4.0, 3.0, 2.0, 1.0, 0.0];
        let labels = vec![1i8, 1, -1, 1, -1, -1];
        let q = qvalues(&scores, &labels, reported());
        let mut previous = 0.0;
        for &i in &descending(&scores) {
            assert!(q[i] + 1e-12 >= previous, "q-values must be non-decreasing");
            previous = q[i];
        }
    }

    /// The complete-null failure mode: a leading run of targets must never be
    /// reported at q = 0 just because no decoy happened to outrank it.
    #[test]
    fn leading_target_run_never_receives_zero_qvalue() {
        // Ten targets above every decoy, exactly balanced overall.
        let mut scores = Vec::new();
        let mut labels = Vec::new();
        for k in 0..10 {
            scores.push(100.0 - k as f64);
            labels.push(1i8);
        }
        for k in 0..30 {
            scores.push(10.0 - k as f64);
            labels.push(if k % 2 == 0 { 1 } else { -1 });
        }
        let q = qvalues(&scores, &labels, reported());
        for (i, value) in q.iter().enumerate() {
            assert!(
                *value > 0.0,
                "row {i} received q = 0 with no decoy evidence"
            );
        }
        // The best achievable estimate is the safeguard decoy over all targets
        // that could be accepted with it.
        assert!(q[0] >= 1.0 / labels.iter().filter(|&&l| l > 0).count() as f64 - 1e-12);
    }

    #[test]
    fn a_list_without_decoys_still_has_positive_qvalues() {
        let scores = vec![3.0, 2.0, 1.0];
        let labels = vec![1i8, 1, 1];
        let q = qvalues(&scores, &labels, reported());
        assert!(q.iter().all(|&v| v > 0.0));
        // One safeguard decoy over three targets.
        assert!((q[0] - 1.0 / 3.0).abs() < 1e-12, "q[0] = {}", q[0]);
    }

    #[test]
    fn qvalues_match_the_documented_closed_form() {
        let scores = vec![5.0, 4.0, 3.0, 2.0];
        let labels = vec![1i8, -1, 1, -1];
        let q = qvalues(&scores, &labels, reported());
        // rank 1: D=1(+0 seen), T=1 -> 1/1 ; rank 2: D=2, T=1 -> 2/1 -> 1
        // rank 3: D=2, T=2 -> 1 ; rank 4: D=3, T=2 -> 1
        assert!((q[0] - 1.0).abs() < 1e-12);
        let scores = vec![5.0, 4.0, 3.0, 2.0, 1.0, 0.0];
        let labels = vec![1i8, 1, 1, 1, 1, -1];
        let q = qvalues(&scores, &labels, reported());
        // Best set: all five targets with the safeguard decoy -> 1/5.
        assert!((q[0] - 0.2).abs() < 1e-12, "q[0] = {}", q[0]);
    }

    /// Every member of an exact tie must receive the same q-value, whatever order
    /// the rows arrive in.
    #[test]
    fn tied_scores_receive_identical_qvalues_in_any_row_order() {
        let scores = vec![2.0, 2.0, 2.0, 2.0, 1.0, 1.0];
        let labels = vec![1i8, -1, 1, -1, 1, -1];
        let q = qvalues(&scores, &labels, reported());
        for i in 0..4 {
            assert!(
                (q[i] - q[0]).abs() < 1e-12,
                "tie member {i} has q = {}, expected {}",
                q[i],
                q[0]
            );
        }
        assert!((q[4] - q[5]).abs() < 1e-12);

        // Reversing the rows must not move any q-value.
        let reversed_scores: Vec<f64> = scores.iter().rev().copied().collect();
        let reversed_labels: Vec<i8> = labels.iter().rev().copied().collect();
        let reversed_q = qvalues(&reversed_scores, &reversed_labels, reported());
        for i in 0..scores.len() {
            assert!(
                (reversed_q[scores.len() - 1 - i] - q[i]).abs() < 1e-12,
                "row {i} q changed under permutation"
            );
        }
    }

    #[test]
    fn all_identical_scores_collapse_to_one_qvalue() {
        let scores = vec![1.0; 8];
        let labels = vec![1i8, -1, 1, -1, 1, -1, 1, -1];
        let q = qvalues(&scores, &labels, reported());
        assert!(q.iter().all(|&v| (v - q[0]).abs() < 1e-12));
        // 5 decoys (4 observed + safeguard) over 4 targets, capped at 1.
        assert!((q[0] - 1.0).abs() < 1e-12);
    }

    /// A target-heavy list must not have its FDP shrunk by a count ratio.  The
    /// old estimator multiplied by `D_total / T_total`, which is anti-conservative
    /// exactly when decoys are scarce.
    #[test]
    fn target_decoy_imbalance_is_not_rescaled_by_a_count_ratio() {
        let mut scores = Vec::new();
        let mut labels = Vec::new();
        for k in 0..100 {
            scores.push(-(k as f64));
            labels.push(1i8);
        }
        for k in 0..10 {
            scores.push(-50.0 - k as f64 * 0.5);
            labels.push(-1i8);
        }
        let q = qvalues(&scores, &labels, reported());
        // Worst rank: 10 observed decoys + safeguard over 100 targets = 0.11.
        let worst = q
            .iter()
            .copied()
            .fold(f64::NEG_INFINITY, |acc, value| acc.max(value));
        assert!(
            (worst - 0.11).abs() < 1e-12,
            "worst q = {worst}, expected 11/100"
        );
    }

    #[test]
    fn the_opportunity_ratio_scales_the_estimate() {
        let scores = vec![5.0, 4.0, 3.0, 2.0, 1.0, 0.0];
        let labels = vec![1i8, 1, 1, 1, 1, -1];
        let balanced = qvalues(&scores, &labels, Tdc::reported(0.5));
        // Two decoys per target: an incorrect target wins with probability 1/3,
        // so each decoy stands for half as many incorrect targets.
        let two_decoys = qvalues(&scores, &labels, Tdc::reported(1.0 / 3.0));
        assert!((two_decoys[0] - balanced[0] * 0.5).abs() < 1e-12);
    }

    #[test]
    fn dropping_the_safeguard_is_confined_to_the_training_estimator() {
        let scores = vec![5.0, 4.0, 3.0, 2.0];
        let labels = vec![1i8, 1, 1, -1];
        assert!(qvalues(&scores, &labels, Tdc::reported(0.5))[0] > 0.0);
        assert_eq!(qvalues(&scores, &labels, Tdc::training(0.5))[0], 0.0);
    }

    #[test]
    fn empty_and_single_row_inputs_are_handled() {
        assert!(qvalues(&[], &[], reported()).is_empty());
        assert_eq!(qvalues(&[1.0], &[1], reported()), vec![1.0]);
        assert_eq!(qvalues(&[1.0], &[-1], reported()), vec![1.0]);
    }

    #[test]
    fn repeated_calls_are_deterministic() {
        let scores = vec![3.0, 1.0, 3.0, 2.0, 2.0, 0.0, -1.0];
        let labels = vec![1i8, -1, -1, 1, 1, -1, 1];
        let first = qvalues_and_peps(&scores, &labels, reported());
        for _ in 0..8 {
            assert_eq!(qvalues_and_peps(&scores, &labels, reported()), first);
        }
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
            for tdc in [Tdc::reported(0.5), Tdc::training(0.5), Tdc::reported(0.25)] {
                let q = qvalues(&scores, &labels, tdc);
                for threshold in [0.0, 0.01, 0.5, 1.0, 1.01, f64::NAN] {
                    let expected = q
                        .iter()
                        .zip(&labels)
                        .filter(|(qvalue, label)| **label > 0 && **qvalue < threshold)
                        .count();
                    assert_eq!(
                        target_count_at_fdr(&scores, &labels, tdc, threshold),
                        expected,
                        "scores={scores:?}, labels={labels:?}, tdc={tdc:?}, threshold={threshold}"
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
            let fresh = target_count_at_fdr(&scores, &labels, reported(), 0.5);
            let reused = target_count_at_fdr_into(&scores, &labels, reported(), 0.5, &mut order);
            assert_eq!(reused, fresh);
        }
    }

    #[test]
    fn reversed_integer_ranks_match_the_negated_score_scan() {
        let mut state = 0x243f_6a88_85a3_08d3u64;
        for len in 0..200 {
            let mut scores = Vec::with_capacity(len);
            let mut labels = Vec::with_capacity(len);
            for index in 0..len {
                state = state
                    .wrapping_mul(6_364_136_223_846_793_005)
                    .wrapping_add(1_442_695_040_888_963_407);
                scores.push(((state >> 60) as i8 - 8) as f64);
                labels.push(if index % 3 == 0 { -1i8 } else { 1i8 });
            }
            let order = descending(&scores);
            let mut ranks = vec![0u32; len];
            let mut rank = 0u32;
            for position in 0..order.len() {
                if position > 0 && scores[order[position]] != scores[order[position - 1]] {
                    rank += 1;
                }
                ranks[order[position]] = rank;
            }
            let negated: Vec<f64> = scores.iter().map(|value| -value).collect();
            for threshold in [0.01, 0.1, 0.5, 1.0] {
                assert_eq!(
                    target_count_at_reversed_ranks_into(
                        &ranks,
                        &labels,
                        reported(),
                        threshold,
                        &mut Vec::new()
                    ),
                    target_count_at_fdr(&negated, &labels, reported(), threshold),
                    "len={len}, threshold={threshold}"
                );
            }
        }
    }

    #[test]
    fn qvalue_workspace_overwrites_reused_buffers() {
        let mut workspace = QValueWorkspace::default();
        let mut reused = Vec::new();
        for (scores, labels, tdc) in [
            (vec![3.0, 2.0, 1.0], vec![1, -1, 1], Tdc::reported(0.5)),
            (vec![5.0], vec![-1], Tdc::reported(0.25)),
            (Vec::new(), Vec::new(), Tdc::reported(0.5)),
            (
                vec![1.0, 4.0, 2.0, 3.0],
                vec![-1, 1, -1, 1],
                Tdc::training(0.5),
            ),
        ] {
            let fresh = qvalues(&scores, &labels, tdc);
            qvalues_into(&scores, &labels, tdc, &mut workspace, &mut reused);
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
                let q = qvalues(&scores, &labels, reported());
                target_mask_at_fdr_into(
                    &scores,
                    &labels,
                    reported(),
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

    // ----- PEP invariants -----------------------------------------------------

    fn mixed_case() -> (Vec<f64>, Vec<i8>) {
        let mut scores = Vec::new();
        let mut labels = Vec::new();
        let mut state = 0x9E37_79B9_7F4A_7C15u64;
        for k in 0..400 {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            let noise = ((state >> 40) as f64 / (1u64 << 24) as f64) - 0.5;
            let target = k % 2 == 0;
            // Targets carry signal in the leading part of the list.
            scores.push(if target {
                100.0 - k as f64 * 0.2 + noise
            } else {
                40.0 - k as f64 * 0.2 + noise
            });
            labels.push(if target { 1i8 } else { -1i8 });
        }
        (scores, labels)
    }

    #[test]
    fn peps_are_bounded_and_strictly_positive() {
        let (scores, labels) = mixed_case();
        for value in peps(&scores, &labels, reported()) {
            assert!(
                value > 0.0 && value <= 1.0,
                "PEP out of the open-lower unit range: {value}"
            );
        }
    }

    /// The old estimator printed exactly zero for every PSM above the first decoy.
    #[test]
    fn a_leading_target_run_never_receives_zero_pep() {
        let mut scores = Vec::new();
        let mut labels = Vec::new();
        for k in 0..500 {
            scores.push(1000.0 - k as f64);
            labels.push(1i8);
        }
        for k in 0..500 {
            scores.push(100.0 - k as f64);
            labels.push(-1i8);
        }
        let pep = peps(&scores, &labels, reported());
        assert!(pep.iter().all(|&value| value > 0.0));
        // Prior mass alone puts a floor of half a false discovery on the list.
        assert!(pep[0] >= 0.5 / 500.0 - 1e-12, "pep[0] = {}", pep[0]);
    }

    #[test]
    fn peps_are_monotone_along_worsening_score() {
        let (scores, labels) = mixed_case();
        let pep = peps(&scores, &labels, reported());
        let mut previous = 0.0;
        for &i in &descending(&scores) {
            assert!(
                pep[i] + 1e-12 >= previous,
                "PEP decreased along worsening score"
            );
            previous = pep[i];
        }
    }

    /// The defining relation: the mean PEP over the top `k` targets is the
    /// q-value at `k`.  Isotonic smoothing and the pseudocount perturb it, so the
    /// test allows the pseudocount's own magnitude plus a small tolerance.
    #[test]
    fn mean_pep_tracks_the_qvalue_it_was_derived_from() {
        let (scores, labels) = mixed_case();
        let (q, pep) = qvalues_and_peps(&scores, &labels, reported());
        let order = descending(&scores);
        let targets: Vec<usize> = order
            .iter()
            .copied()
            .filter(|&i| labels[i] > 0)
            .collect();
        let pseudocount = 0.5 / targets.len() as f64;
        let mut running = 0.0;
        for (k, &i) in targets.iter().enumerate() {
            running += pep[i];
            let mean = running / (k + 1) as f64;
            assert!(
                (mean - q[i]).abs() <= pseudocount + 0.02,
                "k={}, mean PEP={mean}, q={}",
                k + 1,
                q[i]
            );
        }
    }

    #[test]
    fn tied_scores_receive_identical_peps() {
        let scores = vec![2.0, 2.0, 2.0, 2.0, 1.0, 1.0, 0.0, 0.0];
        let labels = vec![1i8, -1, 1, -1, 1, -1, 1, -1];
        let pep = peps(&scores, &labels, reported());
        assert!((pep[0] - pep[2]).abs() < 1e-12);
        assert!((pep[4] - pep[4]).abs() < 1e-12);
        assert!((pep[6] - pep[6]).abs() < 1e-12);
    }

    #[test]
    fn peps_are_invariant_to_input_row_order() {
        let (scores, labels) = mixed_case();
        let pep = peps(&scores, &labels, reported());
        let reversed_scores: Vec<f64> = scores.iter().rev().copied().collect();
        let reversed_labels: Vec<i8> = labels.iter().rev().copied().collect();
        let reversed_pep = peps(&reversed_scores, &reversed_labels, reported());
        for i in 0..scores.len() {
            assert!(
                (reversed_pep[scores.len() - 1 - i] - pep[i]).abs() < 1e-12,
                "row {i} PEP changed under permutation"
            );
        }
    }

    #[test]
    fn a_decoy_only_list_reports_no_confidence() {
        let scores = vec![3.0, 2.0, 1.0];
        let labels = vec![-1i8, -1, -1];
        let pep = peps(&scores, &labels, reported());
        assert!(pep.iter().all(|&value| (value - 1.0).abs() < 1e-12));
    }

    #[test]
    fn peps_handle_empty_and_single_row_inputs() {
        assert!(peps(&[], &[], reported()).is_empty());
        assert_eq!(peps(&[1.0], &[1], reported()), vec![1.0]);
    }
}
