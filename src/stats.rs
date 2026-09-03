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
//! A PEP is the *local* counterpart of the cumulative quantity above: the
//! probability that one particular target is incorrect.  Käll et al. (2008),
//! "Posterior error probabilities and false discovery rates: two sides of the
//! same coin", give the relation between them — the estimated number of false
//! discoveries among a set is the sum of that set's PEPs.  The implication runs
//! from calibrated PEPs to calibrated FDR, **not** the other way round, so a
//! valid cumulative estimate is a starting point for a local one and never a
//! proof of it.
//!
//! The same scan that produces `FDP` also produces the estimated number of
//! incorrect targets at or above a score,
//!
//! ```text
//! F(s) = min(T(s), pi0 * lambda * (D(s) + 1))
//! ```
//!
//! which is exactly `T(s) * FDP(s)`.  `F` is a non-decreasing step function of
//! the score threshold, and the PEP of a target is the increment of `F` it is
//! responsible for: walking the tie groups best-first, each group's increment
//! `F(g) - F(g-1)` is shared equally among the targets in that group.  Groups
//! holding only decoys carry their increment forward to the next group that has
//! targets, so no estimated false discovery is lost.
//!
//! Differencing a step function is high variance — this is why QVALITY fits a
//! smooth nonparametric model instead — so an isotonic (PAVA) fit then enforces
//! the monotonicity a posterior error probability must have in score.  PAVA
//! redistributes mass only inside a pooled block, so the sum over all targets is
//! preserved: **the reported PEPs sum to the reported estimated number of false
//! discoveries.**
//!
//! Two properties follow without any tuned constant.  The leading run of targets
//! above every decoy receives the finite-sample safeguard decoy's `lambda`
//! spread across it, so no finite input produces `PEP = 0`; and the declared
//! opportunity ratio scales PEPs exactly as it scales q-values.
//!
//! # What this estimator does not claim
//!
//! It inherits every assumption of the cumulative estimator above.  If `FDP` is
//! anti-conservative on some data then its increments are too, and no property
//! of the isotonic fit repairs that.  See `validation/` for the measured
//! calibration, which is **not** a validation of these values as posterior
//! probabilities.

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
    #[cfg(feature = "profiling")]
    let old_capacities = (
        workspace.order.capacity(),
        workspace.fdr_at.capacity(),
        q.capacity(),
    );
    // Fail closed.  A non-finite score has no place in a score ordering: it
    // cannot be compared, it cannot join a tie group, and admitting one lets a
    // single row silently move the q-values of every row above it.  The PIN
    // parser already rejects non-finite features, so reaching this means a model
    // diverged.
    assert!(
        scores.iter().all(|value| value.is_finite()),
        "target-decoy estimation requires finite scores; a non-finite score reached the estimator"
    );

    #[cfg(feature = "profiling")]
    let buffer_start = std::time::Instant::now();
    workspace.order.clear();
    workspace.order.extend(0..n);
    workspace.fdr_at.clear();
    workspace.fdr_at.resize(n, 1.0);
    q.clear();
    q.resize(n, 1.0);
    #[cfg(feature = "profiling")]
    crate::profile::record(
        "qvalue",
        "qvalue_buffer_setup",
        buffer_start.elapsed(),
        Some(n as u64),
        None,
    );
    #[cfg(feature = "profiling")]
    crate::profile::allocation_site(
        "stats::qvalue materialization buffers",
        u64::from(workspace.order.capacity() > old_capacities.0)
            + u64::from(workspace.fdr_at.capacity() > old_capacities.1)
            + u64::from(q.capacity() > old_capacities.2),
        (usize::from(workspace.order.capacity() > old_capacities.0)
            * workspace.order.capacity()
            * std::mem::size_of::<usize>()
            + usize::from(workspace.fdr_at.capacity() > old_capacities.1)
                * workspace.fdr_at.capacity()
                * std::mem::size_of::<f64>()
            + usize::from(q.capacity() > old_capacities.2)
                * q.capacity()
                * std::mem::size_of::<f64>()) as u64,
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
    let tie_scan_start = std::time::Instant::now();
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
    #[cfg(feature = "profiling")]
    crate::profile::record(
        "qvalue",
        "qvalue_tie_group_scan",
        tie_scan_start.elapsed(),
        Some(n as u64),
        None,
    );
    // Reverse cumulative minimum: q is the best achievable FDP of any set that
    // still contains this PSM.
    #[cfg(feature = "profiling")]
    let monotonic_start = std::time::Instant::now();
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
        "qvalue_monotonic_scan",
        monotonic_start.elapsed(),
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
    #[cfg(feature = "profiling")]
    let old_capacities = (workspace.order.capacity(), accepted.capacity());
    assert!(
        scores.iter().all(|value| value.is_finite()),
        "target-decoy estimation requires finite scores; a non-finite score reached the estimator"
    );
    #[cfg(feature = "profiling")]
    let buffer_start = std::time::Instant::now();
    workspace.order.clear();
    workspace.order.extend(0..n);
    accepted.clear();
    accepted.resize(n, 0);
    #[cfg(feature = "profiling")]
    crate::profile::record(
        "qvalue",
        "qvalue_buffer_setup",
        buffer_start.elapsed(),
        Some(n as u64),
        None,
    );
    #[cfg(feature = "profiling")]
    crate::profile::allocation_site(
        "stats::training qvalue buffers",
        u64::from(workspace.order.capacity() > old_capacities.0)
            + u64::from(accepted.capacity() > old_capacities.1),
        (usize::from(workspace.order.capacity() > old_capacities.0)
            * workspace.order.capacity()
            * std::mem::size_of::<usize>()
            + usize::from(accepted.capacity() > old_capacities.1) * accepted.capacity())
            as u64,
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
    let tie_scan_start = std::time::Instant::now();
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
    #[cfg(feature = "profiling")]
    crate::profile::record(
        "qvalue",
        "qvalue_tie_group_scan",
        tie_scan_start.elapsed(),
        Some(n as u64),
        None,
    );
    #[cfg(feature = "profiling")]
    let selection_start = std::time::Instant::now();
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
        "qvalue_positive_materialization",
        selection_start.elapsed(),
        Some(n as u64),
        None,
    );
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
    let n = scores.len();
    #[cfg(feature = "profiling")]
    let _qvalues = crate::profile::Scope::with_elements("qvalue", "qvalues_total", n);
    #[cfg(feature = "profiling")]
    let old_order_capacity = order.capacity();
    assert!(
        scores.iter().all(|value| value.is_finite()),
        "target-decoy estimation requires finite scores; a non-finite score reached the estimator"
    );
    #[cfg(feature = "profiling")]
    let buffer_start = std::time::Instant::now();
    order.clear();
    order.extend(0..n);
    #[cfg(feature = "profiling")]
    crate::profile::record(
        "qvalue",
        "qvalue_buffer_setup",
        buffer_start.elapsed(),
        Some(n as u64),
        None,
    );
    #[cfg(feature = "profiling")]
    if order.capacity() > old_order_capacity {
        crate::profile::allocation_site(
            "stats::initial-direction sort buffer",
            1,
            (order.capacity() * std::mem::size_of::<usize>()) as u64,
        );
    }
    #[cfg(feature = "profiling")]
    let sort_start = std::time::Instant::now();
    sort_score_order(order, scores);
    #[cfg(feature = "profiling")]
    crate::profile::record(
        "sort",
        "qvalue_score_order",
        sort_start.elapsed(),
        Some(n as u64),
        None,
    );

    #[cfg(feature = "profiling")]
    let tie_scan_start = std::time::Instant::now();
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
    #[cfg(feature = "profiling")]
    crate::profile::record(
        "qvalue",
        "qvalue_tie_group_scan",
        tie_scan_start.elapsed(),
        Some(n as u64),
        None,
    );
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
    let n = ranks.len();
    #[cfg(feature = "profiling")]
    let _qvalues = crate::profile::Scope::with_elements("qvalue", "qvalues_total", n);
    #[cfg(feature = "profiling")]
    let old_order_capacity = order.capacity();
    #[cfg(feature = "profiling")]
    let buffer_start = std::time::Instant::now();
    order.clear();
    order.extend(0..n);
    #[cfg(feature = "profiling")]
    crate::profile::record(
        "qvalue",
        "qvalue_buffer_setup",
        buffer_start.elapsed(),
        Some(n as u64),
        None,
    );
    #[cfg(feature = "profiling")]
    if order.capacity() > old_order_capacity {
        crate::profile::allocation_site(
            "stats::initial-direction rank-sort buffer",
            1,
            (order.capacity() * std::mem::size_of::<usize>()) as u64,
        );
    }
    #[cfg(feature = "profiling")]
    let sort_start = std::time::Instant::now();
    sort_reversed_rank_order(order, ranks);
    #[cfg(feature = "profiling")]
    crate::profile::record(
        "sort",
        "qvalue_reversed_rank_order",
        sort_start.elapsed(),
        Some(n as u64),
        None,
    );

    #[cfg(feature = "profiling")]
    let tie_scan_start = std::time::Instant::now();
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
    #[cfg(feature = "profiling")]
    crate::profile::record(
        "qvalue",
        "qvalue_tie_group_scan",
        tie_scan_start.elapsed(),
        Some(n as u64),
        None,
    );
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
    peps_from_competition_into(scores, labels, tdc, &mut workspace, &mut pep);
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
/// underflow of a genuinely tiny increment; it is never the value that removes
/// an otherwise exact zero, because the finite-sample safeguard decoy already
/// keeps the leading run strictly positive.
const PEP_FLOOR: f64 = 1e-12;

/// Posterior error probabilities from the same target-decoy scan that produced
/// the q-values, aligned to input order.
///
/// `workspace.order` must hold the descending score order of `scores`.
fn peps_from_competition_into(
    scores: &[f64],
    labels: &[i8],
    tdc: Tdc,
    workspace: &mut QValueWorkspace,
    pep: &mut Vec<f64>,
) {
    let n = workspace.order.len();
    #[cfg(feature = "profiling")]
    let _peps = crate::profile::Scope::with_elements("pep", "pep_from_competition_total", n);
    #[cfg(feature = "profiling")]
    let old_capacities = (
        pep.capacity(),
        workspace.target_ranks.capacity(),
        workspace.target_pep.capacity(),
        workspace.pava_value.capacity(),
        workspace.pava_weight.capacity(),
        workspace.pava_len.capacity(),
    );
    #[cfg(feature = "profiling")]
    let pep_start = std::time::Instant::now();
    pep.clear();
    pep.resize(n, 1.0);
    workspace.target_ranks.clear();
    workspace.target_pep.clear();
    if n == 0 {
        return;
    }

    // Walk the tie groups best-first, accumulating the estimated number of
    // incorrect targets `F` and handing each group's increment to the targets it
    // contains.  A decoy-only group leaves `F` raised and `assigned` untouched,
    // so its increment reaches the next group that has a target to carry it.
    let lambda = tdc.decoy_factor();
    let mut targets = 0.0f64;
    let mut decoys = tdc.initial_decoys();
    let mut assigned = 0.0f64;
    let mut group_first_target = 0usize;
    for rank in 0..n {
        let row = workspace.order[rank];
        if labels[row] > 0 {
            targets += 1.0;
            workspace.target_ranks.push(rank);
            workspace.target_pep.push(0.0);
        } else {
            decoys += 1.0;
        }
        if ends_score_group(&workspace.order, scores, rank) {
            let group_targets = workspace.target_pep.len() - group_first_target;
            if group_targets > 0 {
                // A probability cannot exceed one, so a group of `g` targets can
                // absorb at most `g` further false discoveries.  Applying that
                // bound here rather than clamping the finished curve is what
                // keeps the reported PEPs summing to the reported estimate; the
                // bound only ever binds deep in the tail, where every PEP has
                // already saturated at 1.
                let estimated_false = (tdc.pi0 * lambda * decoys)
                    .min(targets)
                    .min(assigned + group_targets as f64);
                let share = ((estimated_false - assigned) / group_targets as f64).max(0.0);
                for slot in &mut workspace.target_pep[group_first_target..] {
                    *slot = share;
                }
                assigned = estimated_false;
                group_first_target = workspace.target_pep.len();
            }
        }
    }
    let target_count = workspace.target_ranks.len();
    if target_count == 0 {
        #[cfg(feature = "profiling")]
        crate::profile::record(
            "pep",
            "pep_from_competition",
            pep_start.elapsed(),
            Some(n as u64),
            None,
        );
        return;
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

    // Decoys carry no error-rate claim; the column exists so a decoy row can be
    // placed on the same monotone curve as the targets around it.  Each decoy
    // takes the value of the nearest target at or above it, and the leading
    // decoys take the first target's value.
    let mut current = workspace.target_pep[0];
    let mut next_target = 0usize;
    for rank in 0..n {
        let row = workspace.order[rank];
        if labels[row] > 0 {
            current = workspace.target_pep[next_target];
            next_target += 1;
        } else {
            pep[row] = current.clamp(PEP_FLOOR, 1.0);
        }
    }
    #[cfg(feature = "profiling")]
    crate::profile::record(
        "pep",
        "pep_from_competition",
        pep_start.elapsed(),
        Some(n as u64),
        None,
    );
    #[cfg(feature = "profiling")]
    {
        let capacities = [
            pep.capacity(),
            workspace.target_ranks.capacity(),
            workspace.target_pep.capacity(),
            workspace.pava_value.capacity(),
            workspace.pava_weight.capacity(),
            workspace.pava_len.capacity(),
        ];
        let old = [
            old_capacities.0,
            old_capacities.1,
            old_capacities.2,
            old_capacities.3,
            old_capacities.4,
            old_capacities.5,
        ];
        let element_sizes = [
            std::mem::size_of::<f64>(),
            std::mem::size_of::<usize>(),
            std::mem::size_of::<f64>(),
            std::mem::size_of::<f64>(),
            std::mem::size_of::<f64>(),
            std::mem::size_of::<usize>(),
        ];
        crate::profile::allocation_site(
            "stats::PEP materialization/PAVA buffers",
            capacities
                .iter()
                .zip(old)
                .filter(|(new, old)| **new > *old)
                .count() as u64,
            capacities
                .iter()
                .zip(old)
                .zip(element_sizes)
                .filter(|((new, old), _)| **new > *old)
                .map(|((capacity, _), size)| capacity * size)
                .sum::<usize>() as u64,
        );
    }
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

    /// **Frozen requirement.** Every row at the same score shares one rejection
    /// boundary: a threshold cannot admit some members of a tie and not others.
    ///
    /// Three targets and one decoy at score 5, one decoy at score 1. The tie
    /// group is evaluated once, after all four rows have been counted:
    /// `1 * (1 + 1) / 3 = 2/3`. Ending the group after each row instead lets the
    /// leading targets be scored before their tied decoy has been counted, and
    /// the reverse cumulative minimum then hands them 1/3 -- half the estimate,
    /// from the same data.
    #[test]
    fn equal_scores_share_one_rejection_boundary() {
        let scores = vec![5.0, 5.0, 5.0, 5.0, 1.0];
        let labels = vec![1i8, 1, 1, -1, -1];
        let q = qvalues(&scores, &labels, reported());
        for (index, value) in q.iter().take(4).enumerate() {
            assert!(
                (value - 2.0 / 3.0).abs() < 1e-12,
                "row {index} q = {value}, expected 2/3"
            );
        }
        // Same requirement under every permutation of the tie group.
        for rotation in 0..4usize {
            let mut permuted_scores = scores.clone();
            let mut permuted_labels = labels.clone();
            permuted_labels[..4].rotate_left(rotation);
            permuted_scores[..4].rotate_left(rotation);
            let permuted = qvalues(&permuted_scores, &permuted_labels, reported());
            for value in permuted.iter().take(4) {
                assert!(
                    (value - 2.0 / 3.0).abs() < 1e-12,
                    "rotation {rotation}: q = {value}, expected 2/3"
                );
            }
        }
    }

    /// The same requirement for PEPs: tied targets must receive one shared
    /// increment of the estimated false-discovery count, not one each.
    #[test]
    fn equal_scores_share_one_pep_increment() {
        let scores = vec![5.0, 5.0, 5.0, 5.0, 1.0];
        let labels = vec![1i8, 1, 1, -1, -1];
        let pep = peps(&scores, &labels, reported());
        for (index, value) in pep.iter().take(3).enumerate() {
            assert!(
                (value - 2.0 / 3.0).abs() < 1e-12,
                "target {index} PEP = {value}, expected 2/3"
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

    /// **Frozen requirement (independent audit section 5).**  A non-finite score
    /// used to change the q-values of the finite rows above it.  The estimator
    /// now refuses the input instead.
    #[test]
    #[should_panic(expected = "requires finite scores")]
    fn a_non_finite_score_is_refused() {
        let _ = qvalues(&[3.0, 2.0, 1.0, f64::NAN], &[1, -1, 1, 1], reported());
    }

    #[test]
    #[should_panic(expected = "requires finite scores")]
    fn an_infinite_score_is_refused() {
        let _ = qvalues(&[3.0, f64::INFINITY], &[1, -1], reported());
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

    /// A longer list with repeated scores, interleaved labels and a decoy-free
    /// head, so the isotonic fit has real blocks to pool.
    fn long_case() -> (Vec<f64>, Vec<i8>) {
        let mut scores = Vec::new();
        let mut labels = Vec::new();
        let mut state = 0x1234_5678_9abc_def0u64;
        for k in 0..400 {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            let draw = ((state >> 33) % 1000) as f64 / 1000.0;
            // Decoys become common only in the tail.
            let label = if k < 40 || draw > 0.35 { 1i8 } else { -1 };
            scores.push(500.0 - k as f64 * 0.5 - (draw * 0.25));
            labels.push(label);
        }
        (scores, labels)
    }

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

    /// Estimated false discoveries among the top `k` targets, enumerated
    /// independently of the estimator from the raw counts: `lambda * (D + 1)`,
    /// never more than the number of targets accepted so far, and never rising
    /// by more than one per target because a probability cannot exceed one.
    fn estimated_false_discoveries(scores: &[f64], labels: &[i8], tdc: Tdc) -> Vec<f64> {
        let order = descending(scores);
        let lambda = tdc.null_target_win_prob / (1.0 - tdc.null_target_win_prob);
        let mut decoys = if tdc.skip_decoys_plus_one { 0.0 } else { 1.0 };
        let mut targets = 0.0f64;
        let mut carried = 0.0f64;
        let mut per_target = Vec::new();
        let mut rank = 0usize;
        while rank < order.len() {
            let mut end = rank;
            while !ends_score_group(&order, scores, end) {
                end += 1;
            }
            let mut group_targets = 0usize;
            for &row in &order[rank..=end] {
                if labels[row] > 0 {
                    targets += 1.0;
                    group_targets += 1;
                } else {
                    decoys += 1.0;
                }
            }
            if group_targets > 0 {
                carried = (lambda * decoys)
                    .min(targets)
                    .min(carried + group_targets as f64);
                for _ in 0..group_targets {
                    per_target.push(carried);
                }
            }
            rank = end + 1;
        }
        per_target
    }

    /// When the decoys outnumber what the targets can absorb, every PEP
    /// saturates at exactly 1 and the estimate stops growing.
    #[test]
    fn a_decoy_dominated_tail_saturates_at_one() {
        let mut scores = vec![10.0];
        let mut labels = vec![1i8];
        for k in 0..20 {
            scores.push(5.0 - k as f64);
            labels.push(-1);
        }
        scores.push(-100.0);
        labels.push(1);
        let pep = peps(&scores, &labels, reported());
        assert!(
            (pep[0] - 1.0).abs() < 1e-12,
            "leading target PEP {}",
            pep[0]
        );
        assert!(
            (pep[scores.len() - 1] - 1.0).abs() < 1e-12,
            "trailing target PEP {}",
            pep[scores.len() - 1]
        );
    }

    /// **Frozen requirement (Kall et al. 2008, and independent audit M1).**
    ///
    /// The PEPs of the top `k` targets must sum to the estimated number of false
    /// discoveries among them, which is exactly `k` times the raw estimated FDP
    /// at `k`.  This is an identity of the estimator, not an approximation: it
    /// holds to floating-point tolerance at every isotonic block boundary and can
    /// never be broken by adding prior mass on top of the finished curve.
    #[test]
    fn summed_peps_reproduce_the_estimated_false_discovery_count() {
        for (scores, labels) in [mixed_case(), long_case()] {
            let (_, pep) = qvalues_and_peps(&scores, &labels, reported());
            let expected = estimated_false_discoveries(&scores, &labels, reported());
            let order = descending(&scores);
            let targets: Vec<usize> = order.iter().copied().filter(|&i| labels[i] > 0).collect();
            assert_eq!(targets.len(), expected.len());
            // The total is exact; interior partial sums are exact wherever the
            // isotonic fit does not pool across the point.
            let total: f64 = targets.iter().map(|&i| pep[i]).sum();
            assert!(
                (total - expected[expected.len() - 1]).abs() < 1e-9,
                "summed PEP {total} != estimated false discoveries {}",
                expected[expected.len() - 1]
            );
            // A pooled isotonic block can only move mass within itself, so no
            // partial sum may ever exceed the raw count it is fitted to by more
            // than one block's worth of redistribution.
            let mut running = 0.0f64;
            for (k, &i) in targets.iter().enumerate() {
                running += pep[i];
                assert!(
                    running <= expected[k] + 1e-9 || running <= (k + 1) as f64,
                    "partial sum {running} above estimate {} at k={}",
                    expected[k],
                    k + 1
                );
            }
        }
    }

    /// **Frozen hand oracle.**  Five targets above one decoy, `p = 1/2`.
    ///
    /// The only estimated false discovery is the finite-sample safeguard decoy,
    /// so one false discovery is spread over the five targets that outrank
    /// everything: every PEP is exactly 0.2, which is also their q-value.  A
    /// pseudocount added after the fact would push all five to 0.3.
    #[test]
    fn the_leading_run_carries_exactly_the_safeguard_decoy() {
        let scores = vec![5.0, 4.0, 3.0, 2.0, 1.0, 0.0];
        let labels = vec![1i8, 1, 1, 1, 1, -1];
        let (q, pep) = qvalues_and_peps(&scores, &labels, reported());
        for k in 0..5 {
            assert!(
                (pep[k] - 0.2).abs() < 1e-12,
                "target {k} PEP {} is not 0.2",
                pep[k]
            );
            assert!((q[k] - 0.2).abs() < 1e-12);
        }
        let total: f64 = pep[..5].iter().sum();
        assert!((total - 1.0).abs() < 1e-12, "summed PEP {total} is not 1");
    }

    /// **Frozen hand oracle.**  Ten targets, one decoy, ten more targets.
    ///
    /// The safeguard decoy plus the observed one give two estimated false
    /// discoveries over twenty targets: the isotonic fit is flat at 0.1.
    #[test]
    fn a_second_decoy_adds_exactly_one_more_false_discovery() {
        let mut scores: Vec<f64> = Vec::new();
        let mut labels: Vec<i8> = Vec::new();
        for k in 0..10 {
            scores.push(100.0 - k as f64);
            labels.push(1);
        }
        scores.push(50.0);
        labels.push(-1);
        for k in 0..10 {
            scores.push(40.0 - k as f64);
            labels.push(1);
        }
        let pep = peps(&scores, &labels, reported());
        let total: f64 = labels
            .iter()
            .zip(&pep)
            .filter(|(&label, _)| label > 0)
            .map(|(_, &value)| value)
            .sum();
        assert!((total - 2.0).abs() < 1e-12, "summed PEP {total} is not 2");
        for (index, (&label, &value)) in labels.iter().zip(&pep).enumerate() {
            if label > 0 {
                assert!((value - 0.1).abs() < 1e-12, "target {index} PEP {value}");
            }
        }
    }

    /// The estimator must react to the declared opportunity ratio exactly as the
    /// q-value does, because it is differencing the same count.
    #[test]
    fn the_opportunity_ratio_scales_the_peps() {
        let scores = vec![5.0, 4.0, 3.0, 2.0, 1.0, 0.0];
        let labels = vec![1i8, 1, 1, 1, 1, -1];
        let half = peps(&scores, &labels, Tdc::reported(0.5));
        let third = peps(&scores, &labels, Tdc::reported(1.0 / 3.0));
        for k in 0..5 {
            assert!(
                (third[k] - half[k] * 0.5).abs() < 1e-12,
                "p=1/3 PEP {} is not half of p=1/2 PEP {}",
                third[k],
                half[k]
            );
        }
    }

    /// No estimator output may be a copy of a smoothing constant: doubling or
    /// removing prior mass must change the reported probabilities, and this test
    /// pins the exact values so a silent change cannot pass.
    #[test]
    fn leading_pep_is_pinned_by_the_estimator_not_by_a_constant() {
        // 20 targets ahead of every decoy: one safeguard false discovery spread
        // over 20 gives 0.05 exactly.
        let mut scores: Vec<f64> = (0..20).map(|k| 100.0 - k as f64).collect();
        let mut labels: Vec<i8> = vec![1; 20];
        scores.push(0.0);
        labels.push(-1);
        let pep = peps(&scores, &labels, reported());
        for (k, value) in pep.iter().enumerate().take(20) {
            assert!(
                (value - 0.05).abs() < 1e-12,
                "target {k} PEP {value} is not 0.05",
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
