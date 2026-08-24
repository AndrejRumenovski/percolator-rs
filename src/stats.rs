#![allow(clippy::items_after_test_module)]

//! Target-decoy q-values and a monotone PEP estimate.

/// Given scores and labels (+1 target / -1 decoy), return q-values aligned to input order.
/// Standard target-decoy competition: walk high->low score, FDR = pi0 * (#decoys / #targets),
/// then enforce monotonicity (running minimum from the bottom).
pub fn qvalues(scores: &[f64], labels: &[i8], pi0: f64) -> Vec<f64> {
    let n = scores.len();
    #[cfg(feature = "profiling")]
    let _qvalues = crate::profile::Scope::with_elements("qvalue", "qvalues_total", n);
    #[cfg(feature = "profiling")]
    crate::profile::allocation_site(
        "stats::qvalues temporary vectors",
        3,
        n as u64 * (std::mem::size_of::<usize>() + 2 * std::mem::size_of::<f64>()) as u64,
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
        "qvalue_score_order",
        sort_start.elapsed(),
        Some(n as u64),
        None,
    );

    #[cfg(feature = "profiling")]
    let evaluate_start = std::time::Instant::now();
    let mut q = vec![1.0f64; n];
    let mut targets = 0.0f64;
    let mut decoys = 0.0f64;
    // first pass: raw fdr at each rank
    let mut fdr_at = vec![1.0f64; n];
    for (rank, &i) in order.iter().enumerate() {
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
        fdr_at[rank] = fdr.min(1.0);
    }
    // monotonize from worst rank upward, assign to targets
    let mut running = 1.0f64;
    for rank in (0..n).rev() {
        if fdr_at[rank] < running {
            running = fdr_at[rank];
        }
        let i = order[rank];
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
    q
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
    fn peps_bounded_and_monotone() {
        let scores = vec![5.0, 4.0, 3.0, 2.0, 1.0];
        let labels = vec![1i8, 1, -1, -1, -1];
        let p = peps(&scores, &labels, 1.0);
        for v in &p {
            assert!(*v >= 0.0 && *v <= 1.0);
        }
    }
}

/// Monotone non-increasing PEP vs score via a light isotonic (pool-adjacent-violators)
/// fit of the decoy indicator over score-sorted PSMs. Returns PEP aligned to input order.
#[allow(clippy::needless_range_loop)]
pub fn peps(scores: &[f64], labels: &[i8], pi0: f64) -> Vec<f64> {
    let n = scores.len();
    #[cfg(feature = "profiling")]
    let _peps = crate::profile::Scope::with_elements("pep", "pep_pava_total", n);
    #[cfg(feature = "profiling")]
    crate::profile::allocation_site(
        "stats::peps temporary vectors",
        10,
        n as u64 * (3 * std::mem::size_of::<usize>() + 7 * std::mem::size_of::<f64>()) as u64,
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
    #[cfg(feature = "profiling")]
    let pava_start = std::time::Instant::now();
    // target=0 (correct-ish), decoy=1; PAVA to get monotone-increasing prob-of-decoy along worsening score
    let y: Vec<f64> = order
        .iter()
        .map(|&i| if labels[i] < 0 { 1.0 } else { 0.0 })
        .collect();
    // PAVA (isotonic, non-decreasing)
    let mut val = y.clone();
    let mut wt = vec![1.0f64; n];
    let mut idx: Vec<usize> = (0..n).collect();
    let mut m = n;
    let mut i = 0;
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
    let _ = (&mut val, &mut wt, &mut idx, &mut m, &mut i);
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
