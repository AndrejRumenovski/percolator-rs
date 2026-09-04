//! Stage-resolved probe of the reported PEP estimator.
//!
//! Two jobs.  `--selftest` runs hand-computable fixtures against an oracle that
//! is written independently from the module documentation, checking every
//! transformation in isolation.  `--score-file` runs the production estimator
//! over an arbitrary (score,label) list so a simulation harness outside Rust can
//! drive the real code path instead of a reimplementation of it.

#[path = "../src/stats.rs"]
mod stats;

use std::io::{Read, Write};
use stats::{qvalues_and_peps, Tdc};

const EPS: f64 = 1e-11;

fn close(a: f64, b: f64) -> bool {
    (a - b).abs() <= EPS * (1.0 + a.abs().max(b.abs()))
}

fn check(name: &str, actual: &[f64], expected: &[f64]) {
    assert_eq!(actual.len(), expected.len(), "{name}: length");
    for (i, (&a, &e)) in actual.iter().zip(expected).enumerate() {
        assert!(close(a, e), "{name}[{i}]: actual={a:.17e} expected={e:.17e}");
    }
    println!("  ok  {name}");
}

/// Stage 1 of the oracle: per-tie-group increments of
/// `F(g) = min(T_g, pi0*lambda*(D_g+1), assigned + n_g)`, shared inside a group.
/// Returned in best-first target order, before any isotonic step.
fn oracle_increments(scores: &[f64], labels: &[i8], p: f64) -> Vec<f64> {
    let lambda = p / (1.0 - p);
    let mut order: Vec<usize> = (0..scores.len()).collect();
    order.sort_by(|&l, &r| scores[r].total_cmp(&scores[l]));
    let mut out: Vec<f64> = Vec::new();
    let (mut t, mut d, mut assigned) = (0.0f64, 1.0f64, 0.0f64);
    let mut start = 0usize;
    while start < order.len() {
        let s = scores[order[start]];
        let mut end = start;
        while end < order.len() && scores[order[end]] == s {
            end += 1;
        }
        let mut n_g = 0usize;
        for &row in &order[start..end] {
            if labels[row] > 0 {
                t += 1.0;
                n_g += 1;
            } else {
                d += 1.0;
            }
        }
        if n_g > 0 {
            let f = (1.0 * lambda * d).min(t).min(assigned + n_g as f64);
            let share = ((f - assigned) / n_g as f64).max(0.0);
            for _ in 0..n_g {
                out.push(share);
            }
            assigned = f;
        }
        start = end;
    }
    out
}

/// Stage 2 of the oracle: unit-weight isotonic (non-decreasing) regression,
/// written as the slopes of the greatest convex minorant of the cumulative sum,
/// which is a different algorithm from the production pool-adjacent-violators
/// loop and therefore an independent check of it.
fn oracle_gcm_slopes(values: &[f64]) -> Vec<f64> {
    let n = values.len();
    let mut cum = vec![0.0f64; n + 1];
    for i in 0..n {
        cum[i + 1] = cum[i] + values[i];
    }
    // Lower convex hull of the points (i, cum[i]).
    let mut hull: Vec<usize> = Vec::new();
    for i in 0..=n {
        while hull.len() >= 2 {
            let a = hull[hull.len() - 2];
            let b = hull[hull.len() - 1];
            let left = (cum[b] - cum[a]) * ((i - b) as f64);
            let right = (cum[i] - cum[b]) * ((b - a) as f64);
            if left >= right {
                hull.pop();
            } else {
                break;
            }
        }
        hull.push(i);
    }
    let mut out = vec![0.0f64; n];
    for w in hull.windows(2) {
        let (a, b) = (w[0], w[1]);
        let slope = (cum[b] - cum[a]) / ((b - a) as f64);
        for slot in out.iter_mut().take(b).skip(a) {
            *slot = slope;
        }
    }
    out
}

/// Full oracle for reported target PEPs, in best-first target order.
fn oracle_peps(scores: &[f64], labels: &[i8], p: f64) -> Vec<f64> {
    let mut v = oracle_gcm_slopes(&oracle_increments(scores, labels, p));
    for x in &mut v {
        *x = x.clamp(1e-12, 1.0);
    }
    v
}

fn target_peps_best_first(scores: &[f64], labels: &[i8], p: f64) -> Vec<f64> {
    let (_, pep) = qvalues_and_peps(scores, labels, Tdc::reported(p));
    let mut order: Vec<usize> = (0..scores.len()).collect();
    order.sort_by(|&l, &r| scores[r].total_cmp(&scores[l]));
    order
        .into_iter()
        .filter(|&row| labels[row] > 0)
        .map(|row| pep[row])
        .collect()
}

fn fixture(name: &str, scores: &[f64], labels: &[i8], p: f64) {
    let actual = target_peps_best_first(scores, labels, p);
    let expected = oracle_peps(scores, labels, p);
    check(name, &actual, &expected);
}

fn selftest() {
    println!("== stage agreement: production estimator vs independent GCM oracle ==");

    // Complete null: alternating targets and decoys, no signal anywhere.
    let n = 400;
    let scores: Vec<f64> = (0..n).map(|i| -(i as f64)).collect();
    let labels: Vec<i8> = (0..n).map(|i| if i % 2 == 0 { 1 } else { -1 }).collect();
    fixture("complete null, alternating", &scores, &labels, 0.5);

    // Perfect separation: every target above every decoy.
    let scores: Vec<f64> = (0..200).map(|i| 200.0 - i as f64).collect();
    let labels: Vec<i8> = (0..200).map(|i| if i < 100 { 1 } else { -1 }).collect();
    fixture("perfect separation", &scores, &labels, 0.5);

    // Overlapping: deterministic interleave with a target-heavy head.
    let mut scores = Vec::new();
    let mut labels = Vec::new();
    for i in 0..300 {
        scores.push(300.0 - i as f64);
        labels.push(if (i * 7) % 11 < 6 { 1 } else { -1 });
    }
    fixture("overlapping distributions", &scores, &labels, 0.5);

    // Exact ties across the whole list.
    let scores = vec![1.0f64; 50];
    let labels: Vec<i8> = (0..50).map(|i| if i % 3 == 0 { -1 } else { 1 }).collect();
    fixture("all scores tied", &scores, &labels, 0.5);

    // Repeated scores in blocks, mixed labels inside a block.
    let scores: Vec<f64> = (0..120).map(|i| (120 - i) as f64 / 7.0).collect();
    let labels: Vec<i8> = (0..120).map(|i| if (i / 3) % 2 == 0 { 1 } else { -1 }).collect();
    fixture("repeated score blocks", &scores, &labels, 0.5);

    // Small samples.
    for n in 1..=8usize {
        let scores: Vec<f64> = (0..n).map(|i| -(i as f64)).collect();
        let labels: Vec<i8> = (0..n).map(|i| if i % 2 == 0 { 1 } else { -1 }).collect();
        fixture(&format!("small sample n={n}"), &scores, &labels, 0.5);
    }

    // Severe class imbalance both ways.
    let scores: Vec<f64> = (0..1000).map(|i| -(i as f64)).collect();
    let labels: Vec<i8> = (0..1000).map(|i| if i % 100 == 0 { -1 } else { 1 }).collect();
    fixture("imbalance: 1% decoys", &scores, &labels, 0.5);
    let labels: Vec<i8> = (0..1000).map(|i| if i % 100 == 0 { 1 } else { -1 }).collect();
    fixture("imbalance: 1% targets", &scores, &labels, 0.5);

    // Boundary and degenerate score values.
    let scores = vec![0.0, -0.0, f64::MIN_POSITIVE, -f64::MIN_POSITIVE, 1e308, -1e308];
    let labels = vec![1i8, -1, 1, -1, 1, -1];
    fixture("boundary score values", &scores, &labels, 0.5);

    // Non-0.5 opportunity ratios.
    for p in [0.1f64, 0.25, 0.5, 0.75, 0.9] {
        let scores: Vec<f64> = (0..200).map(|i| -(i as f64)).collect();
        let labels: Vec<i8> = (0..200).map(|i| if i % 4 == 0 { -1 } else { 1 }).collect();
        fixture(&format!("opportunity ratio p={p}"), &scores, &labels, p);
    }

    println!("== hand-computable closed forms ==");

    // Five targets above one decoy.  F is 1 at the first group and stays 1, so
    // the whole estimate is the safeguard decoy spread over the leading block.
    let scores = vec![5.0, 4.0, 3.0, 2.0, 1.0, 0.0];
    let labels = vec![1i8, 1, 1, 1, 1, -1];
    check(
        "five targets above one decoy",
        &target_peps_best_first(&scores, &labels, 0.5),
        &[0.2; 5],
    );

    // One target above one decoy: safeguard puts a whole false discovery on it.
    check(
        "single target above single decoy",
        &target_peps_best_first(&[1.0, 0.0], &[1i8, -1], 0.5),
        &[1.0],
    );

    // One safeguard decoy plus one observed decoy above the fourth target.
    // F = 1 at the first target group and 2 at the third, then flat; the
    // trailing decoy is below every target and is never assigned.  Raw
    // increments 1,1,0,0,0 are pooled by the isotonic step into a flat 0.4.
    let scores = vec![5.0, 4.5, 4.0, 3.0, 2.0, 1.0, 0.0];
    let labels = vec![1i8, -1, 1, 1, 1, 1, -1];
    let peps = target_peps_best_first(&scores, &labels, 0.5);
    let total: f64 = peps.iter().sum();
    assert!(close(total, 2.0), "mass {total}");
    check("two false discoveries over five targets", &peps, &[0.4; 5]);

    // The trailing-decoy exclusion, stated on its own: moving that last decoy
    // above the last target raises the mass from 2 to 3.  Raw increments are
    // 1,1,0,0,1; the isotonic step pools only the first four.
    let scores = vec![5.0, 4.5, 4.0, 3.0, 2.0, 1.5, 1.0];
    let labels = vec![1i8, -1, 1, 1, 1, -1, 1];
    let peps = target_peps_best_first(&scores, &labels, 0.5);
    let total: f64 = peps.iter().sum();
    assert!(close(total, 3.0), "mass {total}");
    check("decoy above the last target", &peps, &[0.5, 0.5, 0.5, 0.5, 1.0]);

    println!("== invariants ==");
    let scores: Vec<f64> = (0..5000).map(|i| ((i * 2654435761u64 as usize) % 100003) as f64).collect();
    let labels: Vec<i8> = (0..5000).map(|i| if (i * 37) % 5 == 0 { -1 } else { 1 }).collect();
    let (_, pep) = qvalues_and_peps(&scores, &labels, Tdc::reported(0.5));
    let best_first = target_peps_best_first(&scores, &labels, 0.5);
    assert!(pep.iter().all(|v| v.is_finite() && (0.0..=1.0).contains(v)));
    assert!(best_first.windows(2).all(|w| w[0] <= w[1] + 1e-15), "monotone");
    let mass: f64 = best_first.iter().sum();
    println!("  ok  finite, in [0,1], monotone; target PEP mass = {mass:.6}");

    println!("== reported decoy PEP is not a function of the input ==");
    // Same six rows, only the order of the last two (which share a score) differs.
    // Every target PEP is invariant; the decoy's reported value is not, because the
    // decoy back-fill walks `order` row by row instead of once per tie group.
    let a_scores = vec![5.0, 4.0, 3.0, 2.0, 1.0, 1.0];
    let a_labels = vec![1i8, 1, 1, 1, 1, -1];
    let b_labels = vec![1i8, 1, 1, 1, -1, 1];
    let (_, pa) = qvalues_and_peps(&a_scores, &a_labels, Tdc::reported(0.5));
    let (_, pb) = qvalues_and_peps(&a_scores, &b_labels, Tdc::reported(0.5));
    let decoy_a = pa[5];
    let decoy_b = pb[4];
    let targets_a = target_peps_best_first(&a_scores, &a_labels, 0.5);
    let targets_b = target_peps_best_first(&a_scores, &b_labels, 0.5);
    check("target PEPs are order invariant", &targets_a, &targets_b);
    println!("  !!  decoy PEP: {decoy_a:.6} when listed after the tied target, {decoy_b:.6} when listed before");
    assert!(
        (decoy_a - decoy_b).abs() > 0.5,
        "expected the documented decoy back-fill order dependence"
    );

    println!("SELFTEST PASS");
}

fn run_score_file(path: &str) {
    let mut text = String::new();
    std::fs::File::open(path)
        .expect("open score file")
        .read_to_string(&mut text)
        .expect("read score file");
    let mut scores = Vec::new();
    let mut labels = Vec::new();
    for line in text.lines() {
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        let mut it = line.split_whitespace();
        scores.push(it.next().unwrap().parse::<f64>().unwrap());
        labels.push(it.next().unwrap().parse::<i8>().unwrap());
    }
    let (q, pep) = qvalues_and_peps(&scores, &labels, Tdc::reported(0.5));
    let stdout = std::io::stdout();
    let mut out = std::io::BufWriter::new(stdout.lock());
    for i in 0..q.len() {
        writeln!(out, "{:.17e}\t{:.17e}", q[i], pep[i]).unwrap();
    }
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() >= 3 && args[1] == "--score-file" {
        run_score_file(&args[2]);
    } else {
        selftest();
    }
}
