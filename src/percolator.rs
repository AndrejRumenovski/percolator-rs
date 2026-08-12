//! Semi-supervised Percolator training: 3-fold nested cross-validation around an
//! iterative linear SVM that separates confident targets from decoys.

use crate::pin::Dataset;
use crate::stats;
use crate::svm::{train, Problem};

pub struct Params {
    pub maxiter: usize,       // semi-supervised iterations
    pub test_fdr: f64,        // FDR to pick positive training examples
    pub subset_max_train: usize, // 0 = use all
    pub seed: u64,
    pub max_newton: usize,
    /// Absolute SVM class weights (`C_pos`, `C_neg`), as in the reference.
    /// `None` selects them by cross-validation (the reference's `Cpos=0` behaviour).
    pub c_alpha: Option<f64>,
    pub c_beta: Option<f64>,
    /// Budget for each candidate during the C grid search (abbreviated training).
    pub c_select_maxiter: usize,
    pub c_select_subset: usize,
}

impl Default for Params {
    fn default() -> Self {
        Params {
            maxiter: 10,
            test_fdr: 0.01,
            subset_max_train: 0,
            seed: 1,
            max_newton: 30,
            c_alpha: None,
            c_beta: None,
            c_select_maxiter: 3,
            c_select_subset: 20_000,
        }
    }
}

/// Hyperparameters for one training run: class weights plus the training budget.
#[derive(Clone, Copy)]
struct Hp {
    alpha: f64,
    beta: f64,
    maxiter: usize,
    subset: usize,
}

/// Grid searched when the class weights are not pinned (absolute `C_pos` / `C_neg`).
///
/// The earlier fixed heuristic `C_pos = max(n_neg/n_pos, 1)`, `C_neg = 1` turns out to be
/// the worst corner of this space on real data: when confident targets are scarce the
/// ratio explodes (300x on sparse files) and swamps the decoys. Measured optima across
/// very different files cluster tightly at `C_pos ~ 1`, `C_neg ~ 4-16`, so the grid is
/// centred there.
const C_POS_GRID: [f64; 3] = [0.25, 1.0, 4.0];
const C_NEG_GRID: [f64; 3] = [1.0, 4.0, 16.0];

/// Fallback weights when the grid search is disabled (`--no-select-c`, `--fast`).
pub const C_POS_DEFAULT: f64 = 1.0;
pub const C_NEG_DEFAULT: f64 = 4.0;

struct Rng(u64);
impl Rng {
    fn next(&mut self) -> u64 {
        let mut x = self.0;
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        self.0 = x;
        x
    }
    fn below(&mut self, n: usize) -> usize {
        (self.next() % n as u64) as usize
    }
}

/// Build the normalized design matrix with an appended bias column (=1.0).
/// Returns (matrix row-major, dim) where dim = n_feat + 1.
fn build_matrix(ds: &Dataset) -> (Vec<f64>, usize) {
    let nf = ds.n_feat;
    let n = ds.n_psm;
    let mut mean = vec![0.0f64; nf];
    let mut var = vec![0.0f64; nf];
    for i in 0..n {
        let row = ds.row(i);
        for j in 0..nf {
            mean[j] += row[j];
        }
    }
    for m in mean.iter_mut() {
        *m /= n as f64;
    }
    for i in 0..n {
        let row = ds.row(i);
        for j in 0..nf {
            let d = row[j] - mean[j];
            var[j] += d * d;
        }
    }
    let mut std = vec![1.0f64; nf];
    for j in 0..nf {
        let v = (var[j] / n as f64).sqrt();
        std[j] = if v > 1e-12 { v } else { 1.0 };
    }
    let dim = nf + 1;
    let mut x = vec![0.0f64; n * dim];
    for i in 0..n {
        let row = ds.row(i);
        let base = i * dim;
        for j in 0..nf {
            x[base + j] = (row[j] - mean[j]) / std[j];
        }
        x[base + nf] = 1.0; // bias
    }
    (x, dim)
}

fn score_all(x: &[f64], dim: usize, w: &[f64], rows: &[usize], out: &mut [f64]) {
    for (k, &r) in rows.iter().enumerate() {
        out[k] = crate::simd::dot(&w[..dim], &x[r * dim..(r + 1) * dim]);
    }
}

/// Pick the single best-separating feature (either orientation) as the initial direction.
fn initial_direction(x: &[f64], dim: usize, labels: &[i8], test_fdr: f64) -> Vec<f64> {
    let n = labels.len();
    let pi0 = stats::estimate_pi0(labels);
    let all: Vec<usize> = (0..n).collect();
    let mut best_w = vec![0.0f64; dim];
    best_w[dim - 1] = 0.0;
    let mut best_count = -1i64;
    let mut scores = vec![0.0f64; n];
    for j in 0..dim - 1 {
        for &sign in &[1.0f64, -1.0f64] {
            for i in 0..n {
                scores[i] = sign * x[i * dim + j];
            }
            let q = stats::qvalues(&scores, labels, pi0);
            let count = q
                .iter()
                .zip(labels.iter())
                .filter(|(qi, &l)| l > 0 && **qi < test_fdr)
                .count() as i64;
            if count > best_count {
                best_count = count;
                for v in best_w.iter_mut() {
                    *v = 0.0;
                }
                best_w[j] = sign;
            }
        }
    }
    let _ = all;
    best_w
}

/// Train the semi-supervised SVM on `train_rows`, warm-started from `w0`.
fn train_fold(
    x: &[f64],
    dim: usize,
    labels: &[i8],
    train_rows: &[usize],
    w0: &[f64],
    p: &Params,
    rng: &mut Rng,
    hp: Hp,
) -> Vec<f64> {
    let mut w = w0.to_vec();
    let mut scores = vec![0.0f64; train_rows.len()];
    let sub_labels: Vec<i8> = train_rows.iter().map(|&r| labels[r]).collect();
    let pi0 = stats::estimate_pi0(&sub_labels);

    for _iter in 0..hp.maxiter {
        score_all(x, dim, &w, train_rows, &mut scores);
        let q = stats::qvalues(&scores, &sub_labels, pi0);

        // positives: targets under test_fdr ; negatives: all decoys
        let mut pos: Vec<usize> = Vec::new();
        let mut neg: Vec<usize> = Vec::new();
        for (k, &r) in train_rows.iter().enumerate() {
            if labels[r] > 0 {
                if q[k] < p.test_fdr {
                    pos.push(r);
                }
            } else {
                neg.push(r);
            }
        }
        if pos.is_empty() || neg.is_empty() {
            continue;
        }

        // optional subsampling for training speed
        if hp.subset > 0 {
            let cap = hp.subset;
            subsample(&mut pos, cap / 2, rng);
            subsample(&mut neg, cap - pos.len().min(cap / 2), rng);
        }

        let c_pos = hp.alpha;
        let c_neg = hp.beta;

        let mut rows: Vec<usize> = Vec::with_capacity(pos.len() + neg.len());
        let mut y: Vec<f64> = Vec::with_capacity(pos.len() + neg.len());
        let mut c: Vec<f64> = Vec::with_capacity(pos.len() + neg.len());
        for &r in &pos {
            rows.push(r);
            y.push(1.0);
            c.push(c_pos);
        }
        for &r in &neg {
            rows.push(r);
            y.push(-1.0);
            c.push(c_neg);
        }
        let prob = Problem { x, dim, rows: &rows, y: &y, c: &c };
        train(&prob, &mut w, p.max_newton);
    }
    w
}

fn subsample(v: &mut Vec<usize>, cap: usize, rng: &mut Rng) {
    if cap == 0 || v.len() <= cap {
        return;
    }
    // partial Fisher-Yates
    for i in 0..cap {
        let j = i + rng.below(v.len() - i);
        v.swap(i, j);
    }
    v.truncate(cap);
}

pub struct Output {
    pub score: Vec<f64>,
    pub qval: Vec<f64>,
    pub pep: Vec<f64>,
    /// (alpha, beta) actually used, and whether they were chosen by cross-validation.
    pub c_alpha: f64,
    pub c_beta: f64,
    pub c_selected: bool,
}

/// Full 3-fold pass with the given hyperparameters; returns out-of-fold scores
/// for every PSM (each scored by a model that never saw its fold).
fn cv_scores(
    x: &[f64],
    dim: usize,
    labels: &[i8],
    fold: &[u8],
    w0: &[f64],
    p: &Params,
    hp: Hp,
    seed: u64,
) -> Vec<f64> {
    let n = labels.len();
    // Fresh RNG per pass so candidates are compared under identical subsampling.
    let mut rng = Rng(seed.max(1));
    let mut final_score = vec![0.0f64; n];
    for test in 0..3u8 {
        let train_rows: Vec<usize> = (0..n).filter(|&i| fold[i] != test).collect();
        let test_rows: Vec<usize> = (0..n).filter(|&i| fold[i] == test).collect();
        let w = train_fold(x, dim, labels, &train_rows, w0, p, &mut rng, hp);
        let mut sc = vec![0.0f64; test_rows.len()];
        score_all(x, dim, &w, &test_rows, &mut sc);
        for (k, &r) in test_rows.iter().enumerate() {
            final_score[r] = sc[k];
        }
    }
    final_score
}

/// Number of target PSMs below `test_fdr`, computed on out-of-fold scores.
fn yield_at_fdr(scores: &[f64], labels: &[i8], test_fdr: f64) -> usize {
    let pi0 = stats::estimate_pi0(labels);
    let q = stats::qvalues(scores, labels, pi0);
    q.iter()
        .zip(labels.iter())
        .filter(|(qi, &l)| l > 0 && **qi < test_fdr)
        .count()
}

/// Pick (alpha, beta) by cross-validation: for each candidate, run an abbreviated
/// 3-fold pass and keep the one with the highest out-of-fold yield at `test_fdr`.
/// Mirrors the reference's "Selecting Cpos/Cneg by cross-validation" step.
fn select_c(
    x: &[f64],
    dim: usize,
    labels: &[i8],
    fold: &[u8],
    w0: &[f64],
    p: &Params,
) -> (f64, f64) {
    let mut best = (C_POS_DEFAULT, C_NEG_DEFAULT);
    let mut best_yield = -1i64;
    for &alpha in C_POS_GRID.iter() {
        for &beta in C_NEG_GRID.iter() {
            let hp = Hp {
                alpha,
                beta,
                maxiter: p.c_select_maxiter,
                subset: if p.c_select_subset == 0 { p.subset_max_train } else { p.c_select_subset },
            };
            let sc = cv_scores(x, dim, labels, fold, w0, p, hp, p.seed);
            let y = yield_at_fdr(&sc, labels, p.test_fdr) as i64;
            if y > best_yield {
                best_yield = y;
                best = (alpha, beta);
            }
        }
    }
    best
}

pub fn run(ds: &Dataset, p: &Params) -> Output {
    let n = ds.n_psm;
    let (x, dim) = build_matrix(ds);

    // deterministic 3-fold assignment via seeded shuffle
    let mut rng = Rng(p.seed.max(1));
    let mut idx: Vec<usize> = (0..n).collect();
    for i in (1..n).rev() {
        let j = rng.below(i + 1);
        idx.swap(i, j);
    }
    let mut fold = vec![0u8; n];
    for (rank, &i) in idx.iter().enumerate() {
        fold[i] = (rank % 3) as u8;
    }

    let w0 = initial_direction(&x, dim, &ds.labels, p.test_fdr);

    let selected = p.c_alpha.is_none() || p.c_beta.is_none();
    let (alpha, beta) = if selected {
        select_c(&x, dim, &ds.labels, &fold, &w0, p)
    } else {
        (p.c_alpha.unwrap(), p.c_beta.unwrap())
    };

    let hp = Hp { alpha, beta, maxiter: p.maxiter, subset: p.subset_max_train };
    let final_score = cv_scores(&x, dim, &ds.labels, &fold, &w0, p, hp, p.seed);

    let pi0 = stats::estimate_pi0(&ds.labels);
    let qval = stats::qvalues(&final_score, &ds.labels, pi0);
    let pep = stats::peps(&final_score, &ds.labels, pi0);
    Output { score: final_score, qval, pep, c_alpha: alpha, c_beta: beta, c_selected: selected }
}
