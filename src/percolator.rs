//! Semi-supervised Percolator training: 3-fold nested cross-validation around an
//! iterative fold-local learner that separates confident targets from decoys.

use crate::mlp;
use crate::pin::Dataset;
use crate::stats;
use crate::svm::{train, Problem};
use rayon::prelude::*;
use std::collections::BTreeMap;

pub struct Params {
    pub maxiter: usize,       // semi-supervised iterations
    pub test_fdr: f64,        // FDR to pick positive training examples
    pub subset_max_train: usize, // 0 = use all
    pub seed: u64,
    pub max_newton: usize,
    pub svm_tolerance: f64,
    /// Absolute SVM class weights (`C_pos`, `C_neg`), as in the reference.
    /// `None` selects them by cross-validation (the reference's `Cpos=0` behaviour).
    pub c_alpha: Option<f64>,
    pub c_beta: Option<f64>,
    /// Budget for each candidate during the C grid search (abbreviated training).
    pub c_select_maxiter: usize,
    pub c_select_subset: usize,
    /// Worker threads for the grid search. 1 = fully serial (the default, so that
    /// running many files concurrently does not oversubscribe the machine).
    pub num_threads: usize,
    /// Leakage-free per-outer-fold selection of SVM scale, class weights,
    /// feature count, and solver tolerance.
    pub nested_selection: bool,
    /// Fold-local learner. Both models use the same normalization, folds,
    /// semi-supervised labels, out-of-fold scoring, q-values, and PEPs.
    pub model: Model,
    pub mlp_hidden: usize,
    pub mlp_epochs: usize,
    pub mlp_learning_rate: f64,
    pub mlp_l2: f64,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Model {
    Svm,
    Mlp,
}

impl Model {
    pub fn label(self) -> &'static str {
        match self {
            Model::Svm => "svm",
            Model::Mlp => "mlp",
        }
    }
}

impl Default for Params {
    fn default() -> Self {
        Params {
            maxiter: 10,
            test_fdr: 0.01,
            subset_max_train: 0,
            seed: 1,
            max_newton: 30,
            svm_tolerance: 1e-5,
            c_alpha: None,
            c_beta: None,
            c_select_maxiter: 3,
            c_select_subset: 20_000,
            num_threads: 1,
            nested_selection: false,
            model: Model::Svm,
            mlp_hidden: 8,
            mlp_epochs: 10,
            mlp_learning_rate: 0.02,
            mlp_l2: 0.0,
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
    tolerance: f64,
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
    let rows: Vec<usize> = (0..ds.n_psm).collect();
    build_matrix_fit(ds, &rows)
}

/// Per-feature centering and scaling learned from a training partition.
/// Kept separately so explanatory reports can convert SVM coefficients back to
/// the original PIN units.
struct Normalization {
    mean: Vec<f64>,
    std: Vec<f64>,
}

fn fit_normalization(ds: &Dataset, fit_rows: &[usize]) -> Normalization {
    assert!(!fit_rows.is_empty());
    let nf = ds.n_feat;
    let mut mean = vec![0.0f64; nf];
    let mut var = vec![0.0f64; nf];
    for &i in fit_rows {
        let row = ds.row(i);
        for j in 0..nf {
            mean[j] += row[j];
        }
    }
    for m in &mut mean {
        *m /= fit_rows.len() as f64;
    }
    for &i in fit_rows {
        let row = ds.row(i);
        for j in 0..nf {
            let difference = row[j] - mean[j];
            var[j] += difference * difference;
        }
    }
    let mut std = vec![1.0f64; nf];
    for j in 0..nf {
        let value = (var[j] / fit_rows.len() as f64).sqrt();
        std[j] = if value > 1e-12 { value } else { 1.0 };
    }
    Normalization { mean, std }
}

fn transform_matrix(ds: &Dataset, normalization: &Normalization) -> (Vec<f64>, usize) {
    let nf = ds.n_feat;
    let dim = nf + 1;
    let mut x = vec![0.0f64; ds.n_psm * dim];
    for i in 0..ds.n_psm {
        let row = ds.row(i);
        let base = i * dim;
        for j in 0..nf {
            x[base + j] = (row[j] - normalization.mean[j]) / normalization.std[j];
        }
        x[base + nf] = 1.0;
    }
    (x, dim)
}

/// Normalize from `fit_rows` only, then transform every row. Nested selection
/// uses this to keep outer-test and inner-validation features out of fitting.
fn build_matrix_fit(ds: &Dataset, fit_rows: &[usize]) -> (Vec<f64>, usize) {
    let normalization = fit_normalization(ds, fit_rows);
    transform_matrix(ds, &normalization)
}

fn score_all(x: &[f64], dim: usize, w: &[f64], rows: &[usize], out: &mut [f64]) {
    for (k, &r) in rows.iter().enumerate() {
        out[k] = crate::simd::dot(&w[..dim], &x[r * dim..(r + 1) * dim]);
    }
}

enum FoldModel {
    Svm(Vec<f64>),
    Mlp(mlp::Network),
}

impl FoldModel {
    fn score_rows(&self, x: &[f64], dim: usize, rows: &[usize], out: &mut [f64]) {
        match self {
            FoldModel::Svm(weights) => score_all(x, dim, weights, rows, out),
            FoldModel::Mlp(network) => {
                for (k, &row) in rows.iter().enumerate() {
                    out[k] = network.score(&x[row * dim..(row + 1) * dim]);
                }
            }
        }
    }
}

/// Score held-out rows on a fold-comparable scale fitted from training decoys.
/// This is needed when nested selection chooses different C values per fold.
fn standardized_heldout_scores(
    model: &FoldModel,
    x: &[f64],
    dim: usize,
    labels: &[i8],
    train_rows: &[usize],
    heldout_rows: &[usize],
) -> Vec<f64> {
    let mut reference_rows: Vec<usize> = train_rows
        .iter()
        .copied()
        .filter(|&row| labels[row] < 0)
        .collect();
    if reference_rows.len() < 2 {
        reference_rows = train_rows.to_vec();
    }
    let mut reference_scores = vec![0.0; reference_rows.len()];
    model.score_rows(x, dim, &reference_rows, &mut reference_scores);
    let mean = reference_scores.iter().sum::<f64>() / reference_scores.len() as f64;
    let variance = reference_scores
        .iter()
        .map(|score| {
            let difference = score - mean;
            difference * difference
        })
        .sum::<f64>()
        / reference_scores.len() as f64;
    let standard_deviation = variance.sqrt().max(1e-12);
    let mut heldout_scores = vec![0.0; heldout_rows.len()];
    model.score_rows(x, dim, heldout_rows, &mut heldout_scores);
    for score in &mut heldout_scores {
        *score = (*score - mean) / standard_deviation;
    }
    heldout_scores
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

/// Train the selected semi-supervised learner on `train_rows`, initialized from `w0`.
#[allow(clippy::too_many_arguments)]
fn train_fold(
    x: &[f64],
    dim: usize,
    labels: &[i8],
    train_rows: &[usize],
    w0: &[f64],
    p: &Params,
    rng: &mut Rng,
    hp: Hp,
    model_seed: u64,
    feature_mask: Option<&[bool]>,
) -> FoldModel {
    let mut model = match p.model {
        Model::Svm => FoldModel::Svm(w0.to_vec()),
        Model::Mlp => FoldModel::Mlp(mlp::Network::new(dim, p.mlp_hidden, w0, model_seed)),
    };
    let mut scores = vec![0.0f64; train_rows.len()];
    let sub_labels: Vec<i8> = train_rows.iter().map(|&r| labels[r]).collect();
    let pi0 = stats::estimate_pi0(&sub_labels);

    for _iter in 0..hp.maxiter {
        model.score_rows(x, dim, train_rows, &mut scores);
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
        match &mut model {
            FoldModel::Svm(weights) => {
                let prob = Problem {
                    x,
                    dim,
                    rows: &rows,
                    y: &y,
                    c: &c,
                    feature_mask,
                };
                train(&prob, weights, p.max_newton, hp.tolerance);
            }
            FoldModel::Mlp(network) => network.train(
                x,
                &rows,
                &y,
                &c,
                p.mlp_epochs,
                p.mlp_learning_rate,
                p.mlp_l2,
            ),
        }
    }
    model
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
    pub nested_folds: Vec<FoldSelection>,
}

#[derive(Clone, Debug, PartialEq)]
pub struct FoldSelection {
    pub outer_fold: u8,
    pub c: f64,
    pub positive_weight: f64,
    pub negative_weight: f64,
    pub feature_count: usize,
    pub tolerance: f64,
    pub inner_yield: usize,
}

/// Full 3-fold pass with the given hyperparameters; returns out-of-fold scores
/// for every PSM (each scored by a model that never saw its fold).
#[allow(clippy::too_many_arguments)]
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
    let all_folds: [u8; 3] = [0, 1, 2];

    let per_fold = |&test: &u8| -> (Vec<usize>, Vec<f64>) {
        let train_rows: Vec<usize> = (0..n).filter(|&i| fold[i] != test).collect();
        let test_rows: Vec<usize> = (0..n).filter(|&i| fold[i] == test).collect();
        // Each fold gets its own RNG stream derived from the run seed, so folds never
        // depend on one another's draws — results are identical serial or parallel.
        let fold_seed = seed.max(1) ^ ((test as u64 + 1).wrapping_mul(0x9E37_79B9_7F4A_7C15));
        let mut rng = Rng(fold_seed);
        let model = train_fold(
            x,
            dim,
            labels,
            &train_rows,
            w0,
            p,
            &mut rng,
            hp,
            fold_seed ^ 0xD1B5_4A32_D192_ED03,
            None,
        );
        let mut sc = vec![0.0f64; test_rows.len()];
        model.score_rows(x, dim, &test_rows, &mut sc);
        (test_rows, sc)
    };

    let parts: Vec<(Vec<usize>, Vec<f64>)> = if p.num_threads > 1 {
        all_folds.par_iter().map(per_fold).collect()
    } else {
        all_folds.iter().map(per_fold).collect()
    };

    let mut final_score = vec![0.0f64; n];
    for (test_rows, sc) in parts {
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
    let cands: Vec<(f64, f64)> = C_POS_GRID
        .iter()
        .flat_map(|&a| C_NEG_GRID.iter().map(move |&b| (a, b)))
        .collect();
    let subset = if p.c_select_subset == 0 { p.subset_max_train } else { p.c_select_subset };

    // Every candidate re-seeds from p.seed, so its score is independent of evaluation
    // order — the parallel and serial paths return bit-identical results.
    let eval = |&(alpha, beta): &(f64, f64)| -> usize {
        let hp = Hp {
            alpha,
            beta,
            maxiter: p.c_select_maxiter,
            subset,
            tolerance: p.svm_tolerance,
        };
        let sc = cv_scores(x, dim, labels, fold, w0, p, hp, p.seed);
        yield_at_fdr(&sc, labels, p.test_fdr)
    };

    let yields: Vec<usize> = if p.num_threads > 1 {
        cands.par_iter().map(eval).collect()
    } else {
        cands.iter().map(eval).collect()
    };

    // First maximum wins, so ties resolve by grid order regardless of threading.
    let mut best_i = 0;
    for i in 1..cands.len() {
        if yields[i] > yields[best_i] {
            best_i = i;
        }
    }
    cands.get(best_i).copied().unwrap_or((C_POS_DEFAULT, C_NEG_DEFAULT))
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
        // A private pool keeps thread count under this run's control, so callers that
        // already parallelize across files (the benchmark harness) stay single-threaded.
        match rayon::ThreadPoolBuilder::new().num_threads(p.num_threads).build() {
            Ok(pool) if p.num_threads > 1 => {
                pool.install(|| select_c(&x, dim, &ds.labels, &fold, &w0, p))
            }
            _ => select_c(&x, dim, &ds.labels, &fold, &w0, p),
        }
    } else {
        (p.c_alpha.unwrap(), p.c_beta.unwrap())
    };

    let hp = Hp {
        alpha,
        beta,
        maxiter: p.maxiter,
        subset: p.subset_max_train,
        tolerance: p.svm_tolerance,
    };
    let final_score = cv_scores(&x, dim, &ds.labels, &fold, &w0, p, hp, p.seed);

    let pi0 = stats::estimate_pi0(&ds.labels);
    let qval = stats::qvalues(&final_score, &ds.labels, pi0);
    let pep = stats::peps(&final_score, &ds.labels, pi0);
    Output { score: final_score, qval, pep, c_alpha: alpha, c_beta: beta, c_selected: selected }
}
