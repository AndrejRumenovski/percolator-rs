//! Semi-supervised Percolator training: 3-fold nested cross-validation around an
//! iterative fold-local learner that separates confident targets from decoys.

use crate::mlp;
use crate::pin::Dataset;
use crate::stats;
use crate::svm::{train, Problem};
use rayon::prelude::*;
use std::collections::BTreeMap;

pub struct Params {
    pub maxiter: usize,          // semi-supervised iterations
    pub test_fdr: f64,           // FDR to pick positive training examples
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
    #[cfg(feature = "profiling")]
    let _normalization =
        crate::profile::Scope::with_elements("stage", "normalization_total", ds.n_psm);
    let rows: Vec<usize> = (0..ds.n_psm).collect();
    #[cfg(feature = "profiling")]
    crate::profile::allocation_site(
        "percolator::build_matrix row indices",
        1,
        (rows.capacity() * std::mem::size_of::<usize>()) as u64,
    );
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
    #[cfg(feature = "profiling")]
    let _fit = crate::profile::Scope::with_elements(
        "normalization",
        "fit_mean_and_variance",
        fit_rows.len(),
    );
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
    #[cfg(feature = "profiling")]
    crate::profile::allocation_site(
        "percolator::fit_normalization vectors",
        3,
        (3 * nf * std::mem::size_of::<f64>()) as u64,
    );
    for j in 0..nf {
        let value = (var[j] / fit_rows.len() as f64).sqrt();
        std[j] = if value > 1e-12 { value } else { 1.0 };
    }
    Normalization { mean, std }
}

fn transform_matrix(ds: &Dataset, normalization: &Normalization) -> (Vec<f64>, usize) {
    let nf = ds.n_feat;
    let dim = nf + 1;
    #[cfg(feature = "profiling")]
    let allocation_start = std::time::Instant::now();
    let mut x = vec![0.0f64; ds.n_psm * dim];
    #[cfg(feature = "profiling")]
    {
        crate::profile::record(
            "normalization",
            "matrix_allocation_and_zeroing",
            allocation_start.elapsed(),
            Some(x.len() as u64),
            Some((x.len() * std::mem::size_of::<f64>()) as u64),
        );
        crate::profile::allocation_site(
            "percolator::normalized design matrix",
            1,
            (x.capacity() * std::mem::size_of::<f64>()) as u64,
        );
    }
    #[cfg(feature = "profiling")]
    let transform_start = std::time::Instant::now();
    for i in 0..ds.n_psm {
        let row = ds.row(i);
        let base = i * dim;
        for j in 0..nf {
            x[base + j] = (row[j] - normalization.mean[j]) / normalization.std[j];
        }
        x[base + nf] = 1.0;
    }
    #[cfg(feature = "profiling")]
    crate::profile::record(
        "normalization",
        "matrix_transform",
        transform_start.elapsed(),
        Some(ds.n_psm as u64),
        Some((x.len() * std::mem::size_of::<f64>()) as u64),
    );
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
        #[cfg(feature = "profiling")]
        let _scoring =
            crate::profile::Scope::with_elements("scoring", "model_score_rows", rows.len());
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
    #[cfg(feature = "profiling")]
    let _context = crate::profile::context(Some("initial_direction"), None, None, None);
    #[cfg(feature = "profiling")]
    let _initial =
        crate::profile::Scope::with_elements("stage", "initial_direction_selection", labels.len());
    let n = labels.len();
    let pi0 = stats::estimate_pi0(labels);
    let all: Vec<usize> = (0..n).collect();
    let mut best_w = vec![0.0f64; dim];
    best_w[dim - 1] = 0.0;
    let mut best_count = -1i64;
    let mut scores = vec![0.0f64; n];
    #[cfg(feature = "profiling")]
    crate::profile::allocation_site(
        "percolator::initial_direction buffers",
        3,
        ((dim + n) * std::mem::size_of::<f64>() + n * std::mem::size_of::<usize>()) as u64,
    );
    #[cfg(feature = "profiling")]
    let mut scoring_time = std::time::Duration::ZERO;
    for j in 0..dim - 1 {
        for &sign in &[1.0f64, -1.0f64] {
            #[cfg(feature = "profiling")]
            let scoring_start = std::time::Instant::now();
            for i in 0..n {
                scores[i] = sign * x[i * dim + j];
            }
            #[cfg(feature = "profiling")]
            {
                scoring_time += scoring_start.elapsed();
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
    #[cfg(feature = "profiling")]
    crate::profile::record(
        "scoring",
        "initial_direction_feature_scoring",
        scoring_time,
        Some(((dim - 1) * 2 * n) as u64),
        None,
    );
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
    #[cfg(feature = "profiling")]
    let _fold_training = crate::profile::Scope::with_elements(
        "cross_validation",
        "fold_training_total",
        train_rows.len(),
    );
    #[cfg(feature = "profiling")]
    let setup_start = std::time::Instant::now();
    let mut model = match p.model {
        Model::Svm => FoldModel::Svm(w0.to_vec()),
        Model::Mlp => FoldModel::Mlp(mlp::Network::new(dim, p.mlp_hidden, w0, model_seed)),
    };
    let mut scores = vec![0.0f64; train_rows.len()];
    let sub_labels: Vec<i8> = train_rows.iter().map(|&r| labels[r]).collect();
    let pi0 = stats::estimate_pi0(&sub_labels);
    #[cfg(feature = "profiling")]
    {
        crate::profile::record(
            "cross_validation",
            "fold_training_setup",
            setup_start.elapsed(),
            Some(train_rows.len() as u64),
            None,
        );
        crate::profile::allocation_site(
            "percolator::train_fold persistent buffers",
            3,
            (std::mem::size_of_val(w0)
                + scores.capacity() * std::mem::size_of::<f64>()
                + sub_labels.capacity() * std::mem::size_of::<i8>()) as u64,
        );
    }

    for iteration in 0..hp.maxiter {
        #[cfg(not(feature = "profiling"))]
        let _ = iteration;
        #[cfg(feature = "profiling")]
        let _iteration_context = crate::profile::context(
            Some("semi_supervised_iteration"),
            None,
            Some(iteration),
            None,
        );
        #[cfg(feature = "profiling")]
        let _iteration = crate::profile::Scope::new("semi_supervised_iteration", "iteration_total");
        model.score_rows(x, dim, train_rows, &mut scores);
        let q = stats::qvalues(&scores, &sub_labels, pi0);

        // positives: targets under test_fdr ; negatives: all decoys
        #[cfg(feature = "profiling")]
        let positive_start = std::time::Instant::now();
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
        #[cfg(feature = "profiling")]
        {
            crate::profile::record(
                "semi_supervised_iteration",
                "confident_positive_selection",
                positive_start.elapsed(),
                Some(train_rows.len() as u64),
                None,
            );
            crate::profile::allocation_site(
                "percolator::train_fold positive/negative rows",
                2,
                ((pos.capacity() + neg.capacity()) * std::mem::size_of::<usize>()) as u64,
            );
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

        #[cfg(feature = "profiling")]
        let training_buffer_start = std::time::Instant::now();
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
        #[cfg(feature = "profiling")]
        {
            crate::profile::record(
                "semi_supervised_iteration",
                "training_row_buffer_setup",
                training_buffer_start.elapsed(),
                Some(rows.len() as u64),
                None,
            );
            crate::profile::allocation_site(
                "percolator::train_fold iteration training buffers",
                3,
                (rows.capacity() * std::mem::size_of::<usize>()
                    + (y.capacity() + c.capacity()) * std::mem::size_of::<f64>())
                    as u64,
            );
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
    #[cfg(feature = "profiling")]
    let _cv = crate::profile::Scope::with_elements("stage", "cross_validation_total", labels.len());
    let n = labels.len();
    let all_folds: [u8; 3] = [0, 1, 2];

    let per_fold = |&test: &u8| -> (Vec<usize>, Vec<f64>) {
        #[cfg(feature = "profiling")]
        let _fold_context =
            crate::profile::context(Some("cross_validation_fold"), Some(test), None, None);
        #[cfg(feature = "profiling")]
        let _fold_total = crate::profile::Scope::new("cross_validation", "fold_total");
        #[cfg(feature = "profiling")]
        let setup_start = std::time::Instant::now();
        let train_rows: Vec<usize> = (0..n).filter(|&i| fold[i] != test).collect();
        let test_rows: Vec<usize> = (0..n).filter(|&i| fold[i] == test).collect();
        // Each fold gets its own RNG stream derived from the run seed, so folds never
        // depend on one another's draws — results are identical serial or parallel.
        let fold_seed = seed.max(1) ^ ((test as u64 + 1).wrapping_mul(0x9E37_79B9_7F4A_7C15));
        let mut rng = Rng(fold_seed);
        #[cfg(feature = "profiling")]
        {
            crate::profile::record(
                "cross_validation",
                "fold_setup",
                setup_start.elapsed(),
                Some(n as u64),
                None,
            );
            crate::profile::allocation_site(
                "percolator::cv_scores fold row vectors",
                2,
                ((train_rows.capacity() + test_rows.capacity()) * std::mem::size_of::<usize>())
                    as u64,
            );
        }
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
        #[cfg(feature = "profiling")]
        let _heldout_context =
            crate::profile::context(Some("final_heldout_scoring"), None, None, None);
        let mut sc = vec![0.0f64; test_rows.len()];
        #[cfg(feature = "profiling")]
        crate::profile::allocation_site(
            "percolator::cv_scores heldout score vector",
            1,
            (sc.capacity() * std::mem::size_of::<f64>()) as u64,
        );
        model.score_rows(x, dim, &test_rows, &mut sc);
        (test_rows, sc)
    };

    #[cfg(feature = "profiling")]
    let dispatch_start = std::time::Instant::now();
    let parts: Vec<(Vec<usize>, Vec<f64>)> = if p.num_threads > 1 {
        all_folds.par_iter().map(per_fold).collect()
    } else {
        all_folds.iter().map(per_fold).collect()
    };
    #[cfg(feature = "profiling")]
    crate::profile::record(
        "cross_validation",
        "fold_dispatch_and_join",
        dispatch_start.elapsed(),
        Some(3),
        None,
    );

    #[cfg(feature = "profiling")]
    let merge_start = std::time::Instant::now();
    let mut final_score = vec![0.0f64; n];
    for (test_rows, sc) in parts {
        for (k, &r) in test_rows.iter().enumerate() {
            final_score[r] = sc[k];
        }
    }
    #[cfg(feature = "profiling")]
    {
        crate::profile::record(
            "cross_validation",
            "heldout_score_merge",
            merge_start.elapsed(),
            Some(n as u64),
            None,
        );
        crate::profile::allocation_site(
            "percolator::cv_scores final score vector",
            1,
            (final_score.capacity() * std::mem::size_of::<f64>()) as u64,
        );
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
    let subset = if p.c_select_subset == 0 {
        p.subset_max_train
    } else {
        p.c_select_subset
    };

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
    cands
        .get(best_i)
        .copied()
        .unwrap_or((C_POS_DEFAULT, C_NEG_DEFAULT))
}

#[derive(Clone)]
struct RankedFeature {
    index: usize,
    sign: f64,
    yield_count: usize,
}

#[derive(Clone, Copy)]
struct NestedCandidate {
    c: f64,
    positive_weight: f64,
    negative_weight: f64,
    feature_count: usize,
    tolerance: f64,
}

struct InnerSplit {
    x: Vec<f64>,
    train_rows: Vec<usize>,
    validation_rows: Vec<usize>,
    ranking: Vec<RankedFeature>,
}

fn rank_features(
    x: &[f64],
    dim: usize,
    labels: &[i8],
    rows: &[usize],
    test_fdr: f64,
) -> Vec<RankedFeature> {
    let subset_labels: Vec<i8> = rows.iter().map(|&row| labels[row]).collect();
    let pi0 = stats::estimate_pi0(&subset_labels);
    let mut scores = vec![0.0; rows.len()];
    let mut ranking = Vec::with_capacity(dim.saturating_sub(1));
    for feature in 0..dim - 1 {
        let mut best = RankedFeature {
            index: feature,
            sign: 1.0,
            yield_count: 0,
        };
        for &sign in &[1.0, -1.0] {
            for (k, &row) in rows.iter().enumerate() {
                scores[k] = sign * x[row * dim + feature];
            }
            let qvalues = stats::qvalues(&scores, &subset_labels, pi0);
            let count = qvalues
                .iter()
                .zip(&subset_labels)
                .filter(|(qvalue, label)| **label > 0 && **qvalue < test_fdr)
                .count();
            if count > best.yield_count {
                best.sign = sign;
                best.yield_count = count;
            }
        }
        ranking.push(best);
    }
    ranking.sort_by(|left, right| {
        right
            .yield_count
            .cmp(&left.yield_count)
            .then_with(|| left.index.cmp(&right.index))
    });
    ranking
}

fn feature_mask(dim: usize, ranking: &[RankedFeature], feature_count: usize) -> Vec<bool> {
    let mut mask = vec![false; dim];
    for feature in ranking.iter().take(feature_count) {
        mask[feature.index] = true;
    }
    mask[dim - 1] = true; // bias
    mask
}

fn ranked_initial_direction(dim: usize, ranking: &[RankedFeature]) -> Vec<f64> {
    let mut initial = vec![0.0; dim];
    if let Some(best) = ranking.first() {
        initial[best.index] = best.sign;
    }
    initial
}

fn assign_folds(rows: &[usize], fold_count: u8, seed: u64, total_rows: usize) -> Vec<u8> {
    let mut shuffled = rows.to_vec();
    let mut rng = Rng(seed.max(1));
    for i in (1..shuffled.len()).rev() {
        let j = rng.below(i + 1);
        shuffled.swap(i, j);
    }
    let mut fold = vec![u8::MAX; total_rows];
    for (rank, &row) in shuffled.iter().enumerate() {
        fold[row] = (rank % fold_count as usize) as u8;
    }
    fold
}

/// Keep reports of the same candidate from different engines together.  Otherwise
/// an engine-A report could appear in training while its engine-B counterpart is
/// scored in the held-out fold, which would leak the candidate's label and agreement
/// signal into cross validation.
fn assign_dataset_folds(ds: &Dataset, rows: &[usize], fold_count: u8, seed: u64) -> Vec<u8> {
    if !ds.ensemble {
        return assign_folds(rows, fold_count, seed, ds.n_psm);
    }
    let mut by_candidate: BTreeMap<(i64, i8, String), Vec<usize>> = BTreeMap::new();
    for &row in rows {
        by_candidate
            .entry((ds.scan[row], ds.labels[row], ds.peptide[row].clone()))
            .or_default()
            .push(row);
    }
    let mut candidates: Vec<Vec<usize>> = by_candidate.into_values().collect();
    let mut rng = Rng(seed.max(1));
    for index in (1..candidates.len()).rev() {
        candidates.swap(index, rng.below(index + 1));
    }
    let mut fold = vec![u8::MAX; ds.n_psm];
    for (rank, candidate) in candidates.into_iter().enumerate() {
        for row in candidate {
            fold[row] = (rank % fold_count as usize) as u8;
        }
    }
    fold
}

fn inner_splits(
    ds: &Dataset,
    outer_train: &[usize],
    outer_fold: u8,
    p: &Params,
) -> Vec<InnerSplit> {
    const INNER_FOLDS: u8 = 2;
    let inner_seed = p.seed
        ^ ((outer_fold as u64 + 1).wrapping_mul(0xA24B_AED4_963E_E407))
        ^ 0x9FB2_1C65_1E98_DF25;
    let assignments = assign_dataset_folds(ds, outer_train, INNER_FOLDS, inner_seed);
    (0..INNER_FOLDS)
        .map(|validation_fold| {
            let train_rows: Vec<usize> = outer_train
                .iter()
                .copied()
                .filter(|&row| assignments[row] != validation_fold)
                .collect();
            let validation_rows: Vec<usize> = outer_train
                .iter()
                .copied()
                .filter(|&row| assignments[row] == validation_fold)
                .collect();
            let (x, dim) = build_matrix_fit(ds, &train_rows);
            let ranking = rank_features(&x, dim, &ds.labels, &train_rows, p.test_fdr);
            InnerSplit {
                x,
                train_rows,
                validation_rows,
                ranking,
            }
        })
        .collect()
}

fn evaluate_nested_candidate(
    splits: &[InnerSplit],
    dim: usize,
    labels: &[i8],
    outer_fold: u8,
    candidate: NestedCandidate,
    p: &Params,
) -> usize {
    let mut validation_scores = Vec::new();
    let mut validation_labels = Vec::new();
    for (inner_fold, split) in splits.iter().enumerate() {
        let mask = feature_mask(dim, &split.ranking, candidate.feature_count);
        let initial = ranked_initial_direction(dim, &split.ranking);
        let hp = Hp {
            alpha: candidate.c * candidate.positive_weight,
            beta: candidate.c * candidate.negative_weight,
            // Nested validation evaluates the same training budget used by the
            // final outer model; abbreviated proxy training caused the legacy
            // selector's choices not to transfer.
            maxiter: p.maxiter,
            subset: p.subset_max_train,
            tolerance: candidate.tolerance,
        };
        let fold_seed = p.seed
            ^ ((outer_fold as u64 + 1).wrapping_mul(0xD6E8_FD50_1A4B_8C27))
            ^ ((inner_fold as u64 + 1).wrapping_mul(0x9E37_79B9_7F4A_7C15));
        let mut rng = Rng(fold_seed.max(1));
        let model = train_fold(
            &split.x,
            dim,
            labels,
            &split.train_rows,
            &initial,
            p,
            &mut rng,
            hp,
            fold_seed ^ 0x94D0_49BB_1331_11EB,
            Some(&mask),
        );
        let scores = standardized_heldout_scores(
            &model,
            &split.x,
            dim,
            labels,
            &split.train_rows,
            &split.validation_rows,
        );
        validation_scores.extend(scores);
        validation_labels.extend(split.validation_rows.iter().map(|&row| labels[row]));
    }
    yield_at_fdr(&validation_scores, &validation_labels, p.test_fdr)
}

fn select_stage(
    splits: &[InnerSplit],
    dim: usize,
    labels: &[i8],
    outer_fold: u8,
    candidates: &[NestedCandidate],
    p: &Params,
) -> (NestedCandidate, usize) {
    let mut best = candidates[0];
    let mut best_yield = evaluate_nested_candidate(splits, dim, labels, outer_fold, best, p);
    for &candidate in &candidates[1..] {
        let candidate_yield =
            evaluate_nested_candidate(splits, dim, labels, outer_fold, candidate, p);
        if candidate_yield > best_yield {
            best = candidate;
            best_yield = candidate_yield;
        }
    }
    (best, best_yield)
}

fn select_outer_hyperparameters(
    ds: &Dataset,
    outer_train: &[usize],
    outer_fold: u8,
    p: &Params,
) -> FoldSelection {
    let dim = ds.n_feat + 1;
    let splits = inner_splits(ds, outer_train, outer_fold, p);
    let mut selected = NestedCandidate {
        c: 1.0,
        positive_weight: 1.0,
        negative_weight: 4.0,
        feature_count: ds.n_feat,
        tolerance: p.svm_tolerance,
    };

    // Stage 1 jointly selects identifiable SVM scale and the negative:positive
    // class-weight ratio. Keeping positive weight at one avoids redundant
    // parameterizations where C and both weights can rescale one another.
    let mut scale_weight_candidates = Vec::new();
    for c in [1.0, 0.25, 4.0] {
        for negative_weight in [4.0, 1.0, 16.0] {
            scale_weight_candidates.push(NestedCandidate {
                c,
                negative_weight,
                ..selected
            });
        }
    }
    (selected, _) = select_stage(
        &splits,
        dim,
        &ds.labels,
        outer_fold,
        &scale_weight_candidates,
        p,
    );

    // Stage 2 chooses a training-only univariate ranking cutoff. Ties retain
    // all features because the existing full model is the conservative choice.
    let mut feature_counts = vec![ds.n_feat, ds.n_feat.min(8), ds.n_feat.min(4)];
    feature_counts.dedup();
    let feature_candidates: Vec<NestedCandidate> = feature_counts
        .into_iter()
        .map(|feature_count| NestedCandidate {
            feature_count,
            ..selected
        })
        .collect();
    (selected, _) = select_stage(&splits, dim, &ds.labels, outer_fold, &feature_candidates, p);

    // Stage 3 selects the Newton gradient-norm stopping tolerance.
    let mut tolerances = vec![p.svm_tolerance];
    for tolerance in [1e-3, 1e-5, 1e-7] {
        if !tolerances.contains(&tolerance) {
            tolerances.push(tolerance);
        }
    }
    let tolerance_candidates: Vec<NestedCandidate> = tolerances
        .into_iter()
        .map(|tolerance| NestedCandidate {
            tolerance,
            ..selected
        })
        .collect();
    let (selected, inner_yield) = select_stage(
        &splits,
        dim,
        &ds.labels,
        outer_fold,
        &tolerance_candidates,
        p,
    );

    FoldSelection {
        outer_fold,
        c: selected.c,
        positive_weight: selected.positive_weight,
        negative_weight: selected.negative_weight,
        feature_count: selected.feature_count,
        tolerance: selected.tolerance,
        inner_yield,
    }
}

fn nested_cv_scores(ds: &Dataset, outer_fold: &[u8], p: &Params) -> (Vec<f64>, Vec<FoldSelection>) {
    let outer_folds = [0u8, 1, 2];
    let evaluate_outer = |&test_fold: &u8| {
        let train_rows: Vec<usize> = (0..ds.n_psm)
            .filter(|&row| outer_fold[row] != test_fold)
            .collect();
        let test_rows: Vec<usize> = (0..ds.n_psm)
            .filter(|&row| outer_fold[row] == test_fold)
            .collect();
        let selected = select_outer_hyperparameters(ds, &train_rows, test_fold, p);
        let (x, dim) = build_matrix_fit(ds, &train_rows);
        let ranking = rank_features(&x, dim, &ds.labels, &train_rows, p.test_fdr);
        let mask = feature_mask(dim, &ranking, selected.feature_count);
        let initial = ranked_initial_direction(dim, &ranking);
        let hp = Hp {
            alpha: selected.c * selected.positive_weight,
            beta: selected.c * selected.negative_weight,
            maxiter: p.maxiter,
            subset: p.subset_max_train,
            tolerance: selected.tolerance,
        };
        let fold_seed = p.seed ^ ((test_fold as u64 + 1).wrapping_mul(0x9E37_79B9_7F4A_7C15));
        let mut rng = Rng(fold_seed.max(1));
        let model = train_fold(
            &x,
            dim,
            &ds.labels,
            &train_rows,
            &initial,
            p,
            &mut rng,
            hp,
            fold_seed ^ 0xD1B5_4A32_D192_ED03,
            Some(&mask),
        );
        let scores =
            standardized_heldout_scores(&model, &x, dim, &ds.labels, &train_rows, &test_rows);
        (test_rows, scores, selected)
    };

    let parts: Vec<_> = if p.num_threads > 1 {
        outer_folds.par_iter().map(evaluate_outer).collect()
    } else {
        outer_folds.iter().map(evaluate_outer).collect()
    };
    let mut scores = vec![0.0; ds.n_psm];
    let mut selections = Vec::with_capacity(3);
    for (test_rows, fold_scores, selection) in parts {
        for (k, &row) in test_rows.iter().enumerate() {
            scores[row] = fold_scores[k];
        }
        selections.push(selection);
    }
    selections.sort_by_key(|selection| selection.outer_fold);
    (scores, selections)
}

pub fn run(ds: &Dataset, p: &Params) -> Output {
    let n = ds.n_psm;

    // Deterministic 3-fold assignment; ensemble candidate duplicates stay together.
    #[cfg(feature = "profiling")]
    let fold_start = std::time::Instant::now();
    let all_rows: Vec<usize> = (0..n).collect();
    let fold = assign_dataset_folds(ds, &all_rows, 3, p.seed);
    #[cfg(feature = "profiling")]
    {
        crate::profile::record(
            "stage",
            "fold_creation_and_setup",
            fold_start.elapsed(),
            Some(n as u64),
            None,
        );
        crate::profile::allocation_site(
            "percolator::run fold setup vectors",
            2,
            (all_rows.capacity() * std::mem::size_of::<usize>()
                + fold.capacity() * std::mem::size_of::<u8>()) as u64,
        );
    }

    if p.nested_selection {
        assert_eq!(
            p.model,
            Model::Svm,
            "nested selection currently supports only SVM"
        );
        let nested = || nested_cv_scores(ds, &fold, p);
        let (final_score, nested_folds) = match rayon::ThreadPoolBuilder::new()
            .num_threads(p.num_threads)
            .build()
        {
            Ok(pool) if p.num_threads > 1 => pool.install(nested),
            _ => nested(),
        };
        let pi0 = stats::estimate_pi0(&ds.labels);
        let qval = stats::qvalues(&final_score, &ds.labels, pi0);
        let pep = stats::peps(&final_score, &ds.labels, pi0);
        return Output {
            score: final_score,
            qval,
            pep,
            c_alpha: f64::NAN,
            c_beta: f64::NAN,
            c_selected: true,
            nested_folds,
        };
    }

    let (x, dim) = build_matrix(ds);

    let w0 = initial_direction(&x, dim, &ds.labels, p.test_fdr);

    let selected = p.c_alpha.is_none() || p.c_beta.is_none();
    let (alpha, beta) = if selected {
        // A private pool keeps thread count under this run's control, so callers that
        // already parallelize across files (the benchmark harness) stay single-threaded.
        match rayon::ThreadPoolBuilder::new()
            .num_threads(p.num_threads)
            .build()
        {
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
    #[cfg(feature = "profiling")]
    let _final_context = crate::profile::context(Some("final_psm_statistics"), None, None, None);
    let qval = stats::qvalues(&final_score, &ds.labels, pi0);
    let pep = stats::peps(&final_score, &ds.labels, pi0);
    Output {
        score: final_score,
        qval,
        pep,
        c_alpha: alpha,
        c_beta: beta,
        c_selected: selected,
        nested_folds: Vec::new(),
    }
}

/// One fold-local linear model retained only while producing an explanation.
/// It is deliberately reconstructed from the same seeded folds as `run`, so
/// the report describes the out-of-fold scorer rather than a leaky refit.
struct ExplanationFold {
    weights: Vec<f64>,
    normalization: Normalization,
    test_rows: Vec<usize>,
    active_features: Vec<bool>,
    score_mean: f64,
    score_std: f64,
}

pub struct FeatureStat {
    pub index: usize,
    pub name: String,
    /// Mean coefficient in the original PIN units across the three CV models.
    pub raw_weight: f64,
    pub raw_weight_sd: f64,
    /// Mean signed coefficient after one within-fold standard deviation change.
    pub standardized_effect: f64,
    pub standardized_effect_sd: f64,
    /// Pearson correlation of the raw feature with target=+1 / decoy=-1.
    pub label_correlation: f64,
    pub mean: f64,
    pub std: f64,
    pub selected_folds: usize,
    /// Decrease in accepted target PSMs after deterministic, within-test-fold
    /// permutation. This is conditional on the fitted models (no retraining).
    pub permutation_q01_drop: isize,
    pub permuted_q01: usize,
}

pub struct FeatureReport {
    pub baseline_q01: usize,
    pub features: Vec<FeatureStat>,
}

fn mean_and_sd(values: &[f64]) -> (f64, f64) {
    let mean = values.iter().sum::<f64>() / values.len() as f64;
    let variance = values
        .iter()
        .map(|value| {
            let difference = value - mean;
            difference * difference
        })
        .sum::<f64>()
        / values.len() as f64;
    (mean, variance.sqrt())
}

fn training_score_calibration(
    model: &FoldModel,
    x: &[f64],
    dim: usize,
    labels: &[i8],
    train_rows: &[usize],
) -> (f64, f64) {
    let mut reference_rows: Vec<usize> = train_rows
        .iter()
        .copied()
        .filter(|&row| labels[row] < 0)
        .collect();
    if reference_rows.len() < 2 {
        reference_rows = train_rows.to_vec();
    }
    let mut scores = vec![0.0; reference_rows.len()];
    model.score_rows(x, dim, &reference_rows, &mut scores);
    let (mean, std) = mean_and_sd(&scores);
    (mean, std.max(1e-12))
}

fn outer_fold_assignments(ds: &Dataset, seed: u64) -> Vec<u8> {
    let rows: Vec<usize> = (0..ds.n_psm).collect();
    assign_dataset_folds(ds, &rows, 3, seed)
}

fn explain_fixed_models(
    ds: &Dataset,
    p: &Params,
    output: &Output,
    fold: &[u8],
) -> Vec<ExplanationFold> {
    let all_rows: Vec<usize> = (0..ds.n_psm).collect();
    let normalization = fit_normalization(ds, &all_rows);
    let (x, dim) = transform_matrix(ds, &normalization);
    let initial = initial_direction(&x, dim, &ds.labels, p.test_fdr);
    let hp = Hp {
        alpha: output.c_alpha,
        beta: output.c_beta,
        maxiter: p.maxiter,
        subset: p.subset_max_train,
        tolerance: p.svm_tolerance,
    };
    [0u8, 1, 2]
        .iter()
        .map(|&test_fold| {
            let train_rows: Vec<usize> = (0..ds.n_psm)
                .filter(|&row| fold[row] != test_fold)
                .collect();
            let test_rows: Vec<usize> = (0..ds.n_psm)
                .filter(|&row| fold[row] == test_fold)
                .collect();
            let fold_seed = p.seed ^ ((test_fold as u64 + 1).wrapping_mul(0x9E37_79B9_7F4A_7C15));
            let mut rng = Rng(fold_seed.max(1));
            let model = train_fold(
                &x,
                dim,
                &ds.labels,
                &train_rows,
                &initial,
                p,
                &mut rng,
                hp,
                fold_seed ^ 0xD1B5_4A32_D192_ED03,
                None,
            );
            let weights = match model {
                FoldModel::Svm(weights) => weights,
                FoldModel::Mlp(_) => unreachable!("feature reports require SVM"),
            };
            ExplanationFold {
                weights,
                normalization: Normalization {
                    mean: normalization.mean.clone(),
                    std: normalization.std.clone(),
                },
                test_rows,
                active_features: vec![true; dim],
                score_mean: 0.0,
                score_std: 1.0,
            }
        })
        .collect()
}

fn explain_nested_models(ds: &Dataset, p: &Params, fold: &[u8]) -> Vec<ExplanationFold> {
    [0u8, 1, 2]
        .iter()
        .map(|&test_fold| {
            let train_rows: Vec<usize> = (0..ds.n_psm)
                .filter(|&row| fold[row] != test_fold)
                .collect();
            let test_rows: Vec<usize> = (0..ds.n_psm)
                .filter(|&row| fold[row] == test_fold)
                .collect();
            let selected = select_outer_hyperparameters(ds, &train_rows, test_fold, p);
            let normalization = fit_normalization(ds, &train_rows);
            let (x, dim) = transform_matrix(ds, &normalization);
            let ranking = rank_features(&x, dim, &ds.labels, &train_rows, p.test_fdr);
            let active_features = feature_mask(dim, &ranking, selected.feature_count);
            let initial = ranked_initial_direction(dim, &ranking);
            let hp = Hp {
                alpha: selected.c * selected.positive_weight,
                beta: selected.c * selected.negative_weight,
                maxiter: p.maxiter,
                subset: p.subset_max_train,
                tolerance: selected.tolerance,
            };
            let fold_seed = p.seed ^ ((test_fold as u64 + 1).wrapping_mul(0x9E37_79B9_7F4A_7C15));
            let mut rng = Rng(fold_seed.max(1));
            let model = train_fold(
                &x,
                dim,
                &ds.labels,
                &train_rows,
                &initial,
                p,
                &mut rng,
                hp,
                fold_seed ^ 0xD1B5_4A32_D192_ED03,
                Some(&active_features),
            );
            let (score_mean, score_std) =
                training_score_calibration(&model, &x, dim, &ds.labels, &train_rows);
            let weights = match model {
                FoldModel::Svm(weights) => weights,
                FoldModel::Mlp(_) => unreachable!("feature reports require SVM"),
            };
            ExplanationFold {
                weights,
                normalization,
                test_rows,
                active_features,
                score_mean,
                score_std,
            }
        })
        .collect()
}

fn score_explanation_fold(
    ds: &Dataset,
    fold: &ExplanationFold,
    permuted_feature: Option<usize>,
    seed: u64,
    out: &mut [f64],
) {
    let mut source_rows = fold.test_rows.clone();
    if permuted_feature.is_some() {
        let mut rng = Rng(seed.max(1));
        for i in (1..source_rows.len()).rev() {
            let j = rng.below(i + 1);
            source_rows.swap(i, j);
        }
    }
    for (k, &row) in fold.test_rows.iter().enumerate() {
        let mut score = fold.weights[ds.n_feat];
        for feature in 0..ds.n_feat {
            let value_row = if permuted_feature == Some(feature) {
                source_rows[k]
            } else {
                row
            };
            let value = (ds.row(value_row)[feature] - fold.normalization.mean[feature])
                / fold.normalization.std[feature];
            score += fold.weights[feature] * value;
        }
        out[row] = (score - fold.score_mean) / fold.score_std;
    }
}

fn target_q01(scores: &[f64], labels: &[i8], test_fdr: f64) -> usize {
    let qvalues = stats::qvalues(scores, labels, stats::estimate_pi0(labels));
    qvalues
        .iter()
        .zip(labels)
        .filter(|(qvalue, label)| **label > 0 && **qvalue < test_fdr)
        .count()
}

/// Compute out-of-fold linear-model explanations. The report is intentionally
/// post hoc: it never changes rescoring or model selection, and permutation
/// importance holds the trained models fixed rather than retraining them.
pub fn feature_report(ds: &Dataset, p: &Params, output: &Output) -> FeatureReport {
    assert_eq!(
        p.model,
        Model::Svm,
        "feature reports currently support only SVM"
    );
    let fold = outer_fold_assignments(ds, p.seed);
    let models = if p.nested_selection {
        explain_nested_models(ds, p, &fold)
    } else {
        explain_fixed_models(ds, p, output, &fold)
    };
    let baseline_q01 = target_q01(&output.score, &ds.labels, p.test_fdr);
    let global_normalization = fit_normalization(ds, &(0..ds.n_psm).collect::<Vec<_>>());
    let label_values: Vec<f64> = ds.labels.iter().map(|&label| label as f64).collect();
    let (label_mean, label_std) = mean_and_sd(&label_values);
    let mut features = Vec::with_capacity(ds.n_feat);

    for feature in 0..ds.n_feat {
        let raw_weights: Vec<f64> = models
            .iter()
            .map(|model| model.weights[feature] / model.normalization.std[feature])
            .collect();
        let standardized: Vec<f64> = models.iter().map(|model| model.weights[feature]).collect();
        let (raw_weight, raw_weight_sd) = mean_and_sd(&raw_weights);
        let (standardized_effect, standardized_effect_sd) = mean_and_sd(&standardized);
        let mut covariance = 0.0;
        for row in 0..ds.n_psm {
            covariance += (ds.row(row)[feature] - global_normalization.mean[feature])
                * (ds.labels[row] as f64 - label_mean);
        }
        covariance /= ds.n_psm as f64;
        let label_correlation =
            covariance / (global_normalization.std[feature] * label_std).max(1e-12);

        let mut permuted_scores = vec![0.0; ds.n_psm];
        for (fold_index, model) in models.iter().enumerate() {
            let permutation_seed = p.seed
                ^ ((feature as u64 + 1).wrapping_mul(0xD6E8_FEB8_6659_FD93))
                ^ ((fold_index as u64 + 1).wrapping_mul(0xA24B_AED4_963E_E407));
            score_explanation_fold(
                ds,
                model,
                Some(feature),
                permutation_seed,
                &mut permuted_scores,
            );
        }
        let permuted_q01 = target_q01(&permuted_scores, &ds.labels, p.test_fdr);
        features.push(FeatureStat {
            index: feature,
            name: ds.feature_names[feature].clone(),
            raw_weight,
            raw_weight_sd,
            standardized_effect,
            standardized_effect_sd,
            label_correlation,
            mean: global_normalization.mean[feature],
            std: global_normalization.std[feature],
            selected_folds: models
                .iter()
                .filter(|model| model.active_features[feature])
                .count(),
            permutation_q01_drop: baseline_q01 as isize - permuted_q01 as isize,
            permuted_q01,
        });
    }
    FeatureReport {
        baseline_q01,
        features,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn selection_fixture() -> Dataset {
        let n_psm = 120;
        let n_feat = 4;
        let mut features = Vec::with_capacity(n_psm * n_feat);
        let mut labels = Vec::with_capacity(n_psm);
        for row in 0..n_psm {
            let label = if row % 3 == 0 { -1 } else { 1 };
            labels.push(label);
            let jitter = (row % 11) as f64 / 20.0;
            features.extend_from_slice(&[
                label as f64 * 2.0 + jitter,
                (row % 7) as f64 - 3.0,
                label as f64 * ((row % 5) as f64),
                (row % 13) as f64 / 13.0,
            ]);
        }
        Dataset {
            feature_names: (0..n_feat).map(|index| format!("f{index}")).collect(),
            n_feat,
            n_psm,
            features,
            labels,
            spec_id: (0..n_psm).map(|row| format!("psm{row}")).collect(),
            scan: (0..n_psm as i64).collect(),
            peptide: (0..n_psm).map(|row| format!("K.PEP{row}.R")).collect(),
            proteins: (0..n_psm).map(|row| format!("P{row}")).collect(),
            source: vec![0; n_psm],
            source_names: vec!["selection-fixture.pin".to_string()],
            ensemble: false,
        }
    }

    #[test]
    fn outer_test_changes_cannot_change_nested_selection() {
        let dataset = selection_fixture();
        let params = Params {
            maxiter: 2,
            c_select_maxiter: 2,
            num_threads: 1,
            ..Params::default()
        };
        let all_rows: Vec<usize> = (0..dataset.n_psm).collect();
        let outer_fold = assign_folds(&all_rows, 3, params.seed, dataset.n_psm);
        let outer_train: Vec<usize> = all_rows
            .iter()
            .copied()
            .filter(|&row| outer_fold[row] != 0)
            .collect();
        let expected = select_outer_hyperparameters(&dataset, &outer_train, 0, &params);

        let mut changed_test = selection_fixture();
        for row in all_rows.into_iter().filter(|&row| outer_fold[row] == 0) {
            changed_test.labels[row] *= -1;
            for feature in 0..changed_test.n_feat {
                changed_test.features[row * changed_test.n_feat + feature] +=
                    1_000_000.0 * (feature + 1) as f64;
            }
        }
        let actual = select_outer_hyperparameters(&changed_test, &outer_train, 0, &params);
        assert_eq!(actual, expected);
    }

    #[test]
    fn ensemble_candidate_duplicates_stay_in_the_same_fold() {
        let mut dataset = selection_fixture();
        dataset.ensemble = true;
        // Rows 0 and 1 represent the same target candidate reported by two engines.
        dataset.scan[1] = dataset.scan[0];
        dataset.labels[1] = dataset.labels[0];
        dataset.peptide[1] = dataset.peptide[0].clone();
        let rows: Vec<usize> = (0..dataset.n_psm).collect();
        let fold = assign_dataset_folds(&dataset, &rows, 3, 1);
        assert_eq!(fold[0], fold[1]);
    }
}
