//! Semi-supervised Percolator training: 3-fold nested cross-validation around an
//! iterative fold-local learner that separates confident targets from decoys.

use crate::mlp;
use crate::pin::Dataset;
use crate::stats;
use crate::svm::{train, Problem, Workspace as SvmWorkspace};
use rayon::prelude::*;
use std::collections::BTreeMap;

pub struct Params {
    pub maxiter: usize,          // semi-supervised iterations
    pub test_fdr: f64,           // FDR to pick positive training examples
    /// Probability that an incorrect target outranks its paired decoy under the
    /// null.  0.5 for a 1:1 concatenated target-decoy competition; see
    /// [`crate::stats`] for the supported input contract.
    pub null_target_win_prob: f64,
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
    /// Retention-time alignment inputs when `--rt-features` is on.
    ///
    /// The alignment is label-dependent, so it is refitted inside every outer
    /// training partition rather than once over the whole dataset.
    pub rt: Option<crate::rt::Alignment>,
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
            null_target_win_prob: 0.5,
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
            rt: None,
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

/// Per-feature centering and scaling learned from a training partition.
/// Kept separately so explanatory reports can convert SVM coefficients back to
/// the original PIN units.
struct Normalization {
    mean: Vec<f64>,
    std: Vec<f64>,
}

/// Row features with the fold's retention-time residuals patched in.
#[inline]
fn feature_row<'a>(ds: &'a Dataset, row: usize, rt: Option<&RtColumns<'_>>, scratch: &'a mut [f64]) -> &'a [f64] {
    match rt {
        None => ds.row(row),
        Some(rt) => {
            scratch.copy_from_slice(ds.row(row));
            scratch[rt.first_column] = rt.values[row * 2];
            scratch[rt.first_column + 1] = rt.values[row * 2 + 1];
            scratch
        }
    }
}

/// The two retention-time residual columns computed for one fold.
struct RtColumns<'a> {
    first_column: usize,
    values: &'a [f64],
}

fn fit_normalization(ds: &Dataset, fit_rows: &[usize], rt: Option<&RtColumns<'_>>) -> Normalization {
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
    let mut scratch = vec![0.0f64; nf];
    for &i in fit_rows {
        let row = feature_row(ds, i, rt, &mut scratch);
        for j in 0..nf {
            mean[j] += row[j];
        }
    }
    for m in &mut mean {
        *m /= fit_rows.len() as f64;
    }
    for &i in fit_rows {
        let row = feature_row(ds, i, rt, &mut scratch);
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

fn transform_matrix(
    ds: &Dataset,
    normalization: &Normalization,
    rt: Option<&RtColumns<'_>>,
) -> (Vec<f64>, usize) {
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
    let mut scratch = vec![0.0f64; nf];
    for i in 0..ds.n_psm {
        let row = feature_row(ds, i, rt, &mut scratch);
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
fn build_matrix_fit(ds: &Dataset, fit_rows: &[usize], p: &Params) -> (Vec<f64>, usize) {
    let rt_values = fold_rt_columns(ds, fit_rows, p);
    let rt = rt_columns(p, rt_values.as_deref());
    let normalization = fit_normalization(ds, fit_rows, rt.as_ref());
    transform_matrix(ds, &normalization, rt.as_ref())
}

/// Refit the retention-time alignment inside `fit_rows` and materialize its two
/// residual columns for every row.
fn fold_rt_columns(ds: &Dataset, fit_rows: &[usize], p: &Params) -> Option<Vec<f64>> {
    let alignment = p.rt.as_ref()?;
    let mut values = vec![0.0f64; ds.n_psm * 2];
    alignment.residuals(&ds.labels, &ds.source, fit_rows, &mut values);
    Some(values)
}

fn rt_columns<'a>(p: &Params, values: Option<&'a [f64]>) -> Option<RtColumns<'a>> {
    Some(RtColumns {
        first_column: p.rt.as_ref()?.first_column,
        values: values?,
    })
}

fn score_all(x: &[f64], dim: usize, w: &[f64], rows: &[usize], out: &mut [f64]) {
    if dim == 22 {
        for (k, &r) in rows.iter().enumerate() {
            out[k] = crate::simd::dot_22(&w[..22], &x[r * 22..(r + 1) * 22]);
        }
        return;
    }
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

/// Location and scale of one fold model's null distribution, measured on its own
/// training decoys.
///
/// Independently fitted folds have no reason to share an intercept or a score
/// scale, so pooling their raw margins lets an arbitrary per-fold offset decide
/// the merged ranking -- and therefore the q-values, which depend only on the
/// merged order.  Every fold is instead expressed in standard deviations above
/// its own training-decoy mean, which is a common scale with a meaning: how far
/// a PSM sits from that fold's null.
///
/// The reference anchors instead on the held-out selection boundary and the
/// held-out median decoy.  Training decoys are used here because they keep the
/// transform inside the training partition; the mapping is affine and increasing
/// either way, so it never reorders a fold internally.
fn training_null_calibration(
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
    let mut reference_scores = vec![0.0; reference_rows.len()];
    model.score_rows(x, dim, &reference_rows, &mut reference_scores);
    let (mean, standard_deviation) = mean_and_sd(&reference_scores);
    (mean, standard_deviation.max(1e-12))
}

/// Score held-out rows on the fold-comparable scale of [`training_null_calibration`].
fn standardized_heldout_scores(
    model: &FoldModel,
    x: &[f64],
    dim: usize,
    labels: &[i8],
    train_rows: &[usize],
    heldout_rows: &[usize],
) -> Vec<f64> {
    let (mean, standard_deviation) =
        training_null_calibration(model, x, dim, labels, train_rows);
    let mut heldout_scores = vec![0.0; heldout_rows.len()];
    model.score_rows(x, dim, heldout_rows, &mut heldout_scores);
    for score in &mut heldout_scores {
        *score = (*score - mean) / standard_deviation;
    }
    heldout_scores
}

/// Pick the single best-separating feature (either orientation) as the initial
/// direction, using `rows` only.
///
/// Feature *and* label information both enter this choice, so it must be fitted
/// inside the outer training partition: a direction selected on the full dataset
/// carries every held-out label into the model that later scores it.
fn initial_direction(x: &[f64], dim: usize, labels: &[i8], rows: &[usize], p: &Params) -> Vec<f64> {
    #[cfg(feature = "profiling")]
    let _context = crate::profile::context(Some("initial_direction"), None, None, None);
    #[cfg(feature = "profiling")]
    let _initial =
        crate::profile::Scope::with_elements("stage", "initial_direction_selection", rows.len());
    let n = rows.len();
    // Selecting a direction is a heuristic, not an error-rate claim: the
    // reference likewise drops the finite-sample safeguard here because it is
    // too restrictive on small partitions.
    let tdc = stats::Tdc::training(p.null_target_win_prob);
    let test_fdr = p.test_fdr;
    let fit_labels: Vec<i8> = rows.iter().map(|&row| labels[row]).collect();
    let labels = &fit_labels[..];
    let mut order = Vec::new();
    let mut best_w = vec![0.0f64; dim];
    best_w[dim - 1] = 0.0;
    let mut best_count = -1i64;
    let mut scores = vec![0.0f64; n];
    let mut ranks = vec![0u32; n];
    #[cfg(feature = "profiling")]
    crate::profile::allocation_site(
        "percolator::initial_direction buffers",
        4,
        ((dim + n) * std::mem::size_of::<f64>()
            + n * (std::mem::size_of::<usize>() + std::mem::size_of::<u32>())) as u64,
    );
    #[cfg(feature = "profiling")]
    let mut scoring_time = std::time::Duration::ZERO;
    for j in 0..dim - 1 {
        #[cfg(feature = "profiling")]
        let scoring_start = std::time::Instant::now();
        for (k, &row) in rows.iter().enumerate() {
            scores[k] = x[row * dim + j];
        }
        #[cfg(feature = "profiling")]
        {
            scoring_time += scoring_start.elapsed();
        }
        let count =
            stats::target_count_at_fdr_into(&scores, labels, tdc, test_fdr, &mut order) as i64;
        if count > best_count {
            best_count = count;
            for v in best_w.iter_mut() {
                *v = 0.0;
            }
            best_w[j] = 1.0;
        }

        let count = if n <= u32::MAX as usize && !scores.iter().any(|score| score.is_nan()) {
            let mut rank = 0u32;
            for position in 0..order.len() {
                if position > 0
                    && scores[order[position]].partial_cmp(&scores[order[position - 1]])
                        != Some(std::cmp::Ordering::Equal)
                {
                    rank += 1;
                }
                ranks[order[position]] = rank;
            }
            stats::target_count_at_reversed_ranks_into(&ranks, labels, tdc, test_fdr, &mut order)
                as i64
        } else {
            #[cfg(feature = "profiling")]
            let scoring_start = std::time::Instant::now();
            for (k, &row) in rows.iter().enumerate() {
                scores[k] = -x[row * dim + j];
            }
            #[cfg(feature = "profiling")]
            {
                scoring_time += scoring_start.elapsed();
            }
            stats::target_count_at_fdr_into(&scores, labels, tdc, test_fdr, &mut order) as i64
        };
        if count > best_count {
            best_count = count;
            for v in best_w.iter_mut() {
                *v = 0.0;
            }
            best_w[j] = -1.0;
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
    let tdc = stats::Tdc::training(p.null_target_win_prob);
    let mut qvalue_workspace = stats::QValueWorkspace::default();
    let mut accepted = Vec::new();
    let mut svm_workspace = SvmWorkspace::default();
    let mut packed_x = Vec::new();
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
        stats::target_mask_at_fdr_into(
            &scores,
            &sub_labels,
            tdc,
            p.test_fdr,
            &mut qvalue_workspace,
            &mut accepted,
        );

        // positives: targets under test_fdr ; negatives: all decoys
        #[cfg(feature = "profiling")]
        let positive_start = std::time::Instant::now();
        let mut pos: Vec<usize> = Vec::new();
        let mut neg: Vec<usize> = Vec::new();
        for (k, &r) in train_rows.iter().enumerate() {
            if labels[r] > 0 {
                if accepted[k] != 0 {
                    pos.push(k);
                }
            } else {
                neg.push(k);
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
        let mut initial_scores: Vec<f64> = Vec::with_capacity(pos.len() + neg.len());
        for &k in &pos {
            rows.push(train_rows[k]);
            y.push(1.0);
            c.push(c_pos);
            initial_scores.push(scores[k]);
        }
        for &k in &neg {
            rows.push(train_rows[k]);
            y.push(-1.0);
            c.push(c_neg);
            initial_scores.push(scores[k]);
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
                4,
                (rows.capacity() * std::mem::size_of::<usize>()
                    + (y.capacity() + c.capacity() + initial_scores.capacity())
                        * std::mem::size_of::<f64>()) as u64,
            );
        }
        match &mut model {
            FoldModel::Svm(weights) => {
                packed_x.clear();
                packed_x.reserve(rows.len() * dim);
                for &row in &rows {
                    packed_x.extend_from_slice(&x[row * dim..(row + 1) * dim]);
                }
                let prob = Problem {
                    x: &packed_x,
                    dim,
                    rows: &rows,
                    y: &y,
                    c: &c,
                    packed_rows: true,
                    feature_mask,
                };
                train(
                    &prob,
                    weights,
                    &initial_scores,
                    p.max_newton,
                    hp.tolerance,
                    &mut svm_workspace,
                );
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

/// Everything one outer fold needs, fitted inside its own training partition.
///
/// Both the normalization location/scale and the initial direction are estimated
/// from `train_rows` alone.  Fitting either on the full dataset would let a
/// held-out row shape the model that goes on to score it: the feature
/// distribution through normalization, and the labels through the direction
/// search.
struct FoldSetup {
    train_rows: Vec<usize>,
    test_rows: Vec<usize>,
    x: Vec<f64>,
    dim: usize,
    w0: Vec<f64>,
    /// Per-fold RNG stream derived from the run seed, so folds never depend on
    /// one another's draws and results are identical serial or parallel.
    seed: u64,
}

fn fold_setup(ds: &Dataset, fold: &[u8], test: u8, p: &Params, seed: u64) -> FoldSetup {
    #[cfg(feature = "profiling")]
    let _fold_context = crate::profile::context(Some("cross_validation_fold"), Some(test), None, None);
    #[cfg(feature = "profiling")]
    let setup_start = std::time::Instant::now();
    let n = ds.n_psm;
    let train_rows: Vec<usize> = (0..n).filter(|&i| fold[i] != test).collect();
    let test_rows: Vec<usize> = (0..n).filter(|&i| fold[i] == test).collect();
    assert!(
        !train_rows.is_empty(),
        "outer fold {test} left no training rows"
    );
    let (x, dim) = build_matrix_fit(ds, &train_rows, p);
    let w0 = initial_direction(&x, dim, &ds.labels, &train_rows, p);
    #[cfg(feature = "profiling")]
    crate::profile::record(
        "cross_validation",
        "fold_setup",
        setup_start.elapsed(),
        Some(n as u64),
        None,
    );
    FoldSetup {
        train_rows,
        test_rows,
        x,
        dim,
        w0,
        seed: seed.max(1) ^ ((test as u64 + 1).wrapping_mul(0x9E37_79B9_7F4A_7C15)),
    }
}

impl FoldSetup {
    /// Train on this fold's training partition and score its held-out rows.
    fn train_and_score(&self, ds: &Dataset, p: &Params, hp: Hp) -> Vec<f64> {
        let mut rng = Rng(self.seed);
        let model = train_fold(
            &self.x,
            self.dim,
            &ds.labels,
            &self.train_rows,
            &self.w0,
            p,
            &mut rng,
            hp,
            self.seed ^ 0xD1B5_4A32_D192_ED03,
            None,
        );
        #[cfg(feature = "profiling")]
        let _heldout_context =
            crate::profile::context(Some("final_heldout_scoring"), None, None, None);
        standardized_heldout_scores(
            &model,
            &self.x,
            self.dim,
            &ds.labels,
            &self.train_rows,
            &self.test_rows,
        )
    }
}

const OUTER_FOLDS: [u8; 3] = [0, 1, 2];

/// Full 3-fold pass with the given hyperparameters; returns out-of-fold scores
/// for every PSM (each scored by a model that never saw its fold).
fn cv_scores(ds: &Dataset, fold: &[u8], p: &Params, hp: Hp, seed: u64) -> Vec<f64> {
    #[cfg(feature = "profiling")]
    let _cv = crate::profile::Scope::with_elements("stage", "cross_validation_total", ds.n_psm);

    let per_fold = |&test: &u8| -> (Vec<usize>, Vec<f64>) {
        #[cfg(feature = "profiling")]
        let _fold_total = crate::profile::Scope::new("cross_validation", "fold_total");
        let setup = fold_setup(ds, fold, test, p, seed);
        let scores = setup.train_and_score(ds, p, hp);
        (setup.test_rows, scores)
    };

    #[cfg(feature = "profiling")]
    let dispatch_start = std::time::Instant::now();
    // Folds run one at a time unless the caller asked for fold-level threads:
    // each fold holds its own normalized design matrix, so a parallel pass costs
    // three of them at once.
    let parts: Vec<(Vec<usize>, Vec<f64>)> = if p.num_threads > 1 {
        OUTER_FOLDS.par_iter().map(per_fold).collect()
    } else {
        OUTER_FOLDS.iter().map(per_fold).collect()
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
    let mut final_score = vec![0.0f64; ds.n_psm];
    for (test_rows, scores) in parts {
        for (k, &row) in test_rows.iter().enumerate() {
            final_score[row] = scores[k];
        }
    }
    #[cfg(feature = "profiling")]
    crate::profile::record(
        "cross_validation",
        "heldout_score_merge",
        merge_start.elapsed(),
        Some(ds.n_psm as u64),
        None,
    );
    final_score
}

/// Number of target PSMs below `test_fdr`, computed on out-of-fold scores.
fn yield_at_fdr(scores: &[f64], labels: &[i8], test_fdr: f64, p: &Params) -> usize {
    let q = stats::qvalues(scores, labels, stats::Tdc::reported(p.null_target_win_prob));
    q.iter()
        .zip(labels.iter())
        .filter(|(qi, &l)| l > 0 && **qi < test_fdr)
        .count()
}

/// Pick this outer fold's `(alpha, beta)` inside its own training partition.
///
/// # Why this is nested
///
/// The previous implementation scored the whole `C` grid on the *out-of-fold*
/// predictions that were then reported.  Every held-out row therefore took part
/// in choosing the hyperparameters of the model that scored it: flipping only
/// the labels of one outer fold moved the selected class weights from 4:1 to
/// 0.25:1 and changed every one of that fold's reported scores.  That is
/// selection leakage regardless of how well the folds themselves are isolated,
/// because the selection step sits outside them.
///
/// Selection is therefore nested.  For one outer fold:
///
/// ```text
/// outer training rows
///   -> inner split into training / validation partitions
///   -> each grid candidate trained on inner training, scored on inner validation
///   -> the candidate with the highest pooled inner-validation yield
///   -> final model trained on the whole outer training partition with that C
///   -> the untouched outer held-out fold is scored
/// ```
///
/// The outer held-out fold appears nowhere above, so neither its labels nor its
/// features can reach the selection or the model.
fn select_c_for_fold(
    ds: &Dataset,
    outer_train: &[usize],
    outer_fold: u8,
    p: &Params,
) -> (f64, f64, usize) {
    let dim = ds.n_feat + 1;
    let cands: Vec<(f64, f64)> = C_POS_GRID
        .iter()
        .flat_map(|&a| C_NEG_GRID.iter().map(move |&b| (a, b)))
        .collect();
    let subset = if p.c_select_subset == 0 {
        p.subset_max_train
    } else {
        p.c_select_subset
    };
    let splits = inner_splits(ds, outer_train, outer_fold, p);

    // Every candidate re-seeds from the inner-fold seed, so its score is
    // independent of evaluation order -- the parallel and serial paths return
    // bit-identical results.
    let evaluate = |(index, &(alpha, beta)): (usize, &(f64, f64))| -> (usize, usize) {
        let mut validation_scores: Vec<f64> = Vec::new();
        let mut validation_labels: Vec<i8> = Vec::new();
        for (inner_fold, split) in splits.iter().enumerate() {
            let initial = ranked_initial_direction(dim, &split.ranking);
            let hp = Hp {
                alpha,
                beta,
                maxiter: p.c_select_maxiter,
                subset,
                tolerance: p.svm_tolerance,
            };
            let fold_seed = p.seed
                ^ ((outer_fold as u64 + 1).wrapping_mul(0xD6E8_FD50_1A4B_8C27))
                ^ ((inner_fold as u64 + 1).wrapping_mul(0x9E37_79B9_7F4A_7C15));
            let mut rng = Rng(fold_seed.max(1));
            let model = train_fold(
                &split.x,
                dim,
                &ds.labels,
                &split.train_rows,
                &initial,
                p,
                &mut rng,
                hp,
                fold_seed ^ 0x94D0_49BB_1331_11EB,
                None,
            );
            validation_scores.extend(standardized_heldout_scores(
                &model,
                &split.x,
                dim,
                &ds.labels,
                &split.train_rows,
                &split.validation_rows,
            ));
            validation_labels.extend(split.validation_rows.iter().map(|&row| ds.labels[row]));
        }
        (
            index,
            yield_at_fdr(&validation_scores, &validation_labels, p.test_fdr, p),
        )
    };

    let mut yields = vec![0usize; cands.len()];
    let parts: Vec<(usize, usize)> = if p.num_threads > 1 {
        cands.par_iter().enumerate().map(evaluate).collect()
    } else {
        cands.iter().enumerate().map(evaluate).collect()
    };
    for (index, value) in parts {
        yields[index] = value;
    }

    // First maximum wins, so ties resolve by grid order regardless of threading.
    let mut best_i = 0;
    for i in 1..cands.len() {
        if yields[i] > yields[best_i] {
            best_i = i;
        }
    }
    let (alpha, beta) = cands
        .get(best_i)
        .copied()
        .unwrap_or((C_POS_DEFAULT, C_NEG_DEFAULT));
    (alpha, beta, yields.get(best_i).copied().unwrap_or(0))
}

/// Three-fold pass in which each outer fold first selects its own class weights
/// inside its training partition, then trains and scores its held-out rows.
fn cv_scores_with_selected_c(
    ds: &Dataset,
    fold: &[u8],
    p: &Params,
) -> (Vec<f64>, Vec<FoldSelection>) {
    let per_fold = |&test: &u8| -> (Vec<usize>, Vec<f64>, FoldSelection) {
        let setup = fold_setup(ds, fold, test, p, p.seed);
        let (alpha, beta, inner_yield) = select_c_for_fold(ds, &setup.train_rows, test, p);
        let hp = Hp {
            alpha,
            beta,
            maxiter: p.maxiter,
            subset: p.subset_max_train,
            tolerance: p.svm_tolerance,
        };
        let scores = setup.train_and_score(ds, p, hp);
        (
            setup.test_rows,
            scores,
            FoldSelection {
                outer_fold: test,
                c: 1.0,
                positive_weight: alpha,
                negative_weight: beta,
                feature_count: ds.n_feat,
                tolerance: p.svm_tolerance,
                inner_yield,
            },
        )
    };
    let parts: Vec<(Vec<usize>, Vec<f64>, FoldSelection)> = if p.num_threads > 1 {
        OUTER_FOLDS.par_iter().map(per_fold).collect()
    } else {
        OUTER_FOLDS.iter().map(per_fold).collect()
    };
    let mut final_score = vec![0.0f64; ds.n_psm];
    let mut selections = Vec::with_capacity(parts.len());
    for (test_rows, scores, selection) in parts {
        for (k, &row) in test_rows.iter().enumerate() {
            final_score[row] = scores[k];
        }
        selections.push(selection);
    }
    selections.sort_by_key(|selection| selection.outer_fold);
    (final_score, selections)
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
    p: &Params,
) -> Vec<RankedFeature> {
    let test_fdr = p.test_fdr;
    let subset_labels: Vec<i8> = rows.iter().map(|&row| labels[row]).collect();
    let tdc = stats::Tdc::training(p.null_target_win_prob);
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
            let qvalues = stats::qvalues(&scores, &subset_labels, tdc);
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

/// Split `rows` into `fold_count` folds, keeping every candidate from one
/// spectrum together.
///
/// Candidates from the same spectrum are not independent observations: they come
/// from one measurement, compete with each other, and a target and its decoy
/// counterpart share almost everything except the label.  Splitting them across
/// folds lets a model train on one candidate of a spectrum and then score its
/// sibling, which is how a spectrum's label reaches its own held-out score.  The
/// reference splits by spectrum for the same reason
/// (`Scores::createXvalSetsBySpectrum`).
///
/// Grouping is by `(source, scan)` so that joined files with colliding scan
/// numbers stay separate.  Ensemble input is the exception: there the same
/// spectrum is deliberately reported by several engines under different sources,
/// so it groups by scan alone.
fn assign_dataset_folds(ds: &Dataset, rows: &[usize], fold_count: u8, seed: u64) -> Vec<u8> {
    let mut by_spectrum: BTreeMap<(u32, i64), Vec<usize>> = BTreeMap::new();
    for &row in rows {
        let key = if ds.ensemble {
            (0, ds.scan[row])
        } else {
            (ds.source[row], ds.scan[row])
        };
        by_spectrum.entry(key).or_default().push(row);
    }
    let mut spectra: Vec<Vec<usize>> = by_spectrum.into_values().collect();
    let mut rng = Rng(seed.max(1));
    for index in (1..spectra.len()).rev() {
        spectra.swap(index, rng.below(index + 1));
    }

    // Spectra hold different numbers of candidates, so deal each to whichever
    // fold is currently smallest rather than round-robin; ties go to the lowest
    // fold index, which keeps the assignment deterministic.
    let mut fold = vec![u8::MAX; ds.n_psm];
    let mut sizes = vec![0usize; fold_count as usize];
    for spectrum in spectra {
        let mut target = 0usize;
        for candidate in 1..sizes.len() {
            if sizes[candidate] < sizes[target] {
                target = candidate;
            }
        }
        sizes[target] += spectrum.len();
        for row in spectrum {
            fold[row] = target as u8;
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
            let (x, dim) = build_matrix_fit(ds, &train_rows, p);
            let ranking = rank_features(&x, dim, &ds.labels, &train_rows, p);
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
    yield_at_fdr(&validation_scores, &validation_labels, p.test_fdr, p)
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
        let (x, dim) = build_matrix_fit(ds, &train_rows, p);
        let ranking = rank_features(&x, dim, &ds.labels, &train_rows, p);
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
        let (qval, pep) = stats::qvalues_and_peps(
            &final_score,
            &ds.labels,
            stats::Tdc::reported(p.null_target_win_prob),
        );
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

    let selected = p.c_alpha.is_none() || p.c_beta.is_none();
    if selected {
        // A private pool keeps thread count under this run's control, so callers that
        // already parallelize across files (the benchmark harness) stay single-threaded.
        let select = || cv_scores_with_selected_c(ds, &fold, p);
        let (final_score, selections) = match rayon::ThreadPoolBuilder::new()
            .num_threads(p.num_threads)
            .build()
        {
            Ok(pool) if p.num_threads > 1 => pool.install(select),
            _ => select(),
        };
        let (qval, pep) = stats::qvalues_and_peps(
            &final_score,
            &ds.labels,
            stats::Tdc::reported(p.null_target_win_prob),
        );
        return Output {
            score: final_score,
            qval,
            pep,
            c_alpha: f64::NAN,
            c_beta: f64::NAN,
            c_selected: true,
            nested_folds: selections,
        };
    }
    let (alpha, beta) = (p.c_alpha.unwrap(), p.c_beta.unwrap());

    let hp = Hp {
        alpha,
        beta,
        maxiter: p.maxiter,
        subset: p.subset_max_train,
        tolerance: p.svm_tolerance,
    };
    let final_score = cv_scores(ds, &fold, p, hp, p.seed);

    #[cfg(feature = "profiling")]
    let _final_context = crate::profile::context(Some("final_psm_statistics"), None, None, None);
    let (qval, pep) = stats::qvalues_and_peps(
        &final_score,
        &ds.labels,
        stats::Tdc::reported(p.null_target_win_prob),
    );
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
    let hp = Hp {
        alpha: output.c_alpha,
        beta: output.c_beta,
        maxiter: p.maxiter,
        subset: p.subset_max_train,
        tolerance: p.svm_tolerance,
    };
    OUTER_FOLDS
        .iter()
        .map(|&test_fold| {
            // Rebuilt through the same seeded path as `cv_scores`, so the report
            // describes the model that actually produced the held-out scores.
            let setup = fold_setup(ds, fold, test_fold, p, p.seed);
            let rt_values = fold_rt_columns(ds, &setup.train_rows, p);
            let normalization = fit_normalization(
                ds,
                &setup.train_rows,
                rt_columns(p, rt_values.as_deref()).as_ref(),
            );
            let mut rng = Rng(setup.seed);
            let model = train_fold(
                &setup.x,
                setup.dim,
                &ds.labels,
                &setup.train_rows,
                &setup.w0,
                p,
                &mut rng,
                hp,
                setup.seed ^ 0xD1B5_4A32_D192_ED03,
                None,
            );
            let weights = match model {
                FoldModel::Svm(weights) => weights,
                FoldModel::Mlp(_) => unreachable!("feature reports require SVM"),
            };
            let dim = setup.dim;
            ExplanationFold {
                weights,
                normalization,
                test_rows: setup.test_rows,
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
            let rt_values = fold_rt_columns(ds, &train_rows, p);
            let normalization = fit_normalization(
                ds,
                &train_rows,
                rt_columns(p, rt_values.as_deref()).as_ref(),
            );
            let (x, dim) =
                transform_matrix(ds, &normalization, rt_columns(p, rt_values.as_deref()).as_ref());
            let ranking = rank_features(&x, dim, &ds.labels, &train_rows, p);
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
                training_null_calibration(&model, &x, dim, &ds.labels, &train_rows);
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

fn target_q01(scores: &[f64], labels: &[i8], p: &Params) -> usize {
    let test_fdr = p.test_fdr;
    let qvalues = stats::qvalues(scores, labels, stats::Tdc::reported(p.null_target_win_prob));
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
    let baseline_q01 = target_q01(&output.score, &ds.labels, p);
    // Report-only descriptive statistics over the whole input; nothing here
    // feeds a model, so a fold-local retention-time alignment is not needed.
    let global_normalization = fit_normalization(ds, &(0..ds.n_psm).collect::<Vec<_>>(), None);
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
        let permuted_q01 = target_q01(&permuted_scores, &ds.labels, p);
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
            exp_mass: vec![0.0; n_psm],
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
        let outer_fold = assign_dataset_folds(&dataset, &all_rows, 3, params.seed);
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

    /// Perturb everything about the held-out fold -- features and labels -- and
    /// require that the model trained to score that fold does not move.  This is
    /// the falsifiable form of "leakage-free": if normalization or the initial
    /// direction were still fitted globally, these assertions would fail.
    #[test]
    fn heldout_rows_cannot_influence_the_model_that_scores_them() {
        let dataset = selection_fixture();
        let params = Params {
            maxiter: 3,
            num_threads: 1,
            c_alpha: Some(C_POS_DEFAULT),
            c_beta: Some(C_NEG_DEFAULT),
            ..Params::default()
        };
        let all_rows: Vec<usize> = (0..dataset.n_psm).collect();
        let fold = assign_dataset_folds(&dataset, &all_rows, 3, params.seed);

        let mut corrupted = selection_fixture();
        for row in all_rows.iter().copied().filter(|&row| fold[row] == 0) {
            corrupted.labels[row] *= -1;
            for feature in 0..corrupted.n_feat {
                corrupted.features[row * corrupted.n_feat + feature] =
                    1_000_000.0 * (feature + 1) as f64;
            }
        }

        let baseline = fold_setup(&dataset, &fold, 0, &params, params.seed);
        let perturbed = fold_setup(&corrupted, &fold, 0, &params, params.seed);

        // Repeat with the label-dependent retention-time alignment switched on:
        // it is the one remaining supervised preprocessing step, so it has to be
        // refitted per fold like everything else.
        {
            let mut with_rt = selection_fixture();
            let mut corrupted_rt = corrupted.clone();
            let rt_params = Params {
                rt: crate::rt::augment(&mut with_rt),
                ..Params {
                    maxiter: params.maxiter,
                    num_threads: 1,
                    c_alpha: params.c_alpha,
                    c_beta: params.c_beta,
                    ..Params::default()
                }
            };
            crate::rt::augment(&mut corrupted_rt);
            let clean = fold_setup(&with_rt, &fold, 0, &rt_params, rt_params.seed);
            let dirty = fold_setup(&corrupted_rt, &fold, 0, &rt_params, rt_params.seed);
            assert_eq!(
                dirty.w0, clean.w0,
                "retention-time features let held-out rows move the initial direction"
            );
            let training_rows = |setup: &FoldSetup| -> Vec<f64> {
                setup
                    .train_rows
                    .iter()
                    .flat_map(|&row| {
                        setup.x[row * setup.dim..(row + 1) * setup.dim]
                            .iter()
                            .copied()
                    })
                    .collect()
            };
            assert_eq!(
                training_rows(&dirty),
                training_rows(&clean),
                "retention-time alignment moved with held-out rows"
            );
        }

        assert_eq!(
            perturbed.w0, baseline.w0,
            "initial direction moved with held-out rows"
        );
        let normalized_train = |setup: &FoldSetup| -> Vec<f64> {
            setup
                .train_rows
                .iter()
                .flat_map(|&row| {
                    setup.x[row * setup.dim..(row + 1) * setup.dim]
                        .iter()
                        .copied()
                })
                .collect()
        };
        assert_eq!(
            normalized_train(&perturbed),
            normalized_train(&baseline),
            "training-row normalization moved with held-out rows"
        );

        let hp = Hp {
            alpha: C_POS_DEFAULT,
            beta: C_NEG_DEFAULT,
            maxiter: params.maxiter,
            subset: params.subset_max_train,
            tolerance: params.svm_tolerance,
        };
        let weights = |ds: &Dataset, setup: &FoldSetup| -> Vec<f64> {
            let mut rng = Rng(setup.seed);
            match train_fold(
                &setup.x,
                setup.dim,
                &ds.labels,
                &setup.train_rows,
                &setup.w0,
                &params,
                &mut rng,
                hp,
                setup.seed ^ 0xD1B5_4A32_D192_ED03,
                None,
            ) {
                FoldModel::Svm(w) => w,
                FoldModel::Mlp(_) => unreachable!(),
            }
        };
        assert_eq!(
            weights(&corrupted, &perturbed),
            weights(&dataset, &baseline),
            "fold model weights moved with held-out rows"
        );
    }

    /// Each fold must fit its own normalization; otherwise the "training-only"
    /// claim would be vacuously satisfied by every fold sharing one transform.
    #[test]
    fn folds_fit_their_own_normalization_and_direction() {
        let dataset = selection_fixture();
        let params = Params {
            maxiter: 2,
            num_threads: 1,
            ..Params::default()
        };
        let all_rows: Vec<usize> = (0..dataset.n_psm).collect();
        let fold = assign_dataset_folds(&dataset, &all_rows, 3, params.seed);
        let first = fold_setup(&dataset, &fold, 0, &params, params.seed);
        let second = fold_setup(&dataset, &fold, 1, &params, params.seed);
        assert_ne!(
            first.x, second.x,
            "two folds produced an identical design matrix, so normalization is not fold-local"
        );
        let global = fit_normalization(&dataset, &all_rows, None);
        let fold_local = fit_normalization(&dataset, &first.train_rows, None);
        assert_ne!(
            global.mean, fold_local.mean,
            "fold normalization matched the all-rows fit"
        );
    }

    /// Scoring is out of fold: every reported score comes from a model whose
    /// training partition excluded that row.
    #[test]
    fn every_row_is_scored_by_a_model_that_excluded_it() {
        let dataset = selection_fixture();
        let params = Params {
            maxiter: 2,
            num_threads: 1,
            ..Params::default()
        };
        let all_rows: Vec<usize> = (0..dataset.n_psm).collect();
        let fold = assign_dataset_folds(&dataset, &all_rows, 3, params.seed);
        let mut covered = vec![0u8; dataset.n_psm];
        for &test in OUTER_FOLDS.iter() {
            let setup = fold_setup(&dataset, &fold, test, &params, params.seed);
            for &row in &setup.test_rows {
                covered[row] += 1;
                assert!(!setup.train_rows.contains(&row));
            }
        }
        assert!(covered.iter().all(|&count| count == 1));
    }

    /// Standardization must only move folds relative to each other: within a
    /// fold it is an increasing affine map, so the order cannot change.
    #[test]
    fn fold_standardization_preserves_within_fold_order() {
        let dataset = selection_fixture();
        let params = Params {
            maxiter: 3,
            num_threads: 1,
            c_alpha: Some(C_POS_DEFAULT),
            c_beta: Some(C_NEG_DEFAULT),
            ..Params::default()
        };
        let all_rows: Vec<usize> = (0..dataset.n_psm).collect();
        let fold = assign_dataset_folds(&dataset, &all_rows, 3, params.seed);
        let hp = Hp {
            alpha: C_POS_DEFAULT,
            beta: C_NEG_DEFAULT,
            maxiter: params.maxiter,
            subset: params.subset_max_train,
            tolerance: params.svm_tolerance,
        };
        for &test in OUTER_FOLDS.iter() {
            let setup = fold_setup(&dataset, &fold, test, &params, params.seed);
            let mut rng = Rng(setup.seed);
            let model = train_fold(
                &setup.x,
                setup.dim,
                &dataset.labels,
                &setup.train_rows,
                &setup.w0,
                &params,
                &mut rng,
                hp,
                setup.seed ^ 0xD1B5_4A32_D192_ED03,
                None,
            );
            let mut raw = vec![0.0; setup.test_rows.len()];
            model.score_rows(&setup.x, setup.dim, &setup.test_rows, &mut raw);
            let standardized = standardized_heldout_scores(
                &model,
                &setup.x,
                setup.dim,
                &dataset.labels,
                &setup.train_rows,
                &setup.test_rows,
            );
            let rank = |values: &[f64]| {
                let mut order: Vec<usize> = (0..values.len()).collect();
                order.sort_by(|&a, &b| values[b].total_cmp(&values[a]));
                order
            };
            assert_eq!(rank(&raw), rank(&standardized), "fold {test} reordered");

            // The calibration really is the training-null location and scale.
            let (mean, sd) = training_null_calibration(
                &model,
                &setup.x,
                setup.dim,
                &dataset.labels,
                &setup.train_rows,
            );
            let decoy_rows: Vec<usize> = setup
                .train_rows
                .iter()
                .copied()
                .filter(|&row| dataset.labels[row] < 0)
                .collect();
            let mut decoy_scores = vec![0.0; decoy_rows.len()];
            model.score_rows(&setup.x, setup.dim, &decoy_rows, &mut decoy_scores);
            let centred: Vec<f64> = decoy_scores
                .iter()
                .map(|score| (score - mean) / sd)
                .collect();
            let (centred_mean, centred_sd) = mean_and_sd(&centred);
            assert!(centred_mean.abs() < 1e-9, "training null not centred");
            assert!((centred_sd - 1.0).abs() < 1e-9, "training null not scaled");
        }
    }

    /// Every candidate of one spectrum -- including a target and its decoy
    /// counterpart -- must land in the same fold, or a spectrum trains the model
    /// that scores it.
    #[test]
    fn all_candidates_of_a_spectrum_share_a_fold() {
        let mut dataset = selection_fixture();
        // Give each spectrum four candidates: two targets, two decoys.
        for row in 0..dataset.n_psm {
            dataset.scan[row] = (row / 4) as i64;
            dataset.labels[row] = if row % 4 < 2 { 1 } else { -1 };
        }
        let rows: Vec<usize> = (0..dataset.n_psm).collect();
        let fold = assign_dataset_folds(&dataset, &rows, 3, 1);
        for row in 0..dataset.n_psm {
            let first = (row / 4) * 4;
            assert_eq!(
                fold[row], fold[first],
                "row {row} split away from its spectrum"
            );
        }
        assert!(fold.iter().all(|&f| f < 3), "every row must get a fold");
        // Folds stay usable: greedy smallest-first keeps them within one
        // spectrum's worth of each other.
        let mut sizes = [0usize; 3];
        for &f in &fold {
            sizes[f as usize] += 1;
        }
        let spread = sizes.iter().max().unwrap() - sizes.iter().min().unwrap();
        assert!(spread <= 4, "fold sizes drifted apart: {sizes:?}");
    }

    /// Joined files reuse scan numbers, so grouping must be per source.
    #[test]
    fn joined_files_with_colliding_scans_are_not_merged_into_one_spectrum() {
        let mut dataset = selection_fixture();
        dataset.source_names = vec!["a.pin".to_string(), "b.pin".to_string()];
        for row in 0..dataset.n_psm {
            dataset.source[row] = (row % 2) as u32;
            dataset.scan[row] = (row / 2) as i64;
        }
        let rows: Vec<usize> = (0..dataset.n_psm).collect();
        let fold = assign_dataset_folds(&dataset, &rows, 3, 1);
        assert!(
            (0..dataset.n_psm / 2).any(|pair| fold[pair * 2] != fold[pair * 2 + 1]),
            "same scan number in different files was treated as one spectrum"
        );
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
