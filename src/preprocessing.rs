//! Fold-local feature normalization and retention-time materialization.

use crate::percolator::Params;
use crate::pin::Dataset;

/// Per-feature centering and scaling learned from a training partition.
/// Kept separately so explanatory reports can convert SVM coefficients back to
/// the original PIN units.
pub(crate) struct Normalization {
    pub(crate) mean: Vec<f64>,
    pub(crate) std: Vec<f64>,
}

/// Row features with the fold's retention-time residuals patched in.
#[inline]
fn feature_row<'a>(
    ds: &'a Dataset,
    row: usize,
    rt: Option<&RtColumns<'_>>,
    scratch: &'a mut [f64],
) -> &'a [f64] {
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
pub(crate) struct RtColumns<'a> {
    first_column: usize,
    values: &'a [f64],
}

pub(crate) fn fit_normalization(
    ds: &Dataset,
    fit_rows: &[usize],
    rt: Option<&RtColumns<'_>>,
) -> Normalization {
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

pub(crate) fn transform_matrix(
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
pub(crate) fn build_matrix_fit(ds: &Dataset, fit_rows: &[usize], p: &Params) -> (Vec<f64>, usize) {
    #[cfg(feature = "profiling")]
    let _preprocessing =
        crate::profile::Scope::with_elements("preprocessing", "preprocessing_total", ds.n_psm);
    let rt_values = fold_rt_columns(ds, fit_rows, p);
    let rt = rt_columns(p, rt_values.as_deref());
    #[cfg(feature = "profiling")]
    let _normalization =
        crate::profile::Scope::with_elements("normalization", "normalization_total", ds.n_psm);
    let normalization = fit_normalization(ds, fit_rows, rt.as_ref());
    transform_matrix(ds, &normalization, rt.as_ref())
}

/// Refit the retention-time alignment inside `fit_rows` and materialize its two
/// residual columns for every row.
pub(crate) fn fold_rt_columns(ds: &Dataset, fit_rows: &[usize], p: &Params) -> Option<Vec<f64>> {
    let alignment = p.rt.as_ref()?;
    #[cfg(feature = "profiling")]
    let _rt =
        crate::profile::Scope::with_elements("preprocessing", "rt_fold_preprocessing", ds.n_psm);
    let mut values = vec![0.0f64; ds.n_psm * 2];
    alignment.residuals(&ds.labels, &ds.source, fit_rows, &mut values);
    Some(values)
}

pub(crate) fn rt_columns<'a>(p: &Params, values: Option<&'a [f64]>) -> Option<RtColumns<'a>> {
    Some(RtColumns {
        first_column: p.rt.as_ref()?.first_column,
        values: values?,
    })
}
