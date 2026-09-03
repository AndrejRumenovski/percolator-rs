//! Retention-time features: predict RT from peptide sequence (per-residue coefficient
//! model), align to observed elution (ScanNr proxy) by least squares on targets, and
//! append the residual magnitude as discriminative features. Correct PSMs elute near
//! their predicted RT; random/decoy matches deviate — so |obs - pred| aids separation.

use crate::pin::Dataset;

/// Per-residue retention coefficients (relative hydrophobicity, Meek/Guo-style scale).
/// Arbitrary units; the model is linearly re-aligned to observed RT anyway.
fn coeff(aa: u8) -> f64 {
    match aa {
        b'W' => 11.0,
        b'F' => 10.5,
        b'L' => 9.6,
        b'I' => 8.4,
        b'M' => 5.8,
        b'V' => 5.0,
        b'Y' => 4.0,
        b'A' => 0.5,
        b'T' => 0.4,
        b'P' => 0.2,
        b'C' => 0.8,
        b'G' => 0.0,
        b'S' => -0.1,
        b'Q' => -0.3,
        b'N' => -0.5,
        b'E' => 1.1,
        b'D' => -0.5,
        b'H' => -1.0,
        b'R' => -1.3,
        b'K' => -1.9,
        _ => 0.0,
    }
}

/// Predicted RT = sum of residue coefficients over the (mod-stripped) core sequence.
fn predict(peptide: &str) -> f64 {
    // strip flanks A.PEP.B -> PEP, then keep only A-Z letters (drop mods/brackets/digits)
    let core = {
        let b = peptide.as_bytes();
        let first = peptide.find('.');
        let last = peptide.rfind('.');
        match (first, last) {
            (Some(a), Some(c)) if c > a => &peptide[a + 1..c],
            _ => {
                let _ = b;
                peptide
            }
        }
    };
    let mut s = 0.0;
    for &c in core.as_bytes() {
        if c.is_ascii_uppercase() {
            s += coeff(c);
        }
    }
    s
}

/// Least-squares fit obs ≈ a*pred + b over the given index set. Returns (a, b).
fn fit(pred: &[f64], obs: &[f64], idx: &[usize]) -> (f64, f64) {
    let n = idx.len() as f64;
    if n < 2.0 {
        return (1.0, 0.0);
    }
    let (mut sx, mut sy, mut sxx, mut sxy) = (0.0, 0.0, 0.0, 0.0);
    for &i in idx {
        let x = pred[i];
        let y = obs[i];
        sx += x;
        sy += y;
        sxx += x * x;
        sxy += x * y;
    }
    let denom = n * sxx - sx * sx;
    if denom.abs() < 1e-9 {
        return (0.0, sy / n);
    }
    let a = (n * sxy - sx * sy) / denom;
    let b = (sy - a * sx) / n;
    (a, b)
}

/// Sequence-predicted and observed retention times kept alongside the feature
/// matrix so that the alignment can be refitted inside each fold.
///
/// The alignment `obs ~ a*pred + b` is fitted on *target* PSMs, which makes it
/// label-dependent preprocessing.  Fitting it once over the whole dataset before
/// folds are formed would carry every held-out label into the model that scores
/// it, so the coefficients are not stored here -- only the two label-free
/// series the fold-local fit consumes.
pub struct Alignment {
    pub predicted: Vec<f64>,
    pub observed: Vec<f64>,
    /// Index of the `rt_abs_error` column; `rt_sq_error` follows it.
    pub first_column: usize,
}

impl Alignment {
    /// Fit on the targets of `fit_rows` and write both residual features for
    /// every row into `columns`, laid out as `[abs, sq]` pairs indexed by row.
    ///
    /// One alignment is fitted per input file: `ScanNr` indexes a run's own
    /// acquisition, so scan numbers from different files are not on a common
    /// scale and must not share a fit.
    pub fn residuals(
        &self,
        labels: &[i8],
        source: &[u32],
        fit_rows: &[usize],
        columns: &mut [f64],
    ) {
        let sources = source.iter().copied().max().map_or(0, |m| m as usize + 1);
        let mut by_source: Vec<Vec<usize>> = vec![Vec::new(); sources];
        for &row in fit_rows {
            if labels[row] > 0 {
                by_source[source[row] as usize].push(row);
            }
        }
        let mut fallback: Vec<Vec<usize>> = vec![Vec::new(); sources];
        for &row in fit_rows {
            if by_source[source[row] as usize].is_empty() {
                fallback[source[row] as usize].push(row);
            }
        }
        let coefficients: Vec<(f64, f64)> = (0..sources)
            .map(|s| {
                let rows = if by_source[s].is_empty() {
                    &fallback[s]
                } else {
                    &by_source[s]
                };
                fit(&self.predicted, &self.observed, rows)
            })
            .collect();
        for row in 0..self.predicted.len() {
            let (a, b) = coefficients[source[row] as usize];
            let residual = self.observed[row] - (a * self.predicted[row] + b);
            columns[row * 2] = residual.abs();
            columns[row * 2 + 1] = residual * residual;
        }
    }
}

/// Reserve the two RT residual feature columns and record the label-free series
/// they are computed from.
///
/// The columns are left at zero: their values depend on which rows are allowed
/// to fit the alignment, and that is only known once folds exist.  Observed
/// elution is proxied by `ScanNr`, which is monotone in retention time within a
/// run but is not a retention time; a PIN's own `retentiontime` column is
/// metadata and is not read here.
pub fn augment(ds: &mut Dataset) -> Option<Alignment> {
    let n = ds.n_psm;
    if n == 0 {
        return None;
    }
    let predicted: Vec<f64> = (0..n).map(|i| predict(&ds.peptide[i])).collect();
    let observed: Vec<f64> = ds.scan.iter().map(|&s| s as f64).collect();

    let old_f = ds.n_feat;
    let new_f = old_f + 2;
    let mut feats = vec![0.0f64; n * new_f];
    for i in 0..n {
        let src = &ds.features[i * old_f..i * old_f + old_f];
        feats[i * new_f..i * new_f + old_f].copy_from_slice(src);
    }
    ds.features = feats;
    ds.n_feat = new_f;
    ds.feature_names.push("rt_abs_error".to_string());
    ds.feature_names.push("rt_sq_error".to_string());
    Some(Alignment {
        predicted,
        observed,
        first_column: old_f,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn fixture() -> Dataset {
        let peptides = [
            "K.AAAA.R", "K.WWWW.R", "K.LLLL.R", "K.GGGG.R", "K.KKKK.R", "K.FFFF.R",
        ];
        let n = peptides.len();
        Dataset {
            feature_names: vec!["f".to_string()],
            n_feat: 1,
            n_psm: n,
            features: vec![0.0; n],
            labels: vec![1, 1, 1, -1, -1, -1],
            spec_id: (0..n).map(|i| format!("s{i}")).collect(),
            scan: vec![10, 90, 70, 50, 20, 80],
            exp_mass: vec![0.0; n],
            peptide: peptides.iter().map(|p| p.to_string()).collect(),
            proteins: (0..n).map(|i| format!("P{i}")).collect(),
            source: vec![0; n],
            source_names: vec!["a.pin".to_string()],
            ensemble: false,
        }
    }

    /// The alignment is fitted on targets, so it must only ever see the rows it
    /// is given -- otherwise held-out labels reach the fold that scores them.
    #[test]
    fn the_alignment_depends_only_on_the_rows_it_is_fitted_on() {
        let mut ds = fixture();
        let alignment = augment(&mut ds).unwrap();
        let mut first = vec![0.0; ds.n_psm * 2];
        alignment.residuals(&ds.labels, &ds.source, &[0, 1, 3, 4], &mut first);
        let mut second = vec![0.0; ds.n_psm * 2];
        alignment.residuals(&ds.labels, &ds.source, &[0, 1, 2, 3, 4, 5], &mut second);
        assert_ne!(
            first, second,
            "fitting on more rows must change the alignment, or the fit ignores its argument"
        );
    }

    /// Scan numbers index each run's own acquisition, so joined files must not
    /// share one alignment.
    #[test]
    fn each_input_file_gets_its_own_alignment() {
        let mut ds = fixture();
        ds.source_names.push("b.pin".to_string());
        for row in 0..ds.n_psm {
            ds.source[row] = (row % 2) as u32;
            // Second file's scans live in a completely different range.
            if row % 2 == 1 {
                ds.scan[row] += 100_000;
            }
        }
        let alignment = augment(&mut ds).unwrap();
        let rows: Vec<usize> = (0..ds.n_psm).collect();
        let mut per_source = vec![0.0; ds.n_psm * 2];
        alignment.residuals(&ds.labels, &ds.source, &rows, &mut per_source);
        let mut pooled = vec![0.0; ds.n_psm * 2];
        alignment.residuals(&ds.labels, &vec![0u32; ds.n_psm], &rows, &mut pooled);
        assert_ne!(
            per_source, pooled,
            "a shared alignment across files went unnoticed"
        );
        assert!(per_source.iter().all(|value| value.is_finite()));
    }

    #[test]
    fn residual_columns_are_finite_and_non_negative() {
        let mut ds = fixture();
        let alignment = augment(&mut ds).unwrap();
        assert_eq!(ds.n_feat, 3);
        assert_eq!(alignment.first_column, 1);
        let mut columns = vec![0.0; ds.n_psm * 2];
        let rows: Vec<usize> = (0..ds.n_psm).collect();
        alignment.residuals(&ds.labels, &ds.source, &rows, &mut columns);
        for (index, value) in columns.iter().enumerate() {
            assert!(
                value.is_finite() && *value >= 0.0,
                "column {index} = {value}"
            );
        }
        for row in 0..ds.n_psm {
            let absolute = columns[row * 2];
            assert!((columns[row * 2 + 1] - absolute * absolute).abs() < 1e-9);
        }
    }
}
