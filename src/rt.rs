//! Retention-time features: predict RT from peptide sequence (per-residue coefficient
//! model), align to observed elution (ScanNr proxy) by least squares on targets, and
//! append the residual magnitude as discriminative features. Correct PSMs elute near
//! their predicted RT; random/decoy matches deviate — so |obs - pred| aids separation.

use crate::pin::Dataset;

/// Per-residue retention coefficients (relative hydrophobicity, Meek/Guo-style scale).
/// Arbitrary units; the model is linearly re-aligned to observed RT anyway.
fn coeff(aa: u8) -> f64 {
    match aa {
        b'W' => 11.0, b'F' => 10.5, b'L' => 9.6, b'I' => 8.4, b'M' => 5.8,
        b'V' => 5.0, b'Y' => 4.0, b'A' => 0.5, b'T' => 0.4, b'P' => 0.2,
        b'C' => 0.8, b'G' => 0.0, b'S' => -0.1, b'Q' => -0.3, b'N' => -0.5,
        b'E' => 1.1, b'D' => -0.5, b'H' => -1.0, b'R' => -1.3, b'K' => -1.9,
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

/// Append RT residual features (rt_abs_error, rt_sq_error) to the dataset's feature matrix.
pub fn augment(ds: &mut Dataset) {
    let n = ds.n_psm;
    if n == 0 {
        return;
    }
    let pred: Vec<f64> = (0..n).map(|i| predict(&ds.peptide[i])).collect();
    let obs: Vec<f64> = ds.scan.iter().map(|&s| s as f64).collect();

    // align on target PSMs (the correct bulk dominates the fit)
    let mut tgt_idx: Vec<usize> = (0..n).filter(|&i| ds.labels[i] > 0).collect();
    if tgt_idx.is_empty() {
        tgt_idx = (0..n).collect();
    }
    let (a, b) = fit(&pred, &obs, &tgt_idx);

    let old_f = ds.n_feat;
    let new_f = old_f + 2;
    let mut feats = vec![0.0f64; n * new_f];
    for i in 0..n {
        let src = &ds.features[i * old_f..i * old_f + old_f];
        let dst = &mut feats[i * new_f..i * new_f + new_f];
        dst[..old_f].copy_from_slice(src);
        let resid = obs[i] - (a * pred[i] + b);
        dst[old_f] = resid.abs();
        dst[old_f + 1] = resid * resid;
    }
    ds.features = feats;
    ds.n_feat = new_f;
    ds.feature_names.push("rt_abs_error".to_string());
    ds.feature_names.push("rt_sq_error".to_string());
}
