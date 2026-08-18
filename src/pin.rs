//! Fast zero-copy parser for Percolator INput (.pin) tab-delimited files.
//! Memory-maps the file and parses over the raw byte buffer with fast-float
//! (correctly-rounded, so numeric results are identical to std parsing).
//!
//! Columns: SpecId, Label, ScanNr, [ExpMass, CalcMass,] <features...>, Peptide, Proteins.
//! ExpMass/CalcMass are recognized by name and excluded from the feature matrix.

use memmap2::Mmap;
use std::collections::{BTreeMap, BTreeSet};
use std::fs::File;

pub struct Dataset {
    pub feature_names: Vec<String>,
    pub n_feat: usize,
    pub n_psm: usize,
    pub features: Vec<f64>, // row-major n_psm * n_feat
    pub labels: Vec<i8>,
    pub spec_id: Vec<String>,
    pub scan: Vec<i64>,
    pub peptide: Vec<String>,
    pub proteins: Vec<String>,
    pub source: Vec<u32>,          // index into source_names, per PSM
    pub source_names: Vec<String>, // input file basenames
    pub ensemble: bool,
}

/// Concatenate datasets into one pooled dataset for cross-run joint training.
/// All inputs must share the same feature layout (same n_feat).
pub fn merge(mut parts: Vec<Dataset>) -> Dataset {
    assert!(!parts.is_empty());
    let n_feat = parts[0].n_feat;
    let mut out = parts.remove(0); // its source is already all-0, source_names = [file0]
    for mut p in parts {
        assert_eq!(p.n_feat, n_feat, "cannot join files with differing feature columns");
        let sidx = out.source_names.len() as u32;
        out.source_names.append(&mut p.source_names);
        out.features.append(&mut p.features);
        out.labels.append(&mut p.labels);
        out.spec_id.append(&mut p.spec_id);
        out.scan.append(&mut p.scan);
        out.peptide.append(&mut p.peptide);
        out.proteins.append(&mut p.proteins);
        for _ in 0..p.n_psm {
            out.source.push(sidx);
        }
    }
    out.n_psm = out.labels.len();
    out
}

/// Combine results for the same spectra from distinct search engines.
///
/// Unlike [`merge`], engines need not expose the same PIN features.  Their feature
/// spaces are kept separate (rather than assuming that equal column names have the
/// same scale or meaning), with a one-hot engine indicator and two agreement
/// features appended.  Agreement is deliberately exact: a PSM is supported only
/// when another engine reports the same ScanNr, label, and modified peptide string.
/// This makes the mode appropriate for inputs produced from the same raw run.
pub fn merge_ensemble(parts: Vec<Dataset>, engine_names: Vec<String>) -> Result<Dataset, String> {
    if parts.len() < 2 {
        return Err("--ensemble requires at least two ENGINE=PIN inputs".to_string());
    }
    if parts.len() != engine_names.len() {
        return Err("internal ensemble input mismatch".to_string());
    }
    let mut seen = BTreeSet::new();
    for name in &engine_names {
        if name.is_empty() || !seen.insert(name.clone()) {
            return Err(format!("ensemble engine names must be non-empty and unique (got '{name}')"));
        }
    }

    let total_rows: usize = parts.iter().map(|part| part.n_psm).sum();
    let feature_names: Vec<String> = engine_names
        .iter()
        .map(|name| format!("ensemble_engine={name}"))
        .chain(parts.iter().zip(&engine_names).flat_map(|(part, engine)| {
            part.feature_names
                .iter()
                .map(move |feature| format!("ensemble_engine={engine};feature={feature}"))
        }))
        .chain([
            "ensemble_spectrum_engine_count".to_string(),
            "ensemble_psm_engine_count".to_string(),
        ])
        .collect();
    let n_feat = feature_names.len();
    let mut out = Dataset {
        feature_names,
        n_feat,
        n_psm: total_rows,
        features: Vec::with_capacity(total_rows * n_feat),
        labels: Vec::with_capacity(total_rows),
        spec_id: Vec::with_capacity(total_rows),
        scan: Vec::with_capacity(total_rows),
        peptide: Vec::with_capacity(total_rows),
        proteins: Vec::with_capacity(total_rows),
        source: Vec::with_capacity(total_rows),
        source_names: engine_names.clone(),
        ensemble: true,
    };

    let mut feature_offset = engine_names.len();
    for (engine_idx, part) in parts.into_iter().enumerate() {
        for row in 0..part.n_psm {
            let row_start = out.features.len();
            out.features.resize(row_start + n_feat, 0.0);
            out.features[row_start + engine_idx] = 1.0;
            let source = part.row(row);
            out.features[row_start + feature_offset..row_start + feature_offset + part.n_feat]
                .copy_from_slice(source);
            out.labels.push(part.labels[row]);
            out.spec_id.push(part.spec_id[row].clone());
            out.scan.push(part.scan[row]);
            out.peptide.push(part.peptide[row].clone());
            out.proteins.push(part.proteins[row].clone());
            out.source.push(engine_idx as u32);
        }
        feature_offset += part.n_feat;
    }

    let mut spectrum_engines: BTreeMap<i64, BTreeSet<u32>> = BTreeMap::new();
    let mut psm_engines: BTreeMap<(i64, i8, String), BTreeSet<u32>> = BTreeMap::new();
    for row in 0..out.n_psm {
        spectrum_engines.entry(out.scan[row]).or_default().insert(out.source[row]);
        psm_engines
            .entry((out.scan[row], out.labels[row], out.peptide[row].clone()))
            .or_default()
            .insert(out.source[row]);
    }
    let spectrum_count = n_feat - 2;
    let psm_count = n_feat - 1;
    for row in 0..out.n_psm {
        out.features[row * n_feat + spectrum_count] = spectrum_engines[&out.scan[row]].len() as f64;
        out.features[row * n_feat + psm_count] = psm_engines
            [&(out.scan[row], out.labels[row], out.peptide[row].clone())]
            .len() as f64;
    }
    Ok(out)
}

impl Dataset {
    #[inline]
    pub fn row(&self, i: usize) -> &[f64] {
        &self.features[i * self.n_feat..(i + 1) * self.n_feat]
    }
}

fn is_excluded(name: &[u8]) -> bool {
    name == b"ExpMass" || name == b"CalcMass"
}

#[inline]
fn atoi(b: &[u8]) -> i64 {
    let mut i = 0;
    let mut neg = false;
    if !b.is_empty() && (b[0] == b'-' || b[0] == b'+') {
        neg = b[0] == b'-';
        i = 1;
    }
    let mut v: i64 = 0;
    while i < b.len() {
        let c = b[i];
        if c.is_ascii_digit() {
            v = v * 10 + (c - b'0') as i64;
        }
        i += 1;
    }
    if neg { -v } else { v }
}

#[inline]
fn parse_f64(b: &[u8]) -> f64 {
    fast_float::parse(b).unwrap_or(0.0)
}

/// Iterate tab-delimited field byte-slices of a line into `out` (indices reused).
#[inline]
fn split_fields<'a>(line: &'a [u8], out: &mut Vec<&'a [u8]>) {
    out.clear();
    let mut start = 0;
    for i in 0..line.len() {
        if line[i] == b'\t' {
            out.push(&line[start..i]);
            start = i + 1;
        }
    }
    out.push(&line[start..]);
}

pub fn parse(path: &str) -> std::io::Result<Dataset> {
    let file = File::open(path)?;
    let mmap = unsafe { Mmap::map(&file)? };
    let data: &[u8] = &mmap;

    // line iterator over the byte buffer (handles trailing \r)
    let mut lines = data.split(|&c| c == b'\n').map(|l| {
        if l.last() == Some(&b'\r') { &l[..l.len() - 1] } else { l }
    });

    let header = lines.next().unwrap_or(&[]);
    let mut hfields: Vec<&[u8]> = Vec::new();
    split_fields(header, &mut hfields);
    let idx_label = hfields.iter().position(|c| c.eq_ignore_ascii_case(b"Label")).expect("no Label column");
    let idx_scan = hfields.iter().position(|c| c.eq_ignore_ascii_case(b"ScanNr")).expect("no ScanNr column");
    let idx_pep = hfields.iter().position(|c| c.eq_ignore_ascii_case(b"Peptide")).expect("no Peptide column");

    let mut feat_cols: Vec<usize> = Vec::new();
    let mut feature_names: Vec<String> = Vec::new();
    for j in (idx_scan + 1)..idx_pep {
        if !is_excluded(hfields[j]) {
            feat_cols.push(j);
            feature_names.push(String::from_utf8_lossy(hfields[j]).into_owned());
        }
    }
    let n_feat = feat_cols.len();

    // rough row estimate for pre-allocation
    let approx_rows = data.len() / (header.len().max(1) + 1) + 16;
    let mut ds = Dataset {
        feature_names,
        n_feat,
        n_psm: 0,
        features: Vec::with_capacity(approx_rows * n_feat),
        labels: Vec::with_capacity(approx_rows),
        spec_id: Vec::with_capacity(approx_rows),
        scan: Vec::with_capacity(approx_rows),
        peptide: Vec::with_capacity(approx_rows),
        proteins: Vec::with_capacity(approx_rows),
        source: Vec::new(),
        source_names: vec![std::path::Path::new(path)
            .file_name()
            .map(|s| s.to_string_lossy().into_owned())
            .unwrap_or_else(|| path.to_string())],
        ensemble: false,
    };

    let mut fields: Vec<&[u8]> = Vec::with_capacity(hfields.len());
    let mut first = true;
    for line in lines {
        if line.is_empty() {
            continue;
        }
        if first {
            first = false;
            if line.starts_with(b"DefaultDirection") {
                continue;
            }
        }
        split_fields(line, &mut fields);
        if fields.len() <= idx_pep {
            continue;
        }
        let label: i8 = if atoi(fields[idx_label]) > 0 { 1 } else { -1 };
        ds.spec_id.push(String::from_utf8_lossy(fields[0]).into_owned());
        ds.labels.push(label);
        ds.scan.push(atoi(fields[idx_scan]));
        for &j in &feat_cols {
            ds.features.push(parse_f64(fields[j]));
        }
        ds.peptide.push(String::from_utf8_lossy(fields[idx_pep]).into_owned());
        let prot = if fields.len() > idx_pep + 1 {
            // proteins may contain tabs; rejoin from original line span
            let start_ptr = fields[idx_pep + 1].as_ptr() as usize - line.as_ptr() as usize;
            String::from_utf8_lossy(&line[start_ptr..]).into_owned()
        } else {
            String::new()
        };
        ds.proteins.push(prot);
    }
    ds.n_psm = ds.labels.len();
    ds.source = vec![0u32; ds.n_psm];
    Ok(ds)
}

#[cfg(test)]
mod tests {
    use super::{merge_ensemble, Dataset};

    fn dataset(feature: &str, rows: &[(i64, i8, &str, f64)]) -> Dataset {
        Dataset {
            feature_names: vec![feature.to_string()],
            n_feat: 1,
            n_psm: rows.len(),
            features: rows.iter().map(|row| row.3).collect(),
            labels: rows.iter().map(|row| row.1).collect(),
            spec_id: (0..rows.len()).map(|row| format!("psm{row}")).collect(),
            scan: rows.iter().map(|row| row.0).collect(),
            peptide: rows.iter().map(|row| row.2.to_string()).collect(),
            proteins: vec![String::new(); rows.len()],
            source: vec![0; rows.len()],
            source_names: vec!["input.pin".to_string()],
            ensemble: false,
        }
    }

    #[test]
    fn ensemble_namespaces_features_and_counts_exact_cross_engine_support() {
        let comet = dataset("xcorr", &[(10, 1, "A.PEPTIDE.B", 4.0), (11, 1, "A.OTHER.B", 2.0)]);
        let tide = dataset("exact_p_value", &[(10, 1, "A.PEPTIDE.B", 0.01), (10, -1, "A.DECOY.B", 0.9)]);
        let out = merge_ensemble(vec![comet, tide], vec!["comet".to_string(), "tide".to_string()]).unwrap();

        assert_eq!(out.feature_names, vec![
            "ensemble_engine=comet",
            "ensemble_engine=tide",
            "ensemble_engine=comet;feature=xcorr",
            "ensemble_engine=tide;feature=exact_p_value",
            "ensemble_spectrum_engine_count",
            "ensemble_psm_engine_count",
        ]);
        assert_eq!(out.n_psm, 4);
        assert_eq!(out.row(0), &[1.0, 0.0, 4.0, 0.0, 2.0, 2.0]);
        assert_eq!(out.row(1), &[1.0, 0.0, 2.0, 0.0, 1.0, 1.0]);
        assert_eq!(out.row(2), &[0.0, 1.0, 0.0, 0.01, 2.0, 2.0]);
        // Same spectrum, but a distinct decoy peptide: spectrum support remains 2;
        // exact PSM support correctly remains 1.
        assert_eq!(out.row(3), &[0.0, 1.0, 0.0, 0.9, 2.0, 1.0]);
    }
}
