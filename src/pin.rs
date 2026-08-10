//! Fast zero-copy parser for Percolator INput (.pin) tab-delimited files.
//! Memory-maps the file and parses over the raw byte buffer with fast-float
//! (correctly-rounded, so numeric results are identical to std parsing).
//!
//! Columns: SpecId, Label, ScanNr, [ExpMass, CalcMass,] <features...>, Peptide, Proteins.
//! ExpMass/CalcMass are recognized by name and excluded from the feature matrix.

use memmap2::Mmap;
use std::fs::File;

pub struct Dataset {
    #[allow(dead_code)] // retained for diagnostics / future feature reporting
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
