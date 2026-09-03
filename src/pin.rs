//! Fast zero-copy parser for Percolator INput (.pin) tab-delimited files.
//! Memory-maps the file and parses over the raw byte buffer with fast-float
//! (correctly-rounded, so numeric results are identical to std parsing).
//!
//! Columns: SpecId, Label, ScanNr, [ExpMass, CalcMass,] <features...>, Peptide, Proteins.
//! ExpMass/CalcMass are recognized by name and excluded from the feature matrix.

use memmap2::Mmap;
use std::cmp::Ordering;
use std::collections::{BTreeMap, BTreeSet};
use std::fs::File;

#[derive(Clone)]
pub struct Dataset {
    pub feature_names: Vec<String>,
    pub n_feat: usize,
    pub n_psm: usize,
    pub features: Vec<f64>, // row-major n_psm * n_feat
    pub labels: Vec<i8>,
    pub spec_id: Vec<String>,
    pub scan: Vec<i64>,
    /// Experimental precursor mass, when the PIN supplies an `ExpMass` column.
    ///
    /// Not a feature: it identifies the precursor. Together with `(source, scan)`
    /// it separates the distinct precursors a single scan can produce, which is
    /// the unit spectrum-level target-decoy competition runs over. `0.0` when the
    /// column is absent, which collapses a scan to one precursor.
    pub exp_mass: Vec<f64>,
    pub peptide: Vec<String>,
    pub proteins: Vec<String>,
    pub source: Vec<u32>,          // index into source_names, per PSM
    pub source_names: Vec<String>, // input file basenames
    pub ensemble: bool,
}

impl Dataset {
    /// The precursor a row was matched to: the unit that target-decoy
    /// competition runs over.
    ///
    /// `(source, scan)` separates joined files that reuse scan numbers, and
    /// `ExpMass` separates the distinct precursors one scan can yield -- the
    /// experimental neutral mass differs between charge assignments of the same
    /// peak, so it stands in for the charge the reference also keys on. Ensemble
    /// input drops `source`, because there the same spectrum is deliberately
    /// reported by several engines.
    pub fn spectrum_key(&self, row: usize) -> (u32, i64, u64) {
        let source = if self.ensemble { 0 } else { self.source[row] };
        let mass = self.exp_mass[row];
        // -0.0 and 0.0 are the same precursor.
        let bits = if mass == 0.0 { 0 } else { mass.to_bits() };
        (source, self.scan[row], bits)
    }
}

/// Compare two rows by their complete parsed content, with spectrum and PSM
/// identity first and the target/decoy label only as the final fallback.
///
/// The label fallback is not a target-first or decoy-first competition rule: it
/// is reached only when every label-free identifier and every feature agree.
/// Its sole purpose is to give two otherwise equal records one canonical
/// in-memory order before floating-point accumulation and model fitting.
fn compare_rows(dataset: &Dataset, left: usize, right: usize) -> Ordering {
    dataset.scan[left]
        .cmp(&dataset.scan[right])
        .then_with(|| dataset.exp_mass[left].total_cmp(&dataset.exp_mass[right]))
        .then_with(|| dataset.spec_id[left].cmp(&dataset.spec_id[right]))
        .then_with(|| dataset.peptide[left].cmp(&dataset.peptide[right]))
        .then_with(|| dataset.proteins[left].cmp(&dataset.proteins[right]))
        .then_with(|| {
            for (left, right) in dataset.row(left).iter().zip(dataset.row(right)) {
                let ordering = left.total_cmp(right);
                if ordering != Ordering::Equal {
                    return ordering;
                }
            }
            Ordering::Equal
        })
        .then_with(|| dataset.labels[left].cmp(&dataset.labels[right]))
}

/// Put one parsed source into a content-defined row order.
fn canonicalize_rows(dataset: &mut Dataset) {
    let mut order: Vec<usize> = (0..dataset.n_psm).collect();
    order.sort_unstable_by(|&left, &right| compare_rows(dataset, left, right));

    let old_features = std::mem::take(&mut dataset.features);
    let old_labels = std::mem::take(&mut dataset.labels);
    let old_spec_id = std::mem::take(&mut dataset.spec_id);
    let old_scan = std::mem::take(&mut dataset.scan);
    let old_exp_mass = std::mem::take(&mut dataset.exp_mass);
    let old_peptide = std::mem::take(&mut dataset.peptide);
    let old_proteins = std::mem::take(&mut dataset.proteins);

    dataset.features = Vec::with_capacity(old_features.len());
    dataset.labels = Vec::with_capacity(dataset.n_psm);
    dataset.spec_id = Vec::with_capacity(dataset.n_psm);
    dataset.scan = Vec::with_capacity(dataset.n_psm);
    dataset.exp_mass = Vec::with_capacity(dataset.n_psm);
    dataset.peptide = Vec::with_capacity(dataset.n_psm);
    dataset.proteins = Vec::with_capacity(dataset.n_psm);
    for row in order {
        let start = row * dataset.n_feat;
        dataset
            .features
            .extend_from_slice(&old_features[start..start + dataset.n_feat]);
        dataset.labels.push(old_labels[row]);
        dataset.spec_id.push(old_spec_id[row].clone());
        dataset.scan.push(old_scan[row]);
        dataset.exp_mass.push(old_exp_mass[row]);
        dataset.peptide.push(old_peptide[row].clone());
        dataset.proteins.push(old_proteins[row].clone());
    }
    dataset.source = vec![0; dataset.n_psm];
}

fn compare_parts(left: &Dataset, right: &Dataset) -> Ordering {
    left.source_names[0]
        .cmp(&right.source_names[0])
        .then_with(|| left.n_psm.cmp(&right.n_psm))
        .then_with(|| {
            for row in 0..left.n_psm.min(right.n_psm) {
                let ordering = left.scan[row]
                    .cmp(&right.scan[row])
                    .then_with(|| left.exp_mass[row].total_cmp(&right.exp_mass[row]))
                    .then_with(|| left.spec_id[row].cmp(&right.spec_id[row]))
                    .then_with(|| left.peptide[row].cmp(&right.peptide[row]))
                    .then_with(|| left.proteins[row].cmp(&right.proteins[row]))
                    .then_with(|| {
                        for (a, b) in left.row(row).iter().zip(right.row(row)) {
                            let ordering = a.total_cmp(b);
                            if ordering != Ordering::Equal {
                                return ordering;
                            }
                        }
                        Ordering::Equal
                    })
                    .then_with(|| left.labels[row].cmp(&right.labels[row]));
                if ordering != Ordering::Equal {
                    return ordering;
                }
            }
            Ordering::Equal
        })
}

/// Concatenate datasets into one pooled dataset for cross-run joint training.
/// All inputs must share the same feature layout (same n_feat).
///
/// A joined input is a multiset of named source files and parsed PSM records.
/// Argument positions and row positions carry no scientific information. Parts
/// and their rows are therefore canonicalized before numeric source indices are
/// assigned. This makes those indices stable for fold grouping and seeded tie
/// draws, and gives normalization/training a stable accumulation order.
pub fn merge(mut parts: Vec<Dataset>) -> Dataset {
    assert!(!parts.is_empty());
    for part in &parts {
        assert_eq!(
            part.source_names.len(),
            1,
            "joined input parts must each describe one source"
        );
    }
    for part in &mut parts {
        canonicalize_rows(part);
    }
    parts.sort_by(compare_parts);

    let n_feat = parts[0].n_feat;
    let mut out = parts.remove(0); // its source is already all-0, source_names = [file0]
    for mut p in parts {
        assert_eq!(
            p.n_feat, n_feat,
            "cannot join files with differing feature columns"
        );
        let sidx = out.source_names.len() as u32;
        out.source_names.append(&mut p.source_names);
        out.features.append(&mut p.features);
        out.labels.append(&mut p.labels);
        out.spec_id.append(&mut p.spec_id);
        out.scan.append(&mut p.scan);
        out.exp_mass.append(&mut p.exp_mass);
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
            return Err(format!(
                "ensemble engine names must be non-empty and unique (got '{name}')"
            ));
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
        exp_mass: Vec::with_capacity(total_rows),
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
            out.exp_mass.push(part.exp_mass[row]);
            out.peptide.push(part.peptide[row].clone());
            out.proteins.push(part.proteins[row].clone());
            out.source.push(engine_idx as u32);
        }
        feature_offset += part.n_feat;
    }

    // Cross-engine agreement features.
    //
    // Both keys are deliberately label-free.  Keying the per-candidate count on
    // `(ScanNr, Label, Peptide)` -- the previous behaviour -- built a training
    // feature out of the labels of every row in the file, including the rows
    // that later become held-out.  Whether two engines reported the same
    // candidate for the same spectrum is a property of the searches, and a
    // peptide sequence already determines whether it came from the target or the
    // decoy database, so the label adds nothing except the leak.
    let mut spectrum_engines: BTreeMap<i64, BTreeSet<u32>> = BTreeMap::new();
    let mut psm_engines: BTreeMap<(i64, &str), BTreeSet<u32>> = BTreeMap::new();
    for row in 0..out.n_psm {
        spectrum_engines
            .entry(out.scan[row])
            .or_default()
            .insert(out.source[row]);
        psm_engines
            .entry((out.scan[row], out.peptide[row].as_str()))
            .or_default()
            .insert(out.source[row]);
    }
    let spectrum_count = n_feat - 2;
    let psm_count = n_feat - 1;
    let counts: Vec<(f64, f64)> = (0..out.n_psm)
        .map(|row| {
            (
                spectrum_engines[&out.scan[row]].len() as f64,
                psm_engines[&(out.scan[row], out.peptide[row].as_str())].len() as f64,
            )
        })
        .collect();
    for (row, (spectrum, psm)) in counts.into_iter().enumerate() {
        out.features[row * n_feat + spectrum_count] = spectrum;
        out.features[row * n_feat + psm_count] = psm;
    }
    Ok(out)
}

impl Dataset {
    #[inline]
    pub fn row(&self, i: usize) -> &[f64] {
        &self.features[i * self.n_feat..(i + 1) * self.n_feat]
    }
}

/// Metadata column names the PIN format allows between `Label` and the first
/// feature, matched case-insensitively.
///
/// The reference consumes a *contiguous prefix* of these and starts the feature
/// block at the first unrecognized header
/// (`SetHandler::getOptionalFields`).  percolator-rs previously took every
/// column between `ScanNr` and `Peptide` as a feature except `ExpMass` and
/// `CalcMass`, which meant a Sage `FileName` column was parsed as a number and a
/// raw `retentiontime` was trained on as if it were a search score.
const OPTIONAL_COLUMNS: [&[u8]; 7] = [
    b"ScanNr",
    b"ExpMass",
    b"CalcMass",
    b"rt",
    b"retentiontime",
    b"FileName",
    b"SpectraFile",
];

fn is_optional_column(name: &[u8]) -> bool {
    OPTIONAL_COLUMNS
        .iter()
        .any(|known| name.eq_ignore_ascii_case(known))
}

fn invalid_data(message: String) -> std::io::Error {
    std::io::Error::new(std::io::ErrorKind::InvalidData, message)
}

fn show(field: &[u8]) -> String {
    let text = String::from_utf8_lossy(field);
    if text.chars().count() > 40 {
        format!("{}...", text.chars().take(40).collect::<String>())
    } else {
        text.into_owned()
    }
}

/// Strict base-10 integer. Anything else is a malformed file, not a zero.
#[inline]
fn atoi(b: &[u8]) -> Option<i64> {
    let (negative, digits) = match b.first() {
        Some(b'-') => (true, &b[1..]),
        Some(b'+') => (false, &b[1..]),
        _ => (false, b),
    };
    if digits.is_empty() || !digits.iter().all(u8::is_ascii_digit) {
        return None;
    }
    let mut value: i64 = 0;
    for &c in digits {
        value = value.checked_mul(10)?.checked_add((c - b'0') as i64)?;
    }
    Some(if negative { -value } else { value })
}

/// Strict finite float.
///
/// The parser used to fall back to zero for anything it could not read and to
/// accept NaN and infinity unchecked. A silently zeroed feature is
/// indistinguishable from a real measurement of zero, and a non-finite feature
/// propagates through normalization into scores, q-values and the sort order.
/// Both now stop the run with a located diagnostic.
#[inline]
fn parse_f64(b: &[u8]) -> Option<f64> {
    let value: f64 = fast_float::parse(b).ok()?;
    value.is_finite().then_some(value)
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

#[allow(clippy::needless_range_loop)]
pub fn parse(path: &str) -> std::io::Result<Dataset> {
    #[cfg(feature = "profiling")]
    let parse_start = std::time::Instant::now();
    #[cfg(feature = "profiling")]
    let mmap_start = std::time::Instant::now();
    let file = File::open(path)?;
    let mmap = unsafe { Mmap::map(&file)? };
    let data: &[u8] = &mmap;
    #[cfg(feature = "profiling")]
    crate::profile::record(
        "parser",
        "mmap_setup",
        mmap_start.elapsed(),
        None,
        Some(data.len() as u64),
    );

    // line iterator over the byte buffer (handles trailing \r)
    let mut lines = data.split(|&c| c == b'\n').map(|l| {
        if l.last() == Some(&b'\r') {
            &l[..l.len() - 1]
        } else {
            l
        }
    });

    #[cfg(feature = "profiling")]
    let header_start = std::time::Instant::now();
    let header = lines.next().unwrap_or(&[]);
    let mut hfields: Vec<&[u8]> = Vec::new();
    split_fields(header, &mut hfields);
    let required = |name: &str, bytes: &[u8]| -> std::io::Result<usize> {
        hfields
            .iter()
            .position(|c| c.eq_ignore_ascii_case(bytes))
            .ok_or_else(|| invalid_data(format!("{path}: PIN header has no {name} column")))
    };
    let idx_label = required("Label", b"Label")?;
    let idx_scan = required("ScanNr", b"ScanNr")?;
    let idx_pep = required("Peptide", b"Peptide")?;
    let idx_exp_mass = hfields
        .iter()
        .position(|c| c.eq_ignore_ascii_case(b"ExpMass"));

    let mut feature_start = idx_label + 1;
    while feature_start < idx_pep && is_optional_column(hfields[feature_start]) {
        feature_start += 1;
    }
    if idx_scan >= feature_start {
        return Err(invalid_data(format!(
            "{path}: ScanNr must sit in the metadata columns that follow Label, \
             but a feature column appears before it"
        )));
    }
    let feat_cols: Vec<usize> = (feature_start..idx_pep).collect();
    let feature_names: Vec<String> = feat_cols
        .iter()
        .map(|&j| String::from_utf8_lossy(hfields[j]).into_owned())
        .collect();
    let n_feat = feat_cols.len();
    if n_feat == 0 {
        return Err(invalid_data(format!(
            "{path}: PIN header has no feature columns between the metadata block and Peptide"
        )));
    }
    #[cfg(feature = "profiling")]
    crate::profile::record(
        "parser",
        "header_and_feature_names",
        header_start.elapsed(),
        Some(hfields.len() as u64),
        Some(header.len() as u64),
    );

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
        exp_mass: Vec::with_capacity(approx_rows),
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
    #[cfg(feature = "profiling")]
    let rows_start = std::time::Instant::now();
    #[cfg(feature = "profiling")]
    let mut split_time = std::time::Duration::ZERO;
    #[cfg(feature = "profiling")]
    let mut numeric_time = std::time::Duration::ZERO;
    #[cfg(feature = "profiling")]
    let mut string_copy_time = std::time::Duration::ZERO;
    #[cfg(feature = "profiling")]
    let mut float_fields = 0u64;
    #[cfg(feature = "profiling")]
    let mut copied_string_bytes = 0u64;
    #[cfg(feature = "profiling")]
    let mut string_allocations = 0u64;
    // Header consumed above, so the first row of `lines` is file line 2.
    for (offset, line) in lines.enumerate() {
        let line_number = offset + 2;
        if line.is_empty() {
            continue;
        }
        if first {
            first = false;
            if line.starts_with(b"DefaultDirection") {
                continue;
            }
        }
        #[cfg(feature = "profiling")]
        let split_start = std::time::Instant::now();
        split_fields(line, &mut fields);
        #[cfg(feature = "profiling")]
        {
            split_time += split_start.elapsed();
        }
        if fields.len() <= idx_pep {
            return Err(invalid_data(format!(
                "{path}:{line_number}: row has {} fields but the header needs at least {}",
                fields.len(),
                idx_pep + 1
            )));
        }
        #[cfg(feature = "profiling")]
        let spec_string_start = std::time::Instant::now();
        ds.spec_id
            .push(String::from_utf8_lossy(fields[0]).into_owned());
        #[cfg(feature = "profiling")]
        {
            string_copy_time += spec_string_start.elapsed();
        }
        #[cfg(feature = "profiling")]
        let numeric_start = std::time::Instant::now();
        let raw_label = atoi(fields[idx_label]).ok_or_else(|| {
            invalid_data(format!(
                "{path}:{line_number}: Label is not an integer: '{}'",
                show(fields[idx_label])
            ))
        })?;
        ds.labels.push(if raw_label > 0 { 1 } else { -1 });
        ds.scan.push(atoi(fields[idx_scan]).ok_or_else(|| {
            invalid_data(format!(
                "{path}:{line_number}: ScanNr is not an integer: '{}'",
                show(fields[idx_scan])
            ))
        })?);
        ds.exp_mass.push(match idx_exp_mass {
            Some(column) => parse_f64(fields[column]).ok_or_else(|| {
                invalid_data(format!(
                    "{path}:{line_number}: ExpMass is not a finite number: '{}'",
                    show(fields[column])
                ))
            })?,
            None => 0.0,
        });
        for (feature, &j) in feat_cols.iter().enumerate() {
            ds.features.push(parse_f64(fields[j]).ok_or_else(|| {
                invalid_data(format!(
                    "{path}:{line_number}: feature '{}' is not a finite number: '{}'",
                    ds.feature_names[feature],
                    show(fields[j])
                ))
            })?);
        }
        #[cfg(feature = "profiling")]
        {
            numeric_time += numeric_start.elapsed();
            float_fields += feat_cols.len() as u64;
        }
        #[cfg(feature = "profiling")]
        let string_start = std::time::Instant::now();
        ds.peptide
            .push(String::from_utf8_lossy(fields[idx_pep]).into_owned());
        let prot = if fields.len() > idx_pep + 1 {
            // proteins may contain tabs; rejoin from original line span
            let start_ptr = fields[idx_pep + 1].as_ptr() as usize - line.as_ptr() as usize;
            String::from_utf8_lossy(&line[start_ptr..]).into_owned()
        } else {
            String::new()
        };
        #[cfg(feature = "profiling")]
        {
            string_copy_time += string_start.elapsed();
            copied_string_bytes += (fields[0].len() + fields[idx_pep].len() + prot.len()) as u64;
            string_allocations += if prot.is_empty() { 2 } else { 3 };
        }
        ds.proteins.push(prot);
    }
    ds.n_psm = ds.labels.len();
    ds.source = vec![0u32; ds.n_psm];
    #[cfg(feature = "profiling")]
    {
        let row_time = rows_start.elapsed();
        crate::profile::record(
            "parser",
            "row_loading_total",
            row_time,
            Some(ds.n_psm as u64),
            Some(data.len().saturating_sub(header.len()) as u64),
        );
        crate::profile::record(
            "parser",
            "field_splitting",
            split_time,
            Some(ds.n_psm as u64),
            None,
        );
        crate::profile::record(
            "parser",
            "numeric_and_float_parsing",
            numeric_time,
            Some(float_fields),
            None,
        );
        crate::profile::record(
            "parser",
            "string_allocation_and_copy",
            string_copy_time,
            Some(string_allocations),
            Some(copied_string_bytes),
        );
        crate::profile::allocation_site(
            "pin::parse row strings",
            string_allocations,
            copied_string_bytes,
        );
        let vector_bytes = ds.features.capacity() as u64 * std::mem::size_of::<f64>() as u64
            + ds.labels.capacity() as u64 * std::mem::size_of::<i8>() as u64
            + ds.scan.capacity() as u64 * std::mem::size_of::<i64>() as u64
            + ds.spec_id.capacity() as u64 * std::mem::size_of::<String>() as u64
            + ds.peptide.capacity() as u64 * std::mem::size_of::<String>() as u64
            + ds.proteins.capacity() as u64 * std::mem::size_of::<String>() as u64
            + ds.source.capacity() as u64 * std::mem::size_of::<u32>() as u64;
        crate::profile::allocation_site("pin::parse column vectors", 7, vector_bytes);
        crate::profile::record(
            "parser",
            "pin_parse_total",
            parse_start.elapsed(),
            Some(ds.n_psm as u64),
            Some(data.len() as u64),
        );
    }
    Ok(ds)
}

#[cfg(test)]
mod tests {
    use super::{merge, merge_ensemble, Dataset};

    fn dataset(feature: &str, rows: &[(i64, i8, &str, f64)]) -> Dataset {
        Dataset {
            feature_names: vec![feature.to_string()],
            n_feat: 1,
            n_psm: rows.len(),
            features: rows.iter().map(|row| row.3).collect(),
            labels: rows.iter().map(|row| row.1).collect(),
            spec_id: (0..rows.len()).map(|row| format!("psm{row}")).collect(),
            scan: rows.iter().map(|row| row.0).collect(),
            exp_mass: vec![0.0; rows.len()],
            peptide: rows.iter().map(|row| row.2.to_string()).collect(),
            proteins: vec![String::new(); rows.len()],
            source: vec![0; rows.len()],
            source_names: vec!["input.pin".to_string()],
            ensemble: false,
        }
    }

    fn named_dataset(name: &str, rows: &[(i64, i8, &str, f64)]) -> Dataset {
        let mut result = dataset("score", rows);
        result.source_names = vec![name.to_string()];
        result.spec_id = rows
            .iter()
            .map(|row| format!("{}:{}", name, row.2))
            .collect();
        result
    }

    fn joined_snapshot(dataset: &Dataset) -> Vec<String> {
        (0..dataset.n_psm)
            .map(|row| {
                format!(
                    "{}\t{}\t{}\t{}\t{}\t{}\t{:?}",
                    dataset.source_names[dataset.source[row] as usize],
                    dataset.spec_id[row],
                    dataset.labels[row],
                    dataset.scan[row],
                    dataset.peptide[row],
                    dataset.proteins[row],
                    dataset.row(row),
                )
            })
            .collect()
    }

    /// Joined inputs are a multiset of named scientific records. Neither the
    /// argument order nor the row layout is evidence, so both must materialize
    /// the same internal dataset before folds or floating-point fits are built.
    #[test]
    fn joined_merge_is_canonical_under_file_and_row_permutations() {
        let alpha_rows = [
            (11, 1, "K.ALPHA_TARGET.R", 4.0),
            (11, -1, "K.ALPHA_DECOY.R", 4.0),
            (
                12,
                1,
                "K.ALPHA_NEAR.R",
                f64::from_bits(4.0f64.to_bits() + 1),
            ),
        ];
        let beta_rows = [
            (21, -1, "K.BETA_DECOY.R", -3.0),
            (22, 1, "K.BETA_TARGET.R", 2.0),
        ];
        let canonical = merge(vec![
            named_dataset("alpha.pin", &alpha_rows),
            named_dataset("beta.pin", &beta_rows),
        ]);

        let mut reversed_alpha = alpha_rows;
        reversed_alpha.reverse();
        let mut reversed_beta = beta_rows;
        reversed_beta.reverse();
        let permuted = merge(vec![
            named_dataset("beta.pin", &reversed_beta),
            named_dataset("alpha.pin", &reversed_alpha),
        ]);

        assert_eq!(
            joined_snapshot(&permuted),
            joined_snapshot(&canonical),
            "scientifically equivalent joined permutations produced different internal rows"
        );
    }

    /// **Frozen failure case (independent audit M5).**
    ///
    /// The cross-engine agreement features are built once over every row, before
    /// folds exist.  They must therefore be a function of the search results
    /// alone: flipping any labels, including the labels of rows that later land
    /// in a held-out fold, must leave every feature value untouched.
    #[test]
    fn ensemble_features_do_not_depend_on_any_label() {
        let rows_a = [
            (10, 1, "A.PEPTIDE.B", 4.0),
            (10, -1, "A.YTIDEPEP.B", 1.0),
            (11, 1, "A.OTHER.B", 2.0),
        ];
        let rows_b = [
            (10, 1, "A.PEPTIDE.B", 0.01),
            (10, -1, "A.YTIDEPEP.B", 0.9),
            (12, 1, "A.THIRD.B", 0.5),
        ];
        let clean = merge_ensemble(
            vec![dataset("xcorr", &rows_a), dataset("exact_p_value", &rows_b)],
            vec!["comet".to_string(), "tide".to_string()],
        )
        .unwrap();

        // Every assignment of labels, held fixed in features.
        for mask in 0..(1u32 << 6) {
            let flip = |rows: &[(i64, i8, &'static str, f64)], offset: u32| {
                rows.iter()
                    .enumerate()
                    .map(|(index, row)| {
                        let bit = 1u32 << (offset + index as u32);
                        (
                            row.0,
                            if mask & bit != 0 { -row.1 } else { row.1 },
                            row.2,
                            row.3,
                        )
                    })
                    .collect::<Vec<_>>()
            };
            let flipped = merge_ensemble(
                vec![
                    dataset("xcorr", &flip(&rows_a, 0)),
                    dataset("exact_p_value", &flip(&rows_b, 3)),
                ],
                vec!["comet".to_string(), "tide".to_string()],
            )
            .unwrap();
            assert_eq!(
                flipped.features, clean.features,
                "label mask {mask:#08b} changed an ensemble feature"
            );
        }
    }

    #[test]
    fn ensemble_namespaces_features_and_counts_exact_cross_engine_support() {
        let comet = dataset(
            "xcorr",
            &[(10, 1, "A.PEPTIDE.B", 4.0), (11, 1, "A.OTHER.B", 2.0)],
        );
        let tide = dataset(
            "exact_p_value",
            &[(10, 1, "A.PEPTIDE.B", 0.01), (10, -1, "A.DECOY.B", 0.9)],
        );
        let out = merge_ensemble(
            vec![comet, tide],
            vec!["comet".to_string(), "tide".to_string()],
        )
        .unwrap();

        assert_eq!(
            out.feature_names,
            vec![
                "ensemble_engine=comet",
                "ensemble_engine=tide",
                "ensemble_engine=comet;feature=xcorr",
                "ensemble_engine=tide;feature=exact_p_value",
                "ensemble_spectrum_engine_count",
                "ensemble_psm_engine_count",
            ]
        );
        assert_eq!(out.n_psm, 4);
        assert_eq!(out.row(0), &[1.0, 0.0, 4.0, 0.0, 2.0, 2.0]);
        assert_eq!(out.row(1), &[1.0, 0.0, 2.0, 0.0, 1.0, 1.0]);
        assert_eq!(out.row(2), &[0.0, 1.0, 0.0, 0.01, 2.0, 2.0]);
        // Same spectrum, but a distinct decoy peptide: spectrum support remains 2;
        // exact PSM support correctly remains 1.
        assert_eq!(out.row(3), &[0.0, 1.0, 0.0, 0.9, 2.0, 1.0]);
    }

    // ----- parser fails closed ------------------------------------------------

    fn write_pin(body: &str) -> std::path::PathBuf {
        use std::io::Write;
        let path = std::env::temp_dir().join(format!(
            "percolator-rs-pin-{}-{:p}.pin",
            std::process::id(),
            body.as_ptr()
        ));
        let mut file = std::fs::File::create(&path).unwrap();
        file.write_all(body.as_bytes()).unwrap();
        path
    }

    const GOOD: &str = "SpecId\tLabel\tScanNr\tf1\tf2\tPeptide\tProteins\n\
                        s1\t1\t10\t1.5\t2.5\tK.PEP.R\tP1\n\
                        s2\t-1\t11\t0.5\t1.0\tK.DEP.R\tDECOY_P1\n";

    fn parse_body(body: &str) -> std::io::Result<Dataset> {
        let path = write_pin(body);
        let result = super::parse(path.to_str().unwrap());
        std::fs::remove_file(&path).ok();
        result
    }

    fn rejects(body: &str, needle: &str) {
        match parse_body(body) {
            Ok(_) => panic!("parser accepted malformed input: {needle}"),
            Err(error) => {
                let message = error.to_string();
                assert!(
                    message.contains(needle),
                    "diagnostic '{message}' does not mention '{needle}'"
                );
            }
        }
    }

    #[test]
    fn a_well_formed_pin_still_parses() {
        let ds = parse_body(GOOD).expect("valid PIN should parse");
        assert_eq!(ds.n_psm, 2);
        assert_eq!(ds.n_feat, 2);
        assert_eq!(ds.row(0), &[1.5, 2.5]);
        assert_eq!(ds.labels, vec![1, -1]);
        assert_eq!(ds.scan, vec![10, 11]);
    }

    /// A silently zeroed feature is indistinguishable from a real measurement of
    /// zero, so the run must stop instead.
    #[test]
    fn a_malformed_feature_is_rejected() {
        rejects(
            &GOOD.replace("\t1.5\t", "\tnot-a-number\t"),
            "not a finite number",
        );
    }

    #[test]
    fn a_missing_feature_is_rejected() {
        rejects(&GOOD.replace("\t1.5\t", "\t\t"), "not a finite number");
    }

    #[test]
    fn non_finite_features_are_rejected() {
        for text in ["nan", "NaN", "inf", "-inf", "Infinity"] {
            rejects(
                &GOOD.replace("\t1.5\t", &format!("\t{text}\t")),
                "not a finite number",
            );
        }
    }

    #[test]
    fn the_diagnostic_locates_the_row_and_column() {
        let body = GOOD.replace("\t0.5\t", "\tbad\t");
        let message = match parse_body(&body) {
            Ok(_) => panic!("parser accepted a malformed feature"),
            Err(error) => error.to_string(),
        };
        assert!(message.contains(":3:"), "no line number in '{message}'");
        assert!(message.contains("'f1'"), "no column name in '{message}'");
        assert!(
            message.contains("'bad'"),
            "no offending text in '{message}'"
        );
    }

    #[test]
    fn a_malformed_label_is_rejected() {
        rejects(
            &GOOD.replace("s1\t1\t", "s1\ttarget\t"),
            "Label is not an integer",
        );
    }

    #[test]
    fn a_malformed_scan_number_is_rejected() {
        rejects(
            &GOOD.replace("\t10\t", "\t10a\t"),
            "ScanNr is not an integer",
        );
    }

    #[test]
    fn a_short_row_is_rejected_rather_than_skipped() {
        rejects(&format!("{GOOD}s3\t1\t12\n"), "fields but the header needs");
    }

    #[test]
    fn a_missing_required_column_is_rejected() {
        rejects(
            "SpecId\tScanNr\tf1\tPeptide\tProteins\ns1\t10\t1.0\tK.P.R\tP1\n",
            "no Label column",
        );
    }

    #[test]
    fn metadata_columns_are_excluded_whatever_their_case() {
        let body = "SpecId\tLabel\tScanNr\texpmass\tCALCMASS\tf1\tPeptide\tProteins\n\
                    s1\t1\t10\t500.1\t500.2\t1.5\tK.PEP.R\tP1\n";
        let ds = parse_body(body).expect("valid PIN should parse");
        assert_eq!(ds.feature_names, vec!["f1"]);
        assert_eq!(ds.row(0), &[1.5]);
    }

    /// Sage writes a spectrum filename and a raw retention time in the metadata
    /// block. Neither is a search score; the filename is not even a number.
    #[test]
    fn the_metadata_prefix_stops_at_the_first_feature() {
        let body = "SpecId\tLabel\tScanNr\tExpMass\tCalcMass\tFileName\tretentiontime\trank\tscore\tPeptide\tProteins\n\
                    s1\t1\t10\t500.1\t500.2\trun_a.mzML\t12.5\t1\t3.5\tK.PEP.R\tP1\n";
        let ds = parse_body(body).expect("valid PIN should parse");
        assert_eq!(ds.feature_names, vec!["rank", "score"]);
        assert_eq!(ds.row(0), &[1.0, 3.5]);
    }

    /// A metadata name appearing after the feature block starts is a feature, not
    /// metadata: the prefix is contiguous.
    #[test]
    fn a_metadata_name_after_the_prefix_stays_a_feature() {
        let body = "SpecId\tLabel\tScanNr\tscore\tretentiontime\tPeptide\tProteins\n\
                    s1\t1\t10\t3.5\t12.5\tK.PEP.R\tP1\n";
        let ds = parse_body(body).expect("valid PIN should parse");
        assert_eq!(ds.feature_names, vec!["score", "retentiontime"]);
    }

    #[test]
    fn a_pin_without_feature_columns_is_rejected() {
        rejects(
            "SpecId\tLabel\tScanNr\tExpMass\tPeptide\tProteins\ns1\t1\t10\t500.1\tK.P.R\tP1\n",
            "no feature columns",
        );
    }
}
