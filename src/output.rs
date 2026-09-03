//! Stable tabular output contracts and serialization.

#[cfg(feature = "profiling")]
use crate::profile;
use crate::{percolator, protein};
use std::borrow::Cow;
use std::fs::File;
use std::io::{BufWriter, Write};

pub struct Row<'a> {
    id: Cow<'a, str>,
    score: f64,
    q: f64,
    pep: f64,
    peptide: &'a str,
    proteins: &'a str,
    // `sort_unstable` has size-dependent implementations. Preserve the
    // original 96-byte owned-row layout so equal-score output order remains
    // byte-identical while the text itself is borrowed.
    _sort_layout_padding: [u8; 16],
}

impl<'a> Row<'a> {
    pub fn new(
        id: Cow<'a, str>,
        score: f64,
        q: f64,
        pep: f64,
        peptide: &'a str,
        proteins: &'a str,
    ) -> Self {
        Self {
            id,
            score,
            q,
            pep,
            peptide,
            proteins,
            _sort_layout_padding: [0; 16],
        }
    }

    pub fn q_value(&self) -> f64 {
        self.q
    }

    #[cfg(feature = "profiling")]
    pub fn owned_id_capacity(&self) -> u64 {
        match &self.id {
            Cow::Owned(id) => id.capacity() as u64,
            Cow::Borrowed(_) => 0,
        }
    }
}

#[cfg(feature = "profiling")]
#[derive(Default)]
struct WriteCounters {
    calls: std::sync::atomic::AtomicU64,
    bytes: std::sync::atomic::AtomicU64,
    duration_ns: std::sync::atomic::AtomicU64,
}

#[cfg(feature = "profiling")]
struct ProfiledWriter<W> {
    inner: W,
    counters: std::sync::Arc<WriteCounters>,
}

#[cfg(feature = "profiling")]
impl<W: Write> Write for ProfiledWriter<W> {
    fn write(&mut self, buffer: &[u8]) -> std::io::Result<usize> {
        let start = std::time::Instant::now();
        let result = self.inner.write(buffer);
        let elapsed = start.elapsed().as_nanos().min(u64::MAX as u128) as u64;
        self.counters
            .duration_ns
            .fetch_add(elapsed, std::sync::atomic::Ordering::Relaxed);
        self.counters
            .calls
            .fetch_add(1, std::sync::atomic::Ordering::Relaxed);
        if let Ok(bytes) = result {
            self.counters
                .bytes
                .fetch_add(bytes as u64, std::sync::atomic::Ordering::Relaxed);
        }
        result
    }

    fn flush(&mut self) -> std::io::Result<()> {
        let start = std::time::Instant::now();
        let result = self.inner.flush();
        let elapsed = start.elapsed().as_nanos().min(u64::MAX as u128) as u64;
        self.counters
            .duration_ns
            .fetch_add(elapsed, std::sync::atomic::Ordering::Relaxed);
        self.counters
            .calls
            .fetch_add(1, std::sync::atomic::Ordering::Relaxed);
        result
    }
}

fn write_fixed_6<W: Write>(writer: &mut W, value: f64) -> std::io::Result<()> {
    const SCALE: f64 = 1_000_000.0;
    let scaled = value.abs() * SCALE;
    if scaled.is_finite() && scaled < u64::MAX as f64 - 1.0 {
        let lower = scaled.floor();
        let fraction = scaled - lower;
        // Multiplication can move a value very slightly across the half-way
        // boundary. Let the standard formatter handle every ambiguous case.
        let error_bound = scaled * (2.0 * f64::EPSILON) + 1e-12;
        if (fraction - 0.5).abs() > error_bound {
            let rounded = if fraction < 0.5 { lower } else { lower + 1.0 } as u64;
            let mut buffer = [0u8; 32];
            let mut cursor = buffer.len();
            let mut fraction = rounded % 1_000_000;
            for _ in 0..6 {
                cursor -= 1;
                buffer[cursor] = b'0' + (fraction % 10) as u8;
                fraction /= 10;
            }
            cursor -= 1;
            buffer[cursor] = b'.';
            let mut whole = rounded / 1_000_000;
            loop {
                cursor -= 1;
                buffer[cursor] = b'0' + (whole % 10) as u8;
                whole /= 10;
                if whole == 0 {
                    break;
                }
            }
            if value.is_sign_negative() {
                cursor -= 1;
                buffer[cursor] = b'-';
            }
            return writer.write_all(&buffer[cursor..]);
        }
    }
    write!(writer, "{value:.6}")
}

pub fn write_results(path: &str, mut rows: Vec<Row<'_>>) -> std::io::Result<()> {
    #[cfg(feature = "profiling")]
    let sort_start = std::time::Instant::now();
    rows.sort_unstable_by(|x, y| {
        y.score
            .partial_cmp(&x.score)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    #[cfg(feature = "profiling")]
    profile::record(
        "sort",
        "result_row_score_order",
        sort_start.elapsed(),
        Some(rows.len() as u64),
        None,
    );
    #[cfg(feature = "profiling")]
    let create_start = std::time::Instant::now();
    let f = File::create(path)?;
    #[cfg(feature = "profiling")]
    profile::record(
        "file_io",
        "result_file_create",
        create_start.elapsed(),
        Some(1),
        None,
    );
    #[cfg(feature = "profiling")]
    let counters = std::sync::Arc::new(WriteCounters::default());
    #[cfg(feature = "profiling")]
    let f = ProfiledWriter {
        inner: f,
        counters: std::sync::Arc::clone(&counters),
    };
    #[cfg(feature = "profiling")]
    let serialization_start = std::time::Instant::now();
    #[cfg(feature = "profiling")]
    let row_count = rows.len();
    let mut w = BufWriter::with_capacity(1 << 20, f);
    writeln!(
        w,
        "PSMId\tscore\tq-value\tposterior_error_prob\tpeptide\tproteinIds"
    )?;
    for r in rows {
        w.write_all(r.id.as_bytes())?;
        w.write_all(b"\t")?;
        write_fixed_6(&mut w, r.score)?;
        w.write_all(b"\t")?;
        write_fixed_6(&mut w, r.q)?;
        w.write_all(b"\t")?;
        write_fixed_6(&mut w, r.pep)?;
        w.write_all(b"\t")?;
        w.write_all(r.peptide.as_bytes())?;
        w.write_all(b"\t")?;
        w.write_all(r.proteins.as_bytes())?;
        w.write_all(b"\n")?;
    }
    #[cfg(feature = "profiling")]
    {
        use std::sync::atomic::Ordering;
        w.flush()?;
        drop(w);
        let total = serialization_start.elapsed();
        let write_ns = counters.duration_ns.load(Ordering::Relaxed);
        profile::record(
            "file_io",
            "result_file_write",
            std::time::Duration::from_nanos(write_ns),
            Some(counters.calls.load(Ordering::Relaxed)),
            Some(counters.bytes.load(Ordering::Relaxed)),
        );
        profile::record(
            "serialization",
            "result_format_and_buffer",
            total.saturating_sub(std::time::Duration::from_nanos(write_ns)),
            Some(row_count as u64),
            None,
        );
    }
    Ok(())
}

pub fn write_feature_report(path: &str, report: &percolator::FeatureReport) -> std::io::Result<()> {
    let file = File::create(path)?;
    let mut writer = BufWriter::with_capacity(1 << 16, file);
    writeln!(writer, "# feature_report_version=1")?;
    writeln!(
        writer,
        "# model=linear_svm; coefficients are means across the three out-of-fold models"
    )?;
    writeln!(
        writer,
        "# baseline_target_psms_q<0.01={}",
        report.baseline_q01
    )?;
    writeln!(
        writer,
        "# permutation=deterministic within each held-out fold; models held fixed (no retraining)"
    )?;
    writeln!(
        writer,
        "feature_index\tfeature\traw_weight\traw_weight_fold_sd\tstandardized_effect\tstandardized_effect_fold_sd\tlabel_correlation\tfeature_mean\tfeature_std\tselected_folds\tpermutation_q01_drop\tpermuted_target_psms_q<0.01"
    )?;
    for feature in &report.features {
        writeln!(
            writer,
            "{}\t{}\t{:.8}\t{:.8}\t{:.8}\t{:.8}\t{:.8}\t{:.8}\t{:.8}\t{}\t{}\t{}",
            feature.index,
            feature.name,
            feature.raw_weight,
            feature.raw_weight_sd,
            feature.standardized_effect,
            feature.standardized_effect_sd,
            feature.label_correlation,
            feature.mean,
            feature.std,
            feature.selected_folds,
            feature.permutation_q01_drop,
            feature.permuted_q01,
        )?;
    }
    Ok(())
}

pub fn write_proteins(
    path: &str,
    groups: &[protein::ProtGroup],
    want_decoy: bool,
) -> std::io::Result<()> {
    #[cfg(feature = "profiling")]
    let create_start = std::time::Instant::now();
    let f = File::create(path)?;
    #[cfg(feature = "profiling")]
    profile::record(
        "file_io",
        "protein_file_create",
        create_start.elapsed(),
        Some(1),
        None,
    );
    #[cfg(feature = "profiling")]
    let counters = std::sync::Arc::new(WriteCounters::default());
    #[cfg(feature = "profiling")]
    let f = ProfiledWriter {
        inner: f,
        counters: std::sync::Arc::clone(&counters),
    };
    #[cfg(feature = "profiling")]
    let serialization_start = std::time::Instant::now();
    let mut w = BufWriter::with_capacity(1 << 20, f);
    writeln!(
        w,
        "ProteinGroupId\tq-value\tposterior_error_prob\tscore\tnumPeptides\tproteinIds"
    )?;
    for g in groups
        .iter()
        .filter(|g| g.picked && g.is_decoy == want_decoy)
    {
        // `NA` where the selected method estimates no protein-level posterior.
        // Picked-protein FDR is a cumulative estimate; it has no posterior to
        // report, and the best peptide's PEP is not one.
        let mut pep = String::from("NA");
        if let Some(value) = g.pep {
            pep = format!("{value:.6}");
        }
        writeln!(
            w,
            "{}\t{:.6}\t{}\t{:.6}\t{}\t{}",
            g.proteins.first().map(|s| s.as_str()).unwrap_or(""),
            g.qval,
            pep,
            g.score,
            g.n_peptides,
            g.proteins.join(",")
        )?;
    }
    #[cfg(feature = "profiling")]
    {
        use std::sync::atomic::Ordering;
        w.flush()?;
        drop(w);
        let total = serialization_start.elapsed();
        let write_ns = counters.duration_ns.load(Ordering::Relaxed);
        profile::record(
            "file_io",
            "protein_file_write",
            std::time::Duration::from_nanos(write_ns),
            Some(counters.calls.load(Ordering::Relaxed)),
            Some(counters.bytes.load(Ordering::Relaxed)),
        );
        profile::record(
            "serialization",
            "protein_format_and_buffer",
            total.saturating_sub(std::time::Duration::from_nanos(write_ns)),
            Some(groups.len() as u64),
            None,
        );
    }
    Ok(())
}

#[cfg(test)]
mod output_tests {
    use super::{write_fixed_6, Row};

    fn fast(value: f64) -> String {
        let mut output = Vec::new();
        write_fixed_6(&mut output, value).unwrap();
        String::from_utf8(output).unwrap()
    }

    #[test]
    fn fixed_6_matches_standard_formatter() {
        let edge_cases = [
            0.0,
            -0.0,
            0.00000049,
            -0.00000049,
            0.0000005,
            -0.0000005,
            0.9999995,
            -0.9999995,
            1.23456789,
            -123_456.789_012_3,
            f64::INFINITY,
            f64::NEG_INFINITY,
            f64::NAN,
        ];
        for value in edge_cases {
            assert_eq!(fast(value), format!("{value:.6}"), "value={value:?}");
        }

        let mut state = 0x9e37_79b9_7f4a_7c15u64;
        for _ in 0..100_000 {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            let value = f64::from_bits(state);
            assert_eq!(fast(value), format!("{value:.6}"), "bits={state:#018x}");
        }
    }

    #[test]
    fn borrowed_row_preserves_sort_layout_size() {
        assert_eq!(std::mem::size_of::<Row<'_>>(), 96);
    }
}
