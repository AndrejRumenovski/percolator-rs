//! Reusable composition of rescoring, reported-list, peptide, and protein stages.

use crate::{competition, output, peptide, percolator, pin, protein, protein_bayes, stats};
use std::borrow::Cow;

/// Run the fold-isolated learner while retaining the existing profiling stage.
pub fn rescore(ds: &pin::Dataset, params: &percolator::Params) -> percolator::Output {
    #[cfg(feature = "profiling")]
    let _rescoring = crate::profile::Scope::new("stage", "rescoring");
    percolator::run(ds, params)
}

/// Materialized PSM and peptide report rows plus the identities needed by
/// downstream protein inference.
pub struct Reports<'a> {
    pub target_psms: Vec<output::Row<'a>>,
    pub decoy_psms: Vec<output::Row<'a>>,
    pub target_peptides: Vec<output::Row<'a>>,
    pub decoy_peptides: Vec<output::Row<'a>>,
    pub reported_indices: Vec<usize>,
    pub peptides: peptide::Scores,
    pub target_psms_q01: usize,
    pub target_peptides_q01: usize,
}

fn select_reported_indices(
    ds: &pin::Dataset,
    scores: &[f64],
    psm_competition: bool,
    ensemble: bool,
    seed: u64,
) -> Vec<usize> {
    #[cfg(feature = "profiling")]
    let _selection = crate::profile::Scope::with_elements(
        "competition",
        "psm_competition_and_selection",
        ds.n_psm,
    );
    if psm_competition {
        competition::winner_indices(ds, scores, seed)
    } else if ensemble {
        // Engine reports of the same candidate are training observations, not
        // separate discoveries. Retain the best score per exact candidate.
        let mut best: std::collections::BTreeMap<(i64, i8, String), usize> =
            std::collections::BTreeMap::new();
        for index in 0..ds.n_psm {
            let key = (ds.scan[index], ds.labels[index], ds.peptide[index].clone());
            match best.get(&key) {
                Some(&previous) if scores[previous] >= scores[index] => {}
                _ => {
                    best.insert(key, index);
                }
            }
        }
        best.into_values().collect()
    } else {
        (0..ds.n_psm).collect()
    }
}

/// Apply the reporting policy, recalculate statistics on the surviving PSMs,
/// and derive peptide-level results and output rows.
pub fn build_reports<'a>(
    ds: &'a pin::Dataset,
    rescoring: &percolator::Output,
    params: &percolator::Params,
    psm_competition: bool,
    ensemble: bool,
) -> Reports<'a> {
    #[cfg(feature = "profiling")]
    let psm_context = crate::profile::context(Some("psm_level_processing"), None, None, None);
    #[cfg(feature = "profiling")]
    let psm_processing = crate::profile::Scope::new("stage", "psm_level_processing");

    let reported_indices =
        select_reported_indices(ds, &rescoring.score, psm_competition, ensemble, params.seed);
    // Statistics belong to the list that is actually reported. Whenever
    // competition or ensemble deduplication has removed rows, the q-values
    // computed over the full training list no longer describe it.
    let (reported_qvals, reported_peps) = if reported_indices.len() == ds.n_psm {
        (rescoring.qval.clone(), rescoring.pep.clone())
    } else {
        let reported_scores: Vec<f64> = reported_indices
            .iter()
            .map(|&index| rescoring.score[index])
            .collect();
        let reported_labels: Vec<i8> = reported_indices
            .iter()
            .map(|&index| ds.labels[index])
            .collect();
        stats::qvalues_and_peps(
            &reported_scores,
            &reported_labels,
            stats::Tdc::reported(params.null_target_win_prob),
        )
    };

    let target_capacity = reported_indices
        .iter()
        .filter(|&&index| ds.labels[index] > 0)
        .count();
    let mut target_psms = Vec::with_capacity(target_capacity);
    let mut decoy_psms = Vec::with_capacity(reported_indices.len() - target_capacity);
    #[cfg(feature = "profiling")]
    let mut psm_row_string_bytes = 0u64;
    #[cfg(feature = "profiling")]
    let psm_rows_start = std::time::Instant::now();
    for (output_index, &index) in reported_indices.iter().enumerate() {
        let row = output::Row::new(
            if ensemble {
                Cow::Owned(format!(
                    "{}:{}",
                    ds.source_names[ds.source[index] as usize], ds.spec_id[index]
                ))
            } else {
                Cow::Borrowed(&ds.spec_id[index])
            },
            rescoring.score[index],
            reported_qvals[output_index],
            reported_peps[output_index],
            &ds.peptide[index],
            &ds.proteins[index],
        );
        #[cfg(feature = "profiling")]
        {
            psm_row_string_bytes += row.owned_id_capacity();
        }
        if ds.labels[index] > 0 {
            target_psms.push(row);
        } else {
            decoy_psms.push(row);
        }
    }
    #[cfg(feature = "profiling")]
    {
        crate::profile::record(
            "materialization",
            "psm_row_construction",
            psm_rows_start.elapsed(),
            Some(reported_indices.len() as u64),
            None,
        );
        crate::profile::allocation_site(
            "main::psm output row vectors",
            2,
            ((target_psms.capacity() + decoy_psms.capacity()) * std::mem::size_of::<output::Row>())
                as u64,
        );
        crate::profile::allocation_site(
            "main::psm output row strings",
            u64::from(ensemble) * reported_indices.len() as u64,
            psm_row_string_bytes,
        );
    }
    let target_psms_q01 = target_psms
        .iter()
        .filter(|row| row.q_value() < 0.01)
        .count();
    #[cfg(feature = "profiling")]
    drop(psm_processing);
    #[cfg(feature = "profiling")]
    drop(psm_context);

    #[cfg(feature = "profiling")]
    let peptide_context =
        crate::profile::context(Some("peptide_level_processing"), None, None, None);
    #[cfg(feature = "profiling")]
    let peptide_processing = crate::profile::Scope::new("stage", "peptide_level_processing");
    let peptides = peptide::score(
        ds,
        &reported_indices,
        &rescoring.score,
        params.null_target_win_prob,
    );

    let peptide_target_capacity = peptides
        .indices
        .iter()
        .filter(|&&index| ds.labels[index] > 0)
        .count();
    let mut target_peptides = Vec::with_capacity(peptide_target_capacity);
    let mut decoy_peptides = Vec::with_capacity(peptides.indices.len() - peptide_target_capacity);
    #[cfg(feature = "profiling")]
    let peptide_rows_start = std::time::Instant::now();
    for (peptide_index, &psm_index) in peptides.indices.iter().enumerate() {
        let row = output::Row::new(
            Cow::Borrowed(&ds.spec_id[psm_index]),
            peptides.scores[peptide_index],
            peptides.q_values[peptide_index],
            peptides.peps[peptide_index],
            &ds.peptide[psm_index],
            &ds.proteins[psm_index],
        );
        if ds.labels[psm_index] > 0 {
            target_peptides.push(row);
        } else {
            decoy_peptides.push(row);
        }
    }
    #[cfg(feature = "profiling")]
    {
        crate::profile::record(
            "materialization",
            "peptide_row_construction",
            peptide_rows_start.elapsed(),
            Some(peptides.indices.len() as u64),
            None,
        );
        crate::profile::allocation_site(
            "main::peptide output row vectors",
            2,
            ((target_peptides.capacity() + decoy_peptides.capacity())
                * std::mem::size_of::<output::Row>()) as u64,
        );
    }
    let target_peptides_q01 = target_peptides
        .iter()
        .filter(|row| row.q_value() < 0.01)
        .count();
    #[cfg(feature = "profiling")]
    drop(peptide_processing);
    #[cfg(feature = "profiling")]
    drop(peptide_context);

    Reports {
        target_psms,
        decoy_psms,
        target_peptides,
        decoy_peptides,
        reported_indices,
        peptides,
        target_psms_q01,
        target_peptides_q01,
    }
}

/// Protein inference method selected for a pipeline run.
pub enum ProteinMethod<'a> {
    Picked,
    Bayesian(&'a protein_bayes::Params),
}

/// Protein groups together with comparison counts and optional Bayesian
/// diagnostics needed by the command-line report.
pub struct ProteinResults {
    pub groups: Vec<protein::ProtGroup>,
    pub picked_q01: usize,
    pub classic_q01: usize,
    pub bayesian_diagnostics: Option<protein_bayes::Diagnostics>,
}

/// Infer proteins from the peptide result, preserving complete peptide mapping.
pub fn infer_proteins(
    ds: &pin::Dataset,
    reported_indices: &[usize],
    peptides: &peptide::Scores,
    seed: u64,
    method: ProteinMethod<'_>,
) -> ProteinResults {
    #[cfg(feature = "profiling")]
    let entries_start = std::time::Instant::now();
    let entries = peptide::protein_entries(ds, reported_indices, peptides);
    #[cfg(feature = "profiling")]
    {
        crate::profile::record(
            "protein_inference",
            "protein_entry_materialization",
            entries_start.elapsed(),
            Some(entries.len() as u64),
            None,
        );
        crate::profile::allocation_site(
            "main::protein inference entries",
            (entries.len() + 1) as u64,
            (entries.capacity() * std::mem::size_of::<(f64, f64, String)>()) as u64
                + entries
                    .iter()
                    .map(|entry| entry.2.capacity() as u64)
                    .sum::<u64>(),
        );
    }
    let picked_groups = protein::infer(&entries, seed);
    let picked_q01 = picked_groups
        .iter()
        .filter(|group| !group.is_decoy && group.picked && group.qval < 0.01)
        .count();
    let classic_q01 = protein::classic_target_q01(&picked_groups);
    let (groups, bayesian_diagnostics) = match method {
        ProteinMethod::Picked => (picked_groups, None),
        ProteinMethod::Bayesian(params) => {
            let result = protein_bayes::infer(&entries, params);
            (result.groups, Some(result.diagnostics))
        }
    };
    ProteinResults {
        groups,
        picked_q01,
        classic_q01,
        bayesian_diagnostics,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn dataset() -> pin::Dataset {
        pin::Dataset {
            feature_names: vec!["score".to_string()],
            n_feat: 1,
            n_psm: 3,
            features: vec![0.0; 3],
            labels: vec![1, 1, -1],
            spec_id: vec!["a".to_string(), "b".to_string(), "c".to_string()],
            scan: vec![7, 7, 8],
            exp_mass: vec![500.0; 3],
            peptide: vec![
                "K.SAME.R".to_string(),
                "K.SAME.R".to_string(),
                "K.OTHER.R".to_string(),
            ],
            proteins: vec!["P0".to_string(), "P1".to_string(), "P2".to_string()],
            source: vec![0, 1, 0],
            source_names: vec!["one".to_string(), "two".to_string()],
            ensemble: true,
        }
    }

    #[test]
    fn ensemble_deduplication_keeps_the_best_exact_candidate() {
        let ds = dataset();
        assert_eq!(
            select_reported_indices(&ds, &[1.0, 2.0, 0.5], false, true, 1),
            vec![1, 2]
        );
    }

    #[test]
    fn disabled_competition_on_a_regular_input_keeps_every_row() {
        let ds = dataset();
        assert_eq!(
            select_reported_indices(&ds, &[1.0, 2.0, 0.5], false, false, 1),
            vec![0, 1, 2]
        );
    }
}
