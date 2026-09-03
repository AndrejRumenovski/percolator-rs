//! Peptide identity, representative selection, statistics, and protein mapping.

use crate::{pin, protein, stats};

/// Peptide-level scores and statistics, indexed in canonical input-row order.
pub struct Scores {
    pub indices: Vec<usize>,
    pub scores: Vec<f64>,
    pub q_values: Vec<f64>,
    pub peps: Vec<f64>,
}

/// Return the modified peptide sequence without flanking residues.
pub fn core(sequence: &str) -> &str {
    // Strip flanking residues: A.PEPTIDE.B -> PEPTIDE (keep modifications).
    let bytes = sequence.as_bytes();
    let first = sequence.find('.');
    let last = sequence.rfind('.');
    match (first, last) {
        (Some(a), Some(b)) if b > a => &sequence[a + 1..b],
        _ => {
            let _ = bytes;
            sequence
        }
    }
}

/// Keep the first best-scoring PSM for each `(label, core peptide)` and
/// calculate peptide-level reported-list statistics.
pub fn score(
    ds: &pin::Dataset,
    reported_indices: &[usize],
    psm_scores: &[f64],
    null_target_win_prob: f64,
) -> Scores {
    #[cfg(feature = "profiling")]
    let representative_start = std::time::Instant::now();
    let mut best: ahash::AHashMap<(i8, &str), usize> =
        ahash::AHashMap::with_capacity(reported_indices.len());
    for &index in reported_indices {
        let key = (ds.labels[index], core(&ds.peptide[index]));
        match best.get(&key) {
            Some(&previous) if psm_scores[previous] >= psm_scores[index] => {}
            _ => {
                best.insert(key, index);
            }
        }
    }
    #[cfg(feature = "profiling")]
    crate::profile::record(
        "peptide",
        "peptide_identity_dedup_and_representative",
        representative_start.elapsed(),
        Some(reported_indices.len() as u64),
        None,
    );

    // HashMap iteration is process-randomized. Preserve input order so tied
    // peptide statistics and the loopy-BP message schedule are reproducible.
    let mut indices: Vec<usize> = best.values().copied().collect();
    #[cfg(feature = "profiling")]
    let index_sort = std::time::Instant::now();
    indices.sort_unstable();
    #[cfg(feature = "profiling")]
    crate::profile::record(
        "sort",
        "peptide_input_order",
        index_sort.elapsed(),
        Some(indices.len() as u64),
        None,
    );

    #[cfg(feature = "profiling")]
    let materialization_start = std::time::Instant::now();
    let scores: Vec<f64> = indices.iter().map(|&index| psm_scores[index]).collect();
    let labels: Vec<i8> = indices.iter().map(|&index| ds.labels[index]).collect();
    #[cfg(feature = "profiling")]
    crate::profile::record(
        "peptide",
        "peptide_statistics_materialization",
        materialization_start.elapsed(),
        Some(indices.len() as u64),
        None,
    );
    let (q_values, peps) =
        stats::qvalues_and_peps(&scores, &labels, stats::Tdc::reported(null_target_win_prob));
    Scores {
        indices,
        scores,
        q_values,
        peps,
    }
}

/// Build one protein-inference entry per reported peptide while retaining the
/// complete peptide-to-protein association observed across repeated PSM rows.
///
/// Peptide score and PEP still come from the existing best-PSM selection. The
/// protein mapping is a property of the peptide identity, however, and must not
/// depend on which equal-scoring occurrence happened to be encountered first.
pub fn protein_entries(
    ds: &pin::Dataset,
    reported_indices: &[usize],
    peptides: &Scores,
) -> Vec<(f64, f64, String)> {
    use std::collections::{BTreeMap, BTreeSet};

    debug_assert_eq!(peptides.indices.len(), peptides.scores.len());
    debug_assert_eq!(peptides.indices.len(), peptides.peps.len());

    #[cfg(feature = "profiling")]
    let mapping_union_start = std::time::Instant::now();
    let mut proteins_by_peptide: BTreeMap<(i8, &str), BTreeSet<&str>> = BTreeMap::new();
    for &index in reported_indices {
        let key = (ds.labels[index], core(&ds.peptide[index]));
        let proteins = proteins_by_peptide.entry(key).or_default();
        proteins.extend(protein::split_proteins(&ds.proteins[index]));
    }
    #[cfg(feature = "profiling")]
    crate::profile::record(
        "peptide",
        "peptide_mapping_union",
        mapping_union_start.elapsed(),
        Some(reported_indices.len() as u64),
        None,
    );

    #[cfg(feature = "profiling")]
    let entry_materialization_start = std::time::Instant::now();
    let entries = peptides
        .indices
        .iter()
        .enumerate()
        .map(|(peptide, &index)| {
            let key = (ds.labels[index], core(&ds.peptide[index]));
            let proteins = proteins_by_peptide
                .get(&key)
                .expect("reported peptide is missing its protein associations")
                .iter()
                .copied()
                .collect::<Vec<_>>()
                .join(" ");
            (peptides.scores[peptide], peptides.peps[peptide], proteins)
        })
        .collect();
    #[cfg(feature = "profiling")]
    crate::profile::record(
        "peptide",
        "peptide_mapping_materialization",
        entry_materialization_start.elapsed(),
        Some(peptides.indices.len() as u64),
        None,
    );
    entries
}

#[cfg(test)]
mod tests {
    use super::*;

    fn dataset(labels: Vec<i8>, peptide: Vec<&str>, proteins: Vec<&str>) -> pin::Dataset {
        let n_psm = labels.len();
        pin::Dataset {
            feature_names: vec!["score".to_string()],
            n_feat: 1,
            n_psm,
            features: vec![0.0; n_psm],
            labels,
            spec_id: (0..n_psm).map(|index| format!("row_{index}")).collect(),
            scan: (0..n_psm as i64).collect(),
            exp_mass: vec![500.0; n_psm],
            peptide: peptide.into_iter().map(str::to_string).collect(),
            proteins: proteins.into_iter().map(str::to_string).collect(),
            source: vec![0; n_psm],
            source_names: vec!["fixture.pin".to_string()],
            ensemble: false,
        }
    }

    #[test]
    fn core_identity_removes_flanks_but_retains_modifications() {
        assert_eq!(core("K.PEP[+15.995]TIDE.R"), "PEP[+15.995]TIDE");
        assert_eq!(core("PEPTIDE"), "PEPTIDE");
    }

    #[test]
    fn the_first_exact_best_psm_is_the_representative() {
        let ds = dataset(
            vec![1, 1, 1],
            vec!["K.SAME.R", "A.SAME.B", "K.OTHER.R"],
            vec!["P0", "P1", "P2"],
        );
        let forward = score(&ds, &[0, 1, 2], &[5.0, 5.0, 4.0], 0.5);
        assert_eq!(forward.indices, vec![0, 2]);

        let reversed = score(&ds, &[1, 0, 2], &[5.0, 5.0, 4.0], 0.5);
        assert_eq!(reversed.indices, vec![1, 2]);
    }

    #[test]
    fn target_and_decoy_versions_of_one_peptide_stay_separate() {
        let ds = dataset(
            vec![1, -1],
            vec!["K.SHARED.R", "K.SHARED.R"],
            vec!["TARGET", "DECOY"],
        );
        let result = score(&ds, &[0, 1], &[3.0, 4.0], 0.5);
        assert_eq!(result.indices, vec![0, 1]);
    }

    #[test]
    fn representatives_are_returned_in_input_index_order() {
        let ds = dataset(
            vec![1, 1, 1],
            vec!["K.C.R", "K.A.R", "K.B.R"],
            vec!["P0", "P1", "P2"],
        );
        let result = score(&ds, &[2, 0, 1], &[1.0, 2.0, 3.0], 0.5);
        assert_eq!(result.indices, vec![0, 1, 2]);
    }

    fn two_row_dataset(order: [usize; 2]) -> pin::Dataset {
        let proteins = ["PROT_A", "PROT_B"];
        dataset(
            vec![1, 1],
            vec!["K.AMBIGUOUS.R"; 2],
            order.iter().map(|&index| proteins[index]).collect(),
        )
    }

    /// The two-row minimum for the upstream defect. Protein association is a
    /// set-valued property of peptide identity, independent of which exact- or
    /// near-tied occurrence supplied the peptide score.
    #[test]
    fn complete_peptide_association_survives_representative_and_row_order() {
        for order in [[0, 1], [1, 0]] {
            let ds = two_row_dataset(order);
            for representative in [0, 1] {
                let peptides = Scores {
                    indices: vec![representative],
                    scores: vec![5.0],
                    q_values: vec![0.1],
                    peps: vec![0.25],
                };
                let entries = protein_entries(&ds, &[0, 1], &peptides);
                assert_eq!(entries, vec![(5.0, 0.25, "PROT_A PROT_B".to_string())]);
            }
        }
    }
}
