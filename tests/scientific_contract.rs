//! Regression tests for minimized failures from the scientific audit.

use percolator_rs::{competition, peptide, pin::Dataset, protein, stats};
use std::collections::BTreeMap;

fn tied_candidates(target_copies: usize, decoy_copies: usize) -> Dataset {
    let mut ds = Dataset {
        feature_names: vec!["score".into()],
        n_feat: 1,
        n_psm: 0,
        features: Vec::new(),
        labels: Vec::new(),
        spec_id: Vec::new(),
        scan: Vec::new(),
        exp_mass: Vec::new(),
        peptide: Vec::new(),
        proteins: Vec::new(),
        source: Vec::new(),
        source_names: vec!["fixture.pin".into()],
        ensemble: false,
    };
    for scan in 1..=233 {
        for (label, copies) in [(1, target_copies), (-1, decoy_copies)] {
            for copy in 0..copies {
                ds.features.push(7.0);
                ds.labels.push(label);
                ds.spec_id.push(format!("{scan}_{label}_{copy}"));
                ds.scan.push(scan);
                ds.exp_mass.push(500.0);
                ds.peptide.push(format!("K.PEPTIDE{label}.R"));
                ds.proteins.push(format!("PROTEIN{label}_{copy}"));
                ds.source.push(0);
            }
        }
    }
    ds.n_psm = ds.labels.len();
    ds
}

#[test]
fn duplicate_occurrences_do_not_buy_additional_tie_draws() {
    let reference = tied_candidates(1, 1);
    for (targets, decoys) in [(94, 1), (1, 94), (4, 7)] {
        let duplicated = tied_candidates(targets, decoys);
        for seed in 1..=10 {
            let winners = |ds: &Dataset| -> BTreeMap<_, _> {
                competition::winner_indices(ds, &ds.features, seed)
                    .into_iter()
                    .map(|row| {
                        (
                            ds.scan[row],
                            (ds.labels[row], peptide::core(&ds.peptide[row]).to_owned()),
                        )
                    })
                    .collect()
            };
            assert_eq!(
                winners(&reference),
                winners(&duplicated),
                "{targets}/{decoys}, seed {seed}"
            );
        }
    }
}

#[test]
fn protein_mapping_survives_precursor_competition() {
    let mut ds = tied_candidates(2, 1);
    ds.peptide[0] = "K.SHARED.R".into();
    ds.peptide[1] = "A.SHARED.B".into();
    ds.proteins[0] = "PROT_A".into();
    ds.proteins[1] = "PROT_B".into();
    ds.peptide[2] = "K.SHARED.R".into();
    ds.proteins[2] = "DECOY_ONLY".into();
    let scores = vec![3.0; ds.n_psm];
    for winner in [0, 1] {
        let peptides = peptide::score(&ds, &[winner], &scores, 0.5);
        let entries = peptide::protein_entries(&ds, &[winner], &peptides);
        assert_eq!(entries.len(), 1);
        assert_eq!(entries[0].2, "PROT_A PROT_B");
        assert_eq!(entries[0].0, 3.0);
    }
}

#[test]
fn distinct_protein_member_sets_do_not_share_a_pairing_slot() {
    let entries = vec![
        (10.0, 0.1, "LEFT RIGHT".into()),
        (9.0, 0.1, "LEFT|RIGHT".into()),
        (8.0, 0.1, "DECOY_LEFT DECOY_RIGHT".into()),
    ];
    for arranged in [entries.clone(), entries.into_iter().rev().collect()] {
        let groups = protein::infer(&arranged, 37);
        assert_eq!(groups.len(), 3);
        assert_eq!(groups.iter().filter(|g| !g.is_decoy && g.picked).count(), 2);
        assert!(groups.iter().filter(|g| g.is_decoy).all(|g| !g.picked));
    }
}

#[test]
fn tiny_null_probability_does_not_acquire_artificial_pep_mass() {
    let scores = [5.0, 4.0, 3.0, 2.0, 1.0];
    for probability in [1e-15, 1e-50, 1e-200, 0.5] {
        let (_, pep) = stats::qvalues_and_peps(&scores, &[1; 5], stats::Tdc::reported(probability));
        let expected = probability / (1.0 - probability);
        let observed: f64 = pep.iter().sum();
        assert!(
            (observed / expected - 1.0).abs() < 1e-14,
            "p={probability}: {observed} != {expected}"
        );
        assert!(pep.iter().all(|p| p.is_finite() && *p > 0.0 && *p <= 1.0));
    }
}

#[test]
fn decoy_pep_display_is_invariant_inside_mixed_score_ties() {
    let scores = [
        3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 2.0, 2.0, 1.0, 1.0, 1.0,
    ];
    let labels = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, -1, 1, 1, 1, 1];
    let expected = stats::qvalues_and_peps(&scores, &labels, stats::Tdc::reported(0.5));
    assert_eq!(expected.1[10], expected.1[11]);
    for shift in 0..scores.len() {
        let mut order: Vec<_> = (0..scores.len()).rev().collect();
        order.rotate_left(shift);
        let permuted_scores: Vec<_> = order.iter().map(|&row| scores[row]).collect();
        let permuted_labels: Vec<_> = order.iter().map(|&row| labels[row]).collect();
        let (q, pep) = stats::qvalues_and_peps(
            &permuted_scores,
            &permuted_labels,
            stats::Tdc::reported(0.5),
        );
        for (row, &original) in order.iter().enumerate() {
            assert_eq!(q[row], expected.0[original]);
            assert_eq!(pep[row], expected.1[original]);
        }
    }
}

#[test]
fn exported_probabilities_preserve_threshold_membership_and_tiny_values() {
    use percolator_rs::output;
    use std::borrow::Cow;
    let directory = tempfile::tempdir().unwrap();
    let probabilities = [0.01f64.next_down(), 0.01, 0.01f64.next_up(), 1e-15];
    for (index, &q) in probabilities.iter().enumerate() {
        let psms = directory.path().join(format!("psms-{index}.tsv"));
        output::write_results(
            psms.to_str().unwrap(),
            vec![output::Row::new(
                Cow::Borrowed("psm"),
                3.0,
                q,
                q,
                "K.PEP.R",
                "P1",
            )],
        )
        .unwrap();
        let text = std::fs::read_to_string(psms).unwrap();
        let row: Vec<_> = text.lines().nth(1).unwrap().split('\t').collect();
        assert_eq!(row[2].parse::<f64>().unwrap(), q);
        assert_eq!(row[3].parse::<f64>().unwrap(), q);

        let proteins = directory.path().join(format!("proteins-{index}.tsv"));
        output::write_proteins(
            proteins.to_str().unwrap(),
            &[protein::ProtGroup {
                proteins: vec!["P1".into()],
                score: 3.0,
                qval: q,
                pep: Some(q),
                n_peptides: 1,
                is_decoy: false,
                picked: true,
            }],
            false,
        )
        .unwrap();
        let text = std::fs::read_to_string(proteins).unwrap();
        let row: Vec<_> = text.lines().nth(1).unwrap().split('\t').collect();
        assert_eq!(row[1].parse::<f64>().unwrap(), q);
        assert_eq!(row[2].parse::<f64>().unwrap(), q);
    }
}
