//! percolator-rs command-line orchestration.

mod cli;

use cli::{ensemble_input, parse_args, ProteinInference};
#[cfg(feature = "profiling")]
use percolator_rs::profile;
use percolator_rs::percolator::Model;
use percolator_rs::{competition, output, percolator, pin, protein, protein_bayes, rt, stats};
use std::borrow::Cow;

#[cfg(feature = "profiling")]
#[global_allocator]
static PROFILING_ALLOCATOR: profile::CountingAllocator = profile::CountingAllocator;

fn core_peptide(p: &str) -> &str {
    // strip flanking residues: A.PEPTIDE.B -> PEPTIDE (keep mods)
    let bytes = p.as_bytes();
    let first = p.find('.');
    let last = p.rfind('.');
    match (first, last) {
        (Some(a), Some(b)) if b > a => &p[a + 1..b],
        _ => {
            let _ = bytes;
            p
        }
    }
}

/// Build one protein-inference entry per reported peptide while retaining the
/// complete peptide-to-protein association observed across repeated PSM rows.
///
/// Peptide score and PEP still come from the existing best-PSM selection.  The
/// protein mapping is a property of the peptide identity, however, and must not
/// depend on which equal-scoring occurrence happened to be encountered first.
fn protein_entries(
    ds: &pin::Dataset,
    reported_indices: &[usize],
    peptide_indices: &[usize],
    peptide_scores: &[f64],
    peptide_peps: &[f64],
) -> Vec<(f64, f64, String)> {
    use std::collections::{BTreeMap, BTreeSet};

    debug_assert_eq!(peptide_indices.len(), peptide_scores.len());
    debug_assert_eq!(peptide_indices.len(), peptide_peps.len());

    let mut proteins_by_peptide: BTreeMap<(i8, &str), BTreeSet<&str>> = BTreeMap::new();
    for &index in reported_indices {
        let key = (ds.labels[index], core_peptide(&ds.peptide[index]));
        let proteins = proteins_by_peptide.entry(key).or_default();
        proteins.extend(protein::split_proteins(&ds.proteins[index]));
    }

    peptide_indices
        .iter()
        .enumerate()
        .map(|(peptide, &index)| {
            let key = (ds.labels[index], core_peptide(&ds.peptide[index]));
            let proteins = proteins_by_peptide
                .get(&key)
                .expect("reported peptide is missing its protein associations")
                .iter()
                .copied()
                .collect::<Vec<_>>()
                .join(" ");
            (peptide_scores[peptide], peptide_peps[peptide], proteins)
        })
        .collect()
}

fn main() {
    let mut args = parse_args();
    if args.pins.is_empty() {
        eprintln!("usage: percolator-rs [flags] input.pin [more.pin ...]");
        eprintln!();
        eprintln!("Input contract: a concatenated target-decoy search against a decoy database");
        eprintln!("of the same size as the target database. Spectrum-level target-decoy");
        eprintln!("competition is performed on the rescored values before PSM statistics, so a");
        eprintln!("PIN reporting several candidates per spectrum is handled; --no-psm-competition");
        eprintln!("reports every candidate instead, and its q-values are then not FDR estimates.");
        eprintln!("--null-target-win-prob P (default 0.5) declares the probability that an");
        eprintln!("incorrect target outranks its paired decoy; use 1/(1+k) for k decoys per");
        eprintln!("target. Separate target/decoy searches (mix-max) are not supported.");
        std::process::exit(2);
    }
    if args.pins.len() > 1 && !args.join && !args.ensemble {
        eprintln!("error: multiple inputs require --join (pooled cross-run training), --ensemble (same-run ENGINE=PIN inputs), or separate runs");
        std::process::exit(2);
    }
    if args.ensemble && args.pins.len() < 2 {
        eprintln!("error: --ensemble requires at least two ENGINE=PIN inputs");
        std::process::exit(2);
    }
    if args.ensemble && (args.results_proteins.is_some() || args.decoy_proteins.is_some()) {
        eprintln!("error: protein inference is unavailable with --ensemble; engine-level duplicate evidence needs a dedicated protein model");
        std::process::exit(2);
    }
    #[cfg(feature = "profiling")]
    let profile_session = profile::Session::start(
        args.profile_json.clone(),
        args.profile_cpu.clone(),
        args.profile_allocations,
    )
    .unwrap_or_else(|message| {
        eprintln!("profiling error: {message}");
        std::process::exit(2);
    });
    #[cfg(feature = "profiling")]
    {
        profile::metadata("profile_name", args.profile);
        profile::metadata("seed", args.params.seed);
        profile::metadata("num_threads", args.params.num_threads);
        profile::metadata("maxiter", args.params.maxiter);
        profile::metadata("subset_max_train", args.params.subset_max_train);
        profile::metadata("allocation_counting", args.profile_allocations);
        profile::metadata("input_files", &args.pins);
    }
    eprintln!(
        "profile: {} (model={}, maxiter={}, subset-max-train={}){}{}{}",
        args.profile,
        args.params.model.label(),
        args.params.maxiter,
        if args.params.subset_max_train == 0 {
            "none".to_string()
        } else {
            args.params.subset_max_train.to_string()
        },
        if args.join {
            format!(", join={} files", args.pins.len())
        } else if args.ensemble {
            format!(", ensemble={} engines", args.pins.len())
        } else {
            String::new()
        },
        if args.rt_features {
            ", rt-features"
        } else {
            ""
        },
        if args.params.nested_selection {
            ", nested-selection"
        } else {
            ""
        },
    );
    let t0 = std::time::Instant::now();
    let tp = std::time::Instant::now();
    #[cfg(feature = "profiling")]
    let _input_loading = profile::Scope::new("stage", "input_loading");
    let ensemble_inputs: Vec<(String, String)> = if args.ensemble {
        args.pins
            .iter()
            .map(|input| ensemble_input(input))
            .collect::<Result<_, _>>()
            .unwrap_or_else(|message| {
                eprintln!("error: {message}");
                std::process::exit(2);
            })
    } else {
        args.pins
            .iter()
            .map(|path| (String::new(), path.clone()))
            .collect()
    };
    let mut parts: Vec<pin::Dataset> = Vec::with_capacity(ensemble_inputs.len());
    for (_, path) in &ensemble_inputs {
        parts.push(pin::parse(path).unwrap_or_else(|e| {
            eprintln!("parse error ({path}): {e}");
            std::process::exit(1);
        }));
    }
    let ds = if args.ensemble {
        pin::merge_ensemble(
            parts,
            ensemble_inputs
                .into_iter()
                .map(|(engine, _)| engine)
                .collect(),
        )
        .unwrap_or_else(|message| {
            eprintln!("ensemble error: {message}");
            std::process::exit(2);
        })
    } else if parts.len() == 1 {
        parts.pop().unwrap()
    } else {
        pin::merge(parts)
    };
    let mut ds = ds;
    if args.rt_features {
        // Reserve the residual columns now; the alignment behind them is
        // label-dependent, so it is refitted inside every outer training
        // partition rather than once here.
        args.params.rt = rt::augment(&mut ds);
    }
    let ds = ds;
    #[cfg(feature = "profiling")]
    {
        drop(_input_loading);
        profile::metadata("psms", ds.n_psm);
        profile::metadata("features", ds.n_feat);
        profile::metadata(
            "input_bytes",
            args.pins
                .iter()
                .filter_map(|path| std::fs::metadata(path).ok().map(|metadata| metadata.len()))
                .sum::<u64>(),
        );
    }
    eprintln!("parse: {:.3}s", tp.elapsed().as_secs_f64());
    eprintln!(
        "parsed {} PSMs, {} features ({} targets / {} decoys){}",
        ds.n_psm,
        ds.n_feat,
        ds.labels.iter().filter(|&&l| l > 0).count(),
        ds.labels.iter().filter(|&&l| l < 0).count(),
        if args.join {
            format!(", pooled from {} files", ds.source_names.len())
        } else if args.ensemble {
            format!(", ensemble from {} engines", ds.source_names.len())
        } else {
            String::new()
        },
    );

    #[cfg(feature = "profiling")]
    let _rescoring = profile::Scope::new("stage", "rescoring");
    let out = percolator::run(&ds, &args.params);
    #[cfg(feature = "profiling")]
    drop(_rescoring);
    if out.nested_folds.is_empty() {
        eprintln!(
            "{} class weights: Cpos={:.3} Cneg={:.3}{}{}",
            if args.params.model == Model::Svm {
                "SVM"
            } else {
                "MLP"
            },
            out.c_alpha,
            out.c_beta,
            if out.c_selected {
                " (selected by cross-validation)"
            } else {
                " (fixed)"
            },
            if args.params.model == Model::Mlp {
                format!(
                    "; hidden={}, epochs/iteration={}, learning-rate={}, l2={}",
                    args.params.mlp_hidden,
                    args.params.mlp_epochs,
                    args.params.mlp_learning_rate,
                    args.params.mlp_l2
                )
            } else {
                String::new()
            }
        );
    } else {
        eprintln!("nested SVM selection (outer test folds isolated):");
        for selected in &out.nested_folds {
            eprintln!(
                "  fold {}: C={:.3}, class-weights={:.1}:{:.1}, features={}, tolerance={:.0e}, inner-q01-yield={}",
                selected.outer_fold,
                selected.c,
                selected.positive_weight,
                selected.negative_weight,
                selected.feature_count,
                selected.tolerance,
                selected.inner_yield,
            );
        }
    }

    if let Some(path) = &args.feature_report {
        let report_start = std::time::Instant::now();
        let report = percolator::feature_report(&ds, &args.params, &out);
        output::write_feature_report(path, &report).unwrap_or_else(|error| {
            eprintln!("feature report error ({path}): {error}");
            std::process::exit(1);
        });
        eprintln!(
            "feature report: {} features, baseline target PSMs q<0.01={}, {:.3}s",
            report.features.len(),
            report.baseline_q01,
            report_start.elapsed().as_secs_f64()
        );
    }

    // Engine reports of the same candidate are training observations, not separate
    // discoveries. Retain the best out-of-fold score per exact candidate and
    // recalibrate q-values on that deduplicated candidate set for ensemble output.
    #[cfg(feature = "profiling")]
    let _psm_context = profile::context(Some("psm_level_processing"), None, None, None);
    #[cfg(feature = "profiling")]
    let _psm_processing = profile::Scope::new("stage", "psm_level_processing");
    let reported_indices: Vec<usize> = if args.psm_competition {
        competition::winner_indices(&ds, &out.score, args.params.seed)
    } else if args.ensemble {
        let mut best: std::collections::BTreeMap<(i64, i8, String), usize> =
            std::collections::BTreeMap::new();
        for i in 0..ds.n_psm {
            let key = (ds.scan[i], ds.labels[i], ds.peptide[i].clone());
            match best.get(&key) {
                Some(&previous) if out.score[previous] >= out.score[i] => {}
                _ => {
                    best.insert(key, i);
                }
            }
        }
        best.into_values().collect()
    } else {
        (0..ds.n_psm).collect()
    };
    // Statistics belong to the list that is actually reported. Whenever
    // competition or ensemble deduplication has removed rows, the q-values
    // computed over the full training list no longer describe it.
    let (reported_qvals, reported_peps) = if reported_indices.len() == ds.n_psm {
        (out.qval.clone(), out.pep.clone())
    } else {
        let reported_scores: Vec<f64> = reported_indices.iter().map(|&i| out.score[i]).collect();
        let reported_labels: Vec<i8> = reported_indices.iter().map(|&i| ds.labels[i]).collect();
        stats::qvalues_and_peps(
            &reported_scores,
            &reported_labels,
            stats::Tdc::reported(args.params.null_target_win_prob),
        )
    };

    // PSM-level output
    let target_capacity = reported_indices
        .iter()
        .filter(|&&index| ds.labels[index] > 0)
        .count();
    let mut targets: Vec<output::Row<'_>> = Vec::with_capacity(target_capacity);
    let mut decoys: Vec<output::Row<'_>> =
        Vec::with_capacity(reported_indices.len() - target_capacity);
    #[cfg(feature = "profiling")]
    let mut psm_row_string_bytes = 0u64;
    for (output_index, &i) in reported_indices.iter().enumerate() {
        let r = output::Row::new(
            if args.ensemble {
                Cow::Owned(format!(
                    "{}:{}",
                    ds.source_names[ds.source[i] as usize], ds.spec_id[i]
                ))
            } else {
                Cow::Borrowed(&ds.spec_id[i])
            },
            out.score[i],
            reported_qvals[output_index],
            reported_peps[output_index],
            &ds.peptide[i],
            &ds.proteins[i],
        );
        #[cfg(feature = "profiling")]
        {
            psm_row_string_bytes += r.owned_id_capacity();
        }
        if ds.labels[i] > 0 {
            targets.push(r);
        } else {
            decoys.push(r);
        }
    }
    #[cfg(feature = "profiling")]
    {
        profile::allocation_site(
            "main::psm output row vectors",
            2,
            ((targets.capacity() + decoys.capacity())
                * std::mem::size_of::<output::Row>()) as u64,
        );
        profile::allocation_site(
            "main::psm output row strings",
            u64::from(args.ensemble) * reported_indices.len() as u64,
            psm_row_string_bytes,
        );
    }
    let n_psm_q01 = targets.iter().filter(|r| r.q_value() < 0.01).count();
    #[cfg(feature = "profiling")]
    drop(_psm_processing);
    #[cfg(feature = "profiling")]
    drop(_psm_context);

    // Peptide-level: best-scoring PSM per unique (core) peptide, then re-q-value.
    #[cfg(feature = "profiling")]
    let _peptide_context = profile::context(Some("peptide_level_processing"), None, None, None);
    #[cfg(feature = "profiling")]
    let _peptide_processing = profile::Scope::new("stage", "peptide_level_processing");
    let mut best: ahash::AHashMap<(i8, &str), usize> =
        ahash::AHashMap::with_capacity(reported_indices.len());
    for &i in &reported_indices {
        let key = (ds.labels[i], core_peptide(&ds.peptide[i]));
        match best.get(&key) {
            Some(&j) if out.score[j] >= out.score[i] => {}
            _ => {
                best.insert(key, i);
            }
        }
    }
    // HashMap iteration is process-randomized. Preserve input order so tied
    // peptide statistics and the loopy-BP message schedule are reproducible.
    let mut pep_idx: Vec<usize> = best.values().copied().collect();
    #[cfg(feature = "profiling")]
    let peptide_index_sort = std::time::Instant::now();
    pep_idx.sort_unstable();
    #[cfg(feature = "profiling")]
    profile::record(
        "sort",
        "peptide_input_order",
        peptide_index_sort.elapsed(),
        Some(pep_idx.len() as u64),
        None,
    );
    let pscore: Vec<f64> = pep_idx.iter().map(|&i| out.score[i]).collect();
    let plabel: Vec<i8> = pep_idx.iter().map(|&i| ds.labels[i]).collect();
    let (pq, ppep) = stats::qvalues_and_peps(
        &pscore,
        &plabel,
        stats::Tdc::reported(args.params.null_target_win_prob),
    );

    let peptide_target_capacity = pep_idx
        .iter()
        .filter(|&&index| ds.labels[index] > 0)
        .count();
    let mut ptargets: Vec<output::Row<'_>> = Vec::with_capacity(peptide_target_capacity);
    let mut pdecoys: Vec<output::Row<'_>> =
        Vec::with_capacity(pep_idx.len() - peptide_target_capacity);
    for (k, &i) in pep_idx.iter().enumerate() {
        let r = output::Row::new(
            Cow::Borrowed(&ds.spec_id[i]),
            pscore[k],
            pq[k],
            ppep[k],
            &ds.peptide[i],
            &ds.proteins[i],
        );
        if ds.labels[i] > 0 {
            ptargets.push(r);
        } else {
            pdecoys.push(r);
        }
    }
    #[cfg(feature = "profiling")]
    {
        profile::allocation_site(
            "main::peptide output row vectors",
            2,
            ((ptargets.capacity() + pdecoys.capacity())
                * std::mem::size_of::<output::Row>()) as u64,
        );
    }
    let n_pep_q01 = ptargets.iter().filter(|r| r.q_value() < 0.01).count();
    #[cfg(feature = "profiling")]
    drop(_peptide_processing);
    #[cfg(feature = "profiling")]
    drop(_peptide_context);

    #[cfg(feature = "profiling")]
    let _output_context = profile::context(Some("result_output"), None, None, None);
    #[cfg(feature = "profiling")]
    let _output = profile::Scope::new("stage", "result_output");
    if let Some(p) = &args.results_psms {
        output::write_results(p, targets).unwrap();
    }
    if let Some(p) = &args.decoy_psms {
        output::write_results(p, decoys).unwrap();
    }
    if let Some(p) = &args.results_peptides {
        output::write_results(p, ptargets).unwrap();
    }
    if let Some(p) = &args.decoy_peptides {
        output::write_results(p, pdecoys).unwrap();
    }
    #[cfg(feature = "profiling")]
    drop(_output);
    #[cfg(feature = "profiling")]
    drop(_output_context);

    // Cross-run: per-source yield when pooled (each file's targets scored by the shared model).
    if args.join {
        eprintln!("per-file yield (pooled model, target PSMs q<0.01):");
        for s in 0..ds.source_names.len() as u32 {
            let c = (0..ds.n_psm)
                .filter(|&i| ds.source[i] == s && ds.labels[i] > 0 && out.qval[i] < 0.01)
                .count();
            eprintln!("  [{}] {}", ds.source_names[s as usize], c);
        }
    }

    // Protein inference uses the best score/PEP for each peptide sequence and
    // the union of its protein mappings across all reported PSM occurrences.
    if args.results_proteins.is_some() || args.decoy_proteins.is_some() {
        #[cfg(feature = "profiling")]
        let _protein_inference = profile::Scope::new("stage", "protein_inference_and_output");
        let entries = protein_entries(&ds, &reported_indices, &pep_idx, &pscore, &ppep);
        #[cfg(feature = "profiling")]
        profile::allocation_site(
            "main::protein inference entries",
            (entries.len() + 1) as u64,
            (entries.capacity() * std::mem::size_of::<(f64, f64, String)>()) as u64
                + entries
                    .iter()
                    .map(|entry| entry.2.capacity() as u64)
                    .sum::<u64>(),
        );
        let picked_groups = protein::infer(&entries, args.params.seed);
        let picked_q01 = picked_groups
            .iter()
            .filter(|g| !g.is_decoy && g.picked && g.qval < 0.01)
            .count();
        let n_prot_classic = protein::classic_target_q01(&picked_groups);
        let (groups, method_label) = match args.protein_inference {
            ProteinInference::Picked => (picked_groups, "picked-FDR"),
            ProteinInference::Bayesian => {
                let result = protein_bayes::infer(&entries, &args.protein_bayes);
                eprintln!(
                    "Bayesian protein model: alpha={:.4}, beta={:.4}, gamma={:.4}, peptide-prior={:.4}; components: {} ({} tree-exact, {} loopy); BP iterations: {}, converged: {}",
                    args.protein_bayes.alpha,
                    args.protein_bayes.beta,
                    args.protein_bayes.gamma,
                    args.protein_bayes.peptide_prior,
                    result.diagnostics.components,
                    result.diagnostics.tree_components,
                    result.diagnostics.loopy_components,
                    result.diagnostics.iterations,
                    result.diagnostics.converged,
                );
                (result.groups, "Bayesian")
            }
        };
        let n_prot_q01 = groups
            .iter()
            .filter(|g| !g.is_decoy && g.picked && g.qval < 0.01)
            .count();
        if let Some(p) = &args.results_proteins {
            output::write_proteins(p, &groups, false).unwrap();
        }
        if let Some(p) = &args.decoy_proteins {
            output::write_proteins(p, &groups, true).unwrap();
        }
        if args.protein_inference == ProteinInference::Picked {
            eprintln!(
                "protein groups: {} ({} target, {} decoy); picked entries: {} | target proteins q<0.01: {} (picked-FDR) vs {} (classic)",
                groups.len(),
                groups.iter().filter(|g| !g.is_decoy).count(),
                groups.iter().filter(|g| g.is_decoy).count(),
                groups.iter().filter(|g| g.picked).count(),
                n_prot_q01,
                n_prot_classic
            );
        } else {
            eprintln!(
                "protein groups: {} ({} target, {} decoy); reported entries: {} | target proteins q<0.01: {} ({}) vs {} (picked-FDR) vs {} (classic)",
                groups.len(),
                groups.iter().filter(|g| !g.is_decoy).count(),
                groups.iter().filter(|g| g.is_decoy).count(),
                groups.iter().filter(|g| g.picked).count(),
                n_prot_q01,
                method_label,
                picked_q01,
                n_prot_classic
            );
        }
    }

    eprintln!(
        "target PSMs q<0.01: {} | target peptides q<0.01: {} | {:.2}s",
        n_psm_q01,
        n_pep_q01,
        t0.elapsed().as_secs_f64()
    );
    #[cfg(feature = "profiling")]
    profile_session.finish().unwrap_or_else(|message| {
        eprintln!("profiling error: {message}");
        std::process::exit(1);
    });
}

#[cfg(test)]
mod protein_entry_tests {
    use super::*;

    fn two_row_dataset(order: [usize; 2]) -> pin::Dataset {
        let proteins = ["PROT_A", "PROT_B"];
        pin::Dataset {
            feature_names: vec!["score".to_string()],
            n_feat: 1,
            n_psm: 2,
            features: vec![5.0, f64::from_bits(5.0f64.to_bits() - 1)],
            labels: vec![1, 1],
            spec_id: order.iter().map(|&i| format!("AMB_{i}")).collect(),
            scan: vec![1, 1],
            exp_mass: vec![500.0, 500.0],
            peptide: vec!["K.AMBIGUOUS.R".to_string(); 2],
            proteins: order.iter().map(|&i| proteins[i].to_string()).collect(),
            source: vec![0, 0],
            source_names: vec!["fixture.pin".to_string()],
            ensemble: false,
        }
    }

    /// The two-row minimum for the upstream defect. Protein association is a
    /// set-valued property of peptide identity, independent of which exact- or
    /// near-tied occurrence supplied the peptide score.
    #[test]
    fn complete_peptide_association_survives_representative_and_row_order() {
        for order in [[0, 1], [1, 0]] {
            let ds = two_row_dataset(order);
            for representative in [0, 1] {
                let entries = protein_entries(
                    &ds,
                    &[0, 1],
                    &[representative],
                    &[5.0],
                    &[0.25],
                );
                assert_eq!(
                    entries,
                    vec![(5.0, 0.25, "PROT_A PROT_B".to_string())]
                );
            }
        }
    }
}
