//! percolator-rs command-line orchestration.

mod cli;

use cli::{ensemble_input, parse_args, ProteinInference};
use percolator_rs::percolator::Model;
#[cfg(feature = "profiling")]
use percolator_rs::profile;
use percolator_rs::{output, percolator, pin, pipeline, rt};

#[cfg(feature = "profiling")]
#[global_allocator]
static PROFILING_ALLOCATOR: profile::CountingAllocator = profile::CountingAllocator;

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

    let out = pipeline::rescore(&ds, &args.params);
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

    let pipeline::Reports {
        target_psms,
        decoy_psms,
        target_peptides,
        decoy_peptides,
        reported_indices,
        peptides,
        target_psms_q01,
        target_peptides_q01,
    } = pipeline::build_reports(&ds, &out, &args.params, args.psm_competition, args.ensemble);

    #[cfg(feature = "profiling")]
    let _output_context = profile::context(Some("result_output"), None, None, None);
    #[cfg(feature = "profiling")]
    let _output = profile::Scope::new("stage", "result_output");
    if let Some(p) = &args.results_psms {
        output::write_results(p, target_psms).unwrap();
    }
    if let Some(p) = &args.decoy_psms {
        output::write_results(p, decoy_psms).unwrap();
    }
    if let Some(p) = &args.results_peptides {
        output::write_results(p, target_peptides).unwrap();
    }
    if let Some(p) = &args.decoy_peptides {
        output::write_results(p, decoy_peptides).unwrap();
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
        let method = match args.protein_inference {
            ProteinInference::Picked => pipeline::ProteinMethod::Picked,
            ProteinInference::Bayesian => pipeline::ProteinMethod::Bayesian(&args.protein_bayes),
        };
        let pipeline::ProteinResults {
            groups,
            picked_q01,
            classic_q01,
            bayesian_diagnostics,
        } = pipeline::infer_proteins(&ds, &reported_indices, &peptides, args.params.seed, method);
        if let Some(diagnostics) = bayesian_diagnostics {
            eprintln!(
                "Bayesian protein model: alpha={:.4}, beta={:.4}, gamma={:.4}, peptide-prior={:.4}; components: {} ({} tree-exact, {} loopy); BP iterations: {}, converged: {}",
                args.protein_bayes.alpha,
                args.protein_bayes.beta,
                args.protein_bayes.gamma,
                args.protein_bayes.peptide_prior,
                diagnostics.components,
                diagnostics.tree_components,
                diagnostics.loopy_components,
                diagnostics.iterations,
                diagnostics.converged,
            );
        }
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
                classic_q01
            );
        } else {
            eprintln!(
                "protein groups: {} ({} target, {} decoy); reported entries: {} | target proteins q<0.01: {} ({}) vs {} (picked-FDR) vs {} (classic)",
                groups.len(),
                groups.iter().filter(|g| !g.is_decoy).count(),
                groups.iter().filter(|g| g.is_decoy).count(),
                groups.iter().filter(|g| g.picked).count(),
                n_prot_q01,
                "Bayesian",
                picked_q01,
                classic_q01
            );
        }
    }

    eprintln!(
        "target PSMs q<0.01: {} | target peptides q<0.01: {} | {:.2}s",
        target_psms_q01,
        target_peptides_q01,
        t0.elapsed().as_secs_f64()
    );
    #[cfg(feature = "profiling")]
    profile_session.finish().unwrap_or_else(|message| {
        eprintln!("profiling error: {message}");
        std::process::exit(1);
    });
}
