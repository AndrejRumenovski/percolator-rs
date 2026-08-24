//! percolator-rs — a from-scratch Rust reimplementation of the Percolator
//! semi-supervised PSM rescoring algorithm. CLI mirrors the subset of reference
//! flags used by our benchmark.

mod mlp;
mod percolator;
mod pin;
mod protein;
mod protein_bayes;
mod rt;
mod simd;
mod stats;
mod svm;

use percolator::{Model, Params};
use std::fs::File;
use std::io::{BufWriter, Write};

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

struct Args {
    pins: Vec<String>,
    results_psms: Option<String>,
    decoy_psms: Option<String>,
    results_peptides: Option<String>,
    decoy_peptides: Option<String>,
    results_proteins: Option<String>,
    decoy_proteins: Option<String>,
    feature_report: Option<String>,
    protein_inference: ProteinInference,
    protein_bayes: protein_bayes::Params,
    join: bool,
    ensemble: bool,
    rt_features: bool,
    params: Params,
    profile: &'static str,
}

#[derive(Clone, Copy, PartialEq)]
enum ProteinInference {
    Picked,
    Bayesian,
}

/// Execution profiles: preset (subset_max_train, maxiter, select_c) tuned for a use case.
/// Explicit --subset-max-train / --maxiter / --cpos / --cneg always override the preset.
/// `select_c` enables the per-file SVM class-weight grid search. It is OFF for every profile:
/// measured on PXD032157 it costs ~2.5x wall time and does not beat the fixed default weights
/// (32 files better, 28 worse, 5 tied, aggregate slightly worse). Opt in with --select-c on data where
/// the default Cpos/Cneg may not transfer.
fn preset(name: &str) -> Option<(usize, usize, bool)> {
    match name {
        "fast" => Some((20_000, 5, false)),   // ~quick QA/test pipelines
        "balanced" => Some((40_000, 10, false)), // ~5% yield hit, still fast
        "canonical" => Some((0, 10, false)), // full default sensitivity (0 = no subsetting)
        _ => None,
    }
}

fn parse_args() -> Args {
    let mut a = Args {
        pins: Vec::new(),
        results_psms: None,
        decoy_psms: None,
        results_peptides: None,
        decoy_peptides: None,
        results_proteins: None,
        decoy_proteins: None,
        feature_report: None,
        protein_inference: ProteinInference::Picked,
        protein_bayes: protein_bayes::Params::default(),
        join: false,
        ensemble: false,
        rt_features: false,
        params: Params::default(),
        profile: "canonical",
    };
    let argv: Vec<String> = std::env::args().skip(1).collect();

    // Two-pass so explicit flags win regardless of position relative to the profile flag.
    let mut prof: Option<&'static str> = None;
    let mut maxiter_opt: Option<usize> = None;
    let mut subset_opt: Option<usize> = None;
    let mut alpha_opt: Option<f64> = None;
    let mut beta_opt: Option<f64> = None;
    let mut select_c_opt: Option<bool> = None;

    let mut i = 0;
    while i < argv.len() {
        let s = &argv[i];
        let mut take = || {
            i += 1;
            argv.get(i).cloned().unwrap_or_default()
        };
        match s.as_str() {
            "--results-psms" | "-m" => a.results_psms = Some(take()),
            "--decoy-results-psms" | "-M" => a.decoy_psms = Some(take()),
            "--results-peptides" | "-r" => a.results_peptides = Some(take()),
            "--decoy-results-peptides" | "-B" => a.decoy_peptides = Some(take()),
            "--results-proteins" | "-l" => a.results_proteins = Some(take()),
            "--decoy-results-proteins" | "-L" => a.decoy_proteins = Some(take()),
            "--feature-report" => a.feature_report = Some(take()),
            "--protein-inference" => {
                let v = take();
                a.protein_inference = match v.as_str() {
                    "picked" => ProteinInference::Picked,
                    "bayesian" | "fido" => ProteinInference::Bayesian,
                    _ => {
                        eprintln!("unknown --protein-inference '{v}' (use picked|bayesian)");
                        std::process::exit(2);
                    }
                };
            }
            "--protein-alpha" => a.protein_bayes.alpha = take().parse().unwrap_or(f64::NAN),
            "--protein-beta" => a.protein_bayes.beta = take().parse().unwrap_or(f64::NAN),
            "--protein-gamma" => a.protein_bayes.gamma = take().parse().unwrap_or(f64::NAN),
            "--protein-peptide-prior" => {
                a.protein_bayes.peptide_prior = take().parse().unwrap_or(f64::NAN)
            }
            "--protein-max-iter" => {
                a.protein_bayes.max_iter = take().parse().unwrap_or(0)
            }
            "--rescore-model" | "--model" => {
                let value = take();
                a.params.model = match value.as_str() {
                    "svm" | "linear" => Model::Svm,
                    "mlp" | "neural" => Model::Mlp,
                    _ => {
                        eprintln!("unknown --rescore-model '{value}' (use svm|mlp)");
                        std::process::exit(2);
                    }
                };
            }
            "--mlp-hidden" => a.params.mlp_hidden = take().parse().unwrap_or(0),
            "--mlp-epochs" => a.params.mlp_epochs = take().parse().unwrap_or(0),
            "--mlp-learning-rate" => {
                a.params.mlp_learning_rate = take().parse().unwrap_or(f64::NAN)
            }
            "--mlp-l2" => a.params.mlp_l2 = take().parse().unwrap_or(f64::NAN),
            "--auto-model" | "--nested-select" => a.params.nested_selection = true,
            "--no-auto-model" => a.params.nested_selection = false,
            "--svm-tolerance" => {
                a.params.svm_tolerance = take().parse().unwrap_or(f64::NAN)
            }
            "--join" => a.join = true,
            "--ensemble" => a.ensemble = true,
            "--rt-features" => a.rt_features = true,
            "--seed" => a.params.seed = take().parse().unwrap_or(1),
            "--maxiter" => maxiter_opt = take().parse().ok(),
            "--subset-max-train" | "-N" => subset_opt = take().parse().ok(),
            "--cpos" => alpha_opt = take().parse().ok(),
            "--cneg" => beta_opt = take().parse().ok(),
            "--select-c" => select_c_opt = Some(true),
            "--no-select-c" => select_c_opt = Some(false),
            "--num-threads" => a.params.num_threads = take().parse().unwrap_or(1).max(1),
            "--fast" => prof = Some("fast"),
            "--balanced" => prof = Some("balanced"),
            "--canonical" => prof = Some("canonical"),
            "--profile" => {
                let v = take();
                prof = match v.as_str() {
                    "fast" => Some("fast"),
                    "balanced" => Some("balanced"),
                    "canonical" => Some("canonical"),
                    _ => {
                        eprintln!("unknown --profile '{v}' (use fast|balanced|canonical)");
                        std::process::exit(2);
                    }
                };
            }
            other => {
                if !other.starts_with('-') {
                    a.pins.push(other.to_string());
                }
            }
        }
        i += 1;
    }
    if let Err(message) = a.protein_bayes.validate() {
        eprintln!("invalid Bayesian protein parameter: {message}");
        std::process::exit(2);
    }

    // Resolve: profile sets the baseline, explicit flags override.
    let chosen = prof.unwrap_or("canonical");
    a.profile = chosen;
    let mut select_c = false;
    if let Some((subset, maxiter, sel)) = preset(chosen) {
        a.params.subset_max_train = subset;
        a.params.maxiter = maxiter;
        select_c = sel;
    }
    if let Some(m) = maxiter_opt {
        a.params.maxiter = m;
    }
    if let Some(s) = subset_opt {
        a.params.subset_max_train = s;
    }
    if let Some(s) = select_c_opt {
        select_c = s;
    }
    if a.params.nested_selection && a.params.model != Model::Svm {
        eprintln!("--auto-model currently supports only --rescore-model svm");
        std::process::exit(2);
    }
    if a.ensemble && a.join {
        eprintln!("--ensemble and --join are mutually exclusive");
        std::process::exit(2);
    }
    if a.feature_report.is_some() && a.params.model != Model::Svm {
        eprintln!("--feature-report currently supports only --rescore-model svm");
        std::process::exit(2);
    }
    if a.params.nested_selection && select_c {
        eprintln!("--auto-model and legacy --select-c are mutually exclusive");
        std::process::exit(2);
    }
    if a.params.nested_selection && (alpha_opt.is_some() || beta_opt.is_some()) {
        eprintln!("--auto-model selects class weights; do not combine it with --cpos/--cneg");
        std::process::exit(2);
    }
    // Class weights: pinning either flag pins both (the other takes its default);
    // otherwise --select-c chooses between the per-file grid search (None) and the
    // fixed defaults, which are what every profile uses unless asked otherwise.
    if alpha_opt.is_some() || beta_opt.is_some() {
        a.params.c_alpha = Some(alpha_opt.unwrap_or(percolator::C_POS_DEFAULT));
        a.params.c_beta = Some(beta_opt.unwrap_or(percolator::C_NEG_DEFAULT));
    } else if select_c {
        a.params.c_alpha = None;
        a.params.c_beta = None;
    } else {
        a.params.c_alpha = Some(percolator::C_POS_DEFAULT);
        a.params.c_beta = Some(percolator::C_NEG_DEFAULT);
    }
    if a.params.mlp_hidden == 0 || a.params.mlp_hidden > 256 {
        eprintln!("invalid --mlp-hidden (use 1..256)");
        std::process::exit(2);
    }
    if a.params.mlp_epochs == 0 || a.params.mlp_epochs > 1000 {
        eprintln!("invalid --mlp-epochs (use 1..1000)");
        std::process::exit(2);
    }
    if !a.params.mlp_learning_rate.is_finite() || a.params.mlp_learning_rate <= 0.0 {
        eprintln!("invalid --mlp-learning-rate (must be finite and >0)");
        std::process::exit(2);
    }
    if !a.params.mlp_l2.is_finite() || a.params.mlp_l2 < 0.0 {
        eprintln!("invalid --mlp-l2 (must be finite and >=0)");
        std::process::exit(2);
    }
    if !a.params.svm_tolerance.is_finite() || a.params.svm_tolerance <= 0.0 {
        eprintln!("invalid --svm-tolerance (must be finite and >0)");
        std::process::exit(2);
    }
    a
}

fn ensemble_input(value: &str) -> Result<(String, String), String> {
    let (engine, path) = value
        .split_once('=')
        .ok_or_else(|| format!("invalid ensemble input '{value}'; use ENGINE=PIN"))?;
    if engine.is_empty() || path.is_empty() {
        return Err(format!("invalid ensemble input '{value}'; use ENGINE=PIN"));
    }
    Ok((engine.to_string(), path.to_string()))
}

struct Row {
    id: String,
    score: f64,
    q: f64,
    pep: f64,
    peptide: String,
    proteins: String,
}

fn write_results(path: &str, mut rows: Vec<Row>) -> std::io::Result<()> {
    rows.sort_unstable_by(|x, y| y.score.partial_cmp(&x.score).unwrap_or(std::cmp::Ordering::Equal));
    let f = File::create(path)?;
    let mut w = BufWriter::with_capacity(1 << 20, f);
    writeln!(w, "PSMId\tscore\tq-value\tposterior_error_prob\tpeptide\tproteinIds")?;
    for r in rows {
        writeln!(
            w,
            "{}\t{:.6}\t{:.6}\t{:.6}\t{}\t{}",
            r.id, r.score, r.q, r.pep, r.peptide, r.proteins
        )?;
    }
    Ok(())
}

fn write_feature_report(path: &str, report: &percolator::FeatureReport) -> std::io::Result<()> {
    let file = File::create(path)?;
    let mut writer = BufWriter::with_capacity(1 << 16, file);
    writeln!(writer, "# feature_report_version=1")?;
    writeln!(writer, "# model=linear_svm; coefficients are means across the three out-of-fold models")?;
    writeln!(writer, "# baseline_target_psms_q<0.01={}", report.baseline_q01)?;
    writeln!(writer, "# permutation=deterministic within each held-out fold; models held fixed (no retraining)")?;
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

fn main() {
    let args = parse_args();
    if args.pins.is_empty() {
        eprintln!("usage: percolator-rs [flags] input.pin [more.pin ...]");
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
    eprintln!(
        "profile: {} (model={}, maxiter={}, subset-max-train={}){}{}{}",
        args.profile,
        args.params.model.label(),
        args.params.maxiter,
        if args.params.subset_max_train == 0 { "none".to_string() } else { args.params.subset_max_train.to_string() },
        if args.join { format!(", join={} files", args.pins.len()) }
        else if args.ensemble { format!(", ensemble={} engines", args.pins.len()) }
        else { String::new() },
        if args.rt_features { ", rt-features" } else { "" },
        if args.params.nested_selection { ", nested-selection" } else { "" },
    );
    let t0 = std::time::Instant::now();
    let tp = std::time::Instant::now();
    let ensemble_inputs: Vec<(String, String)> = if args.ensemble {
        args.pins.iter().map(|input| ensemble_input(input)).collect::<Result<_, _>>().unwrap_or_else(|message| {
            eprintln!("error: {message}");
            std::process::exit(2);
        })
    } else {
        args.pins.iter().map(|path| (String::new(), path.clone())).collect()
    };
    let mut parts: Vec<pin::Dataset> = Vec::with_capacity(ensemble_inputs.len());
    for (_, path) in &ensemble_inputs {
        let mut d = pin::parse(path).unwrap_or_else(|e| {
            eprintln!("parse error ({path}): {e}");
            std::process::exit(1);
        });
        if args.rt_features {
            rt::augment(&mut d); // per-file RT alignment, then pool
        }
        parts.push(d);
    }
    let ds = if args.ensemble {
        pin::merge_ensemble(parts, ensemble_inputs.into_iter().map(|(engine, _)| engine).collect())
            .unwrap_or_else(|message| {
                eprintln!("ensemble error: {message}");
                std::process::exit(2);
            })
    } else if parts.len() == 1 {
        parts.pop().unwrap()
    } else {
        pin::merge(parts)
    };
    eprintln!("parse: {:.3}s", tp.elapsed().as_secs_f64());
    eprintln!(
        "parsed {} PSMs, {} features ({} targets / {} decoys){}",
        ds.n_psm,
        ds.n_feat,
        ds.labels.iter().filter(|&&l| l > 0).count(),
        ds.labels.iter().filter(|&&l| l < 0).count(),
        if args.join { format!(", pooled from {} files", ds.source_names.len()) }
        else if args.ensemble { format!(", ensemble from {} engines", ds.source_names.len()) }
        else { String::new() },
    );

    let out = percolator::run(&ds, &args.params);
    if out.nested_folds.is_empty() {
        eprintln!(
            "{} class weights: Cpos={:.3} Cneg={:.3}{}{}",
            if args.params.model == Model::Svm { "SVM" } else { "MLP" },
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
        write_feature_report(path, &report).unwrap_or_else(|error| {
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
    let reported_indices: Vec<usize> = if args.ensemble {
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
    let (reported_qvals, reported_peps) = if args.ensemble {
        let reported_scores: Vec<f64> = reported_indices.iter().map(|&i| out.score[i]).collect();
        let reported_labels: Vec<i8> = reported_indices.iter().map(|&i| ds.labels[i]).collect();
        let reported_pi0 = stats::estimate_pi0(&reported_labels);
        (
            stats::qvalues(&reported_scores, &reported_labels, reported_pi0),
            stats::peps(&reported_scores, &reported_labels, reported_pi0),
        )
    } else {
        (out.qval.clone(), out.pep.clone())
    };

    // PSM-level output
    let mut targets: Vec<Row> = Vec::new();
    let mut decoys: Vec<Row> = Vec::new();
    for (output_index, &i) in reported_indices.iter().enumerate() {
        let r = Row {
            id: if args.ensemble {
                format!("{}:{}", ds.source_names[ds.source[i] as usize], ds.spec_id[i])
            } else {
                ds.spec_id[i].clone()
            },
            score: out.score[i],
            q: reported_qvals[output_index],
            pep: reported_peps[output_index],
            peptide: ds.peptide[i].clone(),
            proteins: ds.proteins[i].clone(),
        };
        if ds.labels[i] > 0 {
            targets.push(r);
        } else {
            decoys.push(r);
        }
    }
    let n_psm_q01 = targets.iter().filter(|r| r.q < 0.01).count();

    // Peptide-level: best-scoring PSM per unique (core) peptide, then re-q-value.
    let mut best: std::collections::HashMap<String, usize> = std::collections::HashMap::new();
    for &i in &reported_indices {
        let key = format!("{}\u{1}{}", ds.labels[i], core_peptide(&ds.peptide[i]));
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
    pep_idx.sort_unstable();
    let pscore: Vec<f64> = pep_idx.iter().map(|&i| out.score[i]).collect();
    let plabel: Vec<i8> = pep_idx.iter().map(|&i| ds.labels[i]).collect();
    let ppi0 = stats::estimate_pi0(&plabel);
    let pq = stats::qvalues(&pscore, &plabel, ppi0);
    let ppep = stats::peps(&pscore, &plabel, ppi0);

    let mut ptargets: Vec<Row> = Vec::new();
    let mut pdecoys: Vec<Row> = Vec::new();
    for (k, &i) in pep_idx.iter().enumerate() {
        let r = Row {
            id: ds.spec_id[i].clone(),
            score: pscore[k],
            q: pq[k],
            pep: ppep[k],
            peptide: ds.peptide[i].clone(),
            proteins: ds.proteins[i].clone(),
        };
        if ds.labels[i] > 0 {
            ptargets.push(r);
        } else {
            pdecoys.push(r);
        }
    }
    let n_pep_q01 = ptargets.iter().filter(|r| r.q < 0.01).count();

    if let Some(p) = &args.results_psms {
        write_results(p, targets).unwrap();
    }
    if let Some(p) = &args.decoy_psms {
        write_results(p, decoys).unwrap();
    }
    if let Some(p) = &args.results_peptides {
        write_results(p, ptargets).unwrap();
    }
    if let Some(p) = &args.decoy_peptides {
        write_results(p, pdecoys).unwrap();
    }

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

    // Protein inference uses the best PSM for each peptide sequence.
    if args.results_proteins.is_some() || args.decoy_proteins.is_some() {
        let entries: Vec<(f64, f64, String)> = pep_idx
            .iter()
            .enumerate()
            .map(|(k, &i)| (pscore[k], ppep[k], ds.proteins[i].clone()))
            .collect();
        let picked_groups = protein::infer(&entries);
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
            write_proteins(p, &groups, false).unwrap();
        }
        if let Some(p) = &args.decoy_proteins {
            write_proteins(p, &groups, true).unwrap();
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
}

fn write_proteins(path: &str, groups: &[protein::ProtGroup], want_decoy: bool) -> std::io::Result<()> {
    let f = File::create(path)?;
    let mut w = BufWriter::with_capacity(1 << 20, f);
    writeln!(w, "ProteinGroupId\tq-value\tposterior_error_prob\tscore\tnumPeptides\tproteinIds")?;
    for g in groups.iter().filter(|g| g.picked && g.is_decoy == want_decoy) {
        writeln!(
            w,
            "{}\t{:.6}\t{:.6}\t{:.6}\t{}\t{}",
            g.proteins.first().map(|s| s.as_str()).unwrap_or(""),
            g.qval,
            g.pep,
            g.score,
            g.n_peptides,
            g.proteins.join(",")
        )?;
    }
    Ok(())
}
