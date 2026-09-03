//! Command-line parsing and validation.

use percolator_rs::percolator::{self, Model, Params};
use percolator_rs::protein_bayes;

pub(crate) struct Args {
    pub(crate) pins: Vec<String>,
    pub(crate) results_psms: Option<String>,
    pub(crate) decoy_psms: Option<String>,
    pub(crate) results_peptides: Option<String>,
    pub(crate) decoy_peptides: Option<String>,
    pub(crate) results_proteins: Option<String>,
    pub(crate) decoy_proteins: Option<String>,
    pub(crate) feature_report: Option<String>,
    pub(crate) protein_inference: ProteinInference,
    pub(crate) protein_bayes: protein_bayes::Params,
    pub(crate) join: bool,
    pub(crate) ensemble: bool,
    /// Perform spectrum-level target-decoy competition on the rescored values
    /// before reporting PSM statistics.
    pub(crate) psm_competition: bool,
    pub(crate) rt_features: bool,
    #[cfg(feature = "profiling")]
    pub(crate) profile_json: Option<String>,
    #[cfg(feature = "profiling")]
    pub(crate) profile_cpu: Option<String>,
    #[cfg(feature = "profiling")]
    pub(crate) profile_allocations: bool,
    pub(crate) params: Params,
    pub(crate) profile: &'static str,
}

#[derive(Clone, Copy, PartialEq)]
pub(crate) enum ProteinInference {
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
        "fast" => Some((20_000, 5, false)), // ~quick QA/test pipelines
        "balanced" => Some((40_000, 10, false)), // ~5% yield hit, still fast
        "canonical" => Some((0, 10, false)), // full default sensitivity (0 = no subsetting)
        _ => None,
    }
}

pub(crate) fn parse_args() -> Args {
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
        #[cfg(feature = "profiling")]
        profile_json: None,
        #[cfg(feature = "profiling")]
        profile_cpu: None,
        #[cfg(feature = "profiling")]
        profile_allocations: false,
        params: Params::default(),
        psm_competition: true,
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
            "--protein-max-iter" => a.protein_bayes.max_iter = take().parse().unwrap_or(0),
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
            "--svm-tolerance" => a.params.svm_tolerance = take().parse().unwrap_or(f64::NAN),
            "--join" => a.join = true,
            "--psm-competition" => a.psm_competition = true,
            "--no-psm-competition" => a.psm_competition = false,
            "--ensemble" => a.ensemble = true,
            "--rt-features" => a.rt_features = true,
            #[cfg(feature = "profiling")]
            "--profile-json" => a.profile_json = Some(take()),
            #[cfg(feature = "profiling")]
            "--profile-cpu" => a.profile_cpu = Some(take()),
            #[cfg(feature = "profiling")]
            "--profile-allocations" => a.profile_allocations = true,
            #[cfg(not(feature = "profiling"))]
            "--profile-json" | "--profile-cpu" | "--profile-allocations" => {
                if s != "--profile-allocations" {
                    let _ = take();
                }
                eprintln!("{s} requires a build with --features profiling");
                std::process::exit(2);
            }
            "--seed" => a.params.seed = take().parse().unwrap_or(1),
            "--null-target-win-prob" => {
                a.params.null_target_win_prob = take().parse().unwrap_or(f64::NAN)
            }
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
    // p = 0 would declare that an incorrect target never outranks its decoy, so
    // every decoy count would convert to zero expected false targets and every
    // q-value and PEP would be exactly zero. That is not a conservative setting,
    // it is a disabled estimator, so the open interval is the contract.
    if !a.params.null_target_win_prob.is_finite()
        || a.params.null_target_win_prob <= 0.0
        || a.params.null_target_win_prob >= 1.0
    {
        eprintln!("invalid --null-target-win-prob (must be finite and in (0, 1))");
        std::process::exit(2);
    }
    a
}

pub(crate) fn ensemble_input(value: &str) -> Result<(String, String), String> {
    let (engine, path) = value
        .split_once('=')
        .ok_or_else(|| format!("invalid ensemble input '{value}'; use ENGINE=PIN"))?;
    if engine.is_empty() || path.is_empty() {
        return Err(format!("invalid ensemble input '{value}'; use ENGINE=PIN"));
    }
    Ok((engine.to_string(), path.to_string()))
}
