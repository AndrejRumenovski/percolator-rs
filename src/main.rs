//! percolator-rs — a from-scratch Rust reimplementation of the Percolator
//! semi-supervised PSM rescoring algorithm. CLI mirrors the subset of reference
//! flags used by our benchmark.

mod percolator;
mod pin;
mod protein;
mod protein_bayes;
mod rt;
mod simd;
mod stats;
mod svm;

use percolator::Params;
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
    protein_inference: ProteinInference,
    protein_bayes: protein_bayes::Params,
    join: bool,
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
/// measured on PXD032157 it costs ~3x wall time and does not beat the fixed default weights
/// (33 files better, 28 worse, aggregate slightly worse). Opt in with --select-c on data where
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
        protein_inference: ProteinInference::Picked,
        protein_bayes: protein_bayes::Params::default(),
        join: false,
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
            "--join" => a.join = true,
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
    a
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

fn main() {
    let args = parse_args();
    if args.pins.is_empty() {
        eprintln!("usage: percolator-rs [flags] input.pin [more.pin ...]");
        std::process::exit(2);
    }
    if args.pins.len() > 1 && !args.join {
        eprintln!("error: multiple inputs require --join (pooled cross-run training), or run them separately");
        std::process::exit(2);
    }
    eprintln!(
        "profile: {} (maxiter={}, subset-max-train={}){}{}",
        args.profile,
        args.params.maxiter,
        if args.params.subset_max_train == 0 { "none".to_string() } else { args.params.subset_max_train.to_string() },
        if args.join { format!(", join={} files", args.pins.len()) } else { String::new() },
        if args.rt_features { ", rt-features" } else { "" },
    );
    let t0 = std::time::Instant::now();
    let tp = std::time::Instant::now();
    let mut parts: Vec<pin::Dataset> = Vec::with_capacity(args.pins.len());
    for path in &args.pins {
        let mut d = pin::parse(path).unwrap_or_else(|e| {
            eprintln!("parse error ({path}): {e}");
            std::process::exit(1);
        });
        if args.rt_features {
            rt::augment(&mut d); // per-file RT alignment, then pool
        }
        parts.push(d);
    }
    let ds = if parts.len() == 1 { parts.pop().unwrap() } else { pin::merge(parts) };
    eprintln!("parse: {:.3}s", tp.elapsed().as_secs_f64());
    eprintln!(
        "parsed {} PSMs, {} features ({} targets / {} decoys){}",
        ds.n_psm,
        ds.n_feat,
        ds.labels.iter().filter(|&&l| l > 0).count(),
        ds.labels.iter().filter(|&&l| l < 0).count(),
        if args.join { format!(", pooled from {} files", ds.source_names.len()) } else { String::new() },
    );

    let out = percolator::run(&ds, &args.params);
    eprintln!(
        "SVM class weights: Cpos={:.3} Cneg={:.3}{}",
        out.c_alpha,
        out.c_beta,
        if out.c_selected { " (selected by cross-validation)" } else { " (fixed)" }
    );

    // PSM-level output
    let mut targets: Vec<Row> = Vec::new();
    let mut decoys: Vec<Row> = Vec::new();
    for i in 0..ds.n_psm {
        let r = Row {
            id: ds.spec_id[i].clone(),
            score: out.score[i],
            q: out.qval[i],
            pep: out.pep[i],
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
    for i in 0..ds.n_psm {
        let key = format!("{}\u{1}{}", ds.labels[i], core_peptide(&ds.peptide[i]));
        match best.get(&key) {
            Some(&j) if out.score[j] >= out.score[i] => {}
            _ => {
                best.insert(key, i);
            }
        }
    }
    let pep_idx: Vec<usize> = best.values().copied().collect();
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
