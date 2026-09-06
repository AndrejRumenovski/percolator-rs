//! Storage orchestration around the existing executable, with independent input runs.
//! No batching/pooling of statistical models is performed here.
use super::{
    cache::Cache,
    download::{hash_file, Budgets, Downloader},
    *,
};
use std::{
    collections::{BTreeMap, BTreeSet},
    fs::{self, File},
    io::{BufRead, BufReader},
    path::Path,
    process::{Command, Stdio},
    sync::atomic::Ordering,
    time::Duration,
};

pub struct RunOptions {
    pub ephemeral: bool,
    pub independent_runs: bool,
    pub batch_size: usize,
    pub pin_retention: Retention,
    pub allow_unverified: bool,
    pub result_bytes_per_input: u64,
    pub analysis_args: Vec<String>,
}
impl Default for RunOptions {
    fn default() -> Self {
        Self {
            ephemeral: false,
            independent_runs: false,
            batch_size: 1,
            pin_retention: Retention::KeepIfPinned,
            allow_unverified: false,
            result_bytes_per_input: 64 * 1024 * 1024,
            analysis_args: vec![],
        }
    }
}
/// Prevent passing filenames, unknown switches, or output locations into the legacy
/// permissive CLI parser. Scientific options are otherwise delegated unchanged.
pub fn validate_analysis_args(args: &[String]) -> Result<()> {
    let flags = [
        "--select-c",
        "--no-select-c",
        "--fast",
        "--balanced",
        "--canonical",
        "--auto-model",
        "--nested-select",
        "--no-auto-model",
        "--psm-competition",
        "--no-psm-competition",
        "--rt-features",
    ];
    let values = [
        "--seed",
        "--maxiter",
        "--subset-max-train",
        "--cpos",
        "--cneg",
        "--num-threads",
        "--profile",
        "--null-target-win-prob",
        "--model",
        "--rescore-model",
        "--mlp-hidden",
        "--mlp-epochs",
        "--mlp-learning-rate",
        "--mlp-l2",
        "--svm-tolerance",
    ];
    let mut i = 0;
    while i < args.len() {
        if flags.contains(&args[i].as_str()) {
            i += 1;
            continue;
        }
        if values.contains(&args[i].as_str()) {
            if args
                .get(i + 1)
                .is_none_or(|s| s.starts_with("--") || s.is_empty())
            {
                return Err(format!("missing analysis value for {}", args[i]).into());
            }
            // The existing parser has historical fallback defaults. Refuse malformed
            // values here instead of recording a request different from execution.
            let v = &args[i + 1];
            let valid = match args[i].as_str() {
                "--profile" => matches!(v.as_str(), "fast" | "balanced" | "canonical"),
                "--model" | "--rescore-model" => matches!(v.as_str(), "svm" | "linear" | "mlp"),
                "--seed" | "--subset-max-train" => v.parse::<u64>().is_ok(),
                "--maxiter" | "--num-threads" | "--mlp-hidden" | "--mlp-epochs" => {
                    v.parse::<usize>().is_ok_and(|x| x > 0)
                }
                _ => v.parse::<f64>().is_ok_and(|x| x.is_finite()),
            };
            if !valid {
                return Err(format!("invalid analysis value for {}: {v:?}", args[i]).into());
            }
            i += 2;
        } else {
            return Err(format!("analysis argument {:?} is unsupported here; inputs/output paths and joint/ensemble models must use the ordinary CLI",args[i]).into());
        }
    }
    Ok(())
}
pub fn validate_run(files: &[RemoteFile], options: &RunOptions) -> Result<()> {
    validate_analysis_args(&options.analysis_args)?;
    if files.len() > 1 && !options.independent_runs {
        return Err("multiple PINs require --independent-runs, explicitly choosing separate models/statistics; use the ordinary CLI for --join/--ensemble".into());
    }
    if options.batch_size == 0 {
        return Err("batch size must be positive".into());
    }
    if options.result_bytes_per_input < 4096 {
        return Err("result budget must be at least 4096 bytes per input".into());
    }
    let mut unique = BTreeSet::new();
    for f in files {
        if !unique.insert(f.object_key()) {
            return Err("selected PINs share a content identity; select each source object once for independent analysis".into());
        }
    }
    for f in files {
        if !f.native_pin() {
            return Err(format!("{:?}: {}", f.filename, f.preparation()).into());
        }
        if f.checksums.is_empty() && !options.allow_unverified {
            return Err(format!("{:?}: repository checksum unavailable; --allow-unverified explicitly accepts local-hash-only provenance",f.filename).into());
        }
    }
    Ok(())
}
/// All output files are bounded independently in the child; the sum of the four
/// per-file ceilings is <= the reserved result budget. No external search/converter
/// is launched and no source/PIN is copied into results.
#[cfg(unix)]
fn limit_child(command: &mut Command, per_file: u64) -> Result<()> {
    use std::os::unix::process::CommandExt;
    let parent = std::process::id();
    unsafe {
        command.pre_exec(move || {
            let limit = libc::rlimit {
                rlim_cur: per_file as libc::rlim_t,
                rlim_max: per_file as libc::rlim_t,
            };
            if libc::setrlimit(libc::RLIMIT_FSIZE, &limit) != 0 {
                return Err(std::io::Error::last_os_error());
            }
            let no_core = libc::rlimit {
                rlim_cur: 0,
                rlim_max: 0,
            };
            if libc::setrlimit(libc::RLIMIT_CORE, &no_core) != 0 {
                return Err(std::io::Error::last_os_error());
            }
            #[cfg(target_os = "linux")]
            {
                if libc::prctl(libc::PR_SET_PDEATHSIG, libc::SIGKILL) != 0 {
                    return Err(std::io::Error::last_os_error());
                }
                if libc::getppid() != parent as libc::pid_t {
                    return Err(std::io::Error::other(
                        "PRIDE parent exited before analysis launch",
                    ));
                }
            }
            Ok(())
        });
    }
    Ok(())
}
#[cfg(not(unix))]
fn limit_child(_: &mut Command, _: u64) -> Result<()> {
    Err("bounded PRIDE analysis currently requires Unix RLIMIT_FSIZE; metadata/download/cache operations are available".into())
}

pub fn run(
    cache: &mut Cache,
    m: &mut Manifest,
    files: &[RemoteFile],
    downloader: &mut Downloader,
    budgets: &Budgets,
    options: &RunOptions,
    executable: &Path,
) -> Result<()> {
    validate_run(files, options)?;
    let total_results = total(files.iter().map(|_| options.result_bytes_per_input))?;
    let release = options.ephemeral && options.pin_retention != Retention::Keep;
    let plan = download::plan(
        cache,
        m,
        files,
        budgets,
        release,
        options.batch_size,
        total_results,
    )?;
    cache.evict(&plan.expected_evictions, false)?;
    m.selected_files = files.iter().map(|f| f.id.clone()).collect();
    cache.save_manifest(m)?;
    for batch in files.chunks(options.batch_size) {
        let protected: BTreeSet<_> = batch.iter().map(RemoteFile::object_key).collect();
        let mut keys = Vec::new();
        for f in batch {
            let key = downloader.fetch(cache, m, f, &protected)?;
            let o = cache.index.objects.get_mut(&key).unwrap();
            // A KEEP reference is never weakened by another workflow.
            if o.retention != Retention::Keep {
                o.retention = options.pin_retention;
            }
            if options.pin_retention == Retention::UntilResultVerified {
                o.result_verified = false;
            }
            cache.save_index()?;
            keys.push(key);
        }
        for (f, key) in batch.iter().zip(&keys) {
            process_one(
                cache,
                m,
                &f.id,
                key,
                &downloader.cancelled,
                budgets,
                options,
                executable,
            )?;
        }

        if options.ephemeral {
            let evict: Vec<_> = keys
                .into_iter()
                .filter(|k| cache.evictable(&cache.index.objects[k], false))
                .collect();
            cache.evict(&evict, false)?;
            // evict updates all on-disk manifests; reload to preserve availability state.
            *m = cache.load_manifest(&m.accession)?;
        }
    }
    Ok(())
}
#[allow(clippy::too_many_arguments)]
fn process_one(
    cache: &mut Cache,
    m: &mut Manifest,
    input_id: &str,
    key: &str,
    cancelled: &std::sync::Arc<std::sync::atomic::AtomicBool>,
    budgets: &Budgets,
    options: &RunOptions,
    executable: &Path,
) -> Result<()> {
    let exe_hash = hash_file(executable)?.sha256;
    let protected = m
        .local_files
        .values()
        .map(|f| f.object_key.clone())
        .collect();
    let evictions = cache.eviction_plan(
        options.result_bytes_per_input,
        total([options.result_bytes_per_input, budgets.safety])?,
        &protected,
    )?;
    cache.evict(&evictions, false)?;
    let timestamp = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)?
        .as_nanos();
    let id = format!("{}-{timestamp}", m.accession);
    let dir = cache.path(&format!("results/{id}"))?;
    fs::create_dir(&dir)?;
    let input = cache.path(&cache.index.objects[key].relative_path)?;
    let input_sha = cache.index.objects[key].local_sha256.clone();
    let experiment = Experiment {
        id: id.clone(),
        input_ids: vec![input_id.to_owned()],
        state: "running".into(),
        error: None,
        started_unix_seconds: now(),
        completed_unix_seconds: None,
        percolator_rs_version: env!("CARGO_PKG_VERSION").into(),
        percolator_rs_commit: option_env!("PERCOLATOR_RS_BUILD_COMMIT").map(str::to_owned),
        executable_sha256: exe_hash.clone(),
        parameters: options.analysis_args.clone(),
        ephemeral: options.ephemeral,
        pin_retention: options.pin_retention,
        result_hashes: BTreeMap::new(),
        lineage: vec![Lineage {
            id: format!("{id}:input-pin"),
            inputs: vec![input_id.to_owned()],
            output_sha256: input_sha,
            kind: if key.starts_with("prepared-") {
                "prepared_pin"
            } else {
                "downloaded_pin"
            }
            .into(),
            tool: if key.starts_with("prepared-") {
                "external preparation; see manifest lineage"
            } else {
                "PRIDE HTTPS download"
            }
            .into(),
            tool_version: None,
            parameters: vec![],
            protein_database: None,
            database_sha256: None,
            decoy_generation: None,
        }],
    };
    m.experiments.push(experiment);
    cache.save_manifest(m)?;
    let outcome = (|| -> Result<BTreeMap<String, String>> {
        if cancelled.load(Ordering::Relaxed) {
            return Err("analysis interrupted".into());
        }
        if cache::available_space(&cache.root)?
            < total([options.result_bytes_per_input, budgets.safety])?
        {
            return Err("insufficient free space for bounded analysis outputs".into());
        }
        // Validate using the exact existing reader. It checks labels/features;
        // database design cannot be established from a file extension.
        drop(crate::pin::parse(
            input.to_str().ok_or("non-UTF8 PIN path")?,
        )?);
        cache.index.objects.get_mut(key).unwrap().pin_validated = true;
        let local = cache.index.objects[key].local();
        m.local_files.insert(input_id.to_owned(), local);
        cache.save_index()?;
        cache.save_manifest(m)?;
        let names = [
            ("--results-psms", "target.psms.tsv"),
            ("--decoy-results-psms", "decoy.psms.tsv"),
            ("--results-peptides", "target.peptides.tsv"),
            ("--decoy-results-peptides", "decoy.peptides.tsv"),
        ];
        let mut staging = BTreeMap::new();
        for (_, name) in names {
            use sha2::{Digest, Sha256};
            let temp_key = format!("remote-{:x}", Sha256::digest(format!("{id}/{name}")));
            let relative_path = format!("tmp/{temp_key}.part");
            cache.index.objects.insert(
                temp_key.clone(),
                cache::Object {
                    key: temp_key.clone(),
                    relative_path: relative_path.clone(),
                    bytes: options.result_bytes_per_input / 4,
                    state: State::Partial,
                    local_sha256: None,
                    verification: vec![],
                    projects: BTreeSet::from([m.accession.clone()]),
                    last_used_unix_seconds: now(),
                    retention: Retention::Evict,
                    result_verified: false,
                    reproducible: true,
                    etag: None,
                    pin_validated: false,
                },
            );
            staging.insert(name, (temp_key, cache.path(&relative_path)?));
        }
        cache.save_index()?;
        let mut command = Command::new(executable);
        command
            .current_dir(&cache.root)
            .args(&options.analysis_args)
            .arg(&input)
            .stdin(Stdio::null())
            .stdout(Stdio::null())
            .stderr(Stdio::inherit());
        for (flag, name) in names {
            command.arg(flag).arg(&staging[name].1);
        }
        limit_child(&mut command, options.result_bytes_per_input / 4)?;
        let mut child = command.spawn()?;
        let status = loop {
            if let Some(s) = child.try_wait()? {
                break s;
            }
            if cancelled.load(Ordering::Relaxed) {
                let _ = child.kill();
                let _ = child.wait();
                return Err("analysis interrupted".into());
            }
            std::thread::sleep(Duration::from_millis(100));
        };
        if !status.success() {
            return Err(format!("percolator-rs failed ({status}); outputs incomplete; source retained. A signal may indicate the configured result-size ceiling.").into());
        }
        let mut hashes = BTreeMap::new();
        for (_, name) in names {
            let path = &staging[name].1;
            validate_result(path)?;
            let hash = hash_file(path)?;
            if hash.bytes > options.result_bytes_per_input / 4 {
                return Err("analysis output exceeded its reservation".into());
            }
            File::open(path)?.sync_all()?;
            hashes.insert(format!("results/{id}/{name}"), hash.sha256);
        }
        // Publish only after all four outputs have passed validation.
        for (_, name) in names {
            let (temp_key, path) = &staging[name];
            fs::rename(path, dir.join(name))?;
            cache.index.objects.get_mut(temp_key).unwrap().state = State::Evicted;
        }
        cache.save_index()?;
        File::open(&dir)?.sync_all()?;
        Ok(hashes)
    })();
    let e = m.experiments.last_mut().unwrap();
    match outcome {
        Ok(hashes) => {
            for (path, sha) in &hashes {
                e.lineage.push(Lineage {
                    id: path.clone(),
                    inputs: vec![format!("{id}:input-pin")],
                    output_sha256: Some(sha.clone()),
                    kind: "percolator_result".into(),
                    tool: "percolator-rs".into(),
                    tool_version: Some(env!("CARGO_PKG_VERSION").into()),
                    parameters: options.analysis_args.clone(),
                    protein_database: None,
                    database_sha256: None,
                    decoy_generation: None,
                });
            }
            e.state = "verified".into();
            e.completed_unix_seconds = Some(now());
            e.result_hashes = hashes;
            cache.save_manifest(m)?;
            cache.index.objects.get_mut(key).unwrap().result_verified = true;
            cache.save_index()?;
            eprintln!("verified results: {}", dir.display());
        }
        Err(error) => {
            e.state = "failed".into();
            e.error = Some(error.to_string());
            cache.save_manifest(m)?;
            return Err(error);
        }
    }
    Ok(())
}

pub fn run_prepared(
    cache: &mut Cache,
    m: &mut Manifest,
    id: &str,
    cancelled: &std::sync::Arc<std::sync::atomic::AtomicBool>,
    budgets: &Budgets,
    options: &RunOptions,
    executable: &Path,
) -> Result<()> {
    let plan = super::prepare::run_plan(cache, m, id, budgets, options.result_bytes_per_input)?;
    validate_analysis_args(&options.analysis_args)?;
    cache.evict(&plan.expected_evictions, false)?;
    let artifact = m
        .prepared_pins
        .get(id)
        .ok_or("unknown prepared PIN")?
        .clone();
    let o = cache
        .index
        .objects
        .get(&artifact.object_key)
        .ok_or("prepared PIN evicted; regenerate using retained lineage")?;
    let path = cache.path(&o.relative_path)?;
    let hashes = hash_file(&path)?;
    if hashes.sha256 != artifact.sha256 || hashes.bytes != artifact.bytes {
        cache
            .index
            .objects
            .get_mut(&artifact.object_key)
            .unwrap()
            .state = State::Corrupt;
        cache.save_index()?;
        return Err("prepared PIN is corrupt; regenerate using retained lineage".into());
    }
    let o = cache.index.objects.get_mut(&artifact.object_key).unwrap();
    if o.retention != Retention::Keep {
        o.retention = options.pin_retention;
    }
    if o.retention == Retention::UntilResultVerified {
        o.result_verified = false;
    }
    m.prepared_pins.get_mut(id).unwrap().retention = o.retention;
    cache.save_index()?;
    cache.save_manifest(m)?;
    process_one(
        cache,
        m,
        id,
        &artifact.object_key,
        cancelled,
        budgets,
        options,
        executable,
    )?;
    if options.ephemeral && cache.evictable(&cache.index.objects[&artifact.object_key], false) {
        cache.evict(&[artifact.object_key], false)?;
        *m = cache.load_manifest(&m.accession)?;
    }
    Ok(())
}

fn validate_result(path: &Path) -> Result<()> {
    let mut rows = BufReader::new(File::open(path)?).lines();
    if rows.next().transpose()?.as_deref()
        != Some("PSMId\tscore\tq-value\tposterior_error_prob\tpeptide\tproteinIds")
    {
        return Err("missing/invalid final result header".into());
    }
    for row in rows {
        let row = row?;
        let cols: Vec<_> = row.split('\t').collect();
        if cols.len() < 6 {
            return Err("truncated final result row".into());
        }
        for (i, c) in cols.iter().enumerate().take(4).skip(1) {
            let n: f64 = c.parse()?;
            if !n.is_finite() || (i >= 2 && !(0.0..=1.0).contains(&n)) {
                return Err("invalid numeric final result value".into());
            }
        }
    }
    Ok(())
}
