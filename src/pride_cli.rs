//! PRIDE CLI dispatch, deliberately separate from the scientific CLI contract.
use clap::{Args, Parser, Subcommand, ValueEnum};
use percolator_rs::pride::{
    self,
    cache::Cache,
    client::PrideClient,
    download::{self, Budgets, Downloader},
    workflow::{self, RunOptions},
    *,
};
use std::{
    collections::BTreeSet,
    path::PathBuf,
    sync::{atomic::AtomicBool, Arc},
};

#[derive(Parser)]
#[command(
    name = "percolator-rs pride",
    about = "Public PRIDE datasets in a bounded, disposable working cache"
)]
struct Cli {
    #[arg(long, global = true, env = "PERCOLATOR_RS_PRIDE_CACHE")]
    cache_dir: Option<PathBuf>,
    #[arg(long,global=true,env="PERCOLATOR_RS_PRIDE_CACHE_LIMIT",default_value="50GB",value_parser=pride::bytes)]
    cache_limit: u64,
    #[arg(
        long,
        global = true,
        help = "Show a plan without modifying files, metadata, or cache state"
    )]
    dry_run: bool,
    #[command(subcommand)]
    command: Action,
}
#[derive(Subcommand)]
enum Action {
    /// Inspect metadata, storage costs and preparation options; no source downloads.
    Info(Inspect),
    /// Show the complete indexed inventory plus checksum-table entries.
    Files {
        #[command(flatten)]
        inspect: Inspect,
        #[command(flatten)]
        selection: Selection,
    },
    /// Print the versioned JSON manifest, including retained provenance.
    Manifest(Inspect),
    /// Discover public projects with the official v3 keyword/filter interface.
    Search {
        #[arg(default_value = "")]
        keyword: String,
        #[arg(long)]
        filter: Option<String>,
        #[arg(long, default_value_t = 0)]
        page: u32,
        #[arg(long, default_value_t = 20)]
        page_size: u32,
    },
    /// Plan or execute a bounded verified download (execution requires --yes).
    Fetch {
        #[command(flatten)]
        inspect: Inspect,
        #[command(flatten)]
        selection: Selection,
        #[command(flatten)]
        budget: BudgetArgs,
        #[arg(long)]
        yes: bool,
    },
    /// Import a validated external PIN with a JSON preparation/search recipe.
    Prepare {
        #[command(flatten)]
        inspect: Inspect,
        #[arg(long)]
        pin: PathBuf,
        #[arg(long)]
        recipe: PathBuf,
        #[arg(long, value_enum, default_value = "keep-if-pinned")]
        retention: PinRetention,
        #[command(flatten)]
        budget: BudgetArgs,
        #[arg(long)]
        yes: bool,
    },
    /// Run existing percolator-rs analysis on selected uncompressed PINs.
    Run {
        #[command(flatten)]
        inspect: Inspect,
        #[command(flatten)]
        selection: Selection,
        #[command(flatten)]
        budget: BudgetArgs,
        #[arg(
            long,
            help = "Use an imported prepared PIN lineage ID instead of a remote selection"
        )]
        prepared: Option<String>,
        #[arg(long)]
        yes: bool,
        #[arg(long)]
        ephemeral: bool,
        #[arg(
            long,
            help = "Explicitly choose separate per-file models/statistics for multiple PINs"
        )]
        independent_runs: bool,
        #[arg(long, default_value_t = 1)]
        batch_size: usize,
        #[arg(long, value_enum, default_value = "keep-if-pinned")]
        pin_retention: PinRetention,
        #[arg(
            long,
            help = "Permit analysis with local SHA-256 only when PRIDE publishes no checksum"
        )]
        allow_unverified: bool,
        #[arg(long,default_value="64MiB",value_parser=pride::bytes)]
        max_results_per_input: u64,
        #[arg(
            last = true,
            help = "Existing scientific CLI options, e.g. -- --profile fast --seed 1"
        )]
        analysis_args: Vec<String>,
    },
    /// Inspect, pin and reclaim managed data while preserving provenance and results.
    Cache {
        #[command(subcommand)]
        command: CacheAction,
    },
}
#[derive(Args)]
struct Inspect {
    accession: Pxd,
    #[arg(
        long,
        help = "Refresh official metadata while preserving local processing history"
    )]
    refresh: bool,
    #[arg(long)]
    json: bool,
}
#[derive(Args, Default)]
struct Selection {
    #[arg(
        long = "file",
        help = "Exact inventory ID or filename; repeat to select multiple files"
    )]
    files: Vec<String>,
    #[arg(
        long = "format",
        help = "PIN, mzML, mzIdentML, mzTab, RAW, FASTA, etc.; repeat for alternatives"
    )]
    formats: Vec<String>,
    #[arg(
        long = "category",
        help = "Official category, processed, or search-engine-output"
    )]
    categories: Vec<String>,
    #[arg(long, help = "Explicitly select the entire project inventory")]
    all: bool,
}
#[derive(Args)]
struct BudgetArgs {
    #[arg(long,default_value="1GB",value_parser=pride::bytes,help="Maximum bytes transferred in this operation, including retries")]
    max_download: u64,
    #[arg(long,value_parser=pride::bytes)]
    max_working_space: Option<u64>,
    #[arg(long,default_value="1GB",value_parser=pride::bytes)]
    safety_margin: u64,
}
impl BudgetArgs {
    fn budgets(&self) -> Budgets {
        Budgets {
            max_download: self.max_download,
            max_working_space: self.max_working_space,
            safety: self.safety_margin,
        }
    }
}
#[derive(Clone, Copy, ValueEnum)]
enum PinRetention {
    Keep,
    Evict,
    KeepIfPinned,
    UntilResultVerified,
}
impl From<PinRetention> for Retention {
    fn from(p: PinRetention) -> Self {
        match p {
            PinRetention::Keep => Self::Keep,
            PinRetention::Evict => Self::Evict,
            PinRetention::KeepIfPinned => Self::KeepIfPinned,
            PinRetention::UntilResultVerified => Self::UntilResultVerified,
        }
    }
}
#[derive(Subcommand)]
enum CacheAction {
    Status,
    Pin {
        accession: Pxd,
    },
    Unpin {
        accession: Pxd,
    },
    /// Remove all disposable unpinned data; --dry-run previews exact objects.
    Prune {
        #[arg(
            long,
            help = "Explicitly reclaim every evictable object, down to zero when nothing is protected"
        )]
        all_evictable: bool,
    },
    /// Override KEEP/unfinished-PIN retention; preserve pins, metadata and final results.
    PurgeData {
        #[arg(
            long,
            help = "Confirm deletion of all unpinned source/prepared/temporary data"
        )]
        yes: bool,
    },
    /// Clean tracked partial downloads left by interrupted operations; respects pins.
    CleanAbandoned,
}
fn print_json(value: &impl serde::Serialize) -> Result<()> {
    println!("{}", serde_json::to_string_pretty(value)?);
    Ok(())
}
fn select(m: &Manifest, s: &Selection, default_all: bool) -> Result<Vec<RemoteFile>> {
    for name in &s.files {
        let matched: Vec<_> = m
            .inventory
            .iter()
            .filter(|f| &f.id == name || &f.filename == name)
            .collect();
        if matched.is_empty() {
            return Err(
                format!("no file matches {name:?}; use `pride files` inventory IDs").into(),
            );
        }
        if matched.len() > 1 {
            return Err(
                format!("ambiguous filename {name:?}; select an exact inventory ID").into(),
            );
        }
    }
    let explicit =
        !s.files.is_empty() || !s.formats.is_empty() || !s.categories.is_empty() || s.all;
    let files = m
        .inventory
        .iter()
        .filter(|f| {
            if !explicit {
                return default_all;
            }
            (s.files.is_empty() || s.files.iter().any(|v| v == &f.id || v == &f.filename))
                && (s.formats.is_empty() || s.formats.iter().any(|v| f.matches(v)))
                && (s.categories.is_empty() || s.categories.iter().any(|v| f.matches(v)))
        })
        .cloned()
        .collect();
    Ok(files)
}
fn metadata(cache: &Cache, i: &Inspect, persist: bool) -> Result<Manifest> {
    let path = cache.path(&format!("manifests/{}.json", i.accession))?;
    let old = if path.exists() {
        Some(cache.load_manifest(&i.accession)?)
    } else {
        None
    };
    if !i.refresh {
        if let Some(m) = old {
            return Ok(m);
        }
    }
    let mut m = PrideClient::new()?.manifest(&i.accession)?;
    if let Some(old) = old {
        m.local_files = old.local_files;
        m.prepared_pins = old.prepared_pins;
        m.preparation_attempts = old.preparation_attempts;
        m.experiments = old.experiments;
        m.lineage = old.lineage;
        m.selected_files = old.selected_files;
        // Historical remote identities must survive PRIDE changing/deleting a record.
        for f in old.inventory {
            if !m
                .inventory
                .iter()
                .any(|n| n.id == f.id && n.object_key() == f.object_key())
            {
                m.inventory_notes.push(format!(
                    "Historical remote record retained for provenance: {}",
                    f.id
                ));
                // History belongs to a separate field; do not present stale records as current selections.
                m.remote_history.push(f);
            }
        }
        m.remote_history.extend(old.remote_history);
    }
    if persist {
        cache.save_manifest(&m)?;
    }
    Ok(m)
}
fn info(m: &Manifest) -> Result<()> {
    println!(
        "{}\n{}\nStatus: {} | submission: {} | retrieved: {} Unix seconds",
        m.accession,
        m.project.title.as_deref().unwrap_or("title unavailable"),
        m.project.status.as_deref().unwrap_or("unknown"),
        m.project.submission_date.as_deref().unwrap_or("unknown"),
        m.retrieved_unix_seconds
    );
    if let Some(d) = &m.project.description {
        println!("{d}");
    }
    for (label, terms) in [
        ("Organisms", &m.project.organisms),
        ("Tissues", &m.project.tissues),
        ("Instruments", &m.project.instruments),
        ("Experiment types", &m.project.experiment_types),
        ("Modifications", &m.project.modifications),
    ] {
        println!(
            "{label}: {}",
            terms
                .as_ref()
                .filter(|v| !v.is_empty())
                .map(|ts| ts
                    .iter()
                    .map(|t| t
                        .name
                        .as_deref()
                        .or(t.value.as_deref())
                        .or(t.accession.as_deref())
                        .unwrap_or("unknown"))
                    .collect::<Vec<_>>()
                    .join(", "))
                .unwrap_or_else(|| "unavailable".into())
        );
    }
    println!(
        "Files: {} ({} indexed; {} supplemental)",
        m.inventory.len(),
        m.indexed_file_count,
        m.inventory.len() - m.indexed_file_count
    );
    for (label, files) in [
        ("Total", m.inventory.iter().collect::<Vec<_>>()),
        (
            "RAW",
            m.inventory.iter().filter(|f| f.matches("raw")).collect(),
        ),
        (
            "Processed",
            m.inventory
                .iter()
                .filter(|f| f.matches("processed"))
                .collect(),
        ),
        (
            "Native PIN candidates",
            m.inventory.iter().filter(|f| f.native_pin()).collect(),
        ),
        (
            "Validated directly compatible PINs",
            m.inventory
                .iter()
                .filter(|f| m.compatibility(f) == Compatibility::DirectlyCompatible)
                .collect(),
        ),
    ] {
        let known = total(
            files
                .iter()
                .filter_map(|f| f.size_bytes.or(f.checksum_table_size)),
        )?;
        let missing = files
            .iter()
            .filter(|f| f.size_bytes.or(f.checksum_table_size).is_none())
            .count();
        println!(
            "{label}: {known} bytes ({:.3} GB), {} files, {missing} sizes unavailable",
            known as f64 / 1e9,
            files.len()
        );
    }
    let kinds: BTreeSet<_> = m.inventory.iter().map(RemoteFile::preparation).collect();
    for option in kinds {
        println!("Preparation: {option}");
    }
    for note in &m.inventory_notes {
        println!("Inventory note: {note}");
    }
    for f in &m.inventory {
        if let (Some(a), Some(b)) = (f.size_bytes, f.checksum_table_size) {
            if a != b {
                println!(
                    "Size conflict: {:?}: API {a}, checksum table {b}; download blocked",
                    f.filename
                );
            }
        }
    }
    Ok(())
}
fn cancelled() -> Result<Arc<AtomicBool>> {
    let c = Arc::new(AtomicBool::new(false));
    let signal = c.clone();
    ctrlc::set_handler(move || {
        signal.store(true, std::sync::atomic::Ordering::Relaxed);
    })?;
    Ok(c)
}

pub fn run() -> Result<()> {
    let cli = Cli::parse_from(
        std::iter::once("percolator-rs pride".to_owned()).chain(std::env::args().skip(2)),
    );
    if let Action::Search {
        keyword,
        filter,
        page,
        page_size,
    } = &cli.command
    {
        return print_json(&PrideClient::new()?.search(
            keyword,
            filter.as_deref(),
            *page,
            *page_size,
        )?);
    }
    let execution = match &cli.command {
        Action::Fetch { yes, .. } | Action::Run { yes, .. } | Action::Prepare { yes, .. } => *yes,
        Action::Cache {
            command: CacheAction::PurgeData { yes },
        } => *yes,
        _ => true,
    };
    let readonly = cli.dry_run
        || !execution
        || matches!(
            cli.command,
            Action::Cache {
                command: CacheAction::Status
            }
        );
    let root = match cli.cache_dir {
        Some(root) => root,
        None => Cache::default_root()?,
    };
    let mut cache = Cache::open(&root, cli.cache_limit, !readonly)?;
    match cli.command {
        Action::Info(i) => {
            let m = metadata(&cache, &i, !readonly)?;
            if i.json {
                print_json(&m)
            } else {
                info(&m)
            }
        }
        Action::Manifest(i) => print_json(&metadata(&cache, &i, !readonly)?),
        Action::Files {
            inspect: i,
            selection,
        } => {
            let m = metadata(&cache, &i, !readonly)?;
            let files = select(&m, &selection, true)?;
            if i.json {
                return print_json(&files);
            }
            info(&m)?;
            for f in files {
                println!(
                    "{}\t{:?}\tcategory={}\tformat={}\tbytes={}\t{:?}\tchecksums={}",
                    f.id,
                    f.filename,
                    f.category.as_deref().unwrap_or("unknown"),
                    f.format_name(),
                    f.size_bytes
                        .or(f.checksum_table_size)
                        .map(|x| x.to_string())
                        .unwrap_or_else(|| "unknown".into()),
                    m.compatibility(&f),
                    f.checksums.len()
                );
            }
            Ok(())
        }
        Action::Fetch {
            inspect: i,
            selection,
            budget,
            yes,
        } => {
            let mut m = metadata(&cache, &i, false)?;
            let files = select(&m, &selection, false)?;
            let budgets = budget.budgets();
            let plan = download::plan(&cache, &m, &files, &budgets, false, 1, 0)?;
            print_json(&plan)?;
            if cli.dry_run || !yes {
                eprintln!("Plan only. Use --yes to execute this selection within its budgets.");
                return Ok(());
            }
            cache.save_manifest(&m)?;
            cache.evict(&plan.expected_evictions, false)?;
            m.selected_files = files.iter().map(|f| f.id.clone()).collect();
            cache.save_manifest(&m)?;
            let protected = files.iter().map(RemoteFile::object_key).collect();
            let mut downloader = Downloader::new(&budgets, cancelled()?)?;
            for f in &files {
                downloader.fetch(&mut cache, &mut m, f, &protected)?;
            }
            print_json(&cache.status()?)
        }
        Action::Prepare {
            inspect: i,
            pin,
            recipe,
            retention,
            budget,
            yes,
        } => {
            let mut m = metadata(&cache, &i, false)?;
            let recipe = pride::cache::read_json(&recipe)?;
            let budgets = budget.budgets();
            print_json(&pride::prepare::import_plan(
                &cache, &m, &pin, &recipe, &budgets,
            )?)?;
            if cli.dry_run || !yes {
                eprintln!("Plan only. Use --yes to import this PIN and recipe.");
                return Ok(());
            }
            cache.save_manifest(&m)?;
            let id = pride::prepare::import(
                &mut cache,
                &mut m,
                &pin,
                recipe,
                &budgets,
                retention.into(),
            )?;
            println!("Prepared PIN ID: {id}");
            print_json(&cache.status()?)
        }
        Action::Run {
            inspect: i,
            selection,
            budget,
            yes,
            ephemeral,
            independent_runs,
            batch_size,
            pin_retention,
            allow_unverified,
            max_results_per_input,
            analysis_args,
            prepared,
        } => {
            let mut m = metadata(&cache, &i, false)?;
            let files = select(&m, &selection, false)?;
            let budgets = budget.budgets();
            let options = RunOptions {
                ephemeral,
                independent_runs,
                batch_size,
                pin_retention: pin_retention.into(),
                allow_unverified,
                result_bytes_per_input: max_results_per_input,
                analysis_args,
            };
            if let Some(id) = prepared {
                if selection.all
                    || !selection.files.is_empty()
                    || !selection.formats.is_empty()
                    || !selection.categories.is_empty()
                {
                    return Err("--prepared cannot be combined with a remote file selection".into());
                }
                workflow::validate_analysis_args(&options.analysis_args)?;
                if options.result_bytes_per_input < 4096 {
                    return Err("result budget must be at least 4096 bytes per input".into());
                }
                print_json(&pride::prepare::run_plan(
                    &cache,
                    &m,
                    &id,
                    &budgets,
                    options.result_bytes_per_input,
                )?)?;
                if cli.dry_run || !yes {
                    eprintln!("Plan only. Use --yes to execute prepared PIN analysis.");
                    return Ok(());
                }
                cache.save_manifest(&m)?;
                workflow::run_prepared(
                    &mut cache,
                    &mut m,
                    &id,
                    &cancelled()?,
                    &budgets,
                    &options,
                    &std::env::current_exe()?,
                )?;
                return print_json(&cache.status()?);
            }
            workflow::validate_run(&files, &options)?;
            let result_bytes = total(files.iter().map(|_| max_results_per_input))?;
            print_json(&download::plan(
                &cache,
                &m,
                &files,
                &budgets,
                ephemeral && options.pin_retention != Retention::Keep,
                batch_size,
                result_bytes,
            )?)?;
            if cli.dry_run || !yes {
                eprintln!(
                    "Plan only. Use --yes to execute; each PIN uses its own model and statistics."
                );
                return Ok(());
            }
            cache.save_manifest(&m)?;
            let mut downloader = Downloader::new(&budgets, cancelled()?)?;
            workflow::run(
                &mut cache,
                &mut m,
                &files,
                &mut downloader,
                &budgets,
                &options,
                &std::env::current_exe()?,
            )?;
            print_json(&cache.status()?)
        }
        Action::Cache { command } => match command {
            CacheAction::Status => print_json(&cache.status()?),
            CacheAction::Pin { accession } => {
                if !cli.dry_run {
                    cache.pin(&accession, true)?;
                }
                println!(
                    "{} {accession}",
                    if cli.dry_run { "Would pin" } else { "Pinned" }
                );
                Ok(())
            }
            CacheAction::Unpin { accession } => {
                if !cli.dry_run {
                    cache.pin(&accession, false)?;
                }
                println!(
                    "{} {accession}",
                    if cli.dry_run {
                        "Would unpin"
                    } else {
                        "Unpinned"
                    }
                );
                Ok(())
            }
            CacheAction::Prune { .. } => print_json(&cache.prune(false, cli.dry_run, false)?),
            CacheAction::CleanAbandoned => print_json(&cache.prune(false, cli.dry_run, true)?),
            CacheAction::PurgeData { yes } => {
                let report = cache.prune(true, cli.dry_run || !yes, false)?;
                print_json(&report)?;
                if !yes {
                    eprintln!("Preview only: --yes confirms removal of all unpinned large data, including KEEP artifacts. Metadata and final results remain.");
                }
                Ok(())
            }
        },
        Action::Search { .. } => unreachable!(),
    }
}
