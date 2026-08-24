//! Feature-gated runtime profiling support.
//!
//! This module is compiled only with `--features profiling`. It deliberately
//! records measurements without changing the numerical path of the scorer.

use pprof::protos::Message;
use serde::Serialize;
use serde_json::Value;
use std::alloc::{GlobalAlloc, Layout, System};
use std::cell::RefCell;
use std::collections::{BTreeMap, BTreeSet};
use std::fs::File;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};
use std::sync::{Mutex, OnceLock};
use std::time::{Duration, Instant};

static ENABLED: AtomicBool = AtomicBool::new(false);
static COUNT_ALLOCATIONS: AtomicBool = AtomicBool::new(false);
static ALLOC_CALLS: AtomicU64 = AtomicU64::new(0);
static ALLOCATED_BYTES: AtomicU64 = AtomicU64::new(0);
static DEALLOC_CALLS: AtomicU64 = AtomicU64::new(0);
static DEALLOCATED_BYTES: AtomicU64 = AtomicU64::new(0);

/// Counts allocator traffic in profiling builds. The production allocator is
/// untouched because this type is only installed under the profiling feature.
pub struct CountingAllocator;

unsafe impl GlobalAlloc for CountingAllocator {
    unsafe fn alloc(&self, layout: Layout) -> *mut u8 {
        if COUNT_ALLOCATIONS.load(Ordering::Relaxed) {
            ALLOC_CALLS.fetch_add(1, Ordering::Relaxed);
            ALLOCATED_BYTES.fetch_add(layout.size() as u64, Ordering::Relaxed);
        }
        unsafe { System.alloc(layout) }
    }

    unsafe fn alloc_zeroed(&self, layout: Layout) -> *mut u8 {
        if COUNT_ALLOCATIONS.load(Ordering::Relaxed) {
            ALLOC_CALLS.fetch_add(1, Ordering::Relaxed);
            ALLOCATED_BYTES.fetch_add(layout.size() as u64, Ordering::Relaxed);
        }
        unsafe { System.alloc_zeroed(layout) }
    }

    unsafe fn dealloc(&self, ptr: *mut u8, layout: Layout) {
        if COUNT_ALLOCATIONS.load(Ordering::Relaxed) {
            DEALLOC_CALLS.fetch_add(1, Ordering::Relaxed);
            DEALLOCATED_BYTES.fetch_add(layout.size() as u64, Ordering::Relaxed);
        }
        unsafe { System.dealloc(ptr, layout) };
    }

    unsafe fn realloc(&self, ptr: *mut u8, layout: Layout, new_size: usize) -> *mut u8 {
        if COUNT_ALLOCATIONS.load(Ordering::Relaxed) {
            ALLOC_CALLS.fetch_add(1, Ordering::Relaxed);
            ALLOCATED_BYTES.fetch_add(new_size as u64, Ordering::Relaxed);
            DEALLOC_CALLS.fetch_add(1, Ordering::Relaxed);
            DEALLOCATED_BYTES.fetch_add(layout.size() as u64, Ordering::Relaxed);
        }
        unsafe { System.realloc(ptr, layout, new_size) }
    }
}

#[derive(Clone, Default)]
struct Context {
    phase: Option<&'static str>,
    fold: Option<u8>,
    iteration: Option<usize>,
    newton_iteration: Option<usize>,
}

thread_local! {
    static CONTEXT: RefCell<Context> = RefCell::new(Context::default());
}

#[derive(Clone, Serialize)]
pub struct Event {
    category: &'static str,
    name: &'static str,
    duration_ns: u64,
    phase: Option<&'static str>,
    fold: Option<u8>,
    iteration: Option<usize>,
    newton_iteration: Option<usize>,
    thread: String,
    elements: Option<u64>,
    bytes: Option<u64>,
}

#[derive(Default)]
struct State {
    start: Option<Instant>,
    events: Vec<Event>,
    allocation_sites: BTreeMap<&'static str, AllocationSite>,
    metadata: BTreeMap<String, Value>,
}

#[derive(Clone, Default, Serialize)]
struct AllocationSite {
    name: &'static str,
    calls: u64,
    bytes: u64,
}

static STATE: OnceLock<Mutex<State>> = OnceLock::new();

fn state() -> &'static Mutex<State> {
    STATE.get_or_init(|| Mutex::new(State::default()))
}

fn allocation_snapshot() -> AllocationTotals {
    AllocationTotals {
        allocation_calls: ALLOC_CALLS.load(Ordering::Relaxed),
        allocated_bytes: ALLOCATED_BYTES.load(Ordering::Relaxed),
        deallocation_calls: DEALLOC_CALLS.load(Ordering::Relaxed),
        deallocated_bytes: DEALLOCATED_BYTES.load(Ordering::Relaxed),
    }
}

#[derive(Clone, Copy, Default, Serialize)]
struct AllocationTotals {
    allocation_calls: u64,
    allocated_bytes: u64,
    deallocation_calls: u64,
    deallocated_bytes: u64,
}

impl AllocationTotals {
    fn difference(self, earlier: Self) -> Self {
        Self {
            allocation_calls: self
                .allocation_calls
                .saturating_sub(earlier.allocation_calls),
            allocated_bytes: self.allocated_bytes.saturating_sub(earlier.allocated_bytes),
            deallocation_calls: self
                .deallocation_calls
                .saturating_sub(earlier.deallocation_calls),
            deallocated_bytes: self
                .deallocated_bytes
                .saturating_sub(earlier.deallocated_bytes),
        }
    }
}

pub struct ContextGuard {
    previous: Context,
}

impl Drop for ContextGuard {
    fn drop(&mut self) {
        CONTEXT.with(|context| *context.borrow_mut() = self.previous.clone());
    }
}

pub fn context(
    phase: Option<&'static str>,
    fold: Option<u8>,
    iteration: Option<usize>,
    newton_iteration: Option<usize>,
) -> ContextGuard {
    let previous = CONTEXT.with(|current| {
        let previous = current.borrow().clone();
        let mut next = previous.clone();
        if phase.is_some() {
            next.phase = phase;
        }
        if fold.is_some() {
            next.fold = fold;
        }
        if iteration.is_some() {
            next.iteration = iteration;
        }
        if newton_iteration.is_some() {
            next.newton_iteration = newton_iteration;
        }
        *current.borrow_mut() = next;
        previous
    });
    ContextGuard { previous }
}

pub struct Scope {
    category: &'static str,
    name: &'static str,
    start: Option<Instant>,
    context: Context,
    elements: Option<u64>,
    bytes: Option<u64>,
}

impl Scope {
    pub fn new(category: &'static str, name: &'static str) -> Self {
        Self::with_size(category, name, None, None)
    }

    pub fn with_elements(category: &'static str, name: &'static str, elements: usize) -> Self {
        Self::with_size(category, name, Some(elements as u64), None)
    }

    pub fn with_size(
        category: &'static str,
        name: &'static str,
        elements: Option<u64>,
        bytes: Option<u64>,
    ) -> Self {
        let enabled = ENABLED.load(Ordering::Relaxed);
        Self {
            category,
            name,
            start: enabled.then(Instant::now),
            context: CONTEXT.with(|context| context.borrow().clone()),
            elements,
            bytes,
        }
    }
}

impl Drop for Scope {
    fn drop(&mut self) {
        if let Some(start) = self.start {
            record_with_context(
                self.category,
                self.name,
                start.elapsed(),
                self.context.clone(),
                self.elements,
                self.bytes,
            );
        }
    }
}

pub fn record(
    category: &'static str,
    name: &'static str,
    duration: Duration,
    elements: Option<u64>,
    bytes: Option<u64>,
) {
    let context = CONTEXT.with(|current| current.borrow().clone());
    record_with_context(category, name, duration, context, elements, bytes);
}

fn record_with_context(
    category: &'static str,
    name: &'static str,
    duration: Duration,
    context: Context,
    elements: Option<u64>,
    bytes: Option<u64>,
) {
    if !ENABLED.load(Ordering::Relaxed) {
        return;
    }
    let event = Event {
        category,
        name,
        duration_ns: duration.as_nanos().min(u64::MAX as u128) as u64,
        phase: context.phase,
        fold: context.fold,
        iteration: context.iteration,
        newton_iteration: context.newton_iteration,
        thread: format!("{:?}", std::thread::current().id()),
        elements,
        bytes,
    };
    state()
        .lock()
        .expect("profiling state poisoned")
        .events
        .push(event);
}

pub fn allocation_site(name: &'static str, calls: u64, bytes: u64) {
    if !ENABLED.load(Ordering::Relaxed) {
        return;
    }
    let mut state = state().lock().expect("profiling state poisoned");
    let site = state
        .allocation_sites
        .entry(name)
        .or_insert(AllocationSite {
            name,
            calls: 0,
            bytes: 0,
        });
    site.calls += calls;
    site.bytes += bytes;
}

pub fn metadata<T: Serialize>(key: &str, value: T) {
    if !ENABLED.load(Ordering::Relaxed) {
        return;
    }
    let value = serde_json::to_value(value).expect("profiling metadata must serialize");
    state()
        .lock()
        .expect("profiling state poisoned")
        .metadata
        .insert(key.to_string(), value);
}

pub struct Session {
    json_path: Option<PathBuf>,
    cpu_prefix: Option<PathBuf>,
    cpu_guard: Option<pprof::ProfilerGuard<'static>>,
    allocations_at_start: AllocationTotals,
    count_allocations: bool,
}

impl Session {
    pub fn start(
        json_path: Option<String>,
        cpu_prefix: Option<String>,
        count_allocations: bool,
    ) -> Result<Self, String> {
        if json_path.is_none() && cpu_prefix.is_none() {
            return Ok(Self {
                json_path: None,
                cpu_prefix: None,
                cpu_guard: None,
                allocations_at_start: allocation_snapshot(),
                count_allocations: false,
            });
        }
        {
            let mut current = state().lock().map_err(|_| "profiling state poisoned")?;
            current.start = None;
            current.events.clear();
            current.allocation_sites.clear();
            current.metadata.clear();
        }
        let cpu_guard = if cpu_prefix.is_some() {
            Some(
                pprof::ProfilerGuardBuilder::default()
                    .frequency(499)
                    .build()
                    .map_err(|error| format!("CPU profiler failed to start: {error}"))?,
            )
        } else {
            None
        };
        {
            let mut current = state().lock().map_err(|_| "profiling state poisoned")?;
            current.start = Some(Instant::now());
        }
        COUNT_ALLOCATIONS.store(count_allocations, Ordering::SeqCst);
        let allocations_at_start = allocation_snapshot();
        ENABLED.store(true, Ordering::SeqCst);
        Ok(Self {
            json_path: json_path.map(PathBuf::from),
            cpu_prefix: cpu_prefix.map(PathBuf::from),
            cpu_guard,
            allocations_at_start,
            count_allocations,
        })
    }

    pub fn finish(mut self) -> Result<(), String> {
        if self.json_path.is_none() && self.cpu_prefix.is_none() {
            return Ok(());
        }
        let elapsed_ns = state()
            .lock()
            .map_err(|_| "profiling state poisoned")?
            .start
            .map(|start| start.elapsed().as_nanos().min(u64::MAX as u128) as u64)
            .unwrap_or(0);
        let allocator = allocation_snapshot().difference(self.allocations_at_start);
        if self.count_allocations {
            COUNT_ALLOCATIONS.store(false, Ordering::SeqCst);
        }

        let cpu = if let Some(guard) = self.cpu_guard.take() {
            let report = guard
                .report()
                .build()
                .map_err(|error| format!("CPU report failed: {error}"))?;
            let summary = write_cpu_artifacts(
                &report,
                self.cpu_prefix
                    .as_deref()
                    .expect("CPU prefix accompanies CPU guard"),
            )?;
            drop(guard);
            Some(summary)
        } else {
            None
        };
        ENABLED.store(false, Ordering::SeqCst);

        if let Some(path) = &self.json_path {
            let current = state().lock().map_err(|_| "profiling state poisoned")?;
            let events = current.events.clone();
            let allocation_sites: Vec<AllocationSite> =
                current.allocation_sites.values().cloned().collect();
            let metadata = current.metadata.clone();
            drop(current);
            let aggregates = aggregate_events(&events);
            let output = ProfileOutput {
                schema_version: 1,
                elapsed_ns,
                metadata,
                allocator,
                allocation_sites,
                aggregates,
                events,
                cpu,
            };
            let file = File::create(path)
                .map_err(|error| format!("cannot create {}: {error}", path.display()))?;
            let mut writer = std::io::BufWriter::with_capacity(1 << 20, file);
            serde_json::to_writer(&mut writer, &output)
                .map_err(|error| format!("cannot serialize {}: {error}", path.display()))?;
            writer
                .flush()
                .map_err(|error| format!("cannot flush {}: {error}", path.display()))?;
        }
        Ok(())
    }
}

#[derive(Serialize)]
struct ProfileOutput {
    schema_version: u32,
    elapsed_ns: u64,
    metadata: BTreeMap<String, Value>,
    allocator: AllocationTotals,
    allocation_sites: Vec<AllocationSite>,
    aggregates: Vec<Aggregate>,
    events: Vec<Event>,
    cpu: Option<CpuSummary>,
}

#[derive(Default, Serialize)]
struct Aggregate {
    category: &'static str,
    name: &'static str,
    calls: u64,
    total_ns: u64,
    min_ns: u64,
    max_ns: u64,
    total_elements: u64,
    total_bytes: u64,
}

fn aggregate_events(events: &[Event]) -> Vec<Aggregate> {
    let mut output: BTreeMap<(&'static str, &'static str), Aggregate> = BTreeMap::new();
    for event in events {
        let aggregate = output
            .entry((event.category, event.name))
            .or_insert(Aggregate {
                category: event.category,
                name: event.name,
                min_ns: u64::MAX,
                ..Aggregate::default()
            });
        aggregate.calls += 1;
        aggregate.total_ns += event.duration_ns;
        aggregate.min_ns = aggregate.min_ns.min(event.duration_ns);
        aggregate.max_ns = aggregate.max_ns.max(event.duration_ns);
        aggregate.total_elements += event.elements.unwrap_or(0);
        aggregate.total_bytes += event.bytes.unwrap_or(0);
    }
    output.into_values().collect()
}

#[derive(Serialize)]
struct CpuSummary {
    frequency_hz: i32,
    duration_ns: u64,
    total_samples: u64,
    leaf_functions: Vec<SymbolSamples>,
    inclusive_functions: Vec<SymbolSamples>,
    protobuf_path: String,
    collapsed_path: String,
    flamegraph_path: String,
}

#[derive(Serialize)]
struct SymbolSamples {
    symbol: String,
    samples: u64,
    percent: f64,
}

fn write_cpu_artifacts(report: &pprof::Report, prefix: &Path) -> Result<CpuSummary, String> {
    let protobuf_path = prefix.with_extension("pb");
    let collapsed_path = prefix.with_extension("collapsed");
    let flamegraph_path = prefix.with_extension("svg");

    let protobuf = report
        .pprof()
        .map_err(|error| format!("cannot build protobuf CPU profile: {error}"))?;
    let mut protobuf_bytes = Vec::new();
    protobuf
        .write_to_vec(&mut protobuf_bytes)
        .map_err(|error| format!("cannot encode protobuf CPU profile: {error}"))?;
    File::create(&protobuf_path)
        .and_then(|mut file| file.write_all(&protobuf_bytes))
        .map_err(|error| format!("cannot write {}: {error}", protobuf_path.display()))?;

    let flamegraph = File::create(&flamegraph_path)
        .map_err(|error| format!("cannot create {}: {error}", flamegraph_path.display()))?;
    report
        .flamegraph(flamegraph)
        .map_err(|error| format!("cannot write flamegraph: {error}"))?;

    let mut collapsed = File::create(&collapsed_path)
        .map_err(|error| format!("cannot create {}: {error}", collapsed_path.display()))?;
    let mut leaf = BTreeMap::<String, u64>::new();
    let mut inclusive = BTreeMap::<String, u64>::new();
    let mut total_samples = 0u64;
    for (frames, &count) in &report.data {
        if count <= 0 {
            continue;
        }
        let count = count as u64;
        total_samples += count;
        let mut stack = Vec::new();
        for frame in frames.frames.iter().rev() {
            for symbol in frame.iter().rev() {
                stack.push(symbol.name().replace(';', ":"));
            }
        }
        if let Some(symbol) = frames
            .frames
            .first()
            .and_then(|frame| frame.first())
            .map(|symbol| symbol.name())
        {
            *leaf.entry(symbol).or_default() += count;
        }
        let mut seen = BTreeSet::new();
        for symbol in &stack {
            if seen.insert(symbol.clone()) {
                *inclusive.entry(symbol.clone()).or_default() += count;
            }
        }
        writeln!(collapsed, "{} {}", stack.join(";"), count)
            .map_err(|error| format!("cannot write {}: {error}", collapsed_path.display()))?;
    }

    Ok(CpuSummary {
        frequency_hz: report.timing.frequency,
        duration_ns: report.timing.duration.as_nanos().min(u64::MAX as u128) as u64,
        total_samples,
        leaf_functions: top_symbols(leaf, total_samples),
        inclusive_functions: top_symbols(inclusive, total_samples),
        protobuf_path: protobuf_path.display().to_string(),
        collapsed_path: collapsed_path.display().to_string(),
        flamegraph_path: flamegraph_path.display().to_string(),
    })
}

fn top_symbols(samples: BTreeMap<String, u64>, total: u64) -> Vec<SymbolSamples> {
    let mut output: Vec<_> = samples
        .into_iter()
        .map(|(symbol, samples)| SymbolSamples {
            symbol,
            samples,
            percent: if total == 0 {
                0.0
            } else {
                100.0 * samples as f64 / total as f64
            },
        })
        .collect();
    output.sort_by(|left, right| {
        right
            .samples
            .cmp(&left.samples)
            .then_with(|| left.symbol.cmp(&right.symbol))
    });
    output.truncate(50);
    output
}
