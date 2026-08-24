# Runtime profile: PXD032157

This is a measurement-only profile of the canonical `percolator-rs` workload. No scoring,
solver, FDR, iteration, fold, or output behavior was optimized. The measurements were made on
2026-08-24 from `d15438da3041a3b379b7a6b50b44b50d7e3a7c7d` plus the feature-gated
instrumentation patch captured with the artifacts.

The host was an AMD Ryzen 5 5600G (6 cores/12 threads) running x86-64 Linux, with Rust 1.97.0
and LLVM 22.1.6. The full dataset contains 65 PIN files, 2,295,401,156 input bytes, and
8,639,746 PSMs. All timing percentages below use the summed in-process elapsed time from the
instrumented full sequential run (57.489 s), unless explicitly labeled as sampled CPU time.
Nested-operation percentages overlap; the top-level stage table does not.

Raw artifacts are under
`$HOME/percolator_rs_out/runtime-profile-20260824-current`. They include per-process JSON,
protobuf CPU profiles, collapsed stacks, an SVG flamegraph for the large-file run and for every
N=4 process, exact output files, build metadata, and the attempted `perf` probe. The compact
machine-readable aggregate is `bench/runtime-profile-results.json`.

## Answer in one paragraph

Rescoring is 64.80% of process time. Inside it, SVM training is 31.49% of total runtime and
active-set/margin evaluation alone is 24.63%; repeated q-value evaluation is 22.44%, of which
sorting is 15.44%. Outside rescoring, PIN parsing is 12.27%, result formatting/buffering is
10.93%, and peptide-level work is 7.35%. Cholesky and the direct solve are not bottlenecks
(0.06% and 0.02% of SVM time). Actual file writes are only 0.88%; formatting the TSV text is
the expensive output operation. These results are independently consistent with the sampled
CPU profile.

## Configurations, scaling, and overhead

| Configuration | Normal wall | Timed wall | Result |
|---|---:|---:|---:|
| largest PIN, `--num-threads 1` | 2.061 s median (3) | 1.997 s median (3) | paired median −3.13% |
| largest PIN, `--num-threads 3` | 1.408 s median (3) | 1.376 s median (3) | paired median −1.82% |
| full PXD032157 sequential | 61.753 s | 59.379 s | paired delta −3.84% |
| full PXD032157, file-level N=4 | 20.330 s | 19.415 s | paired delta −4.50% |
| balanced-order largest-PIN overhead check | 2.173 s median (6) | 2.346 s median (6) | paired median **+0.79%** |
| largest PIN, sampled CPU build | — | 2.124 s | +6.36% versus timed median |
| full N=4, sampled CPU build | — | 21.889 s | +12.75% versus timed run |

Negative deltas are noise/code-layout effects, not speedups. The first matrix always ran the
normal binary first, so a separate six-pair check alternated binary order. Its paired deltas
were −3.65%, −2.64%, +13.43%, +13.55%, −1.39%, and +2.96%, with a median of +0.79%.
Allocation counting was disabled in all timing runs and enabled only in its own large-file run.
CPU sampling used a separate frame-pointer/debug-info build and is not used for wall-time stage
totals.

Three fold threads improve the largest file by 1.46x rather than 3x because parsing, peptide
processing, and output remain serial; measured fold dispatch/join overhead was only 0.62 ms per
run. File-level N=4 improves the full normal run by 3.04x (75.9% parallel efficiency). The N=4
instrumented processes accumulated 73.445 s of process wall time in 19.415 s, or about 94.6%
occupancy of four process slots; contention inflated the summed work relative to sequential,
especially rescoring and output.

## Complete runtime breakdown

| Top-level stage | Time | % total |
|---|---:|---:|
| PIN loading/parsing | 7.128 s | 12.40% |
| Rescoring | 37.253 s | 64.80% |
| PSM-level materialization | 1.259 s | 2.19% |
| Peptide-level processing | 4.224 s | 7.35% |
| Result output | 7.613 s | 13.24% |
| Miscellaneous/unaccounted | 0.011 s | 0.02% |

| Nested operation | Time | % total |
|---|---:|---:|
| SVM training | 18.105 s | 31.49% |
| → active set and margin scoring | 14.158 s | 24.63% |
| → Hessian construction | 3.014 s | 5.24% |
| → gradient | 0.849 s | 1.48% |
| → Cholesky | 0.011 s | 0.02% total / 0.06% SVM |
| → linear solve | 0.004 s | 0.01% total / 0.02% SVM |
| q-value calculation, inclusive | 12.901 s | 22.44% |
| → score-index sorting | 8.879 s | 15.44% |
| → scan and monotonicity | 3.786 s | 6.59% |
| Initial-direction selection, inclusive | 7.607 s | 13.23% |
| Model scoring outside the SVM objective | 2.552 s | 4.44% |
| Confident-positive selection | 0.695 s | 1.21% |
| Normalization | 0.839 s | 1.46% |
| PEP/PAVA, inclusive | 0.964 s | 1.68% |
| Result formatting/buffering | 6.285 s | 10.93% |
| File writes | 0.508 s | 0.88% |

The conditional picked-protein path was measured separately on the largest PIN. It took
350.1 ms (15.25% of that run), including 311.3 ms (13.56%) in
`protein::infer`; protein-specific sorts together took 9.48 ms (0.41%). Protein inference was
not enabled in the standard 65-file benchmark and is therefore absent from its stage total.

## Semi-supervised iterations

Times are means per fold invocation over 65 files (195 fold invocations per iteration). All ten
iterations and all three folds were retained.

| Iteration | Newton steps | Training | Scoring | q-values | Positive selection | Total |
|---:|---:|---:|---:|---:|---:|---:|
| 0 | 2,962 | 39.394 ms | 1.225 ms | 3.006 ms | 0.350 ms | 44.119 ms |
| 1 | 1,427 | 10.942 ms | 1.167 ms | 3.097 ms | 0.354 ms | 15.695 ms |
| 2 | 1,057 | 7.859 ms | 1.168 ms | 3.092 ms | 0.356 ms | 12.609 ms |
| 3 | 885 | 6.617 ms | 1.184 ms | 3.096 ms | 0.359 ms | 11.389 ms |
| 4 | 775 | 5.801 ms | 1.183 ms | 3.074 ms | 0.358 ms | 10.551 ms |
| 5 | 707 | 5.302 ms | 1.176 ms | 3.082 ms | 0.357 ms | 10.057 ms |
| 6 | 651 | 4.844 ms | 1.185 ms | 3.085 ms | 0.358 ms | 9.607 ms |
| 7 | 586 | 4.444 ms | 1.179 ms | 3.081 ms | 0.358 ms | 9.198 ms |
| 8 | 501 | 3.772 ms | 1.186 ms | 3.073 ms | 0.359 ms | 8.528 ms |
| 9 | 445 | 3.870 ms | 1.197 ms | 3.067 ms | 0.357 ms | 8.627 ms |

Iteration 0 accounts for 7.682 s, or 42.4% of all SVM time, because its 2,962 Newton steps are
about twice the count of iteration 1. Scoring, q-values, and positive selection remain nearly
flat across iterations; only solver convergence reduces later-iteration time.

## SVM solver

The 1,950 SVM training calls executed 9,996 Newton iterations. Mean measured Newton-iteration
time was 1.522 ms. `line_search_total` is an inclusive timer containing active-set/margin
rescoring, so its percentage must not be added to the other components.

| Component | Time | Calls | Mean/call | % SVM |
|---|---:|---:|---:|---:|
| Active set and margin scoring | 14.158 s | 10,194 | 1.389 ms | 78.20% |
| Line search, inclusive | 11.315 s | 8,046 | 1.406 ms | 62.50% |
| Hessian construction | 3.014 s | 8,046 | 0.375 ms | 16.65% |
| Gradient | 0.849 s | 9,996 | 0.085 ms | 4.69% |
| Buffer allocation/initialization | 0.037 s | 1,950 | 0.019 ms | 0.21% |
| Cholesky | 0.011 s | 8,046 | 0.0014 ms | 0.06% |
| Linear solve | 0.004 s | 8,046 | 0.0005 ms | 0.02% |
| Convergence logic | 0.0007 s | 9,996 | 0.00007 ms | <0.01% |
| Solver/weight updates | 0.0009 s | 16,290 | <0.0001 ms | <0.01% |

The evidence rules out Cholesky, the direct solve, convergence checks, and buffer initialization
as useful first targets. The solver cost is the repeated full-row objective/active-set pass,
followed by Hessian construction.

## Cross-validation folds

Values are means per file. Training scoring excludes final held-out scoring.

| Fold | Setup | Training | Training scoring | q-values | Held-out scoring | Total |
|---:|---:|---:|---:|---:|---:|---:|
| 0 | 0.650 ms | 137.869 ms | 11.873 ms | 29.564 ms | 1.235 ms | 139.768 ms |
| 1 | 0.688 ms | 139.340 ms | 11.791 ms | 30.229 ms | 1.241 ms | 141.282 ms |
| 2 | 0.688 ms | 144.221 ms | 11.885 ms | 32.465 ms | 1.233 ms | 146.155 ms |

Fold 2 is 4.6% slower than fold 0 on aggregate and was the slowest fold in 35/65 files; fold 1
was slowest in 19 and fold 0 in 11. The median within-file maximum/minimum fold ratio was 1.086
(mean 1.097). This is a small, data-dependent imbalance rather than a scheduling bottleneck.

## Parsing

The parser processed 2.138 GiB at 310.5 MiB/s. `mmap` setup and header handling are negligible;
page access happens during the row loop and is included there.

| Parser operation | Time | % total | Volume |
|---|---:|---:|---:|
| Total PIN parse | 7.051 s | 12.27% | 2,295,401,156 bytes |
| Row loading | 7.047 s | 12.26% | 8,639,746 rows |
| Numeric/float parsing | 1.999 s | 3.48% | 181,434,666 floats |
| String allocation/copy | 1.792 s | 3.12% | 925,703,237 copied bytes |
| Field splitting | 1.573 s | 2.74% | 8,639,746 rows |
| `mmap` setup | 0.0009 s | <0.01% | 65 mappings |
| Header/features | 0.0006 s | <0.01% | 1,820 columns seen |

The remaining 1.683 s in the row loop covers vector pushes/growth, validation/branches, protein
span discovery, page faults, and timer overhead.

## Sorting and q-values

| Sort | Calls | Elements processed | Time | % total |
|---|---:|---:|---:|---:|
| q-value score order | 4,810 | 551,797,661 | 8.879 s | 15.44% |
| result row score order | 260 | 16,133,409 | 0.808 s | 1.40% |
| PEP score order | 130 | 16,133,409 | 0.434 s | 0.75% |
| peptide input order | 65 | 7,493,663 | 0.078 s | 0.14% |

The q-value sorts split into 2,730 initial-direction sorts (3.893 s), 1,950 semi-supervised
iteration sorts (4.545 s), 65 final-PSM sorts (0.237 s), and 65 peptide sorts (0.203 s).
Thus repeated training-iteration sorting alone is 7.91% of total runtime and is significant.
Initial-direction evaluation makes 42 q-value calls per file; semi-supervised training makes 30
per file.

## Allocations

The dedicated largest-PIN allocation run observed 2,799,341 allocation calls and 1,049,560,148
allocated bytes. It completed with byte-identical output and is excluded from timing conclusions.
Site accounting is approximate capacity/owned-string traffic, so sites can overlap allocator
totals and should be used to rank churn, not as a heap-live-size measurement.

| Largest-PIN allocation site | Approx. calls | Approx. bytes | Repetition |
|---|---:|---:|---|
| `stats::qvalues` temporary vectors | 222 | 408,367,728 | 3 vectors × 74 calls |
| PIN column vectors | 7 | 92,971,742 | per file |
| Per-iteration training row/y/C vectors | 90 | 64,334,136 | 3 vectors × 30 fold-iterations |
| Normalized design matrix | 1 | 46,918,960 | per file |
| SVM work buffers | 210 | 43,026,704 | 7 buffers × 30 fits |
| PEP temporary vectors | 20 | 38,964,160 | 10 vectors × 2 calls |
| PSM output row vectors | 2 | 37,748,736 | per file |
| Positive/negative row vectors | 60 | 33,161,216 | 2 vectors × 30 fold-iterations |
| PSM parser strings | 799,755 | 29,644,764 | 3 strings/PSM where nonempty |
| PSM output strings | 799,755 | 29,644,764 | cloned for output |
| Peptide output strings | 661,401 | 24,300,897 | cloned for output |
| Peptide dedup keys | 266,585 | 5,919,801 | one formatted key/PSM |

Across the full run, q-value temporary capacity was approximately 13.24 GB, iteration training
buffers 2.06 GB, SVM work buffers 1.38 GB, and PSM/peptide output strings 1.73 GB. Sampled full
N=4 CPU time also showed allocator internals (`_int_free_chunk` 4.11% inclusive,
`__libc_malloc2` 2.29%, `_int_malloc` 2.15%), confirming that allocation is visible, although
those inclusive percentages overlap and do not establish a single allocator-only wall total.

## CPU hotspots

Linux `perf` was attempted but rejected because `perf_event_paranoid=4` and the process lacked
the required capability. No system setting was changed. The feature instead used `pprof-rs`
SIGPROF sampling at 499 Hz. The large-file profile collected 985 samples; the full N=4 profile
combined 37,062 samples. Percentages are inclusive unless labeled leaf/self.

| Function | Large PIN CPU | Full N=4 CPU | Full N=4 leaf |
|---|---:|---:|---:|
| `percolator::train_fold` | 49.75% | 50.15% | 8.24% |
| `svm::Problem::f_and_active` | 26.40% | 27.91% | 27.89% |
| `stats::qvalues` | 22.23% | 20.24% | 6.96% |
| q-value unstable quicksort | 14.72% | 13.22% | 11.88% |
| `write_results` | 13.20% | 13.51% | — |
| `initial_direction` | 13.20% | 13.13% | 2.88% |
| `FoldModel::score_rows` | 4.16% | 5.02% | 5.01% |
| `memcpy` | 4.06% | 4.45% | 4.45% |

The sampled and event-timed rankings agree. Timer calls themselves accounted for approximately
2% of sampled full-run CPU (`clock_gettime` 1.91% leaf), another reason to use the normal/timed
wall comparison and event totals for magnitude, and samples for independent rank confirmation.

## Determinism and scientific validation

All standard runs completed 65/65 files and produced exactly 107,046 target PSMs and 37,469
target peptides at strict q<0.01. Normal, timed, allocation-counted, and CPU-profiled outputs
were byte-identical: 260/260 full sequential files, 260/260 full N=4 files, and every single-file
comparison. The run retained seed 1, fixed Cpos=1/Cneg=4, ten semi-supervised iterations, three
cross-validation folds, the existing solver convergence path, and all current q-value/PEP logic.

## Ranked optimization opportunities (not implemented)

Percentages are measured wall-time impact. Entries overlap where stated, so they are rankings of
code paths to investigate rather than additive speedup promises.

| Rank | Measured impact and exact path | Scaling / why expensive | Plausible next experiment | Difficulty | Correctness/determinism risk |
|---:|---|---|---|---|---|
| **1** | **24.63% total; 78.20% SVM** — `svm::Problem::f_and_active`, `src/svm.rs:42` | Full margin dot product over training rows on every initial objective and line-search trial; scales with PSMs × features × Newton/line-search evaluations × folds | Measure data-layout/vectorized dot variants or avoid redundant objective passes while preserving the accepted Newton step exactly | Medium | Medium: floating-point order can change margins, active membership, and convergence |
| **2** | **15.44% total** — q-value index sort in `stats::qvalues`, `src/stats.rs:8` | 4,810 sorts over 551.8M elements; scales O(calls × PSMs log PSMs), with 51.2% of sort time inside training iterations | Benchmark deterministic order reuse, specialized numeric ordering, or a stable tie-preserving rank representation | Medium/High | High: ties/order affect q-values and confident positives |
| **3** | **13.23% inclusive; 2.04% feature scan exclusive** — `initial_direction`, `src/percolator.rs:298` | Evaluates both directions of 21 features and performs 42 full q-value calls/file; overlaps rank 2 by 6.78% total | Test pre-ranked feature columns or an equivalent count-at-FDR calculation against byte-identical outputs | High | High: the chosen initial feature controls every fold trajectory |
| **4** | **12.27% total** — `pin::parse`, `src/pin.rs:203` | 2.14 GiB, 8.64M rows, 181.4M floats, and 25.9M strings; scales with bytes, rows, features | Separately prototype lower-overhead field scanning/float conversion and borrowed/deferred strings | Medium | Low/Medium: malformed rows and exact accepted syntax must remain identical |
| **5** | **10.93% total** — `write_results`, `src/main.rs:365` | Per-row generic formatting, especially four fixed-precision floats/text fields; CPU profile confirms 13.51% inclusive | Benchmark a specialized row formatter or larger batched serialization while retaining exact bytes | Medium | Low if golden-byte tests cover every output |
| **6** | **6.59% total** — q-value scan/monotonization in `stats::qvalues`, `src/stats.rs:8` | Rebuilds/reorders temporary vectors and scans every q-value call; scales with PSMs × calls | Profile buffer reuse and fused/in-place scans independently of sorting | Low/Medium | Medium: off-by-one/tie/FDR errors are scientifically consequential |
| **7** | **about 6.0% residual peptide stage** — dedup/materialization at `src/main.rs:773` | Formats a key for every PSM, hashes it, then clones strings/rows; scales with PSMs and unique peptides | Measure indexed/borrowed keys and delayed row materialization | Medium | Medium: peptide canonicalization and deterministic tie choice must not change |
| **8** | **5.24% total; 16.65% SVM** — `Problem::hessian`, `src/svm.rs:99` | Explicit active-row outer products over the upper triangle; scales with active PSMs × features² × Newton steps | Benchmark symmetric packed accumulation or evidence-backed vectorization | Medium | Medium: summation order can change Newton directions/convergence |
| **9** | **4.44% total** — `FoldModel::score_rows`, `src/percolator.rs:245` | Repeated full held/training-row dot products; scales with PSMs × features × iterations/folds | Benchmark matrix traversal/data layout using the same arithmetic order first | Medium | Medium for floating-point/tie effects |
| **10** | **1.68% total** — `stats::peps`, `src/stats.rs:131` | Two full-data sorts plus PAVA buffers; scales with PSMs/peptides × log count | Investigate shared ordering/buffer reuse only after higher-ranked work | Medium | High: PEP calibration and PAVA monotonicity must remain exact |

### First three investigations

1. SVM active-set/margin evaluation: it is the largest exclusive measured hotspot in both wall
   timers and CPU samples.
2. The q-value engine, beginning with its repeated score sort and then its scan: together it is
   22.44% of runtime and dominates both initial direction and every semi-supervised iteration.
3. PIN row parsing: it is the largest non-overlapping stage after those algorithmic costs and has
   clear sub-costs in floats, strings, and field splitting. Output formatting is a close fourth.

Cholesky, direct solving, file writes, fold scheduling, `mmap` setup, header parsing, and
convergence checks should not be prioritized: their measured costs are too small.
