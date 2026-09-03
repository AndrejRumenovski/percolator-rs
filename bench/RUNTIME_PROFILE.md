# Fresh runtime profile: refactored `percolator-rs`

## Decision

**Recommendation: A — do another optimization round, restricted first to exact-order q-value and initial-direction work.**

The current implementation still has a large, coherent opportunity. Q-value work is 40.51% of
sequential process time, including 29.51% in its two index sorts and 10.14% in tie-group/counting
scans. Initial-direction selection contains most of the extra calls and is 29.24% inclusive. This is
not the old hotspot: active-set/margin evaluation fell from 24.63% to 6.94%, output formatting from
10.93% to 0.85%, and peptide processing from 7.35% to 0.67%.

An N=4 runtime below 12 s is aggressive but realistic only as a combined result: it needs a 22.49%
total reduction, equivalent to an optimistic 61.76% reduction of all measured N=4 q-value cost.
Eliminating score-order sorting alone would still predict 12.15 s. A path based on reusing exact
feature orders in initial-direction selection, removing the reverse-rank re-sort, accelerating the
remaining exact q-value sort/scan, and reducing N=4 memory contention is credible enough to measure
in a separate optimization goal. Below 10 s is not currently realistic: it needs 35.41% total, or
97.23% of all N=4 q-value cost in isolation. Even the impossible case of eliminating every q-value
operation predicts only 9.84 s.

This campaign changed no production algorithm, statistical method, output representation, precision,
or concurrency policy. All source changes are behind the existing `profiling` feature; the remaining
changes are acquisition/report tooling.

## 1. Frozen revision and environment

Measurements were taken on 2026-09-03 from base commit
`159b52424a2eb8a3454303ed4d06af008e14ccb8`, plus the feature-gated instrumentation recorded by this
profiling commit. Before instrumentation, `cargo build --release --locked` produced SHA-256
`9f7bbf67e836b2e9da406868e0f3b5499defe1f591d5b98705e440bb1a40abd7`.

The isolated campaign binaries were:

| Build | SHA-256 |
|---|---|
| normal release | `951cd318941030c492e18a954a0db0238ad2f85c6f032411e93b01cde91f9c45` |
| stage-timed | `a0eb9f074f395dc710ad63a892dc7191fd4bb5ffc4c4729e30d7829a29c85c3a` |
| CPU-sampled | `a1df674d301735d1a6ba638bb240fe08199717d0c32d2ea9efde0848b92c7300` |

| Item | Value |
|---|---|
| CPU | AMD Ryzen 5 5600G, 6 cores / 12 threads, 16 MiB L3 |
| Memory | 30 GiB RAM, 8 GiB swap (unused during environment capture) |
| OS | Ubuntu 26.04 LTS; Linux 7.0.0-30-generic x86-64 |
| Rust | rustc 1.97.0 (`2d8144b78`, 2026-07-07), LLVM 22.1.6 |
| Cargo | 1.97.0 |
| Release flags | `opt-level=3`, thin LTO, one codegen unit, `panic=abort` |
| Target flags | `-C target-cpu=x86-64-v3` |
| Input filesystem | `ntfs3` |
| Temporary output filesystem | `tmpfs` (`/tmp`) |
| Sampling | `pprof-rs` SIGPROF at 499 Hz, frame pointers and debug info enabled |
| `perf` | unavailable: `kernel.perf_event_paranoid=4`; the capability probe failed |

The PXD032157 workload is 65 PIN files, 2,295,401,156 bytes, and 8,639,746 PSMs. The largest
file is `9March2015-29-MAGs-pellet-2ndRep-14N-male-02-comet.pin`: 71,971,056 bytes and 266,585
PSMs. Full-corpus files were submitted in a fixed largest/smallest interleave. Thread counts and
binary order were rotated across repetitions. Output hashing ran after the measured interval.

The compact machine-readable aggregate is [`runtime-profile-results.json`](runtime-profile-results.json).
The reproducible collector is [`fresh_runtime_profile.py`](fresh_runtime_profile.py), and the
aggregator is [`runtime_profile_report.py`](runtime_profile_report.py). The acquisition directories
were `/tmp/percolator-fresh-profile-20260903` for allocation/conditional evidence,
`/tmp/percolator-fresh-timings-20260903` for corrected wall timings, and
`/tmp/percolator-corrected-fold-profile` for the corrected fold labels. The first acquisition exposed
that output hashing was inside the harness timer; that harness defect was corrected and all reported
end-to-end medians come from the second acquisition. Its process-internal profiles, allocations, and
CPU samples were unaffected.

## 2. Current baseline medians

| Workload | Runs | Current median | Refactor-audit reference | Change |
|---|---:|---:|---:|---:|
| largest PIN, one thread | 5 | 1.616816 s | ~1.589 s | +1.75% |
| largest PIN, three-thread mode | 5 | 0.890897 s | ~0.883 s | +0.89% |
| 65 files, sequential | 3 | 49.619487 s | ~49.202 s | +0.85% |
| 65 files, file-level N=4 | 3 | 15.482359 s | ~15.154 s | +2.17% |

The current values are authoritative. The small differences from the audit are consistent with
normal host/code-layout variation and the current environment; the module-boundary measurements
below do not identify a matching new cost.

All full-corpus configurations produced the same aggregate outcomes on every repetition: 106,823
target PSMs and 35,886 target peptides at `q < 0.01` across 65 files.

## 3. Instrumentation overhead

Stage timing and CPU sampling used separate isolated binaries. Negative stage-timer deltas are
code-layout/scheduling noise, not speedups.

| Workload | Normal median | Stage-timed median | Paired median delta | CPU-sampled median | CPU-build delta |
|---|---:|---:|---:|---:|---:|
| largest PIN, one thread | 1.616816 s (5) | 1.616927 s (5) | +0.003% | 1.738147 s (3) | +7.50% |
| largest PIN, three-thread mode | 0.890897 s (5) | 0.871085 s (5) | -2.253% | 1.032045 s (3) | +15.84% |
| full sequential | 49.619487 s (3) | 49.434891 s (3) | -0.366% | 56.508479 s (1) | +13.88% |
| full N=4 | 15.482359 s (3) | 15.384554 s (3) | -0.718% | 17.322051 s (3) | +11.88% |

The one-thread paired result shows essentially zero stage-timer overhead. The five paired one-thread
deltas were -1.255%, +1.290%, -1.239%, +0.017%, and +0.003%. Stage proportions use the stage-timed
build; CPU percentages are independent confirmation and are not used as wall-time totals. The CPU
build's 8-16% overhead comes from frame pointers/debug info, 499 Hz sampling, and profiling events.

## 4. Top-level runtime breakdown

The table aggregates three full sequential stage-timed passes (195 processes). Percentages use
145.407848 s of summed in-process elapsed time. The stage rows are exclusive and sum to 100%; nested
tables later in the report intentionally overlap.

| Exclusive stage | Aggregate time | Mean per corpus | % process time |
|---|---:|---:|---:|
| CLI/setup | 0.003 s | 0.001 s | <0.01% |
| PIN loading, parsing, ordinary one-file join | 24.169 s | 8.056 s | 16.62% |
| Rescoring | 116.759 s | 38.920 s | 80.30% |
| PSM competition/statistics/materialization | 1.394 s | 0.465 s | 0.96% |
| Peptide processing/materialization | 0.981 s | 0.327 s | 0.67% |
| Result output | 2.084 s | 0.695 s | 1.43% |
| Miscellaneous/unaccounted | 0.018 s | 0.006 s | 0.01% |

The 49.619 s end-to-end normal median is larger than the 48.469 s mean in-process stage total because
the latter starts after process/CLI startup and excludes harness scheduling/polling across 65 short
processes. Percentages therefore describe scorer process work, not harness overhead.

Conditional paths were measured on workloads that exercise them:

| Path | Workload | Measured cost |
|---|---|---:|
| input joining/canonicalization | two smallest real PINs, 66,094 PSMs | 21.016 ms; 5.90% of that process |
| RT input augmentation | largest PIN with `--rt-features` | 36.601 ms; 2.07% |
| fold-local RT preprocessing | same RT run, three folds | 8.852 ms; 0.50% |
| picked-protein inference and output | real `data/F_3.pin`, 105,560 PSMs | 61.342 ms; 5.31% |

## 5. Rescoring breakdown

Calls and times aggregate the same three sequential corpus profiles. “% rescoring” is relative to
116.759 s. Q-value and initial-direction rows overlap because q-value passes are nested within
initial direction and the semi-supervised iterations.

| Inclusive operation | Calls | Total | Mean/call | % total | % rescoring |
|---|---:|---:|---:|---:|---:|
| all q-value/count/mask passes | 31,005 | 58.910 s | 1.900 ms | 40.51% | 50.45% |
| initial-direction selection | 585 | 42.524 s | 72.690 ms | 29.24% | 36.42% |
| SVM training | 5,850 | 18.912 s | 3.233 ms | 13.01% | 16.20% |
| model scoring | 7,020 | 8.976 s | 1.279 ms | 6.17% | 7.69% |
| normalization/preprocessing | 585 | 6.932 s | 11.850 ms | 4.77% | 5.94% |
| positive selection | 5,850 | 2.236 s | 0.382 ms | 1.54% | 1.91% |
| fold assignment/setup outside fold work | 195 | 2.042 s | 10.469 ms | 1.40% | 1.75% |
| reported PEP calculation | 585 | 0.759 s | 1.297 ms | 0.52% | 0.65% |

The corrected representative largest-file trace gives the requested per-fold breakdown. Q-value time
includes initial-direction and ten iteration passes, so it overlaps setup/training.

| Fold | Setup | Training | Iteration scoring | All q-values | Held-out scoring | Fold total |
|---:|---:|---:|---:|---:|---:|---:|
| 0 | 173.950 ms | 222.125 ms | 26.616 ms | 195.485 ms | 3.953 ms | 403.959 ms |
| 1 | 172.026 ms | 223.651 ms | 26.171 ms | 194.776 ms | 4.128 ms | 404.524 ms |
| 2 | 171.798 ms | 218.014 ms | 25.819 ms | 195.167 ms | 4.116 ms | 397.930 ms |

The maximum/minimum fold-total ratio is 1.017. There is no representative straggler fold.

## 6. Semi-supervised iteration profile

Each row is the mean per fold invocation over 585 invocations (three corpus profiles × 65 files ×
three folds). Evaluation counts are aggregate counts over those 585 invocations.

| Iteration | Training | Scoring | q-values | Positive selection | Newton | Line search | Active set | Total |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0 | 19.139 ms | 1.363 ms | 3.892 ms | 0.425 ms | 9,000 | 8,523 | 8,523 | 28.758 ms |
| 1 | 3.129 ms | 1.330 ms | 3.578 ms | 0.372 ms | 4,245 | 4,005 | 4,005 | 11.940 ms |
| 2 | 1.966 ms | 1.344 ms | 3.600 ms | 0.372 ms | 3,108 | 2,592 | 2,592 | 8.667 ms |
| 3 | 1.566 ms | 1.338 ms | 3.560 ms | 0.372 ms | 2,598 | 2,049 | 2,049 | 8.202 ms |
| 4 | 1.376 ms | 1.325 ms | 3.591 ms | 0.373 ms | 2,346 | 1,782 | 1,782 | 8.051 ms |
| 5 | 1.244 ms | 1.348 ms | 3.575 ms | 0.387 ms | 2,139 | 1,569 | 1,569 | 7.914 ms |
| 6 | 1.143 ms | 1.328 ms | 3.591 ms | 0.402 ms | 1,950 | 1,374 | 1,374 | 7.822 ms |
| 7 | 0.971 ms | 1.318 ms | 3.603 ms | 0.373 ms | 1,746 | 1,167 | 1,167 | 7.655 ms |
| 8 | 1.061 ms | 1.326 ms | 3.669 ms | 0.373 ms | 1,548 | 1,122 | 1,125 | 7.790 ms |
| 9 | 0.733 ms | 1.339 ms | 3.561 ms | 0.373 ms | 1,380 | 795 | 795 | 7.357 ms |

Iteration 0 is still disproportionate: it is 27.61% of aggregate iteration time, 59.20% of all SVM
training time, and 30.0% of Newton iterations. It is much cheaper than the historical profile,
however: mean training fell from 39.394 to 19.139 ms and total iteration time from 44.119 to
28.758 ms. From iteration 1 onward, q-value sorting is larger than training and remains almost flat.

## 7. SVM solver

The profiles contain 5,850 training calls and 30,060 Newton iterations. Mean recorded Newton-iteration
time is 0.607 ms. `line_search_total`, `line_search_objective_evaluation`, and
`active_set_and_margin_scoring` are the same nested work viewed at different boundaries and must not
be added together.

| Component | Calls | Total | Mean/call | % total | % SVM |
|---|---:|---:|---:|---:|---:|
| line search, inclusive | 24,213 | 10.116 s | 0.418 ms | 6.96% | 53.49% |
| line-search objective evaluation | 24,978 | 10.101 s | 0.404 ms | 6.95% | 53.41% |
| active-set/margin scoring | 24,981 | 10.093 s | 0.404 ms | 6.94% | 53.37% |
| Hessian construction | 24,213 | 7.009 s | 0.289 ms | 4.82% | 37.06% |
| gradient | 30,060 | 1.015 s | 0.034 ms | 0.70% | 5.37% |
| initial objective/active set | 5,850 | 0.496 s | 0.085 ms | 0.34% | 2.62% |
| allocation/buffer initialization | 5,850 | 0.166 s | 0.028 ms | 0.11% | 0.88% |
| Cholesky factorization | 24,213 | 0.0317 s | 0.00131 ms | 0.022% | 0.168% |
| linear solve | 24,213 | 0.0124 s | 0.00051 ms | 0.009% | 0.065% |
| buffer/weight/convergence bookkeeping | 79,251 | 0.0052 s | <0.0001 ms | <0.01% | 0.03% |

Cholesky, the direct solve, convergence, and solver bookkeeping are conclusively not useful targets.
Future SVM work would have to improve the repeated row pass and Hessian construction while preserving
the exact optimization result.

## 8. Q-value and sort profile

The fresh instrumentation found 10,335 q-value/count/mask passes per corpus, not the 33 calls seen by
the stale historical probes. Current fast paths for initial-direction counts and training masks had
not been covered by those probes. This is why the current hotspot ranking is materially different.

| Q-value phase | Calls (3 corpora) | Elements | Total | Mean/call | % total |
|---|---:|---:|---:|---:|---:|
| initial direction | 24,570 | 2,177,215,992 | 36.187 s | 1.473 ms | 24.89% |
| semi-supervised iterations | 5,850 | 518,384,760 | 21.189 s | 3.622 ms | 14.57% |
| final full-list PSM statistics | 195 | 25,919,238 | 1.170 s | 6.001 ms | 0.80% |
| reported/competed PSM statistics | 195 | 5,189,574 | 0.196 s | 1.004 ms | 0.13% |
| peptide statistics | 195 | 4,513,410 | 0.168 s | 0.863 ms | 0.12% |

The protein workload adds two small protein q-value calls over 6,502 elements, totaling 0.093 ms.
The two calls are the picked and classic comparison statistics; picked assignment itself is 0.075 ms.

| Sort | Calls | Elements | Total | Mean size | Mean/call | % total |
|---|---:|---:|---:|---:|---:|---:|
| q-value score order | 18,720 | 1,642,614,978 | 35.436 s | 87,747 | 1.893 ms | 24.37% |
| initial-direction reversed-rank order | 12,285 | 1,088,607,996 | 7.476 s | 88,613 | 0.609 ms | 5.14% |
| result-row score order | 780 | 9,702,984 | 0.369 s | 12,440 | 0.473 ms | 0.25% |
| peptide input order | 195 | 4,513,410 | 0.042 s | 23,146 | 0.217 ms | 0.03% |

| Q-value sub-operation | Calls | Total | % total |
|---|---:|---:|---:|
| score and reverse-rank sorting | 31,005 | 42.912 s | 29.51% |
| tie grouping and cumulative target/decoy counting | 31,005 | 14.751 s | 10.14% |
| allocation/buffer setup | 31,005 | 0.677 s | 0.47% |
| monotonic reverse scan on materialized reported q-values | 585 | 0.044 s | 0.03% |
| positive-mask materialization | 5,850 | 0.010 s | <0.01% |

The validated numeric equality rule for score ties, total-order handling, `+1` reported safeguard,
training-only no-safeguard rule, and strict threshold behavior are constraints. An optimization must
reproduce them exactly. The reverse-rank sort is especially actionable because the opposite
orientation can potentially be scanned from the already computed score order without changing tie
groups; that hypothesis belongs in the next goal, with byte-for-byte and adversarial validation.

## 9. PIN parser

One corpus profile processes 2,295,401,156 bytes, 8,639,746 rows, and 181,434,666 feature floats.
The aggregate three-corpus throughput is 274.3 MiB/s.

| Parser component | Aggregate time | Mean per corpus | % total |
|---|---:|---:|---:|
| parser total | 23.939 s | 7.980 s | 16.46% |
| row loop, inclusive | 23.928 s | 7.976 s | 16.46% |
| feature/integer/ExpMass numeric parsing | 7.841 s | 2.614 s | 5.39% |
| string allocation/copy | 6.113 s | 2.038 s | 4.20% |
| field splitting | 3.419 s | 1.140 s | 2.35% |
| residual row-loop work | 6.555 s | 2.185 s | 4.51% |
| mmap setup | 0.0024 s | 0.0008 s | <0.01% |
| header/feature validation | 0.0016 s | 0.0005 s | <0.01% |

The residual is bounds/validation/error-path checks, vector writes and growth, row-loop branches,
page touching, and profiling clock overhead. One corpus copies 925,703,237 string bytes through
25,919,238 row-string allocations. The annotated final column capacities total about 2.76 GiB across
the 65 processes. Protein mappings are deliberately retained as raw strings during parsing; mapping
tokenization is therefore zero at the parser boundary and is measured later in the protein path.

No standard run incurred a major page fault. Mmap setup is negligible, while CPU samples land in
`pin::parse` (5.29% leaf / 16.48% inclusive), `parse_f64` (3.20%), `split_fields` (1.94%), libc
`malloc` (1.42%), and `memcpy` (7.73% cross-cutting). The parser is primarily CPU/allocation and
copy bound, not storage-I/O or mmap-setup bound. N=4 throughput falls to 262.0 MiB/s and parser share
to 14.66%, consistent with modest shared cache/memory pressure. Historical throughput was
310.5 MiB/s; the current 274.3 MiB/s is 11.7% lower and is authoritative for the repaired code.

## 10. Output profile

One corpus emits 3,234,328 PSM/peptide rows, formats 9,702,984 score/q/PEP numbers, and writes
446,895,839 bytes. Formatting retains the exact six-decimal byte representation.

| Output work | Mean per corpus | % total |
|---|---:|---:|
| output stage, exclusive | 0.695 s | 1.43% |
| numeric formatting + borrowed string copies + buffering | 0.412 s | 0.85% |
| pre-output row sort | 0.123 s | 0.25% |
| actual file writes/flushes | 0.155 s | 0.32% |
| PSM row construction | 0.031 s | 0.06% |
| peptide row construction | 0.025 s | 0.05% |
| file creation | 0.004 s | <0.01% |

There are 260 result buffers per corpus at 1 MiB each (260 MiB of cumulative buffer capacity) and
817 underlying writes/flushes after buffering. Standard non-ensemble rows borrow ID, peptide, and
protein text; measured PSM output-row string allocations are zero. Buffer growth is therefore fixed
and bounded, not incremental.

The 0.85% user-space serialization timer is a strict upper bound on score/q/PEP formatting because it
also includes borrowed string copies into the buffer. Per-field clocks would cost more than these
short operations and materially distort them, so numeric and string-copy wall times were not split
further. Allocation annotations and source-level operation counts provide the separation: numeric
formatting uses a 32-byte stack buffer, borrowed strings allocate nothing, buffer capacity is fixed,
and kernel writes are timed independently. Output is no longer a meaningful optimization target.

## 11. Peptide, protein, join, and RT paths

### Peptide path on PXD032157

| Operation | Mean per corpus | % total |
|---|---:|---:|
| peptide stage | 0.327 s | 0.67% |
| identity/deduplication/representative selection | 0.194 s | 0.40% |
| peptide q-values | 0.056 s | 0.12% |
| peptide PEP | 0.029 s | 0.06% |
| deterministic representative-index sort | 0.014 s | 0.03% |
| output-row construction | 0.025 s | 0.05% |
| statistics-vector materialization | 0.003 s | <0.01% |

### Picked-protein path on `data/F_3.pin`

This real bacterial workload has 105,560 PSMs, 32 features, and 51,075 reported peptide entries.
Its total in-process time is 1.156 s; protein inference/output is 61.342 ms (5.31%).

| Protein operation | Time | % of process |
|---|---:|---:|
| complete peptide-mapping entry construction | 52.916 ms | 4.58% |
| → mapping union/tokenization over reported PSMs | 30.643 ms | 2.65% |
| → deterministic mapping-string materialization | 18.479 ms | 1.60% |
| picked-protein inference, inclusive | 7.019 ms | 0.61% |
| evidence collection | 3.943 ms | 0.34% |
| evidence-set construction/grouping | 0.974 ms | 0.08% |
| group scoring | 0.517 ms | 0.04% |
| canonical group sort | 0.516 ms | 0.04% |
| target/decoy pairing | 0.352 ms | 0.03% |
| final score sort | 0.293 ms | 0.03% |
| pairing-key sort | 0.108 ms | <0.01% |
| picked competition | 0.054 ms | <0.01% |
| protein q-values/assignment | 0.075 ms | <0.01% |
| serialization | 0.451 ms | 0.04% |
| actual writes | 0.070 ms | <0.01% |

The cost is preserving the complete set-valued peptide-to-protein mapping, not q-values or picked
competition. Protein `PEP=NA` behavior and grouping/pairing outputs were verified unchanged.

On the joined workload, canonical row/part ordering costs 21.016 ms. On the RT workload, explicit RT
setup/preprocessing costs 45.453 ms; the extra RT columns also increase ordinary normalization,
initial-direction, sorting, and model-scoring work. Neither conditional feature affects the standard
65-file profile.

## 12. Allocations and memory

Allocation counting was isolated from timing conclusions. Global allocator totals are exact for the
profile interval; named sites are capacity/traffic annotations and are approximate. Cumulative traffic
is not peak resident memory.

| Workload | Allocation calls | Allocated bytes | Deallocation calls | Deallocated bytes | Peak RSS |
|---|---:|---:|---:|---:|---:|
| largest PIN, one thread | 921,869 | 690,881,392 (658.9 MiB) | 119,491 | 556,428,232 | 204.5 MiB |
| full corpus, sequential | 30,066,074 | 22,168,119,374 (20.65 GiB) | 3,984,622 | 17,870,102,874 | 204.9 MiB |
| full corpus, N=4 | 30,066,074 | 22,168,119,374 (20.65 GiB) | 3,984,622 | 17,870,102,874 | 786.5 MiB |

Normal-build median peak RSS is 204.2 MiB on the largest one-thread run, 336.1 MiB in parallel-fold
mode, 203.9 MiB for the full sequential harness, 391.8 MiB at N=2, 781.3 MiB at N=4, and
1,162.0 MiB at N=6. Cumulative allocations are identical between sequential and N=4 because the work
is identical; concurrent live working sets drive the peak. Standard timing runs had zero major faults.

| Largest annotated full-corpus sites | Calls | Approx. bytes |
|---|---:|---:|
| normalized design matrices | 195 | 4.562 GB |
| PIN column-vector capacities | 455 | 2.966 GB |
| iteration training buffers | 7,800 | 2.739 GB |
| positive/negative row vectors | 3,900 | 0.995 GB |
| parser row strings | 25,919,238 | 0.926 GB |
| initial-direction buffers | 780 | 0.346 GB |
| q-value materialization buffers | 585 | 0.285 GB |
| result output buffers | 260 | 0.273 GB |
| SVM work buffers | 1,755 | 0.273 GB |
| PEP/PAVA buffers | 1,170 | 0.269 GB |

The repeated counts identify recreated fold/file/iteration buffers, while `memcpy` at 7.73% of
sequential sampled CPU corroborates significant data movement. Allocation traffic is a secondary,
cross-cutting opportunity and likely contributes to N=4 contention, but the campaign did not change
buffer reuse or data layout.

## 13. CPU sampling hotspots

| Configuration | Samples | score sort | rank sort | `memcpy` | active set | SVM | model score | parser |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| largest, one thread | 2,414 | 21.62% | 4.56% | 8.12% | 7.37% | 7.04% | 6.05% | 4.81% |
| largest, parallel folds | 2,035 | 20.59% | 4.28% | 6.54% | 6.83% | 6.54% | 5.95% | 5.55% |
| full sequential | 24,037 | 21.96% | 4.62% | 7.73% | 7.14% | 6.55% | 6.26% | 5.29% |
| full N=4 | 85,781 | 19.79% | 4.33% | 8.75% | 8.04% | 6.52% | 7.47% | 4.82% |

Those are leaf percentages. Relevant inclusive frames are:

| Inclusive frame | Full sequential | Full N=4 |
|---|---:|---:|
| `pipeline::rescore` / `percolator::run` | 80.24% | 81.58% |
| `cv_scores` fold closure | 77.60% | 78.90% |
| `train_fold` | 43.01% | 44.54% |
| `fold_setup` | 33.29% | 32.89% |
| q-value score quicksort | 23.01% | 20.72% |
| `pin::parse` | 16.48% | 14.85% |
| `target_count_at_fdr_into` | 14.65% | 13.51% |
| `svm::train` | 13.85% | 14.75% |
| reverse-rank target count | 8.95% | 7.96% |

The dominant call stacks are:

1. `main → pipeline::rescore → percolator::run → cv_scores → fold_setup → initial_direction → target_count_at_fdr_into → sort_score_order → slice quicksort`.
2. `... → FoldSetup::train_and_score → train_fold → target_mask_at_fdr_into → sort_score_order`.
3. `... → train_fold → svm::train → f_and_active / Hessian construction`.
4. `main → pin::parse → split_fields / parse_f64 / String copy → malloc/memcpy`.

The sampled profile independently confirms the wall-event ranking. Flamegraphs, protobuf profiles, and
collapsed stacks were generated for all sampled processes during acquisition; the aggregate retains
the symbol/sample totals so flamegraph generation is not required to reproduce the conclusions.

## 14. Parallel scaling

| Largest PIN mode | Median | Speedup | Efficiency |
|---:|---:|---:|---:|
| 1 | 1.616816 s | 1.000× | 100.0% |
| 2 | 0.890906 s | 1.815× | 90.7% |
| 3 | 0.890897 s | 1.815× | 60.5% |
| 6 | 0.890598 s | 1.815× | 30.3% |

For fixed-C CV, `num_threads == 1` selects serial folds and every value greater than one selects the
global Rayon pool; it is a serial/parallel switch rather than a pool-size limit. There are only three
outer folds. This explains why modes 2, 3, and 6 are indistinguishable. The corrected fold trace also
rules out imbalance. Parsing, final statistics, peptide work, and output stay serial, while three
parallel normalized matrices raise largest-file peak RSS from 204 to 336 MiB.

| Full corpus file concurrency | Median | Speedup | Efficiency |
|---:|---:|---:|---:|
| N=1 | 49.619487 s | 1.000× | 100.0% |
| N=2 | 27.040644 s | 1.835× | 91.7% |
| N=4 | 15.482359 s | 3.205× | 80.1% |
| N=6 | 12.058479 s | 4.115× | 68.6% |

N=4 stage profiles average 56.985 s of summed process work versus 48.469 s sequential, a 17.57%
inflation. The four slots are occupied an average of 3.704 ways (92.6%); the remaining loss combines
tail scheduling with per-process cache/memory/allocator contention. Under N=4, SVM work inflates 29%,
model scoring 41%, normalization 23%, output serialization 33%, q-values 5.7%, and parsing 4.7% in
summed process time. Output remains only 1.55%, there is no shared output file/lock, and major faults
are zero, so disk serialization is not the scaling limit. Peak RSS rises almost linearly with file
concurrency. N=6 is useful comparison evidence but is not a substitute for the requested N=4 target.

## 15. Refactor/module-boundary cost

There is no measured evidence that the modular refactor introduced meaningful overhead:

- ordinary one-input joining is 0.00003% of process time;
- PSM and peptide row construction together are 0.115%;
- peptide statistics-vector materialization is 0.006%;
- standard output rows borrow all strings and allocate zero row strings;
- result sorting/serialization occurs once at the output boundary, not in duplicate;
- the hot SVM path uses static Rust calls and an enum match, not heap dynamic dispatch;
- thin LTO and one codegen unit allow cross-module inlining;
- no repeated canonicalization is visible in the ordinary single-file path.

The 4.77% normalization cost and 4.56 GB of cumulative matrix capacities are algorithmic working-set
costs, not copies caused by moving code between modules. Joined-input canonicalization is real and
measured at 5.90% on its conditional workload, but it implements the validated permutation invariant
and is not paid in the 65 independent-file benchmark. Protein mapping materialization is similarly a
scientific data-preservation cost on the enabled protein path.

## 16. Comparison with the historical profile

| Historical hotspot | Historical | Current | Classification | Evidence |
|---|---:|---:|---|---|
| active-set/margin evaluation | 24.63% | 6.94% | **SHRANK** | previous SVM work succeeded; 72% relative reduction in share |
| all SVM training | 31.49% | 13.01% | **SHRANK** | iteration-0 training halved; factorization still negligible |
| q-value sorting | 15.44% | 29.51% | **GREW / newly exposed** | current fast count/mask paths were absent from stale instrumentation |
| initial direction | 13.23% | 29.24% | **GREW / moved here** | fold-local repaired selection evaluates both directions for every feature |
| PIN parsing | 12.27% | 16.46% | **GREW** | 310.5 → 274.3 MiB/s; stricter repaired parser plus shifted denominator |
| output formatting | 10.93% | 0.85% | **NO LONGER IMPORTANT** | borrowed rows, stack numeric formatting, fixed buffers |
| peptide processing | 7.35% | 0.67% | **NO LONGER IMPORTANT** | hash dedup and compact materialization |
| Hessian construction | 5.24% total | 4.82% total | **UNCHANGED** | similar share; 37.06% of the smaller SVM total |
| actual writes | 0.88% | 0.32% | **SHRANK** | 1 MiB buffering; tmpfs output in both current measurements |

The main cost moved from SVM/output/peptide work into repeated exact q-value ordering and the repaired
fold-local initial direction. Comparisons involving q-values are conservative because the historical
instrumentation did not observe all current fast paths; the current numbers are the trustworthy ones.

## 17. Ranked opportunities — diagnosis only

Rows overlap where noted, so maximum gains are ceilings and are not additive. “Measured share” uses
sequential process time unless it explicitly says N=4 or conditional.

| Rank | Exact path | Measured share / scaling | Plausible direction for a later goal | Max ceiling | Difficulty | Scientific / determinism / regression risk |
|---:|---|---|---|---:|---|---|
| 1 | `percolator::initial_direction` | 29.24%; 28.24% at N=4 | precompute exact per-feature orders once per file, filter by fold, and scan both orientations with identical tie groups | 29.24% | high | medium / medium / medium |
| 2 | `stats::sort_score_order` in all q-value/count/mask calls | 24.37%; 21.52% at N=4 | exact finite-`f64` ordering keys, specialized index sorting, or reuse where scores provably do not change | 24.37% | high | medium / medium / high |
| 3 | cross-cutting N=4 working-set contention | +17.57% process work; 3.205× at N=4 | reduce matrix/buffer traffic and live working sets while retaining N=4 file concurrency | ~14.9% N=4 wall if all inflation vanished | high | low / low / high |
| 4 | `pin::parse` | 16.46%; 14.66% at N=4 | reduce row-string allocations/copies and improve capacity estimates without weakening validation | 16.46% | medium | low / low / medium |
| 5 | semi-supervised `target_mask_at_fdr_into` | 14.57% inclusive; flat across iterations | specialize the threshold-only exact scan and reuse workspaces/order representations | 14.57% | medium | medium / medium / medium |
| 6 | `svm::train` row passes | 13.01%; 14.27% at N=4 | fuse/cache exact active/Hessian row data or improve vectorization without changing convergence | 13.01% | high | medium / medium / high |
| 7 | q-value tie grouping/cumulative counting | 10.14%; 9.24% at N=4 | improve contiguous label/order locality or fuse exact scans | 10.14% | medium | high / medium / medium |
| 8 | `FoldModel::score_rows` | 6.17%; 7.40% at N=4 | improve row locality/batching around the existing vectorized dot product | 6.17% | medium | low / medium / medium |
| 9 | `sort_reversed_rank_order` | 5.14%; 4.56% at N=4 | scan the already built score order in reverse instead of sorting ranks again | 5.14% | low-medium | medium / medium / low |
| 10 | normalization/design-matrix materialization | 4.77%; 4.97% at N=4; 4.56 GB capacities/corpus | reduce three fold-local full-matrix materializations or reuse storage exactly | 4.77%, potentially more through contention | high | medium / medium / high |

The lowest-risk first experiment is rank 9, but it cannot deliver a large end-to-end gain alone. The
round is justified by ranks 1-2 as a combined exact-order problem. Parser and SVM should remain
secondary until the q-value hypothesis is validated. Output, PEP, peptide reporting, Cholesky, solve,
and protein q-values should not be optimized now. On protein-enabled workloads, mapping entry
construction is a separate 4.58% conditional opportunity, but it does not affect the standard target.

## 18. Amdahl ceilings and target feasibility

The following applies the measured N=4 process fractions to the 15.482359 s normal median. “20%”
means removing 20% of that hotspot's cost. It is optimistic: scheduler/tail behavior is assumed
unchanged. Initial direction and sorts are nested inside all q-value work and cannot be combined by
adding their rows.

| N=4 hotspot | N=4 share | 20% less hotspot | 50% less | 75% less | Eliminated |
|---|---:|---:|---:|---:|---:|
| all q-value work | 36.42% | 14.355 s | 12.663 s | 11.253 s | 9.844 s |
| initial direction | 28.24% | 14.608 s | 13.297 s | 12.204 s | 11.111 s |
| score-order sorts | 21.52% | 14.816 s | 13.816 s | 12.983 s | 12.150 s |
| PIN parser | 14.66% | 15.028 s | 14.347 s | 13.780 s | 13.213 s |
| SVM training | 14.27% | 15.040 s | 14.377 s | 13.825 s | 13.272 s |
| tie-group/count scans | 9.24% | 15.196 s | 14.767 s | 14.409 s | 14.052 s |

### Feasibility of N=4 below 12 s

Required total reduction: 3.482 s, or 22.49% (1.290× speedup). No single realistic isolated
optimization reaches it. It would require 61.76% less total q-value work, or a combination such as a
large initial-direction/order reduction plus smaller q-scan, parser/SVM, and N=4 contention gains.
Because initial direction has an exact-order reuse hypothesis and the reverse-rank sort is removable
in principle, this is **plausible but not yet demonstrated**. The next optimization round should use
12 s as a stretch outcome, not as permission to change statistical behavior.

### Feasibility of N=4 below 10 s

Required total reduction: 5.482 s, or 35.41% (1.548× speedup). It would require 97.23% less q-value
cost in isolation. Eliminating every score sort is insufficient, and eliminating every q-value
operation—scientifically impossible—only predicts 9.844 s. A safe route would have to win across
q-values, parsing, SVM/scoring, normalization, and contention simultaneously. On current evidence,
**below 10 s is not a realistic correctness-preserving N=4 target**.

## 19. Validation and final recommendation

Before instrumentation, `python3 refactor/verify_baseline.py` passed release tests, shell gates, exact
frozen outputs, and adversarial behavior. The release suite reported 138 passing tests. The final
profiling-feature build was also checked directly against the frozen artifacts and matched:

- fixed-C serial and parallel outputs;
- selected-C outputs;
- ensemble outputs;
- PSM, peptide, and protein files;
- protein grouping and picked-protein `PEP=NA` output;
- repeated deterministic output;
- joined-input permutation/path-alias invariants and the frozen adversarial summary.

Across the acquisition matrix, 29 normal/stage/CPU/allocation comparisons were byte-identical; the
corrected timing campaign contains 26 such comparisons and no failures. The final source passes
formatting, strict Clippy for all targets/features, release tests, and the baseline verifier.

**Final answer:** the refactored code has not reached uniform diminishing returns. It has one newly
measured, sufficiently large family of work—exact q-value ordering/scanning concentrated in
initial-direction selection—that justifies another optimization round. That round should begin with
an exactness proof and adversarial tests for order reuse/reverse scanning, then benchmark after every
small change. Stop if those experiments cannot reduce q-value time substantially; do not reopen
output, peptide, PEP, protein statistics, or solver-factorization work without new measurements.
