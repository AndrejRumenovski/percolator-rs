# percolator-rs vs C++ Percolator 3.09 — benchmark

> **These measurements predate the 2026-08-25 statistical repair** and describe commit `d83a7ba`,
> whose q-value and PEP estimators, cross-validation isolation and PIN feature selection were all
> subsequently found defective and replaced. They are kept as the record of what was measured then.
> For what the current implementation does and what it has been revalidated against, see
> [`../validation/REPAIR.md`](../validation/REPAIR.md).

Dataset: PXD032157, 65 Comet `.pin` files (2.30 GB). Host: 12-core AMD Ryzen 5 5600G, 32 GB.
Both run per-file single-threaded (`--num-threads 1`), parallelized across files at concurrency N.
percolator-rs uses the **complete default workload** (maxiter 10, full 3-fold CV, no subsetting).

## Q1 — Sub-60 s with the complete default workload (no cut iterations, no reported-q yield loss)?

| implementation | config | wall (65 files) | PSM q<0.01 | peptide q<0.01 |
|---|---|--:|--:|--:|
| C++ reference | default, N=4 | 376.233 s | 103 038 | 35 852 |
| C++ reference | **fast flags** (subset 20k, maxiter 5), N=5 | 59.398 s | 90 395 (−12%) | 30 530 (−15%) |
| **percolator-rs** (current) | default, 1 thread, sequential | 62.0 s | 107 046 | 37 469 |
| **percolator-rs** (current) | **default, N=4** | **20.768 s** | **107 046 (+3.9%)** | **37 469 (+4.5%)** |
| **percolator-rs** (older optimized build) | default, N=6 | 18.6 s | 101 966 | 35 869 |
| percolator-rs (pre-optimization) | default, N=4 | 41.2 s | 102 094 | 35 951 |

Optimizations (see README "Native optimizations"): explicit-Hessian Newton solver (vs matrix-free CG),
mmap + fast-float parsing, vectorized `axpy`, `target-cpu=native` — **~1.8× faster** on the complete default workload.

**Answer: YES.** percolator-rs reaches 20.768 s at N=4 with full iterations and CV. The C++
reference's observed full-setting N=4 run takes 376.233 s; the sub-60 trimmed run drops 12–15% of
identifications.

## Q2 — Peak RSS under identical thread/concurrency limits (N=4, default settings)

| implementation | wall | peak RSS (whole run) |
|---|--:|--:|
| C++ reference (default, N=4) | 376.233 s | **1,564,824 KiB (1.49 GiB)** |
| **percolator-rs (current, N=4)** | **20.768 s** | **914,972 KiB (0.87 GiB)** |

At identical file concurrency, percolator-rs is **18.1× faster** and the observed whole-run peak is
1.71× lower. It reports 3.9% more PSMs and 4.5% more peptides at the same nominal threshold; the
entrapment study, rather than count agreement, evaluates calibration.

## Summary
- **Speed:** 18.1× faster at N=4; sub-60 s on the complete default workload where the observed C++ run is not.
- **Memory:** 1.71× lower observed aggregate peak RSS at N=4.
- **Yield and calibration:** +3.9% PSMs and +4.5% peptides; see the README fidelity notes and
  signal-present entrapment study before treating nominal q-values as equivalent error rates.

Reproduce: `cargo build --release`, then run `bash bench/run_rs.sh canonical 4` and
`bash bench/run_cpp.sh 4`. The trimmed C++ run is `bash bench/fastrun.sh 5`.
