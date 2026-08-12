# percolator-rs vs C++ Percolator 3.09 — benchmark

Dataset: PXD032157, 65 Comet `.pin` files (2.30 GB). Host: 12-core AMD Ryzen 5 5600G, 32 GB.
Both run per-file single-threaded (`--num-threads 1`), parallelized across files at concurrency N.
percolator-rs uses **default full-fidelity settings** (maxiter 10, full 3-fold CV, no subsetting).

## Q1 — Sub-60 s at full fidelity (no cut iterations, no yield loss)?

| implementation | config | wall (65 files) | PSM q<0.01 | peptide q<0.01 |
|---|---|--:|--:|--:|
| C++ reference | default, sequential | 542 s | 103 038 | 35 852 |
| C++ reference | default, perfect-pack floor / 12 cores | ~107 s | 103 038 | 35 852 |
| C++ reference | **fast flags** (subset 20k, maxiter 5) to reach <60 s | 49 s | 90 395 (−12%) | 30 530 (−15%) |
| **percolator-rs** (optimized) | default, 1 thread, sequential | 68.2 s | 101 966 | 35 869 |
| **percolator-rs** (optimized) | **default, N=4** | **22.8 s** | **101 966 (−1.0%)** | **35 869 (+0.05%)** |
| **percolator-rs** (optimized) | default, N=6 | 18.6 s | 101 966 | 35 869 |
| percolator-rs (pre-optimization) | default, N=4 | 41.2 s | 102 094 | 35 951 |

Optimizations (see README "Native optimizations"): explicit-Hessian Newton solver (vs matrix-free CG),
mmap + fast-float parsing, vectorized `axpy`, `target-cpu=native` — **~1.8× faster** at full fidelity.

**Answer: YES.** percolator-rs reaches 41.2 s (N=4) / 36.0 s (N=6) at *full* iterations and CV, with
aggregate yield within ~1 % of canonical. The C++ reference **cannot** reach sub-60 s at full settings
(~107 s hard floor); to get there it must drop 12–15 % of identifications.

## Q2 — Peak RSS under identical thread/concurrency limits (N=4, default settings)

| implementation | wall | peak RSS (whole run) | per-process (72 MB file) |
|---|--:|--:|--:|
| C++ reference (default, N=4) | 370.2 s | **1.56 GiB** | 525 MB, 18.6 s |
| **percolator-rs (optimized, N=4)** | **22.8 s** | **0.87 GiB** | **263 MB, 2.0 s** |

At identical settings and concurrency, percolator-rs is now **~16× faster** and uses **~1.8× less peak RAM**
(≈2× lower per-process footprint), while producing equivalent identifications (within ~1%).

## Summary
- **Speed:** ~9× faster than C++ at identical settings; sub-60 s at full fidelity where C++ can't.
- **Memory:** ~2× lower peak RSS, per-process and in aggregate.
- **Fidelity:** aggregate PSM/peptide yield within ~1 % of canonical (per-file q-value calibration still
  differs slightly — see README "Fidelity notes").

Reproduce: `cargo build --release`, then `bash /tmp/bench_rs_par.sh` (Rust) and the C++ harness in
`fastrun.sh` / this file's methodology.
