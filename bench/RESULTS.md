# Fast reference run — 65× PXD032157 percolator under 60 s

Goal: run all 65 Comet `.pin` files through the reference C++ Percolator 3.09 in **under 60 s wall**
with **minimal peak RAM**. Baseline (sequential, default settings) was **542 s** / 525 MiB peak.

## Why default settings can't hit 60 s
Total work at default settings = **1290 CPU-seconds** → a hard **~107 s floor** on 12 cores even at
perfect packing. Sub-60 s therefore requires cutting per-run work, using percolator's built-in
speed flags (`--subset-max-train`, `--maxiter`).

## Approach
- Per file: `--num-threads 1 --subset-max-train 20000 --maxiter 5` (full 3-fold CV kept; no `--quick-validation`).
- **Single work queue at concurrency N**, files scheduled **size-interleaved** (largest, smallest, 2nd-largest, …)
  so big and small files mix — this keeps the concurrent total file size (and thus peak RSS) steady,
  since percolator RSS ≈ 5× input file size.
- Results written to **local ext4** (`$HOME`), not the NTFS/ntfs-3g "New Volume" (FUSE writes are ~2× slower:
  the same run writing to the external drive took 95.7 s vs 58.6 s to ext4).

## Pareto frontier (compute only, output → /dev/null)
| N | wall | peak RAM |
|---|------|----------|
| 6 | 35.6 s | 1.69 GiB |
| 5 | 40.7 s | 1.36 GiB |
| 4 | 48.8 s | 1.04 GiB |
| 3 | 62.2 s ❌ | 0.74 GiB |

N=3 misses at default trim; trimming to `--maxiter 4` or `--subset-max-train 15000` brings N=3 to
56–58 s at ~0.8 GiB but with a thin margin and noisier peak.

## Verified runs (real output files, ext4, all 65 valid)
| config | wall | peak RAM | valid | target PSMs q<0.01 | target peptides q<0.01 |
|--------|------|----------|-------|--------------------|------------------------|
| **N=4** (min RAM)      | **58.6 s** | **0.88 GiB** | 65/65 | 90 395 | 30 530 |
| **N=5** (safe margin)  | **49.4 s** | 1.19 GiB     | 65/65 | 90 395 | 30 530 |

Baseline canonical (default, sequential): 542 s, 525 MiB, 103 038 PSMs / 35 852 peptides.

**Recommendation:** N=5 for a comfortable margin (49.4 s, ~1.2 GiB); N=4 if minimizing RAM is
paramount (58.6 s, 0.88 GiB). Both are ~9–11× faster than the sequential baseline.

## Accuracy tradeoff
The speed flags cost identifications: **−12% PSMs, −15% peptides** at q<0.01 vs the canonical default run.
These fast outputs are therefore **not canonical** — the canonical reference stays in `reference/PXD032157/`.
For a smaller accuracy hit (−5.6% PSMs) at ~50 s use `--subset-max-train 40000` at N=5 (slightly more RAM).

Run it: `bash bench/fastrun.sh` (writes to `$HOME/percolator_fast_out`).
