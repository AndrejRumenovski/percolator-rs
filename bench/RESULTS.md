# Fast reference run — 65× PXD032157 percolator under 60 s

Goal: run all 65 Comet `.pin` files through the reference C++ Percolator 3.09 in **under 60 s wall**
with minimal peak RAM. A clean 2026-08-24 rerun of the current driver finishes in **59.398 s** at
N=5, peaks at **1,219,140 KiB**, validates all 65 files, and reproduces 90,395 PSMs / 30,530
peptides exactly. The exploratory scaling figures below predate the final reproducibility ledger and
are retained as development history, not as fresh measurements.

## Historical default-work estimate
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

## Historical Pareto exploration (compute only, output → /dev/null)
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
| **N=5** (fresh rerun)  | **59.398 s** | 1.16 GiB   | 65/65 | 90 395 | 30 530 |

Historical baseline canonical (default, sequential): 542 s, 525 MiB, 103 038 PSMs / 35 852
peptides. The clean current canonical N=4 measurement is 376.233 s and 1,564,824 KiB aggregate peak;
see [`RS_VS_CPP.md`](RS_VS_CPP.md).

**Recommendation:** N=5 meets the target on this host, but the 59.398 s repeat leaves little timing
margin; treat the exact yields and 65/65 validation as stable, and the wall budget as host-sensitive.

## Reported-q yield tradeoff
The speed flags cost identifications: **−12% PSMs, −15% peptides** at q<0.01 vs the canonical default run.
These fast outputs are therefore **not canonical** — the canonical reference stays in `reference/PXD032157/`.
The historical exploration also found a smaller reported-q yield reduction (−5.6% PSMs) with
`--subset-max-train 40000`; that configuration was not part of the 2026-08-24 clean rerun.

Run it: `bash bench/fastrun.sh` (writes to `$HOME/percolator_fast_out`).
