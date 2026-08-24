# Canonical runtime optimization: PXD032157

This report records the 2026-08-24 optimization campaign for the canonical Rust
workflow. The objective was to reduce the wall time of all 65 PXD032157 Comet PINs at
file-level concurrency N=4 without changing scientific behavior. The required target was a
median below 15 seconds; below 12 seconds was a stretch target.

## Result

| Build | N=4 wall times | Median | Peak RSS | PSM q<0.01 | Peptide q<0.01 |
|---|---:|---:|---:|---:|---:|
| `ec4f098` baseline | 20.641 / 20.799 / 20.764 s | 20.764 s | not recorded | 107,046 | 37,469 |
| optimized | 12.063 / 12.002 / 12.434 s | **12.063 s** | **794,952 KiB median** | **107,046** | **37,469** |

The optimized median is 8.701 seconds, or **41.9%**, faster. It clears the required
sub-15-second target and finishes 63 milliseconds above the optional sub-12-second stretch
target. The release binary used the repository's portable `x86-64-v3` target and thin LTO;
`target-cpu=native` did not improve the large-file screen and was not adopted.

Every retained full-dataset candidate was checked against the baseline SHA-256 manifest. The
final three trials each reproduced all **260/260** target/decoy PSM/peptide TSV files byte for
byte. Thus the final result preserves not only aggregate yield but output values, row ordering,
ties, fold behavior, solver convergence, and formatting.

Sequential monitoring confirms that the N=4 improvement did not come from shifting work or changing
concurrency. The optimized full-dataset runs took 36.078, 36.278, and 36.265 seconds: a **36.265-second
median**, 41.3% below the documented 61.753-second pre-optimization run. Median peak RSS was 206,708
KiB and all 260 outputs again matched the baseline manifest. The optimized N=4 speedup over its
sequential median is 3.01x, essentially the same scaling as the baseline's 2.97x; the improvement is
therefore reduced per-file computation rather than increased parallelism.

## Measurement method

- Host: AMD Ryzen 5 5600G, 6 cores / 12 threads, Linux, Rust 1.97.0, LLVM 22.1.6.
- Input: 65 PIN files, 2,295,401,156 bytes, 8,639,746 PSMs.
- Command: `RS_BENCH_BIN=target/release/percolator-rs RS_BENCH_OUT=... bash bench/run_rs.sh canonical 4`.
- Sequential control: the same command with concurrency `1`, run three times after the final N=4
  candidate was fixed.
- Each promising change was first screened with five runs of the largest PIN, then measured with
  three full N=4 runs. A change was retained only when it improved the relevant screen or full
  median and preserved exact outputs.
- Candidate runs were isolated: measure, change, test, rebuild, measure, hash-compare, retain or
  revert. Superseded generated outputs were removed after their summaries and hashes were saved.

## Retained changes

The medians below are cumulative, so each row measures the complete retained build at that point.

| Change | Full N=4 median |
|---|---:|
| Baseline | 20.764 s |
| Count accepted targets without materializing initial q-values | 19.696 s |
| Reuse q-value buffers across folds | 19.483 s |
| Exact fixed-width dot product in the active SVM path | 19.149 s |
| Same dot product in whole-dataset scoring | 18.889 s |
| Guarded exact six-decimal formatter | 18.323 s |
| Reuse initial-direction sort storage | 18.166 s |
| Reuse fold-local SVM workspace | 18.012 s |
| Remove proven-redundant q-sort bounds checks | 17.460 s |
| Carry existing margins into SVM initialization | 16.671 s |
| Pack selected SVM training rows contiguously | 15.034 s |
| Borrow peptide deduplication keys | 14.708 s |
| Write constant TSV text directly | 14.521 s |
| Use a target byte mask during training | 14.206 s |
| Borrow output strings while preserving `Row` sort layout | 13.243 s |
| Sort negated initial directions by exact integer ranks | 12.865 s |
| Preallocate output rows exactly | 12.638 s |
| Use a capacity-sized `AHashMap` for peptide lookup | 12.432 s |
| Remove dead PEP buffers | 12.389 s |
| Reuse the packed training matrix allocation | 12.170 s |
| Share the final q-value/PEP score order | **12.063 s** |

The implementation deliberately retains unstable-sort behavior. In particular, borrowing output
strings changed `Row` size and therefore changed tie permutations inside `sort_unstable`; explicit
padding restores the old 96-byte layout and has a unit test. Fixed-width formatting and fixed-size
dot products also have bit-for-bit property tests. The paired q-value/PEP path is tested against the
separate calculations with ties, signed zero, infinities, NaN, and empty input.

## Rejected changes

The following candidates were reverted because they were neutral, slower, or failed exact-output
requirements: active-set branch hoisting; a distinct-score reverse-order shortcut; a fused row
parser; cached score structs whose changed element layout altered tie permutations; extra fold
buffer reuse; pointer-based dot products; `memchr` field splitting; fat LTO; native CPU tuning;
unchecked parser feature access; special formatter paths for zero and one; staging rows before the
buffered writer; and borrowing canonical q-value/PEP vectors. One unpadded borrowed-row version was
also rejected before the corrected layout-preserving version was retained.

## Final instrumented profile

The final profiling build ran all 65 files at N=4 in 12.187 seconds. Its 65 process runtimes sum to
45.798 seconds; nested percentages therefore use that summed in-process time and overlap their
parents.

| Stage or nested operation | Summed time | % process time | Calls |
|---|---:|---:|---:|
| Rescoring | 31.456 s | 68.68% | 65 |
| PIN parsing | 7.185 s | 15.69% | 65 |
| Result output | 4.275 s | 9.33% | 65 |
| Peptide-level processing | 2.486 s | 5.43% | 65 |
| PSM-level materialization | 0.294 s | 0.64% | 65 |
| SVM training | 8.202 s | 17.91% | 1,950 |
| Initial-direction selection | 7.075 s | 15.45% | 65 |
| Q-value calculation | 4.967 s | 10.84% | 2,080 |
| Score-index sorting for q-values | 3.804 s | 8.31% | 2,080 |
| Active-set and margin scoring | 4.766 s | 10.41% | 8,244 |
| Whole-model row scoring | 3.958 s | 8.64% | 2,145 |
| Result formatting/buffering | 2.764 s | 6.04% | 260 |
| Hessian construction | 2.521 s | 5.50% | 8,046 |
| PEP/PAVA | 0.380 s | 0.83% | 130 |

The remaining frontier is dominated by repeated score-index sorting, initial-direction feature
scans, parser allocation/copy work, and exact TSV number formatting. Further work in those areas
must continue to preserve unstable tie permutations and byte-identical output; changing statistical
thresholds, iteration limits, folds, convergence, or output precision is outside this optimization
contract.
