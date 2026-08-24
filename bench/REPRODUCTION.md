# Reproduction ledger

This ledger records the complete rerun performed on 2026-08-24 from Git commit
`59b68bb3791f7bf8b4c1dbe2453150045d956f96`. The host was Linux 7.0.0-30 on an AMD Ryzen 5
5600G with 12 hardware threads. Rust was 1.97.0; the reference executable reported C++ Percolator
3.09.0. Bulk inputs and generated outputs remain outside Git under
`$HOME/percolator_rs_out`.

All portable gates passed: 19 Rust unit tests, 10 integration tests, exact SVM, MLP,
nested-selection, feature-report, ensemble, and picked-protein shell regressions, plus four PrEST
report tests and the Sage-normalizer order-independence test. The self-hosted 65-file gate also
passed exactly at 107,046 PSMs and 37,469 peptides.

## Fresh study results

| Study | Fresh result | Status / artifact |
|---|---|---|
| PXD032157 Rust canonical, N=4 | 20.768 s; 914,972 KiB; 107,046 PSMs; 37,469 peptides | Repaired driver rerun, 65/65 valid; `$HOME/percolator_rs_out/canonical` |
| PXD032157 Rust balanced, N=4 | 20.4 s; 0.87 GiB; 106,817 / 37,526 | 65/65 valid; `$HOME/percolator_rs_out/balanced` |
| PXD032157 Rust fast, N=4 | 14.6 s; 0.88 GiB; 105,237 / 36,772 | 65/65 valid; `$HOME/percolator_rs_out/fast` |
| PXD032157 C++ canonical auto-input, N=4 | 376.233 s; 1,564,824 KiB; 103,038 / 35,852 | Repaired driver rerun, 65/65 valid; `$HOME/percolator_rs_out/cpp-canonical` |
| PXD032157 C++ trimmed, N=5 | 59.398 s; 1,219,140 KiB; 90,395 / 30,530 | Repaired driver rerun, 65/65 valid; `$HOME/percolator_fast_out/PXD032157` |
| Manifest runner, sequential | Rust 61.972 s, 269,368 KiB, 107,046 / 37,469; C++ 1,034.110 s, 497,112 KiB, 102,781 / 35,765 | 65/65 each, no failures; explicit C++ `--search-input concatenated`; versioned JSON and comparison under `$HOME/percolator_rs_out/benchmark-dataset-reproduction-20260824/PXD032157` |
| Four compact datasets | Rust/C++ PSMs: Tide 29,264/27,617; MSFragger 1,554/1,388; Sage 26,624/25,795; yeast 1,126/1,147 | All source/generated checksums passed; `bench/multidataset/recorded-results.tsv` |
| Pure null, three default files | Select-C/fixed false targets: 0/0, 5/5, 3/1 | Conservative at reported q<0.01; `bench/null-calibration-results.tsv` |
| Retention-time features, three deterministic files | PSM deltas −4.81%, +3.98%, +2.00% | Three repeats; pinned inputs and hashes; `bench/advanced-feature-results.tsv` |
| Joint training, four smallest files | 1,524 → 1,606 PSMs (+82, +5.38%); 3/4 improved | Pinned size rule and hashes; `bench/ADVANCED_FEATURES.md` |
| Intra-file threading, largest file | Fixed 2.05 → 1.39 s; select-C 3.93 → 1.85 s | Three-run medians; outputs byte-identical |
| Six-run signal-present entrapment | At q≤0.01: Rust 19,666 at 2.784% adjusted FDP; C++ 19,126 at 2.616% | All source checksums and the 50% amino-acid database rebuild passed; `$HOME/percolator_rs_out/entrapment/report.tsv` |
| SVM versus MLP, PXD032157 | SVM 107,046 / 37,469 in 20.081 s; MLP 105,521 / 36,447 in 119.761 s | `bench/model-comparison-results.tsv` |
| SVM versus MLP, entrapment | SVM 19,666 at 2.784% FDP; MLP 19,382 at 2.600% | Neither validates nominal 1% FDR |
| Fixed versus legacy versus nested, PXD032157 | Fixed 107,046 / 37,469 in 20.056 s; legacy 106,558 / 37,330 in 49.703 s; nested 106,652 / 37,636 in 206.411 s | `bench/automatic-selection-results.tsv` |
| Fixed versus nested, entrapment | Fixed 19,666 at 2.784% FDP; nested 19,556 at 2.675% | Neither validates nominal 1% FDR |
| Picked versus Bayesian protein inference | All five cases converged; the deterministic Sage case is picked 6,062 groups / 1,597 accepted and Bayesian 6,150 / 1,540 | `bench/protein-inference-results.tsv` |
| PrEST protein calibration | Selected α=0.1, β=0.0001, γ=0.001; test aggregate accepted/absent: picked 613/38, fixed Bayesian 824/118, selected Bayesian 685/14 | Pipeline key `8196aa5ed3caadf4`; `bench/protein-calibration-results.tsv` |
| Real single-organism protein check | 1,410 picked-FDR groups versus 1,369 classic | Picked adds 41 groups at q<0.01 |

Counts in the first five rows use strict reported q<0.01. The compact, entrapment, and PrEST
drivers use their documented q≤0.01 convention. Timing is workload- and host-sensitive; exact
identification counts, checksums, convergence, file counts, and empty failure tables are the primary
reproduction checks.

## Commands

Run portable and self-hosted gates:

```bash
cargo test --release
bash tests/regression.sh
bash tests/model_regression.sh
bash tests/selection_regression.sh
bash tests/feature_report.sh
bash tests/ensemble_regression.sh
bash tests/protein_regression.sh
python3 -m unittest -v bench/multidataset/test_normalize_sage_pin.py
(cd bench/protein_calibration && python3 -m unittest -v test_report.py)
bash bench/regression.sh
```

Run the benchmark studies after installing their documented external dependencies and data:

```bash
REPEATS=3 bash bench/multidataset/run.sh
REPEATS=3 bash bench/protein_inference.sh
bash bench/null_calibration.sh 3 42
bash bench/advanced_features.sh
bash bench/entrapment/run.sh
bash bench/model_entrapment.sh
bash bench/selection_entrapment.sh
REPEATS=3 bash bench/protein_calibration/run.sh
bash bench/model_comparison.sh
bash bench/selection_comparison.sh
bash bench/run_rs.sh canonical 4
bash bench/run_rs.sh balanced 4
bash bench/run_rs.sh fast 4
bash bench/run_cpp.sh 4
bash bench/fastrun.sh 5
bash bench/protein_real.sh
```

The manifest-driven run uses a new output directory each time:

```bash
cargo run --quiet --bin validate-benchmark-manifests -- bench/datasets.toml
PERCOLATOR_BENCH_DATA="$PWD/data" \
LD_LIBRARY_PATH="$HOME/opt/perc-libs:${LD_LIBRARY_PATH:-}" \
cargo run --release --bin benchmark-dataset -- \
  --dataset PXD032157 \
  --output "$HOME/percolator_rs_out/benchmark-dataset-$(date +%Y%m%d)" \
  --percolator "$HOME/opt/percolator-root/usr/bin/percolator"
```

## Reproducibility repairs found during the rerun

Sage 0.14.6/0.14.7 writes identical records in nondeterministic order and derives `SpecId` from that
order. Both Sage normalizers now sort complete normalized records and assign stable sequential IDs;
two independent searches canonicalize to SHA-256
`6b6773606f3640eb86655485df35d94d7ad61779d5505dc6fd889a9d4c11c375`. A regression test protects
this invariant.

Crux can leave a generic `comet-out/comet.pin` beside the six explicitly named entrapment search
directories. The model and selection drivers now exclude this non-sample artifact while retaining
their exact-six guard. Finally, RSS samplers running under `set -euo pipefail` now tolerate the
expected initial interval in which `ps -C` finds no worker process.
