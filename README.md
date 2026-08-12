# percolator-rs

[![CI](https://github.com/AndrejRumenovski/percolator-rs/actions/workflows/ci.yml/badge.svg)](https://github.com/AndrejRumenovski/percolator-rs/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

A from-scratch Rust reimplementation of the [Percolator](https://github.com/percolator/percolator)
semi-supervised PSM rescoring algorithm, built to benchmark against the reference C++ Percolator 3.09.

## Algorithm
Faithful to the Percolator method:
1. **PIN parse** — streaming tab-delimited reader (`ExpMass`/`CalcMass` excluded from features by name).
2. **Feature normalization** — per-feature z-score, bias column appended.
3. **Initial direction** — best single feature (either orientation) by target-decoy yield at q<0.01.
4. **Semi-supervised training** — iterate (default 10×): score → target-decoy q-values → pick confident
   targets (q<0.01) as positives, all decoys as negatives → retrain a class-weighted **L2-regularized
   squared-hinge linear SVM** (primal truncated-Newton / Newton-CG solver, the same L2-loss family as the
   reference L2-SVM-MFN).
   The class weights `Cpos`/`Cneg` are **selected per run by cross-validation** over a small grid
   (the reference's `Cpos=0` behaviour): each candidate gets an abbreviated 3-fold pass and the one
   with the highest out-of-fold yield at `q<0.01` wins. Pin them with `--cpos/--cneg`, or skip the
   search with `--no-select-c`.
5. **3-fold nested cross-validation** — each PSM scored by a model trained without its fold (no overfitting
   of the FDR estimate).
6. **Target-decoy q-values** + monotone (PAVA isotonic) PEP; PSM- and peptide-level output in the reference
   tab format.

## Build & run
```
cargo build --release
./target/release/percolator-rs --seed 1 \
  --results-psms t.psms --decoy-results-psms d.psms \
  --results-peptides t.pep --decoy-results-peptides d.pep \
  input.pin
```
Flags mirror the reference subset: `--seed`, `--maxiter`, `--subset-max-train`, `--num-threads`.
`--num-threads` (default 1) parallelizes the class-weight grid search and the 3 CV folds within a
single file — results are bit-identical at any thread count. It defaults to 1 so that harnesses which
already run many files concurrently don't oversubscribe; use it to speed up **single-file** runs
(~2.3x at `--num-threads 9`).
SVM class weights: `--cpos F` / `--cneg F` pin them, `--select-c` / `--no-select-c` force the grid
search on or off (on for `--canonical` and `--balanced`, off for `--fast`).

### Execution profiles
Presets so you don't have to memorize flag combinations. Pass one of `--fast` / `--balanced` /
`--canonical` (or `--profile <name>`); default is `--canonical`. **Explicit `--maxiter` /
`--subset-max-train` always override the preset** (e.g. `--fast --maxiter 8`).

| profile | expands to | when to use |
|---|---|---|
| `--fast` | `--subset-max-train 20000 --maxiter 5` | quick QA / test pipelines |
| `--balanced` | `--subset-max-train 40000 --maxiter 10` | fast with near-full yield |
| `--canonical` (default) | full defaults (maxiter 10, no subsetting) | max sensitivity for publication/production |

Measured across all 65 PXD032157 files at N=4 concurrency (percolator-rs):

| profile | wall | peak RAM | PSM q<0.01 | peptide q<0.01 |
|---|--:|--:|--:|--:|
| `--canonical` | 43.4 s | 0.79 GiB | 102 094 | 35 951 |
| `--balanced` | 42.9 s | 0.81 GiB | 102 269 | 35 867 |
| `--fast` | **21.5 s** | 0.84 GiB | 97 270 (−4.7%) | 33 588 (−6.6%) |

Note: because percolator-rs's canonical mode is already fast, `--balanced` lands essentially at
canonical speed *and* yield here; `--fast` roughly halves the time for a ~5–7 % yield cost (still far
better than the C++ fast config's −12 %/−15 %, since percolator-rs keeps full 3-fold CV — only the SVM
training-set size is capped).

## Benchmark vs C++ Percolator 3.09 (PXD032157, 65 files, 12-core Ryzen 5 5600G)

**Q1 — Can Rust hit sub-60 s at full fidelity (full iterations, no yield loss)?  YES.**

| implementation | settings | wall (65 files) | yield (PSM / peptide q<0.01) |
|---|---|---|---|
| C++ reference | default, sequential | 542 s | 103 038 / 35 852 (canonical) |
| C++ reference | default — floor on 12 cores | ~107 s (can't reach 60 s) | 103 038 / 35 852 |
| C++ reference | **fast flags** to reach 60 s | 49 s | 90 395 / 30 530 (**−12% / −15%**) |
| **percolator-rs** (optimized) | **default full fidelity, N=4** | **22.8 s** | **101 966 / 35 869 (−1.0% / +0.05%)** |
| **percolator-rs** (optimized) | default full fidelity, N=6 | 18.6 s | 101 966 / 35 869 |
| percolator-rs (pre-optimization) | default, N=4 | 41.2 s | 102 094 / 35 951 |

percolator-rs reaches sub-60 s **without** cutting iterations and **without** the 12–15 % yield loss the
C++ implementation needs to get there — aggregate identifications match the canonical run within ~1 %.

**Q2 — Peak RSS under identical concurrency (N=4).** percolator-rs peaks at **0.73 GiB** vs the C++
reference's much larger footprint (per process ~half: 263 MB peak for percolator-rs vs 377–525 MB for C++).
See `bench/RS_VS_CPP.md` for the full table.

## Advanced biological features

### Protein inference — picked-protein FDR (`--results-proteins` / `--decoy-results-proteins`)
Graph-based inference (`src/protein.rs`): builds the peptide↔protein graph, collapses proteins
sharing peptides into indistinguishable **groups** (union-find), scores each group by its best
peptide, and computes protein-group q-values / PEPs by **picked-protein FDR** (Savitski et al.
2015): each target group is paired with its decoy counterpart (matched on decoy-stripped member
names), only the higher-scoring of the pair is kept, and q-values are computed over that picked
list. This removes double-counting and is provably ≥ as sensitive as classic protein TDA — an
invariant enforced by a unit test (`picked_never_less_sensitive_than_classic`). The run log prints
both `q<0.01 (picked-FDR) vs (classic)` counts. Output columns:
`ProteinGroupId, q-value, posterior_error_prob, score, numPeptides, proteinIds`.

> **Honest caveat on this dataset:** PXD032157 is *metaproteomics* — a huge protein DB with ≈1
> peptide per protein, so the decoy:target *protein-group* ratio is ~1:1 and protein-level FDR is
> inherently near-unachievable at q<0.01 (single files yield 0–4 confident proteins; picked ≈
> classic because target/decoy protein sets barely overlap, leaving little to pick between). This
> is a property of the data, not the method — picked FDR shows its benefit on standard
> single-organism data where a protein has many peptides and its target/decoy versions compete.
> The score uses the best-peptide SVM discriminant (continuous; avoids −ln(PEP=0) saturation ties
> that otherwise scramble the ranking). Still not the full Bayesian Fido (α/β/γ marginalization).

### Retention-time features (`--rt-features`)
`src/rt.rs` predicts RT from peptide sequence (per-residue coefficient model), aligns it to
observed elution (ScanNr proxy) by least-squares on targets, and appends `rt_abs_error` /
`rt_sq_error` as two extra features. Correct PSMs elute near their predicted RT; random/decoy
matches deviate. **Measured effect:** +11 % / +5 % PSMs on two files, −9 % on another — the
framework works and *can* boost separation, but the coarse ScanNr proxy makes it inconsistent.
Plugging in true retention times + a stronger predictor (Elude/DeepLC-style) is the path to
reliable gains.

### Cross-run / file-group joint training (`--join file1.pin file2.pin …`)
Pools PSMs from several runs and trains **one shared model** (3-fold CV over the pool), then
scores each run — small files borrow statistical power from the group. Prints per-file yield.
**Measured:** 4 small files, 1400 → 1426 target PSMs at q<0.01 (+1.9 %), 3 of 4 improved
(the strongest run gives a little back — the expected shared-model regularization trade).

## Native optimizations
Making the Rust build faster *without* subset flags or yield loss:

| optimization | what | effect |
|---|---|---|
| **Explicit-Hessian Newton solver** | `dim` is small (~22), so form the 22×22 Hessian `H = I + 2·ΣCₖxₖxₖᵀ` once per step and Cholesky-solve `H d = −g`, instead of matrix-free CG (many Hessian-vector passes over ~100k samples) | main win: single big file 2.79 s → 2.02 s; full N=4 40.1 s → **22.8 s** |
| **mmap + fast-float parsing** | memory-map the file (`memmap2`) and parse floats over the raw byte buffer with `fast-float` (correctly-rounded → identical values) | parse 0.29 s → **0.16 s** (1.8×), yield bit-identical |
| **Vectorized `axpy`** | Hessian outer-product accumulation uses 4-wide `wide::f64x4` (exact, elementwise) | feeds the solver's hot loop |
| **`target-cpu=native`** | `.cargo/config.toml` | lets the backend use AVX2/AVX-512 |

Net: **~1.8× faster end-to-end** (40.1 → 22.8 s at N=4) at full fidelity.

Notes on what *didn't* help and why (measured, not assumed):
- **SIMD dot-product** over length-22 vectors gave no speedup (lane-load + horizontal-reduce overhead ≈ the scalar cost) and its reassociation perturbed borderline q<0.01 counts, so `dot` is kept as an exact sequential sum. The vectorization payoff is in `axpy` (exact), not `dot`.
- The Cholesky Newton step is *more* accurate than the old truncated CG, which shifted aggregate yield by 128 PSMs (102 094 → 101 966, −0.13%); still within ±1% and re-recorded as the new reference.

## CI / regression gates
Performance and accuracy are locked in so refactors can't silently drift:

- **`tests/regression.sh`** (hosted CI, no data needed) — runs percolator-rs on the committed
  deterministic fixture `tests/fixtures/sample.pin` and asserts q<0.01 yield within **±1%** of the
  recorded reference (`tests/expected.env`) plus wall-time and peak-RSS budgets. Exits non-zero on drift.
- **`bench/regression.sh`** (self-hosted / nightly, needs the PXD032157 data) — full 65-file `--canonical`
  run at N=4: asserts 65/65 valid, aggregate PSM & peptide q<0.01 within ±1% of recorded (102 094 / 35 951),
  wall < 60 s, peak RSS < 1.5 GiB.
- **`.github/workflows/ci.yml`** — on push/PR: `cargo build --release` → `cargo test` → `tests/regression.sh`.
  A manual `workflow_dispatch` job (self-hosted runner labelled `percolator-data`) also runs the C++ budget
  smoke test (`bench/fastrun.sh`, <60 s) and the percolator-rs full gate.

Recorded references live in `tests/expected.env`; update them intentionally when a change is meant to move
the numbers. percolator-rs is seed-deterministic, so fixture yields are exact run-to-run.

## Fidelity notes
Aggregate yield matches within ~1 %, but per-file PSM counts differ (percolator-rs's q-values are slightly
more permissive on some files, less on others) because the q-value/PEP uses a simpler pi0 than percolator's
Storey estimate and does not do per-spectrum best-PSM competition. Closing that gap (Storey pi0, per-scan
competition, C-parameter cross-validation) is the next fidelity step; it does not change the speed/memory story.
