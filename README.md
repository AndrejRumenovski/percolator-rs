# percolator-rs

[![CI](https://github.com/AndrejRumenovski/percolator-rs/actions/workflows/ci.yml/badge.svg)](https://github.com/AndrejRumenovski/percolator-rs/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

A from-scratch Rust reimplementation of the [Percolator](https://github.com/percolator/percolator)
semi-supervised PSM rescoring algorithm, built to benchmark against the reference C++ Percolator 3.09.

## Results

Benchmarked against **C++ Percolator 3.09** on five search configurations spanning mosquito,
human, bacterial, and yeast samples; Comet, Tide, MSFragger, and Sage inputs; a timsTOF Pro and
Orbitrap-family instruments; and search databases from 4,647 to 139,191 target proteins. Across the
four compact extension cases, percolator-rs is **7.3–14.6x faster**, using **37–67%** of C++ peak
RSS. Reported-q yield is not uniformly higher: the PSM delta ranges from **−1.8% to +12.0%**.
See the complete, reproducible [multi-dataset benchmark](bench/MULTI_DATASET.md).
The dated command/result matrix for the full study suite is in
[`bench/REPRODUCTION.md`](bench/REPRODUCTION.md).

An experimental `--rescore-model mlp` path runs a deterministic one-hidden-layer neural model
through the same folds and FDR procedures as the default SVM. It does **not** improve aggregate
yield: on PXD032157 it reports 1.42% fewer PSMs and 2.73% fewer peptides while taking 5.96x longer;
four independent extension cases are also slightly lower in aggregate. The SVM remains the default.
See the [small-MLP benchmark](bench/DEEP_LEARNING.md).

The large-scale headline remains PXD032157 — 65 Comet `.pin` files, 2.3 GB — on a 12-core Ryzen 5
5600G. Identical inputs, matched settings, same machine.

| | C++ Percolator 3.09 | percolator-rs |
|---|---|---|
| **Wall clock** — 65 files, 4 processes | 376.2 s | **12.1 s** (31.2x faster) |
| **PSMs** at reported q < 0.01 | 103 038 | **107 046** (+3.9%) |
| **Peptides** at reported q < 0.01 | 35 852 | **37 469** (+4.5%) |
| **Peak memory** — 4 processes | 1.49 GiB | **0.76 GiB** |

Full iterations and full 3-fold cross-validation — no training-set reduction, so the speedup is not
bought by doing less work. A sequential percolator-rs run takes 62.0 s, while the reference takes
376.2 s even with four-file concurrency. To get under 60 s here, the reference must enable speed
flags that cost it 12–15% of its identifications.

Results are bit-deterministic under a fixed seed and guarded by CI regression gates. A pure-null
control is strongly conservative, but a six-run foreign-proteome entrapment experiment finds that
signal-present q-values from **both** implementations are anti-conservative on this search; see
[FDR calibration](#fdr-calibration).

## Algorithm
Faithful to the Percolator method:
1. **PIN parse** — streaming tab-delimited reader (`ExpMass`/`CalcMass` excluded from features by name).
2. **Feature normalization** — per-feature z-score, bias column appended.
3. **Initial direction** — best single feature (either orientation) by target-decoy yield at q<0.01.
4. **Semi-supervised training** — iterate (default 10×): score → target-decoy q-values → pick confident
   targets (q<0.01) as positives, all decoys as negatives → retrain the fold-local learner. The
   default is a class-weighted **L2-regularized squared-hinge linear SVM** (primal truncated-Newton /
   Newton-CG solver, the same L2-loss family as the reference L2-SVM-MFN). An experimental small MLP
   is available with `--rescore-model mlp`.
   The class weights are **absolute** (`Cpos=1`, `Cneg=4`), not derived from the target/decoy class
   balance. This is the single largest accuracy win here: the obvious heuristic
   `Cpos = max(n_neg/n_pos, 1)` explodes (~300x) when confident targets are scarce and swamps the
   decoys — measured, it is the *worst* corner of the weight space. Pin the weights with
   `--cpos/--cneg`, or search them per file with `--select-c` (see [Fidelity notes](#fidelity-notes)).
   For leakage-free selection of C, class ratio, feature count, and solver tolerance, use the
   experimental nested `--auto-model` path.
5. **3-fold nested cross-validation** — each PSM scored by a model trained without its fold (no overfitting
   of the FDR estimate).
6. **Target-decoy q-values** + monotone (PAVA isotonic) PEP; PSM- and peptide-level output in the reference
   tab format.

## FDR calibration
`bench/null_calibration.sh` runs a pure-null experiment: keep only the decoy rows of a real `.pin`
and randomly relabel half as targets, so **every** reported identification is false by construction.
A calibrated method must report ~0 at `q<0.01`. Measured on PXD032157: **0–5 false IDs out of
8k–49k randomly relabelled null targets** (at most ~0.01% against a nominal 1%). Turning on the
class-weight grid search (`--select-c`) changes the three counts from 0/5/1 to 0/5/3 — neither
setting buys yield by loosening the FDR. Machine-readable results are in
[`bench/null-calibration-results.tsv`](bench/null-calibration-results.tsv).

The stronger signal-present check is `bench/entrapment/run.sh`: six deposited mzML runs are
re-searched against the native database plus an equally sized foreign plant proteome. At reported
q≤0.01, percolator-rs accepts **19,666 PSMs at an entrapment-estimated 2.78% FDP**; C++ 3.09 accepts
**19,126 at 2.62%**. Thus the +2.82% nominal-q yield lead on these runs is real as a count but is
*not* validated at an actual 1% FDR. Both implementations share the larger calibration failure,
with Rust slightly more anti-conservative at this cutoff. Full design, uncertainty intervals, and
all thresholds: [`bench/ENTRAPMENT.md`](bench/ENTRAPMENT.md).

## Build & run
```
cargo build --release
./target/release/percolator-rs --seed 1 \
  --results-psms t.psms --decoy-results-psms d.psms \
  --results-peptides t.pep --decoy-results-peptides d.pep \
  input.pin
```
Flags mirror the reference subset: `--seed`, `--maxiter`, `--subset-max-train`, `--num-threads`.
`--rescore-model svm|mlp` selects the fold-local learner; SVM is the default. The MLP architecture
and negative yield result are documented in [`bench/DEEP_LEARNING.md`](bench/DEEP_LEARNING.md).
`--num-threads` (default 1) parallelizes the 3 CV folds within a single file — and the class-weight
grid too, when `--select-c` is on. Results are bit-identical at any thread count. It defaults to 1 so
harnesses that already run many files concurrently don't oversubscribe; use it for **single-file**
runs. On the largest PXD032157 PIN, three-run medians are 2.05 s → 1.39 s at `--num-threads 3`,
and 3.93 s → 1.85 s at `--num-threads 9` with `--select-c`; outputs are byte-identical. See the
[`advanced-feature benchmark`](bench/ADVANCED_FEATURES.md).
Learner class weights: `--cpos F` / `--cneg F` pin them; `--select-c` opts into the per-file grid
search, which is **off by default for every profile** (it costs ~2.5x wall time without beating the
fixed defaults on this dataset — see [Fidelity notes](#fidelity-notes)).
`--auto-model` instead performs true per-outer-fold nested validation of SVM C, class ratio,
training-only feature subsets, and Newton tolerance. It avoids test-fold contamination but is 10.3x
slower, with 0.37% fewer PSMs and 0.45% more peptides on PXD032157; see the complete
[automatic-selection evaluation](bench/AUTOMATIC_SELECTION.md).

### Feature importance report

For a scientifically inspectable linear-SVM run, add `--feature-report features.tsv`. The report
contains each PIN feature's mean out-of-fold coefficient in raw units, its signed standardized
effect, across-fold coefficient standard deviation, raw feature/target-decoy correlation, and a
deterministic permutation importance: the change in accepted target PSMs at the configured FDR when
that feature is shuffled within each held-out fold. `selected_folds` shows how often a feature was
kept by `--auto-model` feature selection.

The report describes the actual three out-of-fold SVMs: normalization, feature masks, and score
scaling are reconstructed from each model's training partition, and no held-out PSM is used to fit a
coefficient. Permutation importance holds those fitted models fixed, so it measures model reliance,
not the gain from retraining without a feature. It is deterministic for a fixed `--seed` and is
currently intentionally unavailable for the nonlinear MLP.

### Execution profiles
Presets so you don't have to memorize flag combinations. Pass one of `--fast` / `--balanced` /
`--canonical` (or `--profile <name>`); default is `--canonical`. **Explicit `--maxiter` /
`--subset-max-train` always override the preset** (e.g. `--fast --maxiter 8`).

| profile | expands to | when to use |
|---|---|---|
| `--fast` | `--subset-max-train 20000 --maxiter 5` | quick QA / test pipelines |
| `--balanced` | `--subset-max-train 40000 --maxiter 10` | fast with near-full yield |
| `--canonical` (default) | full defaults (maxiter 10, no subsetting) | max sensitivity for publication/production |

Measured across all 65 PXD032157 files at N=4 concurrency (percolator-rs). The canonical row is the
median of three portable `x86-64-v3` runs from the 2026-08-24 optimization campaign; the noncanonical
profiles have not yet been rerun with these common-path optimizations.

| profile | wall | peak RAM | PSM q<0.01 | peptide q<0.01 |
|---|--:|--:|--:|--:|
| `--canonical` (default) | **12.1 s** | **0.76 GiB** | 107 046 | 37 469 |
| `--balanced` | 20.4 s | 0.87 GiB | 106 817 (−0.2%) | 37 526 (+0.2%) |
| `--fast` | **14.6 s** | 0.88 GiB | 105 237 (−1.7%) | 36 772 (−1.9%) |

The balanced and fast rows are earlier measurements and remain useful for their yield tradeoffs, but
their absolute wall times should not be compared with the newly optimized canonical row until they
are rerun. These measurements include writing the result files to local ext4. See
[`bench/OPTIMIZATION.md`](bench/OPTIMIZATION.md) for the complete candidate ledger, exact-output
checks, and final runtime profile.

## Benchmark vs C++ Percolator 3.09 (PXD032157, 65 files, 12-core Ryzen 5 5600G)

**Q1 — Can Rust hit sub-60 s at full fidelity (full iterations, no yield loss)?  YES.**

| implementation | settings | wall (65 files) | yield (PSM / peptide q<0.01) |
|---|---|---|---|
| C++ reference | default, N=4 | 376.2 s | 103 038 / 35 852 (canonical) |
| C++ reference | **fast flags**, N=5 | 59.4 s | 90 395 / 30 530 (**−12% / −15%**) |
| **percolator-rs** | default full fidelity, sequential | **62.0 s** | **107 046 / 37 469 (+3.9% / +4.5%)** |
| **percolator-rs** | **default full fidelity, N=4** | **12.1 s** | **107 046 / 37 469** |
| percolator-rs | `--select-c` per-file weight search, N=4 | 49.7 s | 106 558 / 37 330 |
| percolator-rs | `--auto-model` nested selection, N=4 | 206.4 s | 106 652 / 37 636 |

percolator-rs reaches sub-60 s **without** cutting iterations and **without** the 12–15 % yield loss the
C++ implementation needs to get there — and it identifies ~4% *more* than the canonical reference run.
A single percolator-rs process (62.0 s) finishes far ahead of the reference's observed N=4 run.

**Q2 — Peak RSS under identical concurrency (N=4).** percolator-rs peaks at **0.87 GiB** vs the C++
reference's **1.49 GiB**.
See `bench/RS_VS_CPP.md` for the full table.

## Advanced biological features

### Experimental neural rescoring

`--rescore-model mlp` uses an eight-unit tanh hidden layer with a linear residual initialized to
the same best feature as the SVM. Fold assignment, semi-supervised labels, out-of-fold scoring,
q-values, PEPs, and peptide rollup are shared with the SVM path. On the current 65-file and
multi-dataset benchmarks the added nonlinearity is slower and slightly reduces aggregate yield, so
it remains experimental and peptide-sequence embeddings are deferred. See the full
[SVM-versus-MLP evaluation](bench/DEEP_LEARNING.md).

### Protein inference — picked FDR and Bayesian marginalization

`--results-proteins` / `--decoy-results-proteins` now support two inference methods:

- `--protein-inference picked` (default) performs the existing best-peptide picked target-decoy
  competition.
- `--protein-inference bayesian` uses peptide PEPs in an α/β/γ noisy-OR model, clusters proteins
  with identical peptide connectivity, and marginalizes protein presence with sum-product belief
  propagation. Tree components are exact; cyclic components use deterministic damped loopy BP.

The Bayesian defaults are α=0.1, β=0.01, γ=0.5, and peptide prior 0.1. All are configurable. Its
q-values are cumulative expected posterior error, whereas picked q-values come from empirical
target-decoy competition. See the complete model, CLI options, limitations, and five-input
[picked-versus-Bayesian benchmark](bench/PROTEIN_INFERENCE.md). Output columns remain
`ProteinGroupId, q-value, posterior_error_prob, score, numPeptides, proteinIds`.

On a local 105,560-PSM single-organism bacterial search (`data/F_3.pin`), picked mode produces
**1,410 protein groups at q<0.01 vs 1,369 with classic TDA** (+41, +3.0%; seed 1). The optional
`bench/protein_real.sh` gate reproduces the comparison when that uncommitted dataset is present;
the synthetic gate remains the portable hosted-CI check.

> **Honest caveat on this dataset:** PXD032157 is *metaproteomics* — a huge protein DB with ≈1
> peptide per protein, so the decoy:target *protein-group* ratio is ~1:1 and protein-level FDR is
> inherently near-unachievable at q<0.01 (single files yield 0–4 confident proteins; picked ≈
> classic because target/decoy protein sets barely overlap, leaving little to pick between). The
> single-organism result above demonstrates the expected picked-FDR benefit when proteins have many
> peptides and target/decoy versions compete. Bayesian group posteriors behave very differently on
> this redundant graph, as the dedicated benchmark documents; neither list is validated there.
> The picked score uses the best-peptide SVM discriminant. Bayesian mode instead models all distinct
> peptide evidence and shared-peptide ambiguity; it does not assume the picked grouping or ranking.

Protein-level calibration is evaluated separately on the
[PXD008425 PrEST homology standard](bench/PROTEIN_CALIBRATION.md), which has explicit present/absent
protein pairs and 1,000 entrapment proteins. Replicate 1 is reserved for Bayesian α/β/γ selection;
replicates 2 and 3 remain held out for validation and final testing. Selection chooses α=0.1,
β=0.0001, γ=0.001 and sharply improves fixed-default calibration, but does **not** validate nominal
1% protein FDR across every vial; all methods also report false proteins in the held-out blank.

### Retention-time features (`--rt-features`)
`src/rt.rs` predicts RT from peptide sequence (per-residue coefficient model), aligns it to
observed elution (ScanNr proxy) by least-squares on targets, and appends `rt_abs_error` /
`rt_sq_error` as two extra features. Correct PSMs elute near their predicted RT; random/decoy
matches deviate. On the three lexicographically first PXD032157 files, the measured PSM effects are
**−4.81%, +3.98%, and +2.00%**. The framework can boost separation, but the coarse ScanNr proxy
makes it inconsistent. The pinned inputs and complete protocol are in the
[`advanced-feature benchmark`](bench/ADVANCED_FEATURES.md).
Plugging in true retention times + a stronger predictor (Elude/DeepLC-style) is the path to
reliable gains.

### Cross-run / file-group joint training (`--join file1.pin file2.pin …`)
Pools PSMs from several runs and trains **one shared model** (3-fold CV over the pool), then
scores each run — small files borrow statistical power from the group. Prints per-file yield.
On the four smallest PXD032157 PINs by byte size, the measured aggregate is **1,524 → 1,606** target
PSMs at q<0.01 (+5.38%); three of four improve and one gives back 16 PSMs. Input hashes and per-file
results are in the [`advanced-feature benchmark`](bench/ADVANCED_FEATURES.md).

### Experimental search-engine ensemble (`--ensemble ENGINE=PIN …`)

For results from multiple search engines on the **same raw run**, use one named PIN per engine:

```bash
./target/release/percolator-rs --seed 1 --ensemble \
  comet=comet.pin msfragger=msfragger.pin tide=tide.pin \
  --results-psms ensemble.target.psms --decoy-results-psms ensemble.decoy.psms
```

This is not a cosmetic concatenation. Engine feature spaces are namespaced, so `xcorr` and an
MSFragger score never get treated as interchangeable merely because their columns happen to share a
name. The pooled model also receives an engine indicator, the number of engines that returned each
`ScanNr`, and the number that returned that exact `(ScanNr, label, modified peptide)` PSM. Thus a
model can learn from engine-specific error patterns and reproducible cross-engine agreement while
every PSM is still scored out of fold. Reports of the same exact candidate are kept in the same CV
fold; only the best-scoring report is emitted, and q-values are recalculated over those unique
candidates, so agreement cannot manufacture extra discoveries.

Inputs must use compatible target/decoy databases and refer to the same run: the agreement features
key spectra by `ScanNr`. Do not mix separate raw files in one ensemble invocation, because scan
numbers can collide. Output PSM IDs are prefixed with the engine name to remain unique. Protein
inference is intentionally unavailable in ensemble mode until duplicate engine evidence has a
separately calibrated protein-level treatment. As with every new scoring regime, assess target-decoy
and entrapment calibration on an independent dataset before interpreting increased q-value yield as
an accuracy gain.

## Native optimizations
Making the Rust build faster *without* subset flags or yield loss:

| optimization | what | effect |
|---|---|---|
| **Explicit-Hessian Newton solver** | `dim` is small (~22), so form the 22×22 Hessian `H = I + 2·ΣCₖxₖxₖᵀ` once per step and Cholesky-solve `H d = −g`, instead of matrix-free CG (many Hessian-vector passes over ~100k samples) | main win: single big file 2.79 s → 2.02 s; full N=4 40.1 s → **22.8 s** |
| **mmap + fast-float parsing** | memory-map the file (`memmap2`) and parse floats over the raw byte buffer with `fast-float` (correctly-rounded → identical values) | parse 0.29 s → **0.16 s** (1.8×), yield bit-identical |
| **Vectorized `axpy`** | Hessian outer-product accumulation uses 4-wide `wide::f64x4` (exact, elementwise) | feeds the solver's hot loop |
| **`target-cpu`** | `.cargo/config.toml` | `x86-64-v3` baseline (AVX2/FMA) so binaries stay portable; `RUSTFLAGS="-C target-cpu=native"` for benchmark builds |

Net: **~1.8× faster end-to-end** (40.1 → 22.8 s at N=4) at full fidelity. These figures are from the
optimization work itself, measured *before* the class-weight change; the current default runs the same
workload in 19.4 s at N=4.

Notes on what *didn't* help and why (measured, not assumed):
- **SIMD dot-product** over length-22 vectors gave no speedup (lane-load + horizontal-reduce overhead ≈ the scalar cost) and its reassociation perturbed borderline q<0.01 counts, so `dot` is kept as an exact sequential sum. The vectorization payoff is in `axpy` (exact), not `dot`.
- The Cholesky Newton step is *more* accurate than the old truncated CG, which shifted aggregate yield by 128 PSMs (102 094 → 101 966, −0.13%) against the baselines of the time; re-recorded then. Both figures predate the class-weight fix that took the current default to 107 046.

## CI / regression gates
Performance and accuracy are locked in so refactors can't silently drift:

- **`tests/regression.sh`** (hosted CI, no data needed) — runs percolator-rs on the committed
  deterministic fixture `tests/fixtures/sample.pin` and asserts q<0.01 yield within **±1%** of the
  recorded reference (`tests/expected.env`) plus wall-time and peak-RSS budgets. Exits non-zero on drift.
- **`tests/model_regression.sh`** — exercises the optional MLP and requires serial and three-thread
  output to be byte-identical, with its development-fixture yield recorded explicitly.
- **`tests/selection_regression.sh`** — requires nested SVM choices and outputs to be byte-identical
  across serial and three-thread execution; a unit test separately proves outer-test mutations
  cannot change that fold's selected hyperparameters.
- **`bench/regression.sh`** (self-hosted / nightly, needs the PXD032157 data) — full 65-file `--canonical`
  run at N=4: asserts 65/65 valid, aggregate PSM & peptide q<0.01 within ±1% of recorded (107 046 / 37 469),
  wall < 45 s, peak RSS < 1.5 GiB.
- **`.github/workflows/ci.yml`** — on push/PR: `cargo build --release` → `cargo test` → `tests/regression.sh`.
  A manual `workflow_dispatch` job (self-hosted runner labelled `percolator-data`) also runs the C++ budget
  smoke test (`bench/fastrun.sh`, <60 s) and the percolator-rs full gate.

Recorded references live in `tests/expected.env`; update them intentionally when a change is meant to move
the numbers. percolator-rs is seed-deterministic, so fixture yields are exact run-to-run.

## Fidelity notes
percolator-rs identifies **more** than the C++ reference at the same reported threshold (+3.9 %
PSMs, +4.5 % peptides). The pure-null control rules out hallucinated signal on pure noise, but the
signal-present entrapment experiment above shows that it does not establish exact q-value
calibration: both implementations exceed the nominal FDR on the entrapment search.

Two standard explanations for the original ~1 % deficit were tested and **disproved by measurement**:

- *Per-spectrum best-PSM competition* — the reference does **not** do it. Every scan in these inputs
  carries 5 Comet ranks, and C++ emits all 92 989 target rows unchanged.
- *Storey pi0* — the reference logs `pi0 = 1` on this data (peptide level `0.999168`), landing in the
  same place as the simple decoys/targets estimate.

The actual cause was the SVM class weighting (Algorithm step 4). Switching to absolute `Cpos`/`Cneg`
instead of scaling by the target/decoy balance closed the gap and reversed it, and collapsed the
per-file spread: files below the reference went 36 → 11, worst single-file deficit −507 → −82. The
original "1 % gap" was never a bias — it was two-sided variance that happened to nearly cancel.

The per-file grid search (`--select-c`) was built on top of that fix and, on this dataset, does
**not** pay for itself: 2.5x the wall time, and a coin flip per file (better on 32, worse on 28,
tied on 5, and marginally worse in aggregate) because candidates are ranked by an abbreviated proxy
run. It remains
available for data where the fixed defaults may not transfer. Those defaults were themselves chosen
by measurement on this dataset — treat them as a well-tested starting point, not a universal
constant.

The newer `--auto-model` path removes the legacy selector's evaluation leakage: normalization,
initialization, feature ranking, hyperparameter choice, and fitting all occur inside each outer
training partition, and fold-specific margins are standardized from training decoys before pooling.
It finishes 394 PSMs below but 167 peptides above fixed defaults while costing 10.3x more, then loses
both metrics on independent extension cases. All 195 outer models keep the existing solver tolerance
and 194 keep all features, so the added flexibility is not justified as a default here. Full design
and held-out results: [`bench/AUTOMATIC_SELECTION.md`](bench/AUTOMATIC_SELECTION.md).
