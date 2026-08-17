# percolator-rs

[![CI](https://github.com/AndrejRumenovski/percolator-rs/actions/workflows/ci.yml/badge.svg)](https://github.com/AndrejRumenovski/percolator-rs/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

A from-scratch Rust reimplementation of the [Percolator](https://github.com/percolator/percolator)
semi-supervised PSM rescoring algorithm, built to benchmark against the reference C++ Percolator 3.09.

## Results

Benchmarked against the reference **C++ Percolator 3.09** on PRIDE dataset PXD032157 — 65 Comet
`.pin` files, 2.3 GB — on a 12-core Ryzen 5 5600G. Identical inputs, identical settings, same machine.

| | C++ Percolator 3.09 | percolator-rs |
|---|---|---|
| **Wall clock** — 65 files, one process | 542 s | **54.9 s** (9.9x faster) |
| **Wall clock** — 65 files, 4 processes | 370 s | **19.4 s** (19x faster) |
| **PSMs** at reported q < 0.01 | 103 038 | **107 046** (+3.9%) |
| **Peptides** at reported q < 0.01 | 35 852 | **37 469** (+4.5%) |
| **Peak memory** — 4 processes | 1.56 GiB | **0.85 GiB** |

Full iterations and full 3-fold cross-validation — no training-set reduction, so the speedup is not
bought by doing less work. **One percolator-rs process finishes faster (54.9 s) than the C++
reference manages using all 12 cores** (~107 s floor); to get under 60 s the reference must enable
speed flags that cost it 12-15% of its identifications.

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
   targets (q<0.01) as positives, all decoys as negatives → retrain a class-weighted **L2-regularized
   squared-hinge linear SVM** (primal truncated-Newton / Newton-CG solver, the same L2-loss family as the
   reference L2-SVM-MFN).
   The class weights are **absolute** (`Cpos=1`, `Cneg=4`), not derived from the target/decoy class
   balance. This is the single largest accuracy win here: the obvious heuristic
   `Cpos = max(n_neg/n_pos, 1)` explodes (~300x) when confident targets are scarce and swamps the
   decoys — measured, it is the *worst* corner of the weight space. Pin the weights with
   `--cpos/--cneg`, or search them per file with `--select-c` (see [Fidelity notes](#fidelity-notes)).
5. **3-fold nested cross-validation** — each PSM scored by a model trained without its fold (no overfitting
   of the FDR estimate).
6. **Target-decoy q-values** + monotone (PAVA isotonic) PEP; PSM- and peptide-level output in the reference
   tab format.

## FDR calibration
`bench/null_calibration.sh` runs a pure-null experiment: keep only the decoy rows of a real `.pin`
and randomly relabel half as targets, so **every** reported identification is false by construction.
A calibrated method must report ~0 at `q<0.01`. Measured on PXD032157: **0–6 false IDs out of
22k–60k null targets** (~0.01% against a nominal 1%), and turning on the class-weight grid search
(`--select-c`) shifts that by at most ±2 — neither setting buys yield by loosening the FDR.

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
`--num-threads` (default 1) parallelizes the 3 CV folds within a single file — and the class-weight
grid too, when `--select-c` is on. Results are bit-identical at any thread count. It defaults to 1 so
harnesses that already run many files concurrently don't oversubscribe; use it for **single-file**
runs (1.74 s → 1.10 s at `--num-threads 3`, saturating there since there are only 3 folds; 3.48 s
→ 1.56 s at `--num-threads 9` with `--select-c`).
SVM class weights: `--cpos F` / `--cneg F` pin them; `--select-c` opts into the per-file grid
search, which is **off by default for every profile** (it costs ~3x wall time without beating the
fixed defaults on this dataset — see [Fidelity notes](#fidelity-notes)).

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
| `--canonical` (default) | 18.2 s | 0.73 GiB | 107 046 | 37 469 |
| `--balanced` | 18.3 s | 0.65 GiB | 106 817 (−0.2%) | 37 526 (+0.2%) |
| `--fast` | **12.5 s** | 0.69 GiB | 105 237 (−1.7%) | 36 772 (−1.9%) |

Note: because percolator-rs's canonical mode is already fast, `--balanced` lands essentially at
canonical speed *and* yield here; `--fast` cuts ~30 % of the time for a ~2 % yield cost (far better
than the C++ fast config's −12 %/−15 %, since percolator-rs keeps full 3-fold CV — only the SVM
training-set size is capped). Measured without writing result files, so these run ~1 s under the
CI gate's figure.

## Benchmark vs C++ Percolator 3.09 (PXD032157, 65 files, 12-core Ryzen 5 5600G)

**Q1 — Can Rust hit sub-60 s at full fidelity (full iterations, no yield loss)?  YES.**

| implementation | settings | wall (65 files) | yield (PSM / peptide q<0.01) |
|---|---|---|---|
| C++ reference | default, sequential | 542 s | 103 038 / 35 852 (canonical) |
| C++ reference | default — floor on 12 cores | ~107 s (can't reach 60 s) | 103 038 / 35 852 |
| C++ reference | **fast flags** to reach 60 s | 49 s | 90 395 / 30 530 (**−12% / −15%**) |
| **percolator-rs** | default full fidelity, sequential | **54.9 s** | **107 046 / 37 469 (+3.9% / +4.5%)** |
| **percolator-rs** | **default full fidelity, N=4** | **19.4 s** | **107 046 / 37 469** |
| percolator-rs | `--select-c` per-file weight search, N=4 | 57.2 s | 106 558 / 37 330 |

percolator-rs reaches sub-60 s **without** cutting iterations and **without** the 12–15 % yield loss the
C++ implementation needs to get there — and it identifies ~4 % *more* than the canonical reference run.
A single percolator-rs process (54.9 s) finishes ahead of the reference's 12-core floor (~107 s).

**Q2 — Peak RSS under identical concurrency (N=4).** percolator-rs peaks at **0.85 GiB** vs the C++
reference's **1.56 GiB** (per process roughly half: 263 MB vs 377–525 MB).
See `bench/RS_VS_CPP.md` for the full table.

## Advanced biological features

### Protein inference — picked-protein FDR (`--results-proteins` / `--decoy-results-proteins`)
Graph-based inference (`src/protein.rs`): builds the peptide↔protein graph, collapses proteins
sharing peptides into indistinguishable **groups** (union-find), scores each group by its best
peptide, and computes protein-group q-values / PEPs by **picked-protein FDR** (Savitski et al.
2015): each target group is paired with its decoy counterpart (matched on decoy-stripped member
names), only the higher-scoring of the pair is kept, and q-values are computed over that picked
list. This removes double-counting and is provably ≥ as sensitive as classic protein TDA — enforced
by a direct unit test, a parser-backed synthetic fixture test in `src/protein.rs`, and a synthetic
end-to-end CLI regression (`tests/protein_regression.sh`). The run log prints both
`q<0.01 (picked-FDR) vs (classic)` counts. Output columns:
`ProteinGroupId, q-value, posterior_error_prob, score, numPeptides, proteinIds`.

On a local 105,560-PSM single-organism bacterial search (`data/F_3.pin`), this produces **1,410
protein groups at q<0.01 vs 1,369 with classic TDA** (+41, +3.0%; seed 1). The optional
`bench/protein_real.sh` gate reproduces the comparison when that uncommitted dataset is present;
the synthetic gate remains the portable hosted-CI check.

> **Honest caveat on this dataset:** PXD032157 is *metaproteomics* — a huge protein DB with ≈1
> peptide per protein, so the decoy:target *protein-group* ratio is ~1:1 and protein-level FDR is
> inherently near-unachievable at q<0.01 (single files yield 0–4 confident proteins; picked ≈
> classic because target/decoy protein sets barely overlap, leaving little to pick between). This
> is a property of that dataset, not the method; the single-organism result above demonstrates the
> expected picked-FDR benefit when proteins have many peptides and target/decoy versions compete.
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
**not** pay for itself: ~3x the wall time, and a coin flip per file (better on 33, worse on 28,
marginally worse in aggregate) because candidates are ranked by an abbreviated proxy run. It remains
available for data where the fixed defaults may not transfer. Those defaults were themselves chosen
by measurement on this dataset — treat them as a well-tested starting point, not a universal
constant.
