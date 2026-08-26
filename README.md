# percolator-rs

[![CI](https://github.com/AndrejRumenovski/percolator-rs/actions/workflows/ci.yml/badge.svg)](https://github.com/AndrejRumenovski/percolator-rs/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

A from-scratch Rust reimplementation of the [Percolator](https://github.com/percolator/percolator)
semi-supervised PSM rescoring algorithm, built to benchmark against the reference C++ Percolator 3.09.

> **Status (2026-08-26):** an independent adversarial audit rejected the previous build, and its
> findings have been repaired. Exact target/decoy score ties are no longer resolved by row order — a
> 200-spectrum fixture went from **200 discoveries or 0 depending on which row came first** to the same
> result under all nine permutations. PEPs are now increments of the estimated false-discovery count
> rather than differences of monotonized q-values plus a constant. Protein groups are proteins with
> identical peptide evidence rather than connected components (largest group on a real bacterial search:
> **13 members → 2**), picked protein ties are drawn fairly, and the picked protein PEP — which was the
> best peptide's PEP under a protein-level name — now reads `NA`. `--select-c` is nested, so **all three
> CV modes pass adversarial leakage testing**.
>
> **What is still not true.** q-values and PEPs remain anti-conservative on signal-present entrapment
> data — 1.8× nominal at q<0.01, unchanged by this repair — so they must not be described as
> calibrated. That excess is now attributed rather than merely reported: about a third of it is present
> with **no rescoring at all**, and the rest tracks the semi-supervised training budget. Protein-level
> output has never been checked against protein ground truth. Full record:
> [`validation/SECOND_REPAIR.md`](validation/SECOND_REPAIR.md); the audit it answers is
> [`validation/INDEPENDENT_AUDIT.md`](validation/INDEPENDENT_AUDIT.md); the earlier failure and repair
> are preserved in [`validation/SCIENTIFIC_VALIDATION.md`](validation/SCIENTIFIC_VALIDATION.md),
> [`validation/IMPLEMENTATION_AUDIT.md`](validation/IMPLEMENTATION_AUDIT.md) and
> [`validation/REPAIR.md`](validation/REPAIR.md).

## Results

The headline claim is agreement, not advantage. Under matched post-processing on four independent
datasets — Tide, Sage, MSFragger and the upstream yeast fixture, five seeds each — percolator-rs and
C++ Percolator 3.09 report PSM counts at q<0.01 that differ by **+0.8, +15.2, −2.8 and +15.4**, with
Jaccard 0.92–0.996 and score rank correlation 0.95–0.999. Before the statistical repair those gaps
were +1,619, +828, +130 and +8.8. Almost the whole historical difference was defects on the
percolator-rs side or a post-processing mismatch in how the comparison was configured; it was never
evidence of better identification.

| Dataset | pre-repair Rust − C++ | repaired Rust − C++ (matched) | Jaccard | score Spearman |
|---|---:|---:|---:|---:|
| PXD007145 Tide | +1,619.0 | **+0.8** | 0.9930 | 0.9988 |
| PXD060954 Sage | +828.2 | **+15.2** | 0.9962 | 0.9955 |
| PXD020243 MSFragger | +130.4 | **−2.8** | 0.9235 | 0.9627 |
| Upstream yeast | +8.8 | **+15.4** | 0.9201 | 0.9498 |

**"Matched" needs a caveat that the earlier text did not give it.** The comparison ran C++ under
`--post-processing-tdc`, on the understanding that this makes the reference compete candidates within
a precursor as percolator-rs does. It does not: that flag governs how q-values are assigned when the
input came from *separate* searches, and on concatenated input 3.09 still reports every candidate.
Tide and Sage report one candidate per precursor, so their comparison is genuinely like-for-like;
MSFragger (mostly 3 per precursor) and yeast (mostly 2) are not, which is exactly where agreement
drops to 0.92. Aligning the lists by running percolator-rs with `--no-psm-competition` moves
MSFragger from 0.877 to 0.914 and yeast from 0.934 to 0.909 at seed 1 — it helps one and hurts the
other, so the residual 0.88–0.93 on those two datasets is a genuine model/solver difference, not a
bookkeeping artefact. Full method, per-seed distributions and the corrected comparison:
[`validation/SECOND_REPAIR.md`](validation/SECOND_REPAIR.md) §8.

On the four compact cases the pre-repair build was **7.3–14.6x faster** than the reference, using
**37–67%** of its peak RSS. Those timings have not been re-measured since the repair, which cost
roughly a third of the throughput on the large benchmark below. See the
[multi-dataset benchmark](bench/MULTI_DATASET.md); the dated command/result matrix for the full study
suite is in [`bench/REPRODUCTION.md`](bench/REPRODUCTION.md).

An experimental `--rescore-model mlp` path runs a deterministic one-hidden-layer neural model
through the same folds and FDR procedures as the default SVM. It does **not** improve aggregate
yield and has not been revalidated since the repair. The SVM remains the default. See the
[small-MLP benchmark](bench/DEEP_LEARNING.md).

The large-scale performance benchmark remains PXD032157 — 65 Comet `.pin` files, 2.3 GB — on a
12-core Ryzen 5 5600G. Inputs and machine are matched; post-processing is not, because the reference
default does not compete candidates within a spectrum and percolator-rs does. The identification
counts below are throughput-run outputs, not a scientific accuracy comparison.

| | C++ Percolator 3.09 | percolator-rs |
|---|---|---|
| **Wall clock** — 65 files, 4 processes | 376.2 s | **15.3 s** (24.5x faster) |
| **Peak memory** — 4 processes | 1.49 GiB | **0.75 GiB** |

Full iterations and full 3-fold cross-validation — no training-set reduction. A sequential
percolator-rs run takes 49.4 s, while the reference takes 376.2 s even with four-file concurrency.
The pre-repair build was faster (12.1 s at N=4, 36.3 s sequential); each of the three folds now fits
its own normalization, initial direction and design matrix, which is where the time went. Peak memory
is unchanged because folds run one at a time by default. This 23–25x figure is valid only for its
named host, dataset, command and comparator, and it is not evidence about correctness.

Identification counts from this pair of runs are deliberately **not** tabulated together. The
recorded reference run auto-detected these PINs as separate-search input and used mix-max, while
percolator-rs uses direct competition, so the two numbers answer different questions. percolator-rs
reports **106,823 target PSMs and 35,886 target peptides** at q<0.01 here; that is its own canonical
baseline (106,795 / 35,866 on the audited build, 107,046 / 37,469 under the rejected methodology), not
a comparison. Per-file yield on this dataset is extremely skewed — 142 to 9,938 target PSMs at q<0.01
across the 65 files — and PXD032157 was used to develop class weights, folds and score scaling, so
these counts are a reproducible development baseline and not a sensitivity estimate.

Results are bit-deterministic under a fixed seed and guarded by regression gates. Repeated
complete-null experiments find no false discovery at any threshold; the signal-present entrapment
experiment still rejects nominal 1% control for both implementations. See
[FDR calibration](#fdr-calibration).

Because exact score ties are now settled by a seeded coin rather than by row order, a dataset with
many exact ties will vary across seeds. That variation is real and was previously hidden inside the
input file's layout, not absent.

## Current algorithm

1. **PIN parse** — streaming tab-delimited reader. Metadata columns are the contiguous prefix after
   `Label` drawn from `ScanNr`, `ExpMass`, `CalcMass`, `rt`/`retentiontime`, `FileName`/`SpectraFile`,
   matched case-insensitively, exactly as the reference does; features start at the first
   unrecognized header. Malformed, missing or non-finite required values **stop the run** with the
   file, line, column name and offending text.
2. **Fold construction** — three folds built from spectrum groups keyed by `(source, ScanNr)`, so
   every candidate of one spectrum, target and decoy alike, trains and is scored together.
3. **Per-fold preprocessing** — z-score location and scale, and the initial direction (best single
   feature in either orientation by training-set yield), are both fitted **inside that fold's
   training partition**. Nothing about a held-out row reaches the model that scores it.
4. **Semi-supervised training** — iterate (default 10x): score → target-decoy q-values → pick
   confident targets (q<0.01) as positives, all decoys as negatives → retrain the fold-local learner.
   The default is a class-weighted **L2-regularized squared-hinge linear SVM** (primal
   truncated-Newton solver, the same L2-loss family as the reference L2-SVM-MFN). Class weights are
   **absolute** (`Cpos=1`, `Cneg=4`); pin them with `--cpos/--cneg`, or search per file with
   `--select-c`, which selects each outer fold's weights inside that fold's own training partition
   using an inner split. `--auto-model` additionally selects C, feature subset and Newton tolerance
   the same way. An experimental small MLP is available with `--rescore-model mlp`.
5. **Fold merging** — each fold's held-out scores are expressed in standard deviations above that
   fold's own training-decoy mean before pooling. Independently fitted models share no intercept or
   score unit, and the merged ranking is the only thing q-values see. (The reference instead anchors
   on the held-out selection boundary and median decoy; training decoys are used here to keep the
   transform inside the training partition.)
6. **Spectrum-level competition** — the best-scoring candidate of each precursor
   (`source`, `ScanNr`, `ExpMass`) is kept and the rest dropped, on the rescored values. It is on by
   default because the q-value estimator assumes each spectrum contributes at most one competition
   winner; `--no-psm-competition` reports every candidate, and its q-values are then not FDR
   estimates. **Exact ties are drawn with a fair coin** keyed on the precursor's own identity and the
   run seed, never by row order: a `k`-way tie holding `t` targets goes to a target with probability
   `t/k`, so permuting the input cannot change a single winner. Note that Percolator 3.09's
   `--post-processing-tdc` does *not* do this on concatenated input — it reports every candidate —
   so the two tools produce different lists on multi-candidate PINs.
7. **Target-decoy q-values** — `min(1, pi0 * (D+1) * lambda / max(1, T))` evaluated once per exact-score
   tie group and shared by every member, then the reverse cumulative minimum. `pi0 = 1`, the
   conservative choice for direct competition and what the reference uses for concatenated input;
   `lambda = p/(1-p)` with `p = --null-target-win-prob` (default 0.5, so 1.0 — use `1/(1+k)` for `k`
   decoys per target). The `+1` is the finite-sample safeguard; without it a leading run of targets
   reports exactly zero estimated FDP however thin the evidence.
8. **Posterior error probabilities** — increments of the same estimated false-discovery count the
   q-value scan produces: `F(s) = min(T(s), pi0 * lambda * (D(s) + 1))`, which is `T(s) * FDP(s)`.
   Each tie group's increment is shared among the targets it holds, decoy-only groups carry theirs
   forward, a group of `g` targets absorbs at most `g` more false discoveries, and PAVA then enforces
   monotonicity in score. The reported PEPs therefore **sum exactly to the reported estimate**, no
   constant is added anywhere, and a PEP of exactly zero is unreachable because the finite-sample
   safeguard decoy is spread across the leading run. Calibrated FDR does not imply calibrated PEPs and
   these are **not** calibrated — see [FDR calibration](#fdr-calibration).
9. **Peptide and protein levels** — best PSM per peptide, then the same estimators. Protein groups
   are proteins with *identical* observed peptide evidence; picked target/decoy competition breaks
   exact ties with the same fair coin. Picked inference estimates no protein-level posterior, so its
   `posterior_error_prob` column reads `NA`. No protein output has been validated against protein
   ground truth.

**Input contract.** A concatenated target-decoy search against a decoy database the same size as the
target database. Separate target/decoy searches (mix-max post-processing) are not supported.

## FDR calibration

**Complete null.** Ten exact-balance relabelings of each of three PXD032157 PINs: only original decoy
rows are kept and half are relabelled target, so features are exchangeable with respect to the label
and every accepted pseudo-target is false by construction. Under a complete null, FDP is 1 whenever
anything is accepted and 0 otherwise, so empirical FDR is the probability of any rejection.

| Threshold | pre-repair | **repaired** | C++ 3.09 |
|---|---:|---:|---:|
| q<0.001 through q<0.10 | 17/30 runs | **0/30 runs** | 0/29 runs |

Empirical complete-null FDR falls from 0.567 to **0.000** at every threshold. No seed was excluded.
Zero rejections in 30 replicates bounds the complete-null FDR above by about 0.12 at 95% confidence —
consistent with control at every threshold, and not a demonstration of control at 0.001. Runner:
[`validation/run_null.py`](validation/run_null.py).

**Signal-present entrapment.** Six deposited mzML runs re-searched against the native database plus
an equally sized foreign plant proteome; pure foreign-proteome assignments are known errors. Five
seeds, at reported q<0.01:

| Arm | accepted | entrapment-adjusted FDP |
|---|---:|---:|
| pre-repair percolator-rs (`d83a7ba`) | 19,542.8 | 2.698% |
| percolator-rs, competition off | 19,226.6 | 2.561% |
| C++ 3.09, no competition | 19,178.0 | 2.521% |
| percolator-rs, audited build (`b38c0db`) | 19,533.0 | 1.816% |
| **percolator-rs (current)** | **19,536.4** | **1.810%** |

**This fails nominal 1% control and the second repair did not change that.** Reported q<0.01
corresponds to an estimated 1.8% FDP here, at every seed, and every threshold from 0.001 to 0.10 is
above nominal (6.15×, 2.41×, 1.81×, 1.37×, 1.17×, 1.09×). Do not read q<0.01 as "1% error".

The excess is now attributed rather than only reported
([`validation/SECOND_REPAIR.md`](validation/SECOND_REPAIR.md) §7). Applying the identical entrapment
accounting to the **raw search score with no rescoring at all** gives 1.21× at nominal 1% and 1.01×
at nominal 5%, so about a third of the 1% excess exists before percolator-rs runs and the excess at
loose thresholds is entirely from rescoring. Raising the training budget walks it up monotonically —
1.21× with no training, 1.40× at one iteration, 1.84× at three, 1.79× at the default ten — which
points at the semi-supervised objective itself: a model fitted to separate the target database from
the decoy database promotes every target-database match, foreign entrapment targets included, and
promotes no decoys. That is a property of the Percolator method, and the reference shows it too.
Separately, the study's own plug-in foreign fraction moves the number by a factor of 1.5 either way.
Design and thresholds: [`bench/ENTRAPMENT.md`](bench/ENTRAPMENT.md).

**Posterior error probabilities.** Same entrapment data, all five seeds, ~98,700 targets each:
exact-zero PEPs **0**, known-false matches at PEP=0 **0**, weighted absolute calibration error
**0.0165–0.0204** (audited build, seed 1: 0.0185). **Rebuilding the estimator did not improve
calibration**, which is the expected result once the optimism is traced to the cumulative estimator it
differences: the per-bin observed/reported ratios (2.2, 2.3, 2.1, 1.7, 1.1) track the cumulative
q-value ratios on the same data. Every populated bin is optimistic, on every one of the six datasets
taken separately. **These are not calibrated posterior error probabilities.** What the repair did fix
is the estimator's semantics: the reported PEPs now sum exactly to the reported estimated
false-discovery count, and no tuned constant appears anywhere in the derivation.

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
For base PIN features, `--auto-model` instead performs per-outer-fold nested validation of SVM C,
class ratio, training-only feature subsets, and Newton tolerance. `--rt-features` refits its
retention-time alignment inside each outer training partition and is covered by the fold-isolation
tests. The base path avoids test-fold contamination but is 10.3x
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
| `--canonical` (default) | full defaults (maxiter 10, no subsetting) | the validated configuration; the only one the repair evidence covers |

Measured across all 65 PXD032157 files at N=4 concurrency (percolator-rs). The canonical row is the
median of three runs on the repaired build; the `--balanced` and `--fast` rows predate both the
2026-08-24 optimization campaign and the statistical repair, and have not been re-measured.

| profile | wall | peak RAM | PSM q<0.01 | peptide q<0.01 |
|---|--:|--:|--:|--:|
| `--canonical` (default), **current** | **15.35 s** | **0.75 GiB** | **106 823** | **35 886** |
| `--canonical` (default), audited build `b38c0db` | 16.1 s | 0.78 GiB | 106 795 | 35 866 |
| `--canonical`, pre-repair | 12.1 s | 0.76 GiB | 107 046 | 37 469 |
| `--balanced`, pre-repair | 20.4 s | 0.87 GiB | 106 817 | 37 526 |
| `--fast`, pre-repair | 14.6 s | 0.88 GiB | 105 237 | 36 772 |

Only the first row describes the current method. The balanced and fast rows are kept for their yield
tradeoffs under the old methodology; their counts are not comparable to the repaired canonical row and
their absolute wall times are not comparable to anything measured after the optimization campaign.
These measurements include writing the result files to local ext4. See
[`bench/OPTIMIZATION.md`](bench/OPTIMIZATION.md) for the complete candidate ledger, exact-output
checks, and final runtime profile.

## Benchmark vs C++ Percolator 3.09 (PXD032157, 65 files, 12-core Ryzen 5 5600G)

**Q1 — Can Rust hit sub-60 s with the complete default training workload (full iterations)? YES.**

Yield columns are each implementation's own reported-q count under its own post-processing. The
reference run auto-detected separate-search input and used mix-max; percolator-rs uses direct
competition. **The columns are not comparable to each other** and are listed only to show that the
Rust speed does not come from doing less work. The nearest thing to a matched comparison lives in
[`validation/SECOND_REPAIR.md`](validation/SECOND_REPAIR.md) §8.

| implementation | settings | wall (65 files) | own reported-q count (PSM / peptide q<0.01) |
|---|---|---|---|
| C++ reference | default, N=4 | 376.2 s | 103 038 / 35 852 |
| C++ reference | **fast flags**, N=5 | 59.4 s | 90 395 / 30 530 (**−12% / −15%** vs its own default) |
| **percolator-rs** | complete default workload, sequential | **49.4 s** | **106 823 / 35 886** |
| **percolator-rs** | **complete default workload, N=4** | **15.3 s** | **106 823 / 35 886** |
| percolator-rs (audited build `b38c0db`) | complete default workload, N=4 | 16.1 s | 106 795 / 35 866 |
| percolator-rs (pre-repair `d83a7ba`) | complete default workload, N=4 | 12.1 s | 107 046 / 37 469 |
| percolator-rs | `--select-c` per-file weight search, N=4 | 49.7 s | not re-measured since the repair |
| percolator-rs | `--auto-model` nested selection, N=4 | 206.4 s | not re-measured since the repair |

percolator-rs reaches sub-60 s **without** cutting iterations and without the 12–15% yield loss the
C++ implementation needs to get there. A single percolator-rs process (49.4 s) still finishes ahead of
the reference's observed N=4 run.

The yield columns are each implementation's own reported-q count under its own post-processing, and
those post-processings differ: Percolator 3.09's `--post-processing-tdc` does not compete candidates
within a precursor on concatenated input, while percolator-rs does. On this five-candidates-per-scan
dataset the two therefore report lists that differ by roughly 5x in length.

**Q2 — Peak RSS under identical concurrency (N=4).** percolator-rs peaks at **0.75 GiB** vs the C++
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

> **Not validated against protein ground truth.** The grouping and the target/decoy competition below
> are correct under their stated definitions and are covered by graph-level and adversarial tests, but
> no protein q-value in this repository has been checked against a known protein composition. Treat
> protein output as a ranked list with an untested error estimate.

**Grouping.** Two proteins form one group only when their *observed peptide sets are identical* — the
standard indistinguishability criterion. Grouping by connected components of the peptide-sharing
graph, which this module previously did, is a much coarser operation: it merges proteins the data
separates, and a chain of shared peptides can collapse a whole component. On a 105,560-PSM
single-organism bacterial search (`data/F_3.pin`, seed 1) the largest reported group went from **13
members to 2**, and groups with more than one member from 31 to 1.

No parsimony is applied. A peptide mapping to several distinguishable groups is evidence for each of
them, and a subset protein is reported separately from its superset because it *is* distinguishable.
That is a stated limitation, not an inference.

**Competition.** `--protein-inference picked` (default) pairs each target group with the decoy group
carrying the same decoy-stripped member names and keeps the higher-scoring of the pair. Exact ties
are drawn with the same fair coin used at the PSM level, keyed on the label-free pairing key: 200
exactly tied pairs previously gave 200 target wins at q<0.01 and now split about evenly with no
discovery, identically under any input order.

**Posterior.** Picked-protein FDR estimates a cumulative error rate over protein groups and no
posterior, so its `posterior_error_prob` column reads **`NA`**. It previously carried the best
peptide's PEP under a protein-level name. `--protein-inference bayesian` uses peptide PEPs in an
α/β/γ noisy-OR model and marginalizes protein presence with sum-product belief propagation (tree
components exact, cyclic components deterministic damped loopy BP); it does produce a protein-level
posterior and fills the column in. That posterior is **experimental and unvalidated**: its inputs are
PEPs that fail calibration, and its defaults (α=0.1, β=0.01, γ=0.5, peptide prior 0.1) come from a
PrEST study that did not establish nominal 1% protein FDR. See the complete model, CLI options and
limitations in [`bench/PROTEIN_INFERENCE.md`](bench/PROTEIN_INFERENCE.md). Output columns remain
`ProteinGroupId, q-value, posterior_error_prob, score, numPeptides, proteinIds`.

On `data/F_3.pin` at seed 1 the corrected grouping reports 1,456 target protein groups at q<0.01
under picked FDR and 1,417 under classic TDA (audited build: 1,409 and 1,368). **The difference is
recorded, not claimed as a benefit**: picking removes the loser of each paired competition, which is a
structural property of the method and not evidence that the resulting list is more accurate. The
optional `bench/protein_real.sh` gate reproduces the comparison when that uncommitted dataset is
present; the synthetic gate remains the portable hosted-CI check.

> **On PXD032157:** it is *metaproteomics* — a huge protein DB with ≈1 peptide per protein, so the
> decoy:target *protein-group* ratio is ~1:1 and protein-level FDR is inherently near-unachievable at
> q<0.01 (single files yield 0–4 confident proteins). Bayesian group posteriors behave very
> differently on this redundant graph, as the dedicated benchmark documents; neither list is validated
> there.

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
`ScanNr`, and the number that returned that exact `(ScanNr, modified peptide)` candidate. Both
agreement features are **label-free**: they were previously keyed on `Label` as well, which built a
whole-dataset training feature out of labels before the folds existed, and a unit test now sweeps
every relabeling of a small ensemble fixture and requires the feature matrix to be byte-identical.
Thus a model can learn from engine-specific error patterns and reproducible cross-engine agreement
while every PSM is still scored out of fold. Reports of the same exact candidate are kept in the same CV
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

Net: **~1.8× faster end-to-end** (40.1 → 22.8 s at N=4) on the complete default workload. These figures are from the
optimization work itself, measured *before* the class-weight change; the current default runs the same
workload in 19.4 s at N=4.

Notes on what *didn't* help and why (measured, not assumed):
- **SIMD dot-product** over length-22 vectors gave no speedup (lane-load + horizontal-reduce overhead ≈ the scalar cost) and its reassociation perturbed borderline q<0.01 counts, so `dot` is kept as an exact sequential sum. The vectorization payoff is in `axpy` (exact), not `dot`.
- The Cholesky Newton step converges more precisely than the old truncated CG, which shifted aggregate yield by 128 PSMs (102 094 → 101 966, −0.13%) against the baselines of the time; re-recorded then. Both figures predate the class-weight fix that took the current default to 107 046.

## CI / regression gates
Performance and recorded output yield are locked in so refactors cannot silently drift:

- **`tests/regression.sh`** (hosted CI, no data needed) — runs percolator-rs on the committed
  deterministic fixture `tests/fixtures/sample.pin` and asserts q<0.01 yield within **±1%** of the
  recorded reference (`tests/expected.env`) plus wall-time and peak-RSS budgets. Exits non-zero on drift.
- **`tests/model_regression.sh`** — exercises the optional MLP and requires serial and three-thread
  output to be byte-identical, with its development-fixture yield recorded explicitly.
- **`tests/selection_regression.sh`** — requires nested SVM choices and outputs to be byte-identical
  across serial and three-thread execution; a unit test separately proves outer-test mutations
  cannot change that fold's selected hyperparameters.
- **`bench/regression.sh`** (self-hosted / nightly, needs the PXD032157 data) — full 65-file `--canonical`
  run at N=4: asserts 65/65 valid, aggregate PSM & peptide q<0.01 within ±1% of recorded
  (**106 823 / 35 886**), wall < 45 s, peak RSS < 1.5 GiB.
- **`.github/workflows/ci.yml`** — on push/PR: `cargo build --release` → `cargo test` → `tests/regression.sh`.
  A manual `workflow_dispatch` job (self-hosted runner labelled `percolator-data`) also runs the C++ budget
  smoke test (`bench/fastrun.sh`, <60 s) and the percolator-rs full gate.

Recorded references live in `tests/expected.env`; update them intentionally when a change is meant to move
the numbers. percolator-rs is seed-deterministic, so fixture yields are exact run-to-run. They were all
re-recorded on the repaired method — the pre-repair values are kept in that file as history, not as an
acceptance criterion, because the methodology that produced them failed validation. Peptide gates on the
12,000-row fixture assert q<0.05 rather than q<0.01: with 38 target peptides above every decoy the
corrected estimator's best possible statement is 1/38 = 0.026, so a q<0.01 peptide gate would assert 0
and pass for any broken build.

These gates lock **implementation determinism**, not scientific validity. The evidence for the latter
is in [`validation/SECOND_REPAIR.md`](validation/SECOND_REPAIR.md), and it is partial. The suite is
itself checked by `validation/mutation_test.py`, which reintroduces twelve known scientific defects —
input-order-dependent ties, a missing finite-sample decoy, broken tie grouping, an added PEP constant,
a raised PEP floor, leaked normalization, a leaked initial direction, leaked C selection, a
label-keyed ensemble feature, target-favouring protein ties, connected-component protein grouping, and
a peptide PEP emitted as a protein PEP — and requires the suite to fail on each. All twelve are
caught.

## Fidelity notes

Under matched post-processing percolator-rs and C++ 3.09 now agree to within ±15 PSMs at q<0.01 on
four independent datasets (see [Results](#results)). The historical "+3.9% PSMs, +4.5% peptides"
against the auto-mode reference was never a like-for-like comparison — C++ auto-detected those PINs
as separate-search input and used mix-max — and the part of the gap that survived a forced-
concatenated rerun turned out to be defects in percolator-rs, not sensitivity. Identification counts
are not accuracy in either direction.

Two standard explanations for the original ~1% deficit were tested by measurement. One of them was
read the wrong way round:

- *Per-spectrum best-PSM competition* — the observation was right: the reference does not compete by
  default, and every scan in these inputs carries 5 Comet ranks, so C++ emits all 92,989 target rows
  unchanged. The conclusion drawn from it — that percolator-rs should not compete either — was wrong.
  Five candidates from one scan are not five independent hypotheses, and target-decoy competition
  assumes each spectrum contributes at most a competition winner. Not competing is why *both*
  implementations were anti-conservative on the entrapment study; competing moves both to the same,
  better place. percolator-rs now competes by default. **Correction to an earlier claim here:** the
  reference does *not* offer the same behaviour under `--post-processing-tdc`. That flag changes how
  q-values are assigned when the input came from separate searches; on concatenated input 3.09 still
  reports every candidate. Verified directly on PXD032157 file 0001 with the flag set: 92,989 C++
  target rows against 19,082 competed percolator-rs rows. Any "matched post-processing" comparison on
  a multi-candidate PIN therefore requires competing the reference output in a separate step.
- *Storey pi0* — confirmed: the reference logs `pi0 = 1` on this data. percolator-rs fixes `pi0` at 1
  for direct competition, for the same reason the reference does.

SVM class weighting (algorithm step 4) was the largest measured cause of the historical count
reversal. Switching to absolute `Cpos`/`Cneg` instead of the previous balance heuristic changed the
reported-yield gap and reduced its per-file spread. Without ground truth it is not an accuracy
result.

The per-file grid search (`--select-c`) does **not** pay for itself on this dataset: 2.5x the wall
time for a coin flip per file, because candidates are ranked by an abbreviated proxy run. It is now
**nested** — each outer fold selects its weights on an inner split of its own training partition —
so the selection bias that used to make its reported q-values optimistic is gone: flipping every
label in outer fold 0 changed all 200 of that fold's scores and moved its weights from 4:1 to 0.25:1
before, and now changes nothing. It remains off by default on cost grounds, not correctness grounds.

`--auto-model` performs nested selection inside each outer training partition. Its recorded
comparison against fixed defaults predates the repair and has not been re-measured; the numbers in
[`bench/AUTOMATIC_SELECTION.md`](bench/AUTOMATIC_SELECTION.md) describe the rejected methodology.

**Not revalidated on truth data.** `--rescore-model mlp`, `--auto-model`, `--join`, and the Bayesian
protein mode. `--select-c`, `--ensemble` and the default path now pass adversarial leakage testing
([`tests/cv_leakage.rs`](tests/cv_leakage.rs), `validation/adversarial_cv.py`), and `--rt-features`
is covered by the fold-isolation tests, but none of these has a truth-based accuracy study. They
inherit the corrected
estimators and the corrected fold isolation, which is necessary but not sufficient. Protein-level
output in particular still carries the PrEST calibration failures recorded in
[`bench/PROTEIN_CALIBRATION.md`](bench/PROTEIN_CALIBRATION.md).
