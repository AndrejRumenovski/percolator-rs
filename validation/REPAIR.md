# Statistical repair of percolator-rs

Repair date: 2026-08-25
Repaired implementation: commit `1348b0f` (frozen for every experiment below)
Pre-repair implementation: commit `d83a7ba`, audited in
[`IMPLEMENTATION_AUDIT.md`](IMPLEMENTATION_AUDIT.md) and rejected in
[`SCIENTIFIC_VALIDATION.md`](SCIENTIFIC_VALIDATION.md)
Reference source: Percolator `rel-3-09`, commit `7238ac50`
Reference executable: Percolator 3.09.0, build date 2026-05-21
Threshold convention: strict `q < threshold`

This document records what was verified, what was changed, and what the predeclared
experiments returned afterwards. The failed pre-repair experiments are preserved unchanged; nothing
in this repair rewrites or deletes them.

## 1. Verified root causes

Each defect was traced to a specific operation in the pre-repair source and checked against
Percolator 3.09 before anything was changed. The audit's proposed explanations were treated as
hypotheses; two of them turned out to be incomplete, and one entirely new cause was found that the
audit had not identified.

### 1.1 No finite-sample decoy safeguard — `stats::qvalues_into`

Pre-repair, the scan started at zero decoys, so any leading run of targets received an estimated
FDP of exactly zero:

```rust
let mut targets = 0.0f64;
let mut decoys  = 0.0f64;          // <- starts at zero
...
let fdr = if targets > 0.0 { (pi0 * decoys) / targets } else { 1.0 };
```

The reference starts its decoy count at one and drops it only for the two training heuristics
(`PosteriorEstimator::getQValues`, `n_z_ge_w = 1` unless `skipDecoysPlusOne`).

**Verification on the preserved failed runs.** Across the 30 complete-null replicates, the
pre-repair build accepted 41 pseudo-targets in 17 runs. Every one of the 41 had a reported q-value
of exactly `0.0`. In 26 of the 30 replicates, whether the run made a discovery is predicted exactly
by whether the top-ranked row happened to be a pseudo-target — which is what a missing safeguard
predicts under an exactly balanced null, and which also explains why the count is near 15 of 30.
The remaining four replicates reached a small `D/T` ratio deeper in the list instead.

### 1.2 No tie grouping — `stats::qvalues_into`

The FDP was evaluated once per *row*, in `sort_unstable` order, so equal-scoring rows could receive
different q-values and the result depended on an arbitrary permutation. The reference evaluates once
per tie group (`myPair->first != (myPair+1)->first`) and assigns the value to every member. A
q-value must be a function of the score threshold alone.

### 1.3 Invalid count-ratio `pi0` — `stats::estimate_pi0`

`estimate_pi0` returned `min(1, D_total / T_total)` and was applied as a `pi0` multiplier, giving
approximately `(D/T)^2` at the end of a list. It is neither of the two things it resembles:

- It is not Storey's `pi0`. The reference estimates that from TDC p-values
  (`PosteriorEstimator::estimatePi0`) and, decisively, only *uses* it for mix-max:
  `Caller.cpp` sets `allScores.setUsePi0(useMixMax_)`, so a concatenated search runs at `pi0 = 1`.
- It is not an opportunity correction either, and it points the wrong way. The reference's
  correction is `decoyFactor = p / (1 - p)` with `p = nullTargetWinProb` (default 0.5, hence 1.0),
  which for a target-heavy list would be greater than one, not less.

Note that the audit described the reference as applying `targetDecoySizeRatio_ = T/D` in the
q-value. That is not correct: `targetDecoySizeRatio_` appears only in `CrossValidation::
initializeGridSearch`, where it scales the SVM class-weight grid. The q-value carries no size ratio
at all.

### 1.4 Unsupported PEP transform — `stats::peps_from_order`

Pre-repair PEPs were `min(1, 2 * pi0 * p)` over a PAVA fit `p` of the pooled decoy indicator. A
fitted local decoy fraction `p` corresponds to local false-target odds `p / (1 - p)`; `2p` is a
linearisation valid only near `p = 0.5`, and the `pi0` factor is the invalid count ratio of §1.3.
More importantly the fit is unregularised, so a leading level set containing no decoy is fitted at
exactly zero. That is the mechanism behind 9,155 PSMs printed at PEP = 0, 87 of them known-false
foreign-proteome matches.

The reference derives PEPs from target q-values instead (`Scores::calcPep`,
`InferPEP::q_to_pep`), using the identity of Käll et al. (2008): `q_k` is the mean PEP over the top
`k` targets, so `raw PEP_k = k*q_k - (k-1)*q_{k-1}`. It then adds a Bayesian pseudocount — half a
false discovery spread over the list — before a monotone fit, with the comment that this "prevents
exact-zero PEPs in the leading tail".

### 1.5 Transductive normalization — `percolator::build_matrix`

`build_matrix(ds)` fitted z-score location and scale on every row before any fold was trained.

### 1.6 Supervised initial direction on all labels — `percolator::initial_direction`

One direction was chosen by scanning every feature against **every label**, including all held-out
labels, and was then handed to all three fold models. The reference selects a direction inside each
outer training set (`SanityCheck::calcInitDirection` on `trainset[fold]`).

### 1.7 Raw fold margins pooled — `percolator::cv_scores`

Three independently fitted models have no shared intercept or score unit; their raw margins were
concatenated and ranked together. The merged ranking is the only thing q-values see. The reference
rescales each fold before merging (`Scores::normalizeScores`, boundary to 0 and median decoy to −1).

### 1.8 Rows split across folds within one spectrum — `percolator::assign_dataset_folds`

Folds were assigned per row outside ensemble mode, so several candidates of one spectrum — including
a target and the decoy it competed against — could land in different folds. The reference splits by
spectrum (`Scores::createXvalSetsBySpectrum`).

### 1.9 Parser coercion and non-finite acceptance — `pin::parse`

`fast_float::parse(...).unwrap_or(0.0)` turned any unreadable feature into zero, `atoi` skipped
non-digits, short rows were dropped silently, and NaN/infinity passed through into normalization,
the SVM and the sort order.

### 1.10 Wrong feature block, and a column of zeros — `pin::parse` *(not in the audit)*

Running the strict parser over the evaluated datasets found that Sage PINs carry a `FileName`
column, and that every row of it had been parsed as `0.0` and handed to the SVM as a feature.
percolator-rs took every column between `ScanNr` and `Peptide` as a feature except `ExpMass` and
`CalcMass` (matched case-sensitively). The reference reads a contiguous prefix of metadata columns
after `Label` — scannr, expmass, calcmass, rt/retentiontime, filename/spectrafile, matched
case-insensitively — and starts features at the first unrecognized header
(`SetHandler::getOptionalFields`). percolator-rs was therefore also training on a raw retention time
in minutes as though it were a search score, on every MSFragger- and Sage-style input.

### 1.11 Top-N input breaks the competition assumption *(not in the audit; the largest single cause)*

The entrapment study found *both* implementations anti-conservative at nominal 1% — 2.70% for Rust
and 2.52% for C++ — which no defect on one side can explain. Counting candidates per spectrum in the
study's own inputs gives the reason:

| input | rows | spectra | candidates per spectrum |
|---|---:|---:|---:|
| entrapment `comet.pin` | 46,804 | 9,380 | 5 |
| PXD032157 `comet.pin` | 34,680 | 6,937 | 5 |
| Tide `hogrebe.pin` | 55,398 | 55,398 | 1 |

Five candidates from one scan are not five independent hypotheses. Target-decoy competition assumes
each spectrum contributes at most the winner of a competition against the decoy database; that is
the assumption which makes one observed decoy stand for one incorrect target. Neither implementation
performed the competition: percolator-rs had no such step, and C++ was invoked with
`--search-input concatenated`, which in `Caller.cpp` clears `targetDecoyCompetition_`.

### 1.12 Supervised retention-time alignment — `rt::augment`

`--rt-features` fitted `obs ~ a*pred + b` on all target PSMs before folds existed, so every held-out
label contributed to a feature the fold model then trained on.

## 2. Corrections

Every correction is a separate commit with its own regression tests. Yield and runtime were treated
as outcomes, never as objectives.

| # | Commit | Correction |
|---|---|---|
| 1 | `e6bda81` | Tie-grouped TDC q-values with the finite-sample safeguard and a declared opportunity ratio; count-ratio `pi0` removed; PEPs derived from q-values with a Bayesian pseudocount |
| 2 | `d640fd6` | Normalization and initial direction fitted inside each outer training partition |
| 3 | `9a3b1e7` | Fold scores standardized to their own training-decoy null before merging |
| 4 | `49e5be2` | Folds built from spectrum groups |
| 5 | `8a1fd52` | Parser fails closed; metadata columns follow the reference contract |
| 6 | `50ea639` | Spectrum-level target-decoy competition on the rescored values, on by default |
| 7 | `63be826` | Retention-time alignment refitted inside each fold, per input file |
| 8 | `46db52f` | Declared opportunity ratio threaded through the nested paths |
| 9 | `1348b0f` | Regression gates re-recorded on the corrected method |

### 2.1 q-values

Before, for scores sorted high to low, reverse-cumulative-minimum of

```
min(1, min(1, D_total/T_total) * D_running / T_running)
```

After, evaluated once per exact-score tie group and shared by every member of it:

```
min(1, pi0 * (D + 1) * lambda / max(1, T)),    lambda = p / (1 - p)
```

with `pi0 = 1` and `p = --null-target-win-prob` (default 0.5, so `lambda = 1`). The `+1` is dropped
only by `Tdc::training`, which drives the semi-supervised positive selection and the initial
direction search — both selection heuristics that make no error-rate claim, and both places where
the reference likewise sets `skipDecoysPlusOne = true` "because the decoys+1 in the FDR estimates is
too restrictive for small datasets".

`pi0` is fixed at 1 because that is the conservative and correct choice for direct target-decoy
competition, and because it is what the reference itself uses for concatenated input. Unequal
target/decoy opportunity is handled explicitly instead, by declaring it: `--null-target-win-prob`
is `1/(1+k)` for `k` decoys per target. Separate target/decoy searches (mix-max) are documented as
unsupported rather than silently mis-estimated.

Sorting uses a deterministic total order with NaN placed last, instead of a comparator that returned
`Equal` for NaN and let it take an arbitrary rank.

### 2.2 PEPs

Derived from the corrected q-values through the Käll et al. (2008) identity, over targets ranked
best-first:

```
raw PEP_k = k*q_k - (k-1)*q_{k-1}
PEP       = PAVA( raw PEP + 0.5/n_targets )   clipped to (0, 1]
```

The pseudocount is prior mass — half a false discovery spread uniformly over the target list — not a
clamp applied to a zero after the fact. Combined with the q-value safeguard it makes an exact-zero
PEP unreachable: `q_1 >= lambda/T_total > 0`, so every raw PEP is already positive before the prior
is added. Decoys carry no error-rate claim and are interpolated in q between the targets that
bracket them, holding the end values outside that range.

### 2.3 Cross-validation isolation

`FoldSetup` fits normalization from its training rows alone, searches the initial direction over its
training rows alone, and refits the retention-time alignment inside its training rows (per input
file, since `ScanNr` indexes a single run's acquisition). Held-out scores are expressed in standard
deviations above that fold's own training-decoy mean before merging.

The reference instead anchors each fold on its held-out selection boundary and held-out median
decoy. Training decoys are used here because they keep the whole transform inside the training
partition; both anchorings are increasing affine maps, so neither can reorder a fold internally, and
the choice affects only how folds interleave.

Folds are built from spectrum groups keyed by `(source, scan)` — by scan alone for ensemble input,
where the same spectrum is deliberately reported by several engines. Groups are dealt to the
currently smallest fold rather than round-robin, since they hold different numbers of candidates.

`--select-c` remains selection-biased and non-nested: it still scores candidates on the same
out-of-fold predictions it later reports. That is now stated in the code and the documentation
rather than left implicit.

### 2.4 Input contract

percolator-rs now states what it requires and enforces what it can:

- a concatenated target-decoy search against a decoy database the same size as the target database;
- `--null-target-win-prob P` declares the null probability that an incorrect target outranks its
  paired decoy;
- spectrum-level competition is performed on the rescored values before PSM statistics, so a PIN
  reporting several candidates per spectrum is handled rather than silently mis-estimated;
- `--no-psm-competition` restores the every-candidate list, whose q-values are then not FDR
  estimates;
- malformed, missing or non-finite required values stop the run with a located diagnostic;
- separate target/decoy (mix-max) input is unsupported.

## 3. Tests added

The tests encode methodological invariants rather than recorded outputs. Where a test is meant to
catch a specific defect, the pre-repair behaviour was restored temporarily and the test confirmed to
fail; those checks are noted below.

**Target-decoy q-values** (`src/stats.rs`) — leading-run positivity; zero-decoy positivity; bounded
in [0,1]; monotone in score; the documented closed form; tie members share a q-value; invariance to
row permutation; all-identical scores collapse to one value; target/decoy imbalance is not rescaled
by a count ratio; the opportunity ratio scales the estimate; the safeguard is dropped only by the
training estimator; empty and single-row inputs; repeated calls are deterministic; the fast
count/mask paths agree with materialized q-values including on NaN and `-0.0` inputs.

**PEPs** (`src/stats.rs`) — strictly positive and bounded; a 500-target leading run never reaches
zero; monotone along worsening score; the mean PEP over the top `k` targets tracks `q_k`; tie
members share a PEP; invariance to row permutation; decoy-only and empty inputs.

**Cross-validation isolation** (`src/percolator.rs`) — corrupting every held-out feature and
flipping every held-out label leaves the fold's initial direction, training-row normalization and
fitted weights bit-identical, both with and without `--rt-features`; folds do not share one design
matrix; every row is scored by a model whose training partition excluded it; standardization
preserves within-fold order and really does centre and scale that fold's training decoys.
*Verified to fail against the pre-repair all-rows fit.*

**Fold construction** (`src/percolator.rs`) — every candidate of a spectrum shares a fold; fold
sizes stay within one spectrum of each other; joined files reusing scan numbers are not fused.

**Competition** (`src/main.rs`) — one winner per precursor; distinct precursors of one scan compete
separately; joined files compete separately; ties resolve by input order and are stable across
repeats; already-competed input passes through unchanged; a PIN without `ExpMass` competes per scan.

**Parser** (`src/pin.rs`) — malformed, missing and non-finite features rejected; the diagnostic
names the line, column and offending text; malformed labels and scan numbers rejected; short rows
rejected rather than skipped; a missing required column rejected; the metadata prefix stops at the
first feature; a metadata name appearing after the feature block stays a feature; a PIN with no
feature columns rejected; a well-formed PIN still parses.

**Retention time** (`src/rt.rs`) — the alignment depends only on the rows it is fitted on; each
input file gets its own alignment; residual columns are finite, non-negative and consistent.
*Verified to fail against an all-rows fit.*

Existing tests were not weakened. Three carried constants produced by the rejected estimator and
were re-recorded, with the reason stated at the assertion: the picked-protein fixture (classic FDR
reported 81 groups at q<0.01 only through a decoy-free q=0 run; 1/81 = 0.0123 does not clear 1%),
and the peptide-level gates on the 12,000-row fixture, which move from q<0.01 to q<0.05 because the
corrected estimator's best possible statement there is 1/38 = 0.026 and a q<0.01 gate would assert 0
and pass for any broken build.

## 4. Complete-null revalidation

The predeclared experiment was rerun unchanged: the same three PXD032157 PINs, the same ten
relabeling seeds (1001–1010), the same exact-balance construction, the same runner
([`run_null.py`](run_null.py)), the same six thresholds. The constructed null PINs are byte-identical
to the pre-repair study — every `(input, relabel_seed)` pair has the same SHA-256 — so this is a
paired comparison, not a re-derivation. The repaired build ran from a clean worktree at `1348b0f`.

Under a complete null every accepted pseudo-target is false by construction, so FDP is 1 whenever
there is any rejection and 0 otherwise; empirical FDR is the probability of any rejection across
replicates.

| Threshold | pre-repair Rust | **repaired Rust** | C++ 3.09 |
|---|---:|---:|---:|
| q<0.001 | 17/30 | **0/30** | 0/29 |
| q<0.005 | 17/30 | **0/30** | 0/29 |
| q<0.01  | 17/30 | **0/30** | 0/29 |
| q<0.02  | 17/30 | **0/30** | 0/29 |
| q<0.05  | 17/30 | **0/30** | 0/29 |
| q<0.10  | 17/30 | **0/30** | 0/29 |

Empirical complete-null FDR falls from **0.567** to **0.000** at every threshold. No seed was
excluded; all 30 repaired replicates completed. C++ reproduces its previous result exactly, including
the one replicate that fails closed with no initial direction, so 29 of its 30 are usable.

Zero rejections in 30 replicates bounds the complete-null FDR above by roughly 0.12 at 95%
confidence (Wilson). That is consistent with control at every tested threshold and does not
demonstrate control at 0.001; a null experiment can only reject calibration, never confirm it at a
resolution finer than its replicate count.

Manifest: `$HOME/percolator_rs_out/scientific-validation/null-repaired-20260825/manifest.json`.

## 5. Signal-present entrapment revalidation

Rerun unchanged on the same six deposited runs, the same five seeds, the same estimated entrapment
fraction (0.50000027) and the same runner. Two extra arms were added *after* seeing that competition
matters, and are labelled as the separate experiments they are: the repaired build with
`--no-psm-competition`, and C++ with `--post-processing-tdc`. The predeclared arms were not altered.

### At reported q<0.01, five seeds

| Arm | accepted | adjusted FDP mean | median | sd | range |
|---|---:|---:|---:|---:|---|
| pre-repair Rust | 19,542.8 | 2.698% | 2.681% | 0.051 | 2.660–2.784% |
| repaired Rust, competition **off** | 19,226.6 | 2.561% | — | 0.070 | — |
| C++ 3.09, `--search-input concatenated` | 19,178.0 | 2.521% | 2.526% | 0.072 | 2.432–2.615% |
| **repaired Rust, competition on (default)** | **19,533.0** | **1.816%** | 1.838% | 0.059 | 1.730–1.884% |

The repaired default accepts essentially the same number of PSMs as the rejected implementation
(19,533 against 19,542.8) at **two-thirds of the false discovery proportion**. This is not a yield
trade.

The two added arms separate the causes. The estimator and cross-validation repairs alone move Rust
from 2.698% to 2.561%, which lands it on C++'s 2.521% — that is, they removed percolator-rs's
*excess* over the reference. The remaining 0.75 percentage points come from spectrum-level
competition, which neither implementation was performing.

### Full threshold curve, five seeds

| Reported q | pre-repair Rust | repaired Rust (default) | C++ (no competition) |
|---:|---:|---:|---:|
| 0.001 | 1.174% | **0.615%** | 1.274% |
| 0.005 | 1.991% | **1.207%** | 1.896% |
| 0.010 | 2.698% | **1.816%** | 2.521% |
| 0.020 | 3.909% | **2.753%** | 3.616% |
| 0.050 | 7.209% | **5.848%** | 6.695% |
| 0.100 | 12.604% | **10.873%** | 11.672% |

The repaired estimator is conservative at q<0.001 (0.615% against a nominal 0.1%… still above
nominal, but now the *only* threshold where it is below both other arms) and runs at roughly 1.1–1.8×
nominal across the rest of the curve, against 2.4–2.7× before. **It is still anti-conservative from
q<0.005 upward and must not be described as calibrated.**

### Why did both implementations fail nominal control?

The hypothesis was that the inputs violate the target-decoy competition assumption, not that either
estimator was uniquely broken. It was tested rather than assumed, and it holds: applying competition
moves the repaired Rust arm from 2.561% to 1.816%, and applying it to C++ through
`--post-processing-tdc` moves C++ to the same place. Competition, not implementation, was the
dominant shared cause.

Rival explanations that were considered and are *not* resolved by this evidence, and which plausibly
account for the residual gap between 1.8% and 1.0%:

- **The entrapment correction itself.** `estimated_false = pure_entrapment / f`, with `f` the
  fraction of accepted non-mixed decoys that are pure entrapment. At tight thresholds `f` is a
  plug-in from a handful of decoys, so the adjusted FDP inherits that sampling noise and any bias in
  the assumption that a false target lands in entrapment space at the same rate a decoy does.
- **Dependence between accepted PSMs.** Peptides shared across spectra and proteins are not
  independent hypotheses, which the interval arithmetic here does not model.
- **The search space itself.** Adding an equal-residue foreign proteome changes the search, so the
  measured FDP describes the entrapment search, not the native-only search a user would run.

None of these were adjusted to move the number toward 1%. The experiment was not altered after the
results were seen.

Manifests: `entrapment-repaired-20260825`, `entrapment-nocompete-20260825` and
`entrapment-cppcompeted-20260825` under `$HOME/percolator_rs_out/scientific-validation/`.

## 6. PEP calibration revalidation

Entrapment-adjusted aggregate calibration through the unchanged
[`pep_entrapment.py`](pep_entrapment.py). The first three rows use the identical 435,261-target
denominator, so they are directly comparable; the repaired default competes, which reduces the list
to one row per precursor and is reported separately rather than compared across denominators.

| Method | targets | PEP exactly 0 | known-false at PEP=0 | weighted abs. calibration error | signed observed − reported |
|---|---:|---:|---:|---:|---:|
| pre-repair Rust | 435,261 | **9,155** | **87** | 0.0634 | +0.0619 |
| repaired Rust, competition off | 435,261 | **0** | **0** | **0.0171** | +0.0171 |
| C++ 3.09 | 435,261 | 0 | 0 | 0.0178 | +0.0160 |
| repaired Rust, default (competed) | 98,680 | **0** | **0** | 0.0185 | +0.0185 |

The pathological zeros are gone, and they are gone because the estimator changed, not because a
floor was applied: with the finite-sample safeguard the smallest achievable q-value is
`lambda / T_total > 0`, so every raw PEP is already positive before the 0.5-false-discovery prior is
added. On the matched denominator the repaired PEPs are marginally better calibrated than the
reference's.

Per-bin calibration for the repaired default (full tables in the `pep-entrapment-*.tsv` artifacts):

| reported PEP | PSMs | mean reported | adjusted observed | error |
|---|---:|---:|---:|---:|
| [0, 1e-12) | 0 | — | — | — |
| [1e-4, 1e-3) | 10,926 | 0.00039 | 0.00384 | +0.0035 |
| [1e-3, 0.005) | 2,978 | 0.00374 | 0.00985 | +0.0061 |
| [0.005, 0.01) | 1,399 | 0.00717 | 0.01654 | +0.0094 |
| [0.01, 0.02) | 1,703 | 0.01411 | 0.03014 | +0.0160 |
| [0.02, 0.05) | 1,326 | 0.03548 | 0.06021 | +0.0247 |
| [0.05, 0.10) | 1,514 | 0.07202 | 0.08134 | +0.0093 |
| [0.10, 0.20) | 1,570 | 0.15863 | 0.19091 | +0.0323 |
| [0.20, 0.50) | 3,796 | 0.34328 | 0.35957 | +0.0163 |
| [0.50, 1.0] | 73,468 | 0.91355 | 0.93493 | +0.0214 |

Every bin is anti-conservative in the same direction, by one to three percentage points. That is a
large improvement on a signed gap of +0.062 and on 87 known-false matches printed at zero, and it is
**not** a demonstration that the PEPs are calibrated. Bins remain sparse at the confident end, and
the observed column is a bin-level entrapment estimate, not individual truth.

## 7. Edge cases

Regenerated through the unchanged [`run_edge_cases.py`](run_edge_cases.py).

| Case | pre-repair | repaired |
|---|---|---|
| 6-row input | exit 0; 2 exact-zero q and PEP | exit 0; **0 exact zeros** |
| 100 targets / 10 decoys | exit 0; 3 exact-zero q and PEP | exit 0; **0 exact zeros** |
| 2,000 duplicate PSMs | exit 0; 54 exact zeros; **equal printed scores got different q** | exit 0; 0 exact zeros; **ties consistent** |
| all features tied | exit 0; **equal printed scores got different q** | exit 0; **ties consistent** |
| malformed feature | exit 0, silently coerced to zero | **exit 1**, `line 2: feature 'lnrSp' is not a finite number: 'not-a-number'` |
| missing feature | exit 0, silently coerced to zero | **exit 1**, located diagnostic |
| NaN feature | exit 0, **non-finite output** | **exit 1**, located diagnostic |
| feature near 1e308 | exit 0; 16 exact zeros | exit 0; 0 exact zeros |
| unusual protein mapping | exit 0; 20 exact zeros | exit 0; 0 exact zeros |

Manifest: `edge-cases-repaired-20260825/manifest.json`.

## 8. Multi-seed stability and cross-dataset behaviour

Five deterministic seeds (1–5) on all four compact datasets, through the unchanged
[`run_multiseed.py`](run_multiseed.py). No per-dataset setting was tuned. The whole study was run
twice — once before and once after a `None`-handling fix in the aggregator, described in §11 — and
every per-run count was identical.

### PSM counts at q<0.01, mean over five seeds

| Dataset | pre-repair Rust | repaired Rust | C++ 3.09 | pre-repair Rust − C++ | repaired Rust − C++ |
|---|---:|---:|---:|---:|---:|
| PXD007145 Tide | 29,252.0 | 27,633.8 | 27,633.0 | +1,619.0 | **+0.8** |
| PXD060954 Sage | 26,617.6 | 25,796.6 | 25,789.4 | +828.2 | **+7.2** |
| PXD020243 MSFragger | 1,541.6 | 1,351.4 | 1,411.2 | +130.4 | **−59.8** |
| Upstream yeast | 1,120.6 | 1,151.4 | 1,111.8 | +8.8 | **+39.6** |

Seed-to-seed standard deviation for the repaired build is 19.5 (Tide), 11.5 (Sage), 25.9
(MSFragger) and 20.5 (yeast) PSMs — comparable to C++'s 27.9 / 4.0 / 38.1 / 34.3 and to the
pre-repair build's own spread. Stability was not traded away.

The Tide and Sage yield gaps — the two largest and most seed-stable differences in the whole
pre-repair study, and the ones most often cited as a percolator-rs advantage — **were statistical
artifacts.** They disappear almost entirely. That is a direct, dataset-level confirmation of what
the complete-null experiment implied.

MSFragger and yeast move in the opposite direction and need the post-processing caveat below.

### Reference agreement before and after

| Dataset | Jaccard at q<0.01 | | score Spearman | | q Spearman | |
|---|---:|---:|---:|---:|---:|---:|
| | before | after | before | after | before | after |
| Tide | 0.9446 | **0.9930** | 0.9986 | 0.9988 | — | 0.9966 |
| Sage | 0.9688 | **0.9959** | 0.9965 | 0.9955 | — | 0.9728 |
| MSFragger | 0.8936 | 0.8825 | 0.8920 | **0.9602** | — | 0.9548 |
| yeast | 0.9257 | 0.9166 | 0.9072 | **0.9424** | — | 0.9423 |

Score rank correlation improves everywhere it was weak. MSFragger — the dataset with the largest
pre-repair ranking divergence — improves from 0.892 to 0.960, consistent with it being the dataset
where percolator-rs had been training on a raw `retentiontime` column that the reference treats as
metadata (§1.10).

### The comparison is no longer matched by default

Counting candidates per precursor explains the remaining Jaccard gaps:

| Dataset | rows | precursors | candidates per precursor |
|---|---:|---:|---|
| Tide | 55,398 | 55,398 | 1 |
| Sage | 35,093 | 35,084 | 1 (9 precursors have 2) |
| yeast | 19,674 | 9,921 | 2 (a target and a decoy) |
| MSFragger | 9,475 | 3,389 | mostly 3 |

On Tide and Sage, competition is a no-op and the comparison stays matched. On yeast and MSFragger
the repaired percolator-rs competes and C++ under `--search-input concatenated` does not, so those
two rows compare different post-processing rather than different implementations.

A separate arm was therefore run with C++ under `--post-processing-tdc`. This is a **new
experiment**, designed after seeing the competition result, and is reported as one:

| Dataset | repaired Rust − C++ (unmatched) | repaired Rust − C++ (matched) | Jaccard (matched) | PEP Spearman (matched) |
|---|---:|---:|---:|---:|
| Tide | +0.8 | +0.8 | 0.9930 | 0.9825 |
| Sage | +7.2 | +15.2 | **0.9962** | 0.9408 |
| MSFragger | −59.8 | **−2.8** | **0.9235** | 0.9492 |
| yeast | +39.6 | **+15.4** | **0.9201** | **0.9325** (from 0.6991) |

Under matched post-processing the repaired percolator-rs and C++ 3.09 agree to within **±15 PSMs on
every one of the four datasets**, against a pre-repair spread of +8.8 to +1,619. Nearly all of the
historical Rust/C++ difference was implementation defects; the residual was a post-processing
mismatch.

Note the direction of the argument: C++ is not a correctness oracle, and this is not offered as
evidence that percolator-rs is now correct. It is evidence that the earlier disagreements had
identifiable causes, all of which were on the percolator-rs side or in how the comparison was
configured.

Manifests: `multiseed-repaired-20260825` and `multiseed-cppcompeted-20260825`.

### C++ under matched competition, entrapment

The same `--post-processing-tdc` arm was run on the entrapment study. Five seeds, at reported
q<0.01:

| Arm | accepted | adjusted FDP |
|---|---:|---:|
| repaired Rust, competition on | 19,533.0 | 1.816% |
| C++ 3.09, `--post-processing-tdc` | 19,515.4 | 1.729% |

Both implementations land in the same place once both compete, and both sit at roughly 1.8× nominal.
That is the clearest available statement of what remains: the residual anti-conservatism on this
experiment is **shared with the reference** and is a property of the data, the search design or the
entrapment correction, not of percolator-rs. It is not fixed, and it is not claimed to be.

## 9. New canonical baseline

PXD032157, all 65 Comet PINs, canonical profile, seed 1. The historical numbers are labelled as
pre-validation and are **not** an acceptance criterion: they were produced by methodology that failed
validation.

| | pre-repair (historical) | repaired | change |
|---|---:|---:|---:|
| Target PSMs at q<0.01 | 107,046 | **106,795** | −251 (−0.23%) |
| Target peptides at q<0.01 | 37,469 | **35,866** | −1,603 (−4.28%) |

The PSM figure barely moves because two corrections very nearly cancel. The finite-sample safeguard
and the fold-score rescaling remove calls; spectrum-level competition adds them back by clearing
losing candidates out of the head of the ranking. On one PXD032157 file the competition effect alone
is +82 PSMs at q<0.01 (185 without, 267 with). The peptide figure falls further because the peptide
list is an order of magnitude shorter, so the safeguard is proportionally larger there.

The right way to read this is *not* "the repair was nearly free". It is that the pre-repair count
happened to be near the corrected one on this dataset while resting on an estimator that rejected a
complete null in 17 of 30 replicates. On Tide and Sage the same repair removes 1,619 and 828 PSMs.
Yield is not evidence either way.

## 10. What each correction cost or gained

A single-factor decomposition across every dataset was not attempted; the corrections interact and
several of them change the ranking rather than only the threshold. What can be stated from the runs
above:

| Correction | Measured effect |
|---|---|
| Safeguard + ties + no count ratio | complete-null rejection 17/30 → 0/30; Tide/Sage yield gap against C++ essentially eliminated |
| PEP from q-values + prior mass | exact-zero PEPs 9,155 → 0; known-false at PEP=0 87 → 0; weighted calibration error 0.0634 → 0.0171 |
| Fold-local preprocessing | isolation tests now fail against the previous all-rows fit; yield effect not separated from the rest |
| Fold-score standardization | fixture 127 → 106 PSMs at q<0.01, the size of the cross-fold scale mismatch |
| Spectrum-grouped folds | fixture 106 → 117 |
| Metadata column contract | MSFragger score Spearman against C++ 0.892 → 0.960 |
| Spectrum-level competition | entrapment adjusted FDP 2.561% → 1.816% at q<0.01, with more PSMs accepted |

## 11. Deviations from a strictly frozen protocol

Recorded rather than smoothed over:

1. **A first pass of the null and entrapment experiments was discarded.** They were launched against
   `target/release/percolator-rs` while later commits were still being built into that path, so the
   binary changed under them. Both output directories were deleted and both were rerun against a
   frozen copy of the `1348b0f` build, whose SHA-256 is recorded next to it. No result from the
   discarded pass is reported.
2. **`run_multiseed.py` had a `None`-handling bug in its aggregator.** With the corrected estimator,
   the yeast fixture accepts nothing at q<0.001 in either implementation, so Jaccard is 0/0 — a case
   the pre-repair estimator never reached because it could report q=0. `summary()` now reports
   undefined replicates as such instead of crashing, and the whole study was rerun; every per-run
   count was identical between the two passes.
3. **Two arms were added after seeing results**, and are labelled as new experiments throughout:
   the repaired build with `--no-psm-competition`, and C++ with `--post-processing-tdc`. The
   predeclared arms were rerun unchanged and are reported unchanged. The `--post-processing-tdc`
   arm needs a wrapper because `--search-input concatenated` and `-Y` together give *no* competition
   — in `Caller.cpp` the `-I` branch runs second and clears `targetDecoyCompetition_`. The wrapper is
   preserved next to the frozen binary.

## 12. Performance rebenchmark

Measured only after the scientific work was finished and frozen; nothing was optimized during or in
response to this. Same host and input as the pre-repair campaign: AMD Ryzen 5 5600G, 6 cores /
12 threads, 65 PXD032157 PINs, 2,295,401,156 bytes, canonical profile, seed 1, three trials each,
machine otherwise idle. Command: `RS_BENCH_BIN=<frozen> bash bench/run_rs.sh canonical N`.

| Measurement | pre-repair | repaired | change |
|---|---:|---:|---:|
| N=4 median wall | 12.063 s | **16.104 s** | 1.34× slower |
| N=4 median peak RSS | 794,952 KiB | **801,940 KiB** | +0.9% |
| Sequential median wall | 36.265 s | **51.061 s** | 1.41× slower |
| Sequential median peak RSS | 206,708 KiB | **209,332 KiB** | +1.3% |
| N=4 speedup over sequential | 3.01× | **3.17×** | — |
| Largest file, 1 thread | not recorded | **1.59 s / 209,056 KiB** | — |
| Largest file, 3 threads | not recorded | **0.88 s / 344,204 KiB** | — |

N=4 trials: 16.104 / 16.092 / 16.185 s. Sequential trials: 50.970 / 51.105 / 51.061 s. Largest-file
screens are the median of five.

The cost is where it should be. Each of the three folds now fits its own normalization, searches its
own initial direction and materializes its own design matrix, which is roughly three times the
preprocessing the single shared matrix used to do. Peak memory is essentially unchanged at the
default because folds run one at a time and each matrix is dropped before the next is built — the
largest-file screens show the difference directly: 209 MB at one thread against 344 MB at three,
where all three matrices are alive at once.

The old 12.063 s figure remains a valid engineering measurement of the old methodology. It is not a
baseline the corrected algorithm failed to hold, because it is not the same algorithm.

Against the recorded C++ default run of 376.2 s at four-file concurrency, the repaired build is
23.4× faster rather than 31.2×. That comparison is between different post-processing on both sides
and is a throughput observation, not a scientific result.

## 13. What still fails, and what remains uncertain

**Still failing.**

1. **Signal-present calibration.** The repaired estimator is anti-conservative from q<0.005 upward
   on the entrapment study: 1.82% adjusted FDP at nominal 1%, against 2.70% before. Better is not
   calibrated. The reference lands in the same place under matched competition (1.73%), so the
   residual is shared, but it is not thereby excused.
2. **PEP calibration.** Every bin is anti-conservative by one to three percentage points. The
   pathological zeros are gone and the weighted error fell from 0.063 to 0.017, which is marginally
   better than the reference on the same denominator — and it is still a measurable gap.
3. **Protein-level inference.** Not revalidated here. Picked-protein q-values inherit the corrected
   estimator and its safeguard, and Bayesian inference now consumes PEPs that are no longer
   unsupported, but the PrEST protein-FDR failures in
   [`../bench/PROTEIN_CALIBRATION.md`](../bench/PROTEIN_CALIBRATION.md) were not re-measured. Protein
   output must still be treated as failing validation.

**Uncertain.**

1. **No complete PSM ground truth exists here.** Entrapment gives partial truth with an estimated
   search-space correction; the null gives complete truth only under a hypothesis that erases all
   signal. Neither measures recall or false negatives.
2. **Zero rejections in 30 null replicates** bounds complete-null FDR at roughly 0.12 above with 95%
   confidence. It is consistent with control at every threshold and demonstrates control at none of
   them below that resolution.
3. **PXD032157 remains development data.** It has been used to choose class weights, to benchmark,
   to build the null, and now to check the repair. It cannot serve as an untouched generalization
   set.
4. **Calibration is tested on one signal-present search design.** Tide, Sage, MSFragger and yeast
   establish agreement and seed stability, not calibration. DIA and other acquisition regimes remain
   untested.
5. **The entrapment correction is a plug-in.** Its effective foreign fraction comes from a small
   number of accepted decoys at tight thresholds, and it assumes a false target lands in entrapment
   space at the rate a decoy does.
6. **Dependence between PSMs is not modelled** anywhere in these intervals — shared peptides and
   proteins make accepted PSMs correlated.
7. **`--select-c` is still selection-biased**, and `--auto-model`, `--rescore-model mlp`,
   `--ensemble`, `--join` and `--rt-features` were not individually revalidated. They inherit the
   corrected estimators and the corrected fold isolation, which is necessary but not sufficient.
8. **The fold-score anchoring differs from the reference** (training-decoy null against held-out
   boundary and median decoy). Both are defensible and neither was shown superior; the choice was
   made for isolation, not for measured effect.

## 14. Verdicts

Each dimension is one of STRONG EVIDENCE, MODERATE EVIDENCE, WEAK EVIDENCE, FAILED VALIDATION or
NOT YET TESTED.

| Dimension | Pre-repair | Repaired | Basis |
|---|---|---|---|
| **IMPLEMENTATION CORRECTNESS** | FAILED VALIDATION | **MODERATE EVIDENCE** | Estimators now match the documented method and the reference's structure; each defect has a regression test verified to fail against the old behaviour; agreement with C++ under matched post-processing is within ±15 PSMs on four datasets. Not STRONG: no independent reimplementation checks these numbers. |
| **CROSS-VALIDATION ISOLATION** (default path) | FAILED VALIDATION | **STRONG EVIDENCE** | Corrupting every held-out feature and flipping every held-out label leaves the fold's direction, normalization and weights bit-identical, with and without `--rt-features`; folds are spectrum-grouped; the property is falsifiable and the test fails against the previous code. Scoped to the default path: `--select-c` remains non-nested by construction, and `--ensemble` builds a label-keyed agreement feature over all rows (§18). |
| **Q-VALUE VALIDITY** | FAILED VALIDATION | **MODERATE EVIDENCE** | Tie-invariant, bounded, monotone, with the finite-sample safeguard and a declared opportunity ratio; matches the reference formula; 0/30 complete-null rejections. Not STRONG: still anti-conservative on the signal-present study. |
| **FDR CALIBRATION** | FAILED VALIDATION | **WEAK EVIDENCE** | Complete null is clean at every threshold, but entrapment sits at 1.8× nominal at q<0.01 and above. Improved by a third, not calibrated. |
| **PEP VALIDITY** | FAILED VALIDATION | **WEAK EVIDENCE** | Derived from a published identity, strictly positive by construction, zero known-false matches at PEP=0, weighted calibration error 0.063 → 0.017 and slightly better than the reference. Every bin is still anti-conservative. |
| **REFERENCE AGREEMENT** | MODERATE EVIDENCE | **STRONG EVIDENCE** | Under matched post-processing, within ±15 PSMs at q<0.01 on all four compact datasets; Jaccard 0.92–0.996; score Spearman 0.95–0.999. C++ is not an oracle, so this bounds implementation divergence rather than establishing correctness. |
| **CROSS-DATASET GENERALIZATION** | WEAK EVIDENCE | **WEAK EVIDENCE** | Five workflows run and agree; calibration is still tested on one signal-present search design and one complete-null construction, both from PXD032157. |
| **REPRODUCIBILITY** | STRONG EVIDENCE | **STRONG EVIDENCE** | Bit-deterministic for a fixed seed; serial and threaded paths byte-identical; the multi-seed study reproduced every per-run count on a full rerun; all experiments recorded with hashes, argv and a frozen binary. |

## 15. What percolator-rs can and cannot claim

**It can claim** that it implements target-decoy competition q-values with the finite-sample
safeguard, exact-tie grouping and a declared opportunity ratio, matching the estimator of
Percolator 3.09; that its posterior error probabilities are derived from those q-values through the
Käll et al. (2008) identity and cannot be zero; that its cross-validation is fold-isolated in a
falsifiable sense, with normalization, initial direction, retention-time alignment and every model
fitted inside the training partition, and with spectra kept whole; that it performs spectrum-level
target-decoy competition, so its statistics apply to a list that satisfies the assumption they
rest on; that it made **no** false discovery in 30 exchangeable-label complete-null replicates at
every threshold from q<0.001 to q<0.10; that it agrees with C++ Percolator 3.09 to within ±15 PSMs
at q<0.01 on four independent datasets under matched post-processing; that it rejects malformed,
missing and non-finite input instead of silently coercing it; and that it is bit-deterministic for a
fixed seed and roughly 23× faster than the recorded reference run on the named host and input.

**It must qualify** any statement about accuracy on signal-present data. Reported q<0.01 corresponds
to an entrapment-estimated 1.8% false discovery proportion on the one search design where that can
be measured — better than the 2.7% it used to be and than the 2.5% the reference reaches without
competition, and matching the 1.7% the reference reaches with it, but not 1%. Reported PEPs are
systematically one to three percentage points optimistic. Both statements come from one PSM search
design; agreement and seed stability generalize across datasets, calibration has not been shown to.

**It cannot claim** that its q-values or PEPs are calibrated; that q<0.01 means 1% empirical error;
that its protein-level output is validated, at either the picked or the Bayesian setting; that
identification yield is evidence of accuracy or sensitivity in either direction; that it is
numerically equivalent to C++ Percolator; that `--select-c` is leakage-free; that `--auto-model`,
`--rescore-model mlp`, `--ensemble`, `--join` or `--rt-features` have been validated; or that any of
this generalizes to DIA or to acquisition regimes not tested here.

**A note on the counts.** On PXD032157 the corrected method reports 106,795 target PSMs at q<0.01
against the historical 107,046. That near-equality is a coincidence of two corrections cancelling,
not a sign that the repair was cosmetic: on Tide and Sage the same corrections remove 1,619 and 828
PSMs that the complete-null experiment shows were statistical artifacts. A lower count with valid
statistics would have been the right outcome; the count happening not to fall is neither evidence
for nor against the repair.

## 16. Artifacts

Under `$HOME/percolator_rs_out/scientific-validation/`:

| Artifact | Contents |
|---|---|
| `frozen-repaired/` | the `1348b0f` binary used for every experiment, its SHA-256 and commit, and the two comparison wrappers |
| `null-repaired-20260825/` | 30 repaired + 30 C++ complete-null runs, hashes, argv, threshold curve |
| `entrapment-repaired-20260825/` | five-seed predeclared entrapment, repaired default vs C++ |
| `entrapment-nocompete-20260825/` | added arm: repaired build with competition off |
| `entrapment-cppcompeted-20260825/` | added arm: C++ with `--post-processing-tdc` |
| `multiseed-repaired-20260825/` | five seeds × four datasets, repaired vs C++ concatenated |
| `multiseed-cppcompeted-20260825/` | added arm: same, against competed C++ |
| `pep-entrapment-*-20260825.{json,tsv}` | PEP calibration for pre-repair, repaired, repaired-nocompete and C++ |
| `edge-cases-repaired-20260825/` | regenerated adversarial fixtures and checks |

The pre-repair artifacts (`null-20260825-v2`, `entrapment-multiseed-20260825`,
`multiseed-20260825-v3`, `edge-cases-20260825`, `qvalue-null-ablation-20260825`) are unchanged and
remain the record of the failure this work responded to.

## 17. README claim audit

Every scientific claim in the README as it stood after the repair, classified. SUPPORTED means the
evidence in this document backs it; SUPPORTED WITH CAVEATS means it holds only with a stated
qualification that the README now carries.

| Claim | Classification | Evidence and wording used |
|---|---|---|
| q-values are tie-invariant, bounded, monotone, with the finite-sample safeguard and a declared opportunity ratio | **SUPPORTED** | Unit invariants; matches `PosteriorEstimator::getQValues` |
| No false discovery in 30 complete-null replicates at every threshold | **SUPPORTED** | `null-repaired-20260825`, 0/30, byte-identical inputs to the failed study |
| Complete-null FDR is controlled | **SUPPORTED WITH CAVEATS** | 0/30 bounds it at ~0.12 at 95%; stated as consistent-with, not demonstrated-at, 0.001 |
| PEPs cannot be exactly zero | **SUPPORTED** | Construction plus 0 observed in 98,680 and 435,261 target PSMs |
| PEPs are calibrated posterior error probabilities | **INSUFFICIENT EVIDENCE** | Every bin anti-conservative by 1–3 points; README says improved, not calibrated |
| q<0.01 means 1% error | **INCORRECT** | Entrapment gives 1.8%; README states this explicitly and says not to read it that way |
| Cross-validation is leakage-free | **SUPPORTED** | Falsifiable isolation tests, verified to fail against the previous code; `--select-c` excluded by name |
| Agreement with C++ 3.09 within ±15 PSMs under matched post-processing | **SUPPORTED** | Four datasets, five seeds, `multiseed-cppcompeted-20260825` |
| percolator-rs identifies more than the reference | **INCORRECT** | Was an artifact; removed. Matched deltas are +0.8, +15.2, −2.8, +15.4 |
| More identifications mean better accuracy or sensitivity | **INCORRECT** | Removed; the README now says counts are not accuracy in either direction |
| Both implementations fail nominal 1% on entrapment | **SUPPORTED** | 1.816% vs 1.729% under matched competition |
| Not competing within a spectrum is correct because the reference does not | **INCORRECT** | The observation was right, the inference wrong; corrected in Fidelity notes |
| pi0 = 1 is right for direct competition | **SUPPORTED** | `Caller.cpp` sets `setUsePi0(useMixMax_)`; reference logs pi0 = 1 on this data |
| Malformed / non-finite input is rejected | **SUPPORTED** | Parser tests and regenerated edge cases |
| Deterministic for a fixed seed | **SUPPORTED** | Byte-identical serial/threaded outputs; multi-seed study reproduced on full rerun |
| 23.4x faster than the reference at N=4 | **SUPPORTED WITH CAVEATS** | One host, one input, unmatched post-processing on both sides; labelled a throughput observation |
| 7.3–14.6x faster on the compact cases | **SUPPORTED WITH CAVEATS** | Pre-repair measurement; README now says so and that it was not re-measured |
| Protein inference is calibrated | **FAILED VALIDATION** | Not revalidated; PrEST failures stand and the README says so |
| `--rescore-model mlp`, `--auto-model`, `--ensemble`, `--join`, `--rt-features` are validated | **NOT YET TESTED** | Listed by name in the README as not revalidated |
| Scientific validity generalizes across datasets | **INSUFFICIENT EVIDENCE** | Calibration tested on one signal-present design and one null construction |

## 18. Skeptical audit of the repair

An attempt to falsify each claim, made after the work was finished.

### 18.1 Is an exact-zero q-value or PEP still reachable?

Scanned every target and decoy PSM row produced by the repaired build across the null, entrapment and
multi-seed studies: **200 files, 2,546,082 rows, zero q-values equal to 0 and zero PEPs equal to 0.**
This is a check of the printed output at six decimals, which is the precision the pathological
pre-repair values were measured at.

### 18.2 Do ties still leak row order into the result, end to end?

The unit tests cover the estimator. The pipeline was checked separately: the yeast PIN was rewritten
with its 19,674 data rows in a different random order and rerun. Every one of the 3,946 output PSMs
has an identical printed score, q-value and PEP, and the q<0.01 count is 665 either way. Fold
assignment is permutation-invariant by construction — spectra are collected into a `BTreeMap` and
shuffled from the seed, not from file order — so this checks the estimator, the competition
tie-break and the reporting path together. It compares printed values; floating-point summation
order inside normalization does change with the permutation, so this bounds the effect below the
printed precision rather than proving bit-equality.

### 18.3 Does the complete-null result survive a construction it was not checked against?

Two additional constructions, declared in [`run_null_variants.py`](run_null_variants.py) before
running and labelled there as new experiments. Ten relabeling seeds (2001–2010) on each of two
PXD032157 PINs, twenty replicates per arm:

| Arm | rows relabelled | pseudo target:decoy | discoveries |
|---|---|---|---:|
| `decoy_balanced` | original decoys | 1:1 | 0/20 at every threshold |
| `target_balanced` | **original targets** | 1:1 | 0/20 at every threshold |
| `decoy_2to1_default_p` | original decoys | **2:1**, `p` left at 0.5 | 0/20 at every threshold |
| `decoy_2to1_declared_p` | original decoys | 2:1, `p` declared 2/3 | 0/20 at every threshold |

The target-derived null is the informative one: those rows are heterogeneous in quality, unlike the
decoy rows the predeclared study uses, and a random label on them is still exchangeable with the
features. It behaves the same.

**The imbalance contrast has no power and should not be read as evidence.** The two 2:1 arms both
report zero because the smallest q-value anywhere in those runs is 0.276 (default `p`) and 0.552
(declared `p`) — an order of magnitude above the loosest threshold tested. A complete null cannot
discriminate an opportunity-ratio error, because the ratio only scales an FDP that is already near
0.5. What the pair does show is that the parameter does the arithmetic it claims: declaring `p = 2/3`
multiplies every q-value by exactly 2.0000, which is `lambda = p/(1-p)` for a 2:1 competition.
Whether the *right* value has been declared for a given search is the user's responsibility, and no
experiment here tests that.

### 18.4 Is there leakage left?

Yes, in two places, both outside the default path and both now named rather than implied:

- **`--select-c`** selects class weights by ranking candidates on the same out-of-fold predictions it
  later reports. This is structural, not a bug, and was not repaired. `--auto-model` is the nested
  alternative.
- **`--ensemble`** builds its `ensemble_psm_engine_count` feature from a map keyed by
  `(ScanNr, Label, Peptide)` over every row, so a label-dependent feature is constructed across the
  whole dataset before folds exist. Ensemble mode was already on the not-revalidated list; this is
  the specific reason it must stay there.

Everything else that reads a label — the peptide rollup, the reported statistics, the per-source
summary — runs after training and is reporting, not preprocessing.

### 18.5 Do the tests encode methodology or the current implementation?

Mostly methodology, with two honest exceptions.

`qvalues_match_the_documented_closed_form` asserts 0.2 and 1.0, both derived by hand from the
documented formula rather than read off a run. The isolation, tie, permutation, positivity and
monotonicity tests assert properties that any correct implementation must have, and the two
leakage tests were confirmed to fail against the previous code.

The exceptions: `dropping_the_safeguard_is_confined_to_the_training_estimator` asserts that
`Tdc::training` *can* report zero, which encodes a deliberate design choice rather than a
methodological necessity; and `a_leading_target_run_never_receives_zero_pep` asserts a floor of
`0.5/n`, which is the specific prior this implementation chose. Both are documented at the assertion.
The yield gates in `tests/` and `bench/` encode current output values by design and prove
determinism, not validity — the README now says so.

### 18.6 Is the improvement seed- or dataset-specific?

Complete null: 0/30 across three source PINs and ten relabelings each, plus 0/20 in each of four
variant arms across two PINs. Entrapment: five seeds, standard deviation 0.059 percentage points, all
five between 1.730% and 1.884%. Multi-seed agreement: five seeds on four datasets, per-dataset PSM
standard deviation 11.5–25.9, comparable to the reference's own 4.0–38.1. The MSFragger arm moved
from Rust reporting +130 PSMs to −2.8 under matched post-processing, which is the largest
single-dataset swing and is explained by the metadata-column correction plus competition.

No seed or dataset was excluded from any table above.

### 18.7 What would falsify the remaining claims?

- An independent signal-present dataset with partial truth on which the repaired q<0.01 exceeds
  roughly 2% adjusted FDP would show the entrapment result does not generalize.
- A complete-null construction with enough leading-target runs to reach q<0.01 — a much larger input,
  or one where the learner has real signal to exploit — would test the safeguard at a resolution
  30 replicates cannot.
- A PEP calibration on a dataset with individual PSM truth would replace the bin-level entrapment
  estimate that every PEP statement here rests on.
