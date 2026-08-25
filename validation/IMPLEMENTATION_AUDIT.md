# Read-only implementation and leakage audit

Audit date: 2026-08-25  
Audited percolator-rs commit: `d83a7ba281f45112b453862820feb94296aaebd3`  
Reference source: C++ Percolator tag `rel-3-09`, commit
`7238ac501494b0f485594f96c77eef662d97785e`  
Reference executable: Percolator 3.09.0, build date 2026-05-21

This report was written before making any methodological change. It audits the production/default
SVM path unless a section explicitly discusses an experimental option. Existing regression tests
are treated as implementation checks, not as evidence of statistical validity.

## Executive finding

The current implementation is recognizably Percolator-style, but the default workflow is **not a
faithful or leakage-free implementation of the reference procedure**, and its PEPs do not implement
the reference 3.09 estimator. Four differences are scientifically material:

1. Default normalization is fitted to all PSM feature vectors before folds are trained.
2. A single initial direction is selected using all target/decoy labels, including every held-out
   fold, and is supplied to all three fold models.
3. Held-out scores from the three default fold models are pooled without the fold-specific score
   normalization used by C++ Percolator. Consequently, cross-fold score scale affects the final
   ranking and q-values.
4. The q-value and PEP estimators differ materially from C++ 3.09. Rust omits the final decoy
   pseudocount and tie grouping; its PEP estimator is an ad hoc transformation of a pooled-label
   PAVA fit and is not justified as a target PSM posterior probability.

The optional `--auto-model` path fixes items 1 and 2 for its outer folds and standardizes held-out
scores from training decoys, but it is a different experimental workflow. It does not repair the
final q-value or PEP estimators. The legacy `--select-c` path additionally selects weights using the
same out-of-fold predictions later reported, so its performance estimate is selection-biased.

These findings invalidate the README statements that the default has “3-fold nested
cross-validation,” that each reported PSM is protected from “overfitting of the FDR estimate,” and
that the implementation is unqualifiedly “Faithful to the Percolator method.” They do not prove
that every identification is wrong. They mean that nominal error control and equivalence require
independent evidence and currently cannot be inferred from the design.

## Pipeline trace: one PSM

For an input row `r` in outer fold 0:

1. `pin::parse` reads its label, scan, feature values, peptide, and proteins. `ExpMass` and
   `CalcMass` are excluded by name. Malformed numeric features are silently converted to zero.
2. `percolator::run` assigns `r` to a deterministic fold after shuffling individual rows. Except in
   ensemble mode, rows from the same spectrum or duplicate candidate are not grouped.
3. The default path fits feature means and standard deviations using **all rows, including `r`**,
   and transforms the entire matrix.
4. It tests every feature in both orientations against **all labels, including the label of `r`**,
   and chooses the direction with the largest target yield below the training q threshold.
5. For the model that scores fold 0, only folds 1 and 2 enter semi-supervised iterations. At each
   iteration, fold-local scores select target positives by the Rust q-value estimator; every decoy
   in folds 1 and 2 is a negative. `r` is not passed to the SVM objective.
6. The trained fold-0 model scores `r`. In the default path this raw margin is pooled directly with
   raw margins from the other two independently fitted models.
7. Global Rust q-values and PEPs are calculated from the pooled out-of-fold score vector. The label
   of `r` and all other labels legitimately enter this final error-estimation step, but the earlier
   use of `r` in normalization and direction selection means the score was not produced under full
   isolation.
8. For peptide output, the best score is selected independently within `(label, core peptide)`;
   target and decoy versions do not directly compete. Q-values and PEPs are recomputed on this
   reduced list.
9. Protein inference consumes the peptide PEPs. Therefore invalid or miscalibrated peptide PEPs
   propagate directly into Bayesian protein inference.

## Data-availability audit

| Stage | Features available | Labels available | Held-out influence in default path |
|---|---|---|---|
| Parsing | All rows | All rows | None yet |
| Feature normalization | All rows | Not used | **Yes: held-out feature distribution sets mean/SD** |
| Fold construction | Row identity; ensemble candidate key only | Ensemble key includes label | Fold assignment itself is seeded; ordinary duplicate/spectrum grouping absent |
| Initial direction | All normalized rows | All rows | **Yes: held-out labels choose feature and sign** |
| Positive selection | Outer-training rows | Outer-training labels | No direct held-out access after the global initialization |
| SVM training | Outer-training rows | Selected targets and all training decoys | No direct held-out rows in objective |
| Legacy C selection | All normalized rows through repeated out-of-fold scores | All labels | **Yes for model selection/evaluation reuse** |
| Fixed-model stopping | Fixed iteration count by default | Training labels | No data-dependent early stopping in canonical profile |
| Held-out prediction | Held-out features | No held-out labels at prediction | Prediction itself is held out; preprocessing/init are not |
| Fold merge | Raw margins from all folds | None | **No fold-specific scale calibration** |
| Final q/PEP | All out-of-fold scores | All labels | Expected for final error estimation, subject to upstream leakage |
| Peptide rollup | All PSM scores and peptide strings | Labels kept separate | Best target and best decoy are selected separately |
| Protein reporting | Peptide score, PEP, mappings | Decoy prefix | Inherits peptide-statistic validity |

## Component audit

### Parsing and feature engineering

- The parser recognizes `Label`, `ScanNr`, and `Peptide` case-insensitively, but feature exclusion
  for `ExpMass` and `CalcMass` is case-sensitive.
- `fast_float::parse(...).unwrap_or(0.0)` silently turns malformed or missing feature values into
  zero. Non-positive or malformed labels become decoys. `atoi` ignores non-digits rather than
  rejecting malformed scans. Short rows are silently skipped.
- No finite-value validation occurs before normalization or training. NaN comparison falls back to
  equality in score sorting, so NaNs can acquire arbitrary ranks.
- Constant features receive scale 1.0, which is numerically safe.
- Optional retention-time features are fitted per complete input file using targets before folds
  are formed. Thus `--rt-features` is supervised preprocessing with direct held-out-label leakage.
- Joint and ensemble modes introduce additional dependency/grouping questions. Ensemble mode keeps
  identical `(scan, label, peptide)` candidates together, but ordinary mode assigns rows
  independently, so duplicate PSMs and multiple candidates from one spectrum may cross folds.

### Normalization

The default calls `build_matrix(ds)`, which fits z-score location and scale on every row. This is
transductive preprocessing. It is weaker leakage than fitting a supervised transform, but it still
violates the requested condition that held-out features cannot influence normalization. The
`--auto-model` outer path fits normalization only on outer-training rows; its inner splits also fit
normalization only on inner-training rows.

C++ 3.09 also normalizes PIN features globally before creating folds. Agreement with that behavior
does not make the validation leakage-free. Separately, C++ normalizes each fold's final SVM score
scale before merging; this Rust default does not.

### Initial direction

The default Rust path chooses one feature and sign on the complete dataset. This is supervised
leakage into every outer test fold. C++ 3.09 instead chooses the initial direction separately from
each outer training set (`SanityCheck::calcInitDirection` calls `trainset[fold]`). The Rust
`--auto-model` path likewise ranks features from outer-training rows and is isolated here.

### Fold construction

- Three folds are balanced by row count after deterministic Fisher-Yates shuffling.
- Folds are not stratified by label or spectrum. Small datasets can therefore have folds missing a
  class, and spectra/candidates can be split across folds.
- Only ensemble mode groups identical candidates, using a key that contains the target/decoy label.
  Target and decoy reports for the same spectrum can therefore still occupy different folds.
- Seeds deterministically change the row permutation and any training subsampling. Fixed full-data
  training has no other stochastic optimizer component.

### Semi-supervised positive selection and iterations

Within a fold, target positives are those whose current training-fold Rust q-value is strictly less
than `test_fdr`; every training decoy is negative. This part does not directly read outer-test rows.
The canonical path performs ten iterations unconditionally. If no positives or no negatives are
available, an iteration silently continues with the previous model rather than failing.

C++ 3.09 uses an initial training FDR on the first iteration, `<=` when materializing positives,
and fold-local q-values with the decoy pseudocount disabled for training. Rust uses the same 0.01
parameter for all iterations and `<`. These are additional, smaller procedural differences.

### SVM objective and class weights

Rust minimizes

`0.5 * ||w||^2 + sum_i C_i * max(0, 1 - y_i w'x_i)^2`,

including the bias in the L2 penalty. This is the same broad L2-loss linear-SVM family as C++
L2-SVM-MFN, but solver, initialization, convergence, and penalty scaling are not demonstrated to be
numerically equivalent.

The canonical fixed weights are `Cpos=1`, `Cneg=4`. C++ 3.09 defaults to a per-fold C grid whose
negative:positive ratios are scaled by the target/decoy input-size ratio. Therefore the canonical
benchmark does not hold model selection or class weighting constant. The repository correctly
identifies class weighting as a major source of yield, but calls it an “accuracy win” without
ground truth; that language is scientifically unsupported.

### Hyperparameters, selection, and stopping

- Canonical fixed hyperparameters are preselected based on PXD032157 benchmark yield, according to
  repository history and README discussion. PXD032157 therefore cannot serve as an untouched test
  set for claims about those defaults.
- Legacy `--select-c` evaluates candidates on the same three out-of-fold predictions used for its
  reported run. It is selection-biased and not nested.
- `--auto-model` performs two inner folds within each outer training partition and keeps outer-test
  features and labels out of model selection. Its staged rather than joint search is legitimate if
  predeclared, but it is not equivalent to C++'s candidate set or selection rule.
- Canonical stopping is a fixed ten semi-supervised iterations plus a Newton tolerance. There is no
  held-out stopping decision. C++ has validation and direction fallback behavior that Rust lacks.

### Held-out scoring and fold merging

Every Rust PSM is scored by an SVM whose loss did not include that row. That necessary property is
present. It is not sufficient for leakage-free CV because preprocessing and initialization differ.

C++ 3.09 linearly rescales each held-out fold so that its score at the selection-FDR boundary is 0
and its median decoy is -1 before merging. Default Rust pools raw SVM margins. Different fold
training sets can produce different intercepts/scales, so raw pooling can mis-rank PSMs across
folds and alter q-values. Only the Rust `--auto-model` path performs a different standardization,
using the mean and SD of training-decoy scores; it is not reference-equivalent.

### PSM q-values and target-decoy competition

Rust sorts all supplied target and decoy PSMs together and computes a running `D/T`, followed by a
reverse cumulative minimum. With `pi0 = estimate_pi0(labels)`, it actually uses
`min(1, pi0 * D/T)`, where `estimate_pi0` is the global count ratio `D_total/T_total` capped at 1.

Problems:

- For an imbalanced list, multiplying by `D_total/T_total` is not a valid correction for unequal
  target/decoy opportunity. If anything, a known unequal-decoy search space requires the reciprocal
  target/decoy opportunity ratio. When targets outnumber decoys, the current formula squares the
  anti-conservative direction: approximately `(D/T)^2` at the end of the list.
- For a balanced concatenated competition, C++ final TDC begins with one decoy (`(D+1)/T`). Rust
  begins with zero. Thus arbitrarily many leading targets can receive q=0 before the first decoy.
- Equal scores are processed in the unspecified order produced by `sort_unstable`. C++ groups all
  members of a tie before assigning the same FDR. Rust q-values can therefore depend on row order
  and can differ within an exact tie.
- Rust does not perform spectrum-level target-decoy competition. Its direct procedure is comparable
  only when the PIN already contains the intended competition winners. The large PXD032157 C++
  headline used auto-detected **mix-max**, while Rust used its direct estimator; those 103,038 and
  107,046 counts are not based on matched post-processing. A later manifest run forced C++
  concatenated behavior, but it still did not add a PSM-level agreement analysis.
- The q-value implementation is bounded and becomes monotone when rows are ordered by the exact
  arbitrary sort permutation. It is not tie-invariant.

These points are sufficient to reject a blanket claim of correct target-decoy q-values. Empirical
calibration must be evaluated separately; the existing entrapment experiment already contradicts
nominal 1% control on its six runs.

### PEP/PAVA

Rust fits PAVA to the binary decoy indicator over the combined target+decoy ranking, then reports
`min(1, 2 * pi0 * fitted_decoy_probability)` as every row's PEP.

For equal target/decoy opportunity, a fitted pooled decoy probability `p` corresponds to the local
target false/target ratio `p/(1-p)`, not `2p`; C++ 3.09's direct TDC helper uses that odds transform.
The Rust `2*pi0*p` formula is only a small-`p` approximation to neither expression generally, and
its count-ratio `pi0` has no established role here. Exact ties are again order-dependent. Sparse
tails can produce exact zero PEPs, which are then consumed as certainty by protein inference.

C++ 3.09 does not use this Rust estimator. Its default derives raw local errors from target
q-values and smooths them against score with a monotone I-spline; `--pava-pep` changes only the
smoother. It also includes a pseudocount. PEP correlation with C++ is therefore a comparison of
different statistical procedures, not a numerical port check.

The Rust PEPs are bounded and monotone under its arbitrary sorted order, but boundedness and
monotonicity do not establish probability calibration. Until validated against legitimate labels,
the PEPs must be classified as **statistically unsupported**.

### Peptide-level statistics

Rust keeps the highest-scoring PSM for each `(label, core peptide)` and recomputes the same q-value
and PEP estimators. Modified peptide syntax is reduced only by removing flanking residues; distinct
modification forms remain distinct. Because target and decoy representatives are selected
separately, this is not explicit target-decoy competition between matched peptide hypotheses.
Every PSM score and statistical issue above propagates to peptide results.

### Protein-level statistics

Picked-protein inference pairs target/decoy groups by decoy-stripped protein names, retains the
higher score, then calls the same no-pseudocount, non-tie-aware q-value function with `pi0=1`.
Unpaired target-only groups always enter as target wins. Its synthetic regression demonstrates an
implementation property, not calibration.

Bayesian protein inference consumes peptide PEPs as evidence. Since PSM/peptide PEP validity has not
been established, Bayesian protein posterior probabilities cannot be validly interpreted even if
belief propagation is numerically correct. Existing PrEST results independently show failed 1%
protein-FDR validation and false groups in blank controls.

## Leakage verdict by requested influence

| Possible held-out influence | Default canonical | `--select-c` | `--auto-model` |
|---|---|---|---|
| Normalization | **Yes** | **Yes** | No at outer/inner fit |
| Feature engineering (`--rt-features`) | **Yes, supervised** | **Yes** | **Yes before nested path** |
| Initial direction | **Yes, labels** | **Yes, labels** | No at outer/inner fit |
| Hyperparameters | Fixed, but developed on headline data | **Yes, evaluation reuse** | No outer-test access |
| Positive selection | Training fold only | Training fold only | Training fold only |
| SVM objective | Training fold only | Training fold only | Training fold only |
| Stopping decisions | Fixed iteration count | Fixed iteration count | Fixed iteration count |
| Model selection | None within run | **Not nested** | Nested |
| Held-out prediction | Model loss excludes row | Model loss excludes row | Isolated |
| Final q/PEP | All labels, as required for reporting | Same | Same |

Verdict: **default and `--select-c` are not genuinely leakage-free**. `--auto-model` is substantially
cleaner for the base PIN features, but the final estimators remain unsupported and optional
feature engineering can reintroduce leakage.

## Existing evidence and what it establishes

- Exact-output regression tests establish determinism for a fixed seed and guard code changes.
  They do not establish correctness or calibration.
- The three-file pure-null experiment is correctly motivated—relabelled decoy rows make every
  accepted pseudo-target false—but one realization per file and evaluation mainly at 1% do not
  characterize FDR control. “Strongly conservative” is too strong without repeated null datasets
  and a threshold curve.
- The six-run entrapment experiment is the strongest current PSM calibration evidence. It finds
  adjusted FDP 2.78% for Rust at reported q<=0.01 and 2.62% for C++; both fail nominal 1%. Its
  Wilson intervals condition on a plug-in decoy-derived entrapment fraction and treat PSMs as
  independent, so the reported intervals understate clustering uncertainty.
- Multi-dataset results show software compatibility and dataset-dependent yield. They are not
  ground-truth or calibration evidence.
- PXD032157 was used repeatedly to choose class weights and benchmark variants. It is development
  data, not independent evidence for a general identification gain.

## Audit verdicts before experiments

| Question | Audit-only verdict |
|---|---|
| Intended Percolator-style learning | **Partial implementation with material departures** |
| Leakage-free nested 3-fold CV | **Failed for default; substantially improved only in experimental auto-model** |
| Target-decoy q-values | **Methodologically invalid for imbalance; incomplete TDC safeguards** |
| PEPs | **No valid statistical justification found; not reference-equivalent** |
| q<0.01 empirical control | **Already failed on the available signal-present entrapment study** |
| Extra Rust yield | **Confounded by class weights, fold-score merging, q estimator, PEP method, and mismatched C++ post-processing** |
| Seed stability/generalization | **Insufficient audit evidence; requires experiments** |

## Required next evidence

The next phases must use matched target-decoy post-processing, preserve every seeded run, compare
PSMs rather than totals, repeat null and entrapment controls over the full threshold grid, and
ablate one methodological difference at a time without replacing the canonical implementation.
Until then, scientific claims must use “reported-q identification yield,” never “accuracy” or
“sensitivity,” and must prominently state the failed entrapment calibration.
