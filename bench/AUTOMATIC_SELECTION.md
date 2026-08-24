# Leakage-free automatic SVM selection

## Outcome

`--auto-model` implements genuine nested validation, but it does **not** improve aggregate
identification yield on the current benchmarks. It remains opt-in; fixed SVM defaults remain the
production default.

On the 65-file PXD032157 benchmark, nested selection reports 106,652 PSMs and 37,636 peptides at
reported q<0.01 versus 107,046 and 37,469 with fixed defaults: -0.37% PSMs but +0.45% distinct
peptides. It wins 33 files, loses 31, and ties one at PSM level, but takes 206.4 seconds versus 20.1
seconds at four-file concurrency (10.3x slower). The older whole-dataset `--select-c` experiment
reports 106,558 PSMs and 37,330 peptides in 49.7 seconds. It wins 32 files, loses 28, and ties five
against fixed weights at PSM level; proper nesting changes the tradeoff but does not produce a
consistent gain.

The four independent extension cases are also lower in aggregate:

| case | fixed PSM / peptide | nested PSM / peptide | nested delta |
|---|---:|---:|---:|
| PXD020243, MSFragger | 1,554 / 1,177 | 1,347 / 1,015 | -207 / -162 |
| PXD060954, Sage | 26,624 / 11,420 | 26,642 / 11,426 | +18 / +6 |
| Hogrebe, Tide | 29,264 / 20,614 | 29,274 / 20,610 | +10 / -4 |
| Percolator yeast | 1,126 / 903 | 1,123 / 878 | -3 / -25 |
| **aggregate** | **58,568 / 34,114** | **58,386 / 33,929** | **-182 / -185** |

On the six-run signal-present entrapment search, nested selection accepts 19,556 PSMs at reported
q<=0.01 versus 19,666 fixed. Adjusted FDP is slightly lower: 2.68% (95% CI 2.42-2.95%)
versus 2.78% (2.53-3.06%). Neither validates nominal 1% FDR.

All 195 PXD032157 outer models retained the existing 1e-5 tolerance; 194 retained all 21 features
and one selected eight. Inner validation mainly varied C and the class-weight ratio, but that
flexibility did not improve aggregate held-out PSM yield. Machine-readable results are in
[`automatic-selection-results.tsv`](automatic-selection-results.tsv).

## Nested design

Every reported PSM is scored by an outer-fold model for which that PSM's fold was excluded from all
of the following:

- feature mean and variance estimation;
- best-feature initialization;
- feature ranking and subset construction;
- hyperparameter scoring and selection;
- semi-supervised positive selection and final model fitting.

Each of the three outer training partitions is split into two inner folds. Candidates are evaluated
with the full final training budget, not the abbreviated proxy used by legacy `--select-c`. After
selection, normalization, feature ranking, and the final model are refit using the complete outer
training partition, then applied once to its untouched outer test fold. Target-decoy q-values and
PEPs are computed only after the three outer held-out score vectors are assembled. Because folds
can choose different C values, held-out margins are standardized using the corresponding training-
decoy mean and variance before pooling; held-out scores never determine this transformation.

Selection is staged to keep the search tractable:

1. Jointly choose SVM regularization scale C from `{0.25, 1, 4}` and negative:positive class-weight
   ratio from `{1, 4, 16}`. Positive weight is fixed at one because scaling C and both class weights
   independently is non-identifiable.
2. Holding that choice fixed, select all features, or the top 8 or top 4 features. Ranking uses
   single-feature target yield on the relevant training partition only; feature identities can
   therefore differ between inner and outer folds.
3. Select Newton gradient-norm tolerance from `{1e-3, 1e-5, 1e-7}`.

Candidates maximize inner held-out target PSM yield at reported q<0.01. First-candidate tie breaking
prefers the current C/weight defaults, all features, and the existing 1e-5 tolerance. Selected
settings are logged separately for each outer fold.

A unit test changes both labels and features in one outer test fold and verifies that its selected
hyperparameters remain identical. The portable regression gate also requires serial and parallel
selection choices and result files to be byte-identical.

## Usage

```bash
cargo build --release
target/release/percolator-rs --canonical --seed 1 --auto-model \
  --results-psms target.psms.tsv input.pin
```

`--auto-model` currently supports only the linear SVM and cannot be combined with `--select-c` or
explicit `--cpos`/`--cneg`. Fixed-mode convergence can be set with `--svm-tolerance`; automatic
mode searches the tolerance grid around that baseline.

Reproduce fixed-versus-nested yield comparisons with:

```bash
bash bench/selection_comparison.sh
SELECTION_BENCH_INPUT=/path/to/pins SELECTION_BENCH_OUT=/path/to/output \
  bash bench/selection_comparison.sh
```

After creating the entrapment searches with `bench/entrapment/run.sh`:

```bash
bash bench/selection_entrapment.sh
```

## Interpretation and limitations

The outer estimate is uncontaminated by model selection, but the search is deliberately staged
rather than a full 81-combination Cartesian grid. Two inner folds also leave selection variance on
small inputs. More repetitions or a one-standard-error rule could stabilize choices, at still
higher cost.

The result does not support enabling automatic selection by default. It trades a small PSM loss for
a small peptide gain on PXD032157, loses both metrics on independent extension cases, and costs an
order of magnitude more. Fixed defaults remain a strong, far cheaper regularizer; per-file
flexibility should not be assumed to improve sensitivity.
