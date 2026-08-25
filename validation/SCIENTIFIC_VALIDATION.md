# Adversarial scientific validation of percolator-rs

Validation date: 2026-08-25  
Audited implementation: `d83a7ba281f45112b453862820feb94296aaebd3`  
Reference source: Percolator `rel-3-09`, commit
`7238ac501494b0f485594f96c77eef662d97785e`  
Reference executable: Percolator 3.09.0, build date 2026-05-21  
Threshold convention in new experiments: strict `q < threshold`

## Executive answer

The canonical/default percolator-rs workflow **fails scientific validation as an estimator of
PSM q-values and PEPs**. It remains a deterministic, fast Percolator-style rescoring program, and
its reported-q identification counts are reproducible. Those facts do not validate its error
probabilities.

The decisive negative control is an exact-balance, exchangeable-label complete null. Every
pseudo-target discovery is false by construction. Rust made at least one false discovery in 17 of
30 independent relabelings at every predeclared threshold from 0.001 to 0.10. Thus its empirical
complete-null FDR was 56.7%, not the nominal 0.1%--10%. The accepted targets all had reported q=0.
C++ made no discoveries in 29 successful runs and failed closed once. This is consistent with the
source audit: Rust omits the final `+1` decoy safeguard, does not group ties, and applies an invalid
count-ratio `pi0` correction.

On a signal-present, six-run native-plus-foreign-proteome entrapment search, five model seeds gave
mean adjusted FDP 2.70% for Rust and 2.52% for C++ at reported q<0.01. Both are anti-conservative;
Rust is modestly worse at this cutoff. Rust also prints PEP=0 for 9,155 target PSMs, including 87
known-false foreign-proteome PSMs. Therefore neither Rust q<0.01 nor Rust PEP can be presented as a
validated 1% error statement.

Rust and C++ agree closely on the Tide and Sage rankings, less closely on the legacy yeast fixture,
and poorly on the MSFragger and a forced-concatenated PXD032157 file. Most q<0.01-exclusive PSMs in
the compact cases are near the other method's cutoff, but some are not, and the overall MSFragger
rank divergence is large. The Rust yield advantage is stable on Tide and Sage, variable on
MSFragger, and reverses on yeast. It is best described as **dataset-dependent reported-q yield
caused by several model and statistical differences**, not increased accuracy or sensitivity.

## 1. Implementation audit

The full read-only audit and one-PSM trace are in
[`IMPLEMENTATION_AUDIT.md`](IMPLEMENTATION_AUDIT.md). The implemented path is recognizably
Percolator-style:

1. parse a PIN and exclude selected metadata columns;
2. normalize features and append a bias;
3. choose a single feature/sign as an initial direction;
4. create three folds;
5. iterate score, q-value-based positive selection, and squared-hinge linear-SVM fitting within
   each outer training partition;
6. score the held-out partition and merge all fold scores;
7. calculate PSM q-values and PAVA-derived PEPs;
8. select best PSMs per peptide and recalculate statistics;
9. optionally perform picked or Bayesian protein inference.

It is not a faithful numerical or statistical port of C++ 3.09. Material differences include
fixed absolute Rust weights (`Cpos=1`, `Cneg=4`) versus C++ fold-local grid selection, one globally
selected Rust initial direction versus a C++ direction selected from each outer training set, raw
Rust fold-margin pooling versus C++ fold-specific score rescaling, and different q-value/PEP
estimators. C++ maps each held-out fold's selection boundary to zero and median decoy to -1 before
merging, a behavior described in the original cross-validation work and present in the inspected
3.09 source. Rust does not do this in the default path.

The Rust squared-hinge objective belongs to the intended broad linear-SVM family, but matching loss
family does not establish solver, penalty, selection, score-scale, or output equivalence.

### PSM-level q-values

For scores sorted high to low, Rust reports the reverse cumulative minimum of

`min(1, min(1, D_total/T_total) * D_running/T_running)`.

This has three critical defects:

- no final-decoy pseudocount, so a leading target run receives q=0;
- no score-tie grouping, so equal-score results can depend on unstable row order;
- a count ratio applied in the wrong direction for target-heavy lists, where it can make an already
  small `D/T` still smaller.

Rust does not itself perform spectrum-level target-decoy competition. Its statistic is defensible
only if its input already represents the intended competition winners and the estimator correctly
accounts for that design. Neither is enforced.

### PEPs

Rust applies PAVA to pooled binary decoy indicators and transforms the fitted probability `p` as
`min(1, 2*pi0*p)`. With equal target/decoy opportunity, the corresponding local false-target odds
would be `p/(1-p)`, not `2p`. C++ 3.09 instead starts from target q-values and uses a monotone
score smoother by default, with a pseudocount. Boundedness and monotonicity of Rust's PEP vector do
not supply a statistical interpretation.

### Peptide and protein propagation

Peptide statistics reuse the same q-value and PEP procedures after selecting the best score within
separate `(label, peptide)` groups. Picked-protein q-values reuse the no-pseudocount, non-tie-aware
q-value function. Bayesian protein inference consumes the unsupported peptide PEPs as probabilities.
The existing PrEST study independently finds anti-conservative protein inference and false groups in
blank controls; see section 8.

## 2. Leakage audit

The default workflow is **not genuinely leakage-free or nested**.

| Operation | Can an outer held-out row influence it? | Finding |
|---|---:|---|
| PIN parsing | no selection yet | malformed/missing numbers silently become zero; NaN is accepted |
| Feature normalization | yes, features | default mean/SD use all rows |
| Initial direction | **yes, features and labels** | selected once on all rows, then supplied to every fold |
| Fold-local positive selection | no direct access | uses outer-training scores/labels only |
| SVM loss | no direct access | the scored row is excluded from the corresponding model objective |
| Default hyperparameters | development reuse | fixed after extensive PXD032157 yield exploration |
| `--select-c` | **yes, evaluation reuse** | candidates are selected on predictions later reported |
| `--auto-model` base features | no | outer/inner preprocessing and selection are training-only |
| `--rt-features` | **yes, supervised** | fitted to all targets before even the nested path |
| Fold merge | no labels | raw fold scales are pooled in default mode |
| Final q/PEP | all scores and labels | expected for error reporting, but upstream scores leaked |

The experimental `--auto-model` path substantially improves isolation for base PIN features, but
it uses a different training-decoy score standardization, does not repair q-values or PEPs, and can
still be contaminated by `--rt-features`. Ordinary mode also fails to group duplicate reports or
all candidates from the same spectrum into one fold.

## 3. Rust versus C++ PSM-level agreement

[`psm_agreement.py`](psm_agreement.py) qualifies rows by result-directory parent, label, PSM ID,
peptide, and protein string. Exact observable duplicate keys are preserved for multiset threshold
counts but excluded from manufactured one-to-one correlations. It computes Pearson and Spearman
correlations, normalized rank differences, threshold intersections, method-exclusive IDs, and
representative disagreements at all six predeclared thresholds.

### Seed-1 compact datasets at q<0.01

| Dataset | Rust / C++ | Intersection | Rust-only / C++-only | Jaccard | Score Spearman | q Spearman | PEP Spearman |
|---|---:|---:|---:|---:|---:|---:|---:|
| PXD007145 Tide | 29,264 / 27,617 | 27,616 | 1,648 / 1 | 0.9437 | 0.9985 | 0.9949 | 0.9764 |
| PXD020243 MSFragger | 1,554 / 1,388 | 1,377 | 177 / 11 | 0.8799 | 0.8828 | 0.8823 | 0.8560 |
| PXD060954 Sage | 26,624 / 25,795 | 25,795 | 829 / 0 | 0.9689 | 0.9970 | 0.9791 | 0.9611 |
| Upstream yeast fixture | 1,126 / 1,147 | 1,106 | 20 / 41 | 0.9477 | 0.8997 | 0.8996 | 0.9066 |

The legacy fixture reuses PSM IDs, so its correlations apply only to unambiguous qualified keys.
The threshold counts retain all rows as multisets.

### Are Rust-only discoveries borderline?

Mostly, but not exclusively. Among unambiguous Rust-only q<0.01 PSMs, the counterpart C++ q-value
median/p95 and the number with C++ q>=0.05 were:

| Dataset | C++ q median | C++ q p95 | C++ q>=0.05 |
|---|---:|---:|---:|
| Tide | 0.01870 | 0.03488 | 3 |
| MSFragger | 0.01603 | 0.04050 | 7 |
| Sage | 0.02526 | 0.05274 | 55 |
| Yeast | 0.02065 | 0.03154 | 0 |

Tide's q<0.01 count gap is almost entirely a boundary shift on a nearly identical ranking. Sage is
also highly rank-correlated, but 55 Rust-only PSMs are not close by the q>=0.05 criterion.
MSFragger's absolute normalized rank-difference p95 is 29.1% over all matching PSMs, so its
disagreement is a substantial ranking difference, not merely threshold arithmetic.

One forced-concatenated PXD032157 file is more severe: score Spearman is 0.471, Rust reports 208
q<0.01 targets, C++ reports zero, and Jaccard is zero. Across the 65-file forced-concatenated
manifest, Rust has more q<0.01 PSMs in 53 files and fewer in 12; three files with zero C++ but
positive Rust contribute 637 Rust calls. The original headline comparison is even less controlled:
C++ auto-detected `separate` and used mix-max while Rust used a direct statistic. It was never a
matched statistical comparison.

Representative disagreement rows, q bands, and full distributions are preserved in
`$HOME/percolator_rs_out/scientific-validation/agreement-*-v2.json` and the v3 multi-seed manifest.

## 4. Multi-seed results

Five seeds (1--5) were predeclared and retained for both implementations on every compact dataset.
C++ was explicitly forced to concatenated input; no per-dataset Rust setting was tuned.

### PSM counts at q<0.01

| Dataset | Rust mean / median / SD / range | C++ mean / median / SD / range | Mean Rust-C++ | Mean Jaccard |
|---|---:|---:|---:|---:|
| Tide | 29,252 / 29,247 / 13.8 / 32 | 27,633 / 27,617 / 28.0 / 60 | +1,619 | 0.9446 |
| MSFragger | 1,541.6 / 1,554 / 22.3 / 46 | 1,411.2 / 1,399 / 38.1 / 94 | +130.4 | 0.8936 |
| Sage | 26,617.6 / 26,620 / 6.2 / 14 | 25,789.4 / 25,790 / 4.0 / 11 | +828.2 | 0.9688 |
| Yeast | 1,120.6 / 1,125 / 11.2 / 29 | 1,111.8 / 1,111 / 34.3 / 88 | +8.8 | 0.9257 |

Tide and Sage show a stable, practically large count gap. MSFragger shows a persistent mean gap but
greater seed sensitivity. The yeast difference changes direction and seed 1 favored C++. Thus the
headline seed is representative for Tide/Sage, within the observed MSFragger distribution, and not
evidence of a universal Rust gain.

The machine-readable successful rerun is
`$HOME/percolator_rs_out/scientific-validation/multiseed-20260825-v3/manifest.json`; it records
software, platform, seeds, exact argv vectors, input and output hashes, and evaluation-script hashes.

## 5. Pure-null calibration

### Construction and interpretation

For each of three PXD032157 PINs, the runner takes only original decoy rows. It assigns an exact half
to pseudo-target and half to pseudo-decoy using ten deterministic relabeling seeds (1001--1010).
Features are therefore exchangeable with respect to the new label and no native target signal
remains. The two classes are exactly balanced. Every accepted pseudo-target is false by
construction.

Under a complete null, FDP is 1 if the program reports at least one pseudo-target and 0 otherwise.
Consequently FDR is the probability of **any** discovery over repeated relabelings; dividing false
discoveries by the very large total number of pseudo-target candidates is the wrong calculation.
Relabelings of one source PIN are not independent biological datasets, so intervals are descriptive
and source-stratified results are retained.

### Results

| Input | Rust runs with any false discovery | Empirical complete-null FDR (Wilson 95%) | Mean / max false targets | C++ successful runs with any |
|---|---:|---:|---:|---:|
| PXD032157 file 1 | 6/10 | 0.60 (0.313--0.832) | 1.8 / 9 | 0/10 |
| PXD032157 file 2 | 6/10 | 0.60 (0.313--0.832) | 1.3 / 4 | 0/10 |
| PXD032157 file 3 | 5/10 | 0.50 (0.237--0.763) | 1.0 / 3 | 0/9 successful; 1 failed closed |
| Aggregate | **17/30** | **0.567 descriptive** | 1.37 / 9 | **0/29** |

These results are identical for thresholds 0.001, 0.005, 0.01, 0.02, 0.05, and 0.10 because every
Rust false discovery in these runs has reported q=0. The observed 56.7% FDR is grossly above every
tested nominal threshold. It cannot be explained as Monte Carlo uncertainty around valid 1%
control.

The manifest and threshold table are
`$HOME/percolator_rs_out/scientific-validation/null-20260825-v2/{manifest.json,calibration.tsv}`.
They include each exact label permutation hash and argv. This experiment fails the Rust q-value
validation.

## 6. Threshold calibration curves

### Signal-present entrapment, five seeds

The native PXD032157 FASTA was augmented with an approximately equal-residue foreign plant
proteome. Pure foreign-proteome target assignments are known errors. The effective foreign search
fraction is estimated among non-mixed decoys at each method/threshold; therefore these are
bin-level adjusted estimates, not individual truth for native assignments.

| Reported q | Rust mean accepted / adjusted FDP | C++ mean accepted / adjusted FDP |
|---:|---:|---:|
| 0.001 | 14,315 / 1.174% | 12,892 / 1.274% |
| 0.005 | 18,126 / 1.991% | 17,599 / 1.896% |
| 0.010 | 19,543 / **2.698%** | 19,178 / **2.521%** |
| 0.020 | 21,098 / 3.909% | 20,673 / 3.616% |
| 0.050 | 23,814 / 7.209% | 23,273 / 6.695% |
| 0.100 | 26,985 / 12.604% | 26,210 / 11.672% |

Both implementations are anti-conservative over the whole curve. At 0.001, C++ has the slightly
larger relative departure; Rust is worse from 0.005 upward. Rust q<0.01 FDP ranges only
2.660%--2.784% across seeds (SD 0.051 percentage points), while C++ ranges 2.432%--2.615% (SD
0.072 points). The calibration failure is stable, not a cherry-picked seed.

The complete-null curve is discontinuous in the most concerning way: Rust's discovery count is
constant from 0.001 to 0.10 because the top calls are assigned q=0. Together the null and
signal-present controls reject calibrated threshold behavior.

## 7. Negative controls

Two defensible controls were used:

1. **Exchangeable relabeling of original decoys** tests the full learner and statistic when labels
   carry no feature signal. It fails as described above. It does not imply biological independence
   between relabelings, so results are stratified by source PIN.
2. **All features tied** makes labels unlearnable and directly stresses tie handling. The program
   exits successfully, but identical printed scores receive different q-values and the q sequence
   is not monotone when sorted by the printed score. Because output scores are rounded, printed-tie
   tests can sometimes overstate internal ties; in this fixture the input feature vectors are
   actually identical, so the broader lack of tie grouping remains a source-confirmed defect.

Arbitrary feature shuffling was not added as a separate claim because the exact-label relabeling is
cleaner: it preserves the complete joint feature distribution while enforcing the null hypothesis.

## 8. Ground-truth evidence

No public dataset in the completed PSM comparison supplies complete, defensible labels for every
candidate PSM. Therefore true positives, false positives, false negatives, recall, and a full
precision-recall curve were **not** fabricated. Target-decoy calls were not renamed ground truth.

The PXD032157 entrapment search supplies **partial PSM truth**: foreign-only assignments are known
false, but a native assignment is not automatically true. It supports an adjusted FDP estimate and
strongly falsifies exact calibration, but cannot measure recall or total false negatives.

PXD008425 PrEST supplies **protein-level present/absent truth**, not complete PSM truth and not a
matched Rust/C++ protein comparison. On held-out test replicate 3 at reported protein q<=0.01,
picked inference has adjusted FDP 3.23%, 4.14%, and 5.66% in A, B, and A+B, and reports 17 absent
PrEST groups in the blank. Fixed Bayesian inference is much worse; calibration-selected Bayesian
parameters improve it but still fail nominal 1% in B and report three false blank groups. This
validates neither Rust protein mode globally. Full results are in
[`../bench/PROTEIN_CALIBRATION.md`](../bench/PROTEIN_CALIBRATION.md).

## 9. Multi-dataset results

The evaluated matrix includes human, mosquito, bacterial, and yeast data; Q Exactive HF, Fusion
Lumos-study, timsTOF Pro, and legacy unknown instrument metadata; Comet, Tide, MSFragger, Sage, and
SEQUEST-style inputs; small fixtures through a 2.3 GB metaproteomics-scale benchmark. Details and
source accessions are in [`../bench/MULTI_DATASET.md`](../bench/MULTI_DATASET.md).

| Case | PSM evidence | Peptide evidence | Calibration evidence | Scientific interpretation |
|---|---|---|---|---|
| PXD032157 Comet | 65-file yield plus one-file PSM agreement | count comparison | complete null and six-run entrapment | Rust q/FDP fails; original C++ post-processing mismatched |
| PXD007145 Tide | five-seed PSM agreement | five-seed counts | none with legitimate truth | stable near-identical ranking, larger Rust reported-q list |
| PXD020243 MSFragger | five-seed PSM agreement | five-seed counts | none with legitimate truth | largest rank divergence; no accuracy conclusion |
| PXD060954 Sage | five-seed PSM agreement | five-seed counts | none with legitimate truth | stable high agreement; some substantial q disagreements |
| Yeast fixture | five-seed PSM agreement | five-seed counts | none with legitimate truth | yield direction changes; PSM IDs partly ambiguous |
| PXD008425 PrEST | protein modes only | upstream evidence | present/absent and blank at protein level | protein-FDR failure; no matched C++ PSM comparison |

Compatibility and seed stability generalize better than calibration. Calibration has been directly
tested on one underlying PSM search design and one protein standard. The compact datasets cannot be
declared calibrated merely because their Rust/C++ rankings agree. DIA behavior remains untested.

## 10. Ablations and causal explanation

No canonical implementation was modified. Available controlled or quasi-controlled evidence is:

| Factor | Evidence | Result | Causal strength |
|---|---|---|---|
| Final decoy pseudocount | fixed printed-score null ablation | `D+1`, with tie grouping, produces no false calls in these preserved score vectors | strong for the zero-q mechanism; not an end-to-end fix |
| Tie grouping | same ablation and tied fixture | grouping reduces row-order artifacts, but `D/T` without `+1` still rejects often | strong implementation evidence |
| Count-ratio `pi0` | source derivation and imbalance fixture | multiplication is anti-conservative for target-heavy lists | strong mathematical evidence; no full retrain |
| Class weighting | repository's historical one-factor measurements | changing the balance heuristic to fixed 1/4 caused the largest historical yield reversal | moderate, because this predates the present validation harness |
| Legacy C selection | fixed versus `--select-c` benchmark | 107,046/37,469 versus 106,558/37,330 on PXD032157 | count effect measured; selection is not nested |
| Nested path | fixed versus `--auto-model` | 106,652/37,636 on PXD032157; lower aggregate extension yield | leakage is improved; q/PEP remain unchanged and invalid |
| MLP learner | same outer/statistical path | fewer reported-q IDs and similar entrapment failure | learner family does not repair estimator calibration |

At q<0.01, a post-hoc fixed-score ablation gives complete-null discovery probabilities of 66.7%
for Rust output scores under rowwise `D/T`, 43.3% under tie-grouped `D/T`, 6.7% under rowwise
`(D+1)/T`, and 0% under tie-grouped `(D+1)/T`. On C++ output scores they are 100%, 37.9%, 0%, and
0%, respectively. The scores are printed to six decimal places, so this is **not an exact internal
ablation**: printed ties and ordering can differ from full precision. It demonstrates the mechanism
on preserved output and motivates a separately validated estimator change; it does not prove that
adding `+1` alone would validate the complete pipeline or signal-present calibration.

The answer to “why more IDs?” is necessarily plural: fixed class weights, different initialization,
unscaled fold margins, different model selection, direct versus mix-max post-processing in the
headline run, no final-decoy pseudocount, no tie grouping, and the count-ratio adjustment all move
rankings or thresholds. The null result proves that at least part of the extra tail can be a
statistical artifact. There is no evidence that the whole difference is false, nor evidence that it
is a sensitivity gain.

## PEP calibration supplement

Entrapment-adjusted calibration bins use the empirical foreign fraction among non-mixed decoys in
the same PEP interval. Across 435,261 target PSMs per method:

| Method | Printed PEP=0 | Known-false foreign PSMs with PEP=0 | Weighted absolute bin gap | Signed observed-minus-reported gap |
|---|---:|---:|---:|---:|
| Rust | 9,155 | **87** | 0.0634 | +0.0619 |
| C++ | 0 | 0 | 0.0357 | +0.0350 |

“PEP=0” means exactly zero after six-decimal output rounding; an internal value may have been a
small positive number. Even with that precision caveat, assigning reported zero to 87 known-false
PSMs and an average positive calibration gap of 6.2 percentage points is incompatible with a
calibrated posterior error probability. C++ is also anti-conservative in this analysis, but less so.

This is aggregate partial-truth calibration. Native PSM correctness remains unknown and sparse
decoy bins are noisy. These limitations weaken numerical precision, not the qualitative Rust PEP
failure. Bayesian protein posteriors inherit the problem when they consume peptide PEPs.

## 11. Edge cases and numerical behavior

[`run_edge_cases.py`](run_edge_cases.py) generated hashed fixtures and recorded exact commands.

| Case | Outcome |
|---|---|
| 6-row input | exits 0; finite/bounded output; two q=0 and PEP=0 calls |
| 100 targets / 10 decoys | exits 0; demonstrates unsupported imbalance correction |
| 2,000 duplicate PSMs | exits 0; 54 exact-zero q/PEP; printed-score tie inconsistency |
| all features tied | exits 0; equal scores can receive different q-values |
| malformed feature | silently accepted and coerced to zero |
| missing feature | silently accepted and coerced to zero |
| NaN feature | exits 0 and emits non-finite score/statistic output |
| feature near `1e308` | exits 0 and happened to remain finite in this fixture |
| unusual/mixed protein mapping | exits 0; syntax tolerated |

Constant features are numerically protected with scale one. That does not compensate for accepting
malformed, missing, or non-finite input silently. The implementation should eventually fail with a
specific diagnostic, but this validation did not repair it.

## 12. Statistical analysis

- **Effect size:** the mean q<0.01 Rust count advantage is +1,619 Tide PSMs, +130 MSFragger,
  +828 Sage, and +9 yeast, but the denominator and ground truth differ, so aggregating them into one
  percentage would be misleading.
- **Variability:** Rust seed SD is small relative to the Tide/Sage count gaps, while MSFragger and
  yeast show material method/seed interaction. Full mean, median, SD, min, max, and range live in
  the manifest rather than only selected summaries here.
- **Paired evidence:** every seed uses the identical PIN for Rust and C++. Overlap and Jaccard are
  more informative than unpaired total counts. There are only four compact datasets, so no
  asymptotic dataset-level significance test is justified.
- **Calibration versus yield:** the five-seed entrapment Rust-C++ FDP difference at 1% is about
  +0.177 percentage points, small compared with their shared 1.5--1.7 point excess over nominal.
  Rust nevertheless reports about 365 more accepted PSMs on average and about 20 more estimated
  false PSMs.
- **Intervals:** null Wilson intervals are binomial descriptions within each source PIN. Entrapment
  intervals condition on a plug-in effective search fraction and ignore correlated PSMs, so they
  are descriptive. These dependencies make stronger p-values inappropriate.
- **Accuracy:** no complete PSM ground truth exists here. Identification yield and correlation are
  not accuracy, precision, recall, or sensitivity.

## 13. Failed experiments and provenance

Failures were retained rather than deleted:

- `multiseed-20260825/manifest.partial.json`: C++ exit 127 because its local shared-library path was
  not supplied. The runner was amended to require and record the library directory; successful
  results are in v3.
- `null-20260825/manifest.partial.json`: first full null run stopped on a C++ no-initial-direction
  case. The runner was changed only to preserve that fail-closed status and continue other
  predeclared replicates, not to count it as a successful zero.
- `null-smoke-20260825`: a two-replicate smoke test, retained but excluded from the 30-run estimate.
- The entrapment runner detected an extra generic `comet-out/comet.pin` and refused to treat it as a
  biological run. The final study uses exactly the six named deposited runs.

Primary machine-readable artifacts under
`$HOME/percolator_rs_out/scientific-validation`:

| Artifact | Contents |
|---|---|
| `multiseed-20260825-v3/manifest.json` | 40 matched compact runs, hashes, commands, aggregate agreement |
| `null-20260825-v2/manifest.json` | 59 successful null runs plus one fail-closed C++ run |
| `null-20260825-v2/calibration.tsv` | threshold curve and Wilson intervals |
| `entrapment-multiseed-20260825/manifest.json` | 10 six-PIN aggregate runs and threshold/overlap curves |
| `qvalue-null-ablation-20260825.json` | fixed printed-score estimator ablations |
| `pep-entrapment-{rust,cpp}-20260825.{json,tsv}` | bin-level PEP calibration |
| `edge-cases-20260825/manifest.json` | generated fixtures, commands, hashes, checks |
| `agreement-*-v2.json` | detailed seed-1 PSM disagreement characterization |

The repository was dirty only because validation scripts and documentation were being added; the
audited Rust binary remained the build of `d83a7ba`. No production methodology was changed during
validation. Public data accessions, input construction, and compatibility transformations are
recorded in [`../bench/DATASETS.md`](../bench/DATASETS.md) and
[`../bench/MULTI_DATASET.md`](../bench/MULTI_DATASET.md).

## 14. Limitations

The experiments can falsify current calibration more strongly than they can identify a replacement
method. The PSM truth is partial, null relabelings are clustered within three source PINs,
entrapment correction uses an estimated effective search fraction, and the component ablation uses
rounded output scores. Agreement datasets without truth cannot validate either method. The findings
therefore support stopping invalid claims now, but a corrected implementation would require a new,
independent validation campaign rather than inheriting validity from this one.

## 15. README claim audit

This table classifies the major claims as they existed before this validation and records the
required wording. The README has been rewritten accordingly.

| Claim | Classification | Evidence and corrected claim |
|---|---|---|
| Faithful Percolator implementation | **INCORRECT** | Percolator-style learner with material initialization, fold-merge, selection, q, and PEP differences |
| Default nested/leakage-free 3-fold CV | **INCORRECT** | loss is out of fold, but global normalization and supervised initialization leak; only base `--auto-model` isolation is substantially nested |
| Correct/calibrated target-decoy q-values | **INCORRECT / FAILED VALIDATION** | invalid imbalance factor, missing pseudocount/ties, 17/30 complete-null rejection rate |
| PEPs are posterior error probabilities | **INCORRECT / FAILED VALIDATION** | unsupported transform; 87 known-false PSMs printed at PEP=0; large calibration gap |
| Pure-null behavior is strongly conservative | **INCORRECT** | old denominator was wrong; complete-null FDR is probability of any call, empirically 56.7% |
| +3.9% PSMs / +4.5% peptides in headline run | **SUPPORTED WITH CAVEATS** as counts | exact reported-q counts, but C++ used mix-max and neither threshold nor model was matched |
| More IDs mean improved accuracy/sensitivity | **MISLEADING / INCORRECT** | no complete PSM truth; null and entrapment show anti-conservative statistical tails |
| Rust always yields more | **INCORRECT** | yeast seed 1 favors C++; seed and dataset interaction exists |
| Results are deterministic for fixed seed | **SUPPORTED** | regression gates and repeated exact runs; this is software reproducibility, not validity |
| Multi-dataset format compatibility | **SUPPORTED** | five distinct workflows run successfully after documented normalization |
| Scientific generalization | **INSUFFICIENT EVIDENCE** | calibration tested on one PSM search design; compact cases lack truth; DIA absent |
| Picked/Bayesian protein inference is calibrated | **INCORRECT / FAILED VALIDATION** | PrEST q<=0.01 FDP failures and false blank groups; Bayesian mode inherits invalid peptide PEPs |
| Speed and RSS gains | **SUPPORTED WITH CAVEATS** | reproducible on the named host/input; irrelevant to scientific correctness and workload-sensitive |

Historical documentation was also corrected where it called the old null study conservative or
called identification loss an accuracy loss. Performance tables remain counts, but surrounding
language now explicitly separates yield from accuracy.

## 16. Remaining uncertainties

1. Complete PSM ground truth, recall, and false-negative counts were unavailable. Entrapment is
   partial truth and its search-space correction is estimated.
2. Only one underlying signal-present PSM search design has direct calibration. Tide, MSFragger,
   Sage, and yeast comparisons establish agreement and variability, not calibration.
3. Five model seeds are adequate to reveal major instability and summarize routine variability,
   but not rare-tail probabilities. The 30 null relabelings are decisive against 1% because the
   effect is enormous, yet remain clustered within three source PINs.
4. Post-hoc q ablations use six-decimal printed scores and cannot replace an end-to-end,
   full-precision, independently validated implementation experiment.
5. No one-factor end-to-end matrix covers initialization, normalization, fold scaling, C grid,
   q-value method, and PEP method across every dataset/seed. Existing evidence is sufficient to
   reject current validity, not to assign an exact percentage of the yield gap to each component.
6. PEP calibration bins estimate false native matches using foreign-decoy fractions; sparse bins and
   PSM dependence limit inferential precision.
7. Reference agreement is not a correctness oracle. C++ itself is anti-conservative on the tested
   entrapment curve.
8. DIA, additional acquisition regimes, and independent signal-present calibration datasets remain
   untested.
9. Parser safety failures were documented, not repaired, because the mandate prohibited changing
   methodology in response to failed validation.

## Final verdicts

| Dimension | Verdict | Conclusion |
|---|---|---|
| **IMPLEMENTATION CORRECTNESS** | **FAILED VALIDATION** | Percolator-style but materially non-equivalent; default CV leaks and statistics are defective |
| **FDR VALIDITY** | **FAILED VALIDATION** | 56.7% complete-null rejection probability and anti-conservative entrapment curve |
| **PEP VALIDITY** | **FAILED VALIDATION** | unsupported formula, known-false zero PEPs, and substantial empirical miscalibration |
| **REPRODUCIBILITY** | **STRONG EVIDENCE** | deterministic seeded behavior and machine-readable, hashed experiment manifests |
| **REFERENCE AGREEMENT** | **MODERATE EVIDENCE** | high on Tide/Sage, weaker on yeast/MSFragger, severe file-level PXD032157 divergence |
| **GENERALIZATION** | **WEAK EVIDENCE** | schema compatibility and some ranking behavior generalize; calibration does not |
| **IDENTIFICATION YIELD** | **MODERATE EVIDENCE** | count differences are real and often seed-stable, but dataset-dependent and not accuracy evidence |

## What percolator-rs can and cannot scientifically claim

percolator-rs can claim that it is a deterministic, high-performance, Percolator-style
semi-supervised rescoring implementation; that it produces reproducible reported-q PSM and peptide
lists; that its rank agreement with C++ 3.09 is high on some workflows and materially lower on
others; and that its reported-q yield is often, but not always, higher under the tested unmatched
model/statistical procedures.

It cannot claim that the default is faithful to C++ Percolator, that its cross-validation is
leakage-free, that its q-values control FDR, that q<0.01 means 1% empirical error, that its PEPs are
validated posterior probabilities, that its extra identifications improve accuracy or sensitivity,
that its protein probabilities are calibrated, or that scientific validity generalizes across
datasets. Those claims require corrected methods followed by new, independent null,
signal-present, PEP, protein, and multi-dataset calibration; the current failed experiments must
remain part of that record.

## Method references

- Käll et al., Percolator cross-validation and score normalization:
  <https://pmc.ncbi.nlm.nih.gov/articles/PMC3489528/>
- Official Percolator source inspected for the 3.09 comparison:
  <https://github.com/percolator/percolator>
- Emery et al., target-decoy competition and the final-decoy `+1` correction:
  <https://pmc.ncbi.nlm.nih.gov/articles/PMC6919216/>
- qvality PEP/q-value methodology: <https://noble.gs.washington.edu/proj/qvality/>
