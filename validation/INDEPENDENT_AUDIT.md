# Independent adversarial scientific audit of `percolator-rs`

Audit date: 2026-08-25  
Repository reviewed: `b38c0db` (`main`)  
Frozen repaired executable: commit `1348b0f4696c407af2dea17b1cf97935c35c1f19`, SHA-256 `d233c18c9c74915be4789e6754af37dfc64a55be484da440f742851651152a2e`  
Reviewer posture: independent falsification, not repair  
Production statistical code modified during this audit: **no**

The machine-readable synopsis is in `validation/INDEPENDENT_AUDIT_RESULTS.json`. The minimal attacks are `validation/independent_stats_probe.rs`, `validation/adversarial_competition.py`, `validation/adversarial_cv.py`, `validation/adversarial_parser.py`, and `validation/adversarial_protein.rs`. Raw prior and repaired study artifacts are under `/home/andrej-rumenovski/percolator_rs_out/scientific-validation`.

## 1. Executive summary

The repair corrected several real defects. The standalone reported-q scan now has the finite-sample `+1`, opportunity ratio, exact numeric score grouping, bounds, and reverse cumulative minimum expected from its stated formula. The canonical fixed-hyperparameter SVM path fits normalization, retention-time alignment, initial direction, positive selection, and weights only on each outer training partition. The 30 predefined post-repair complete-null runs now yield zero discoveries, and all frozen multi-seed outputs reproduce exactly.

The implementation nevertheless **fails scientific validation as a whole**. A valid exact-score-tie input makes PSM target/decoy competition depend entirely on row order. With 200 otherwise identical target/decoy pairs, target-first order reports all 200 targets at `q < 0.01` (minimum q 0.005); reversing each pair reports zero target winners. This violates the fair null-win assumption before the mathematically correct q scan even begins. The existing unit test explicitly requires this input-order behavior.

The repaired PEP procedure removes exact zeros but is not independently justified as a posterior-probability estimator. On entrapment data every populated repaired-PEP bin is optimistic, with weighted absolute and signed errors both 0.0185. Published PEP methodology warns that valid cumulative FDR estimates cannot simply be differentiated into valid PEPs. Two large probability-changing mutations pass all 82 tests.

The optional protein output is not scientifically defensible: connected proteins are called indistinguishable even when their peptide neighborhoods differ; exact target/decoy protein-score ties always select the target; and picked output copies the best peptide PEP into the column named `posterior_error_prob`. These are methodology failures, not merely missing validation.

Cross-validation isolation is strong only for the canonical fixed-C path (and appears correctly nested for `--auto-model`). It is false as an unqualified program-wide claim: `--select-c` uses the same out-of-fold labels both to select hyperparameters and to report predictions, and `--ensemble` constructs a label-keyed global feature before folds exist.

The entrapment curve remains anti-conservative at every predefined threshold: mean adjusted FDP is 0.00615 at nominal 0.001, 0.01816 at nominal 0.01, and 0.10873 at nominal 0.10. This falsifies empirical calibration in this experiment, but does not uniquely identify a Rust coding bug; matched C++ post-processing behaves similarly, and TDC/entrapment/search-space assumptions may be responsible.

Paper-level answer: I would accept claims of deterministic Rust implementation, exact reproducibility, the standalone q-scan formula conditional on valid competition input, canonical fixed-C fold isolation, broad PSM/peptide reference compatibility, and the named throughput measurement. I would require explicit qualification of complete-null improvement, identification counts, cross-dataset agreement, and all statements scoped to the canonical path. I would reject claims of generally valid or calibrated reported q-values/FDR, calibrated PEPs, leakage-free behavior across supported modes, scientifically valid protein q-values/PEPs/grouping, increased accuracy or sensitivity, generalization, or statistical equivalence to C++.

## 2. Critical findings

### CRITICAL

**C1 — PSM competition has a target/decoy row-order attack.** `src/main.rs:331-347` retains the earlier row when scores tie (`score[previous] >= score[i]`). A 200-spectrum fixture with one exactly tied target/decoy candidate per spectrum changes from 200 target winners and 200 discoveries at `q < 0.01` to zero target winners solely by reversing row order. Root cause: deterministic first-row tie resolution is not a fair or conservative target/decoy competition rule. Affected claims: reported q-value validity, FDR control, permutation robustness, and implementation correctness. The behavior is certified, not caught, by `src/main.rs:1271-1277`.

**C2 — Picked-protein exact ties are target-favoring.** `src/protein.rs:583-590` uses target `>=` decoy. The analogous 200-pair attack produces 200 target wins, all at `q < 0.01`, minimum q 0.005. Affected claims: picked-protein q-values, protein FDR, and protein sensitivity.

### MAJOR

**M1 — Repaired PEPs are not calibrated posterior probabilities.** The code differences cumulative q-values, adds an arbitrary `0.5/N` to every target, and applies PAVA (`src/stats.rs:491-512`). In the repaired competed entrapment result, observed adjusted false fractions exceed reported PEP in every populated bin; ECE and signed observed-minus-reported error are both 0.0185144. Removing zeros is an improvement, not validation.

**M2 — Protein grouping merges distinguishable proteins.** `src/protein.rs:373-386` unions any proteins connected by any shared peptide. Proteins A and B supported by shared peptide AB, with A also having a unique peptide, are emitted as one group even though their peptide evidence sets differ. Connected-component grouping is not indistinguishability grouping.

**M3 — Picked output reports peptide PEP as protein PEP.** `src/protein.rs:395-396,420-425,455-462` stores and emits the best peptide's PEP unchanged. A fixture with peptide PEP 0.123456 produces protein `posterior_error_prob` 0.123456. No protein posterior is estimated.

**M4 — The signal-present entrapment curve is anti-conservative.** Across five predefined seeds, nominal q 0.01 gives mean adjusted FDP 0.01816 (range 0.01730–0.01884); every predefined threshold has mean FDP above nominal. This is an empirical calibration failure. Similar C++ behavior prevents attribution to a Rust-only bug but does not rescue validity.

**M5 — Supported CV modes retain label leakage.** `--select-c` changes selected class weights from 4:1 to 0.25:1 and changes all 200 attacked held-out scores when only that held-out fold's labels are flipped. `--ensemble` creates a whole-dataset feature keyed by `(ScanNr, Label, Peptide)` before outer folds. Unqualified “leakage-free” wording is false.

### MODERATE

**O1 — Fold-score calibration is isolated but not reference-equivalent or empirically validated.** Rust maps held-out margins using each training fold's decoy mean and standard deviation. Official Percolator instead uses held-out-fold anchors (selection boundary and median decoy). Rust's map does not reorder within a fold, but it changes cross-fold interleaving and therefore global discoveries.

**O2 — The complete-null study is underpowered for small nominal error rates and misses structured ties.** Zero failures in 30 runs is encouraging, but the two-sided exact 95% upper bound on a per-run failure probability is about 0.116. Random relabeling does not exercise the deterministic target-first tie attack.

**O3 — Input validation is finite-value fail-closed but label-permissive.** NaN/Inf features, textual labels, and malformed scans are rejected. Labels `0`, `2`, and `-2` are accepted and mapped by sign, contrary to a strict PIN `+1/-1` contract.

### MINOR

**N1 — Documentation contains stale or overstrong wording.** The README says `--rt-features` remains outside fold isolation even though code and tests now refit it fold-locally; it calls the repaired configuration “validated”; it describes the PEP transformation as a published identity; and it presents protein picked-yield as an expected benefit despite no truth-based protein validation. `validation/REPAIR.md` calls 0.615% adjusted FDP at nominal 0.1% “conservative,” which is numerically incorrect.

### INFORMATIONAL

**I1 — Several repairs are real and reproducible.** The standalone q scan passes independent hand oracles; default fold-local preprocessing resists held-out corruption; the complete-null result improves from 17/30 to 0/30 runs with discoveries; current and frozen repaired binaries hash identically; and 40 implementation/dataset/seed output sets reproduce byte-for-byte.

## 3. Mathematical implementation audit

The reconstruction below began from current source rather than the repair report.

| Procedure | Intended method | Actual code and formula | Required assumptions | Principal failure modes |
|---|---|---|---|---|
| PSM target-decoy competition | One fair target/decoy winner per precursor | `competition_winners`; key `(source, ScanNr, exact ExpMass bits)`, maximum score, first row wins exact tie | Competition units are correct; null target and decoy are exchangeable; ties are fair/conservative | **Exact ties are label-order biased**; exact mass-bit grouping can split near-equal precursor representations; charge is not separately keyed |
| Reported raw FDP | Finite-sample TDC+ | `FDP_k = min(1, pi0 * lambda * (D_k+1)/max(T_k,1))`, `pi0=1`, `lambda=p/(1-p)` | Declared `p` is the actual null target-win probability; decoys model incorrect targets | Invalid `p` silently falls back to 1; competition/dependence violations invalidate the estimate |
| Training target selection | Heuristic positives for semi-supervised fitting | Same scan but without `+1`; targets below training q threshold selected | No error-rate claim is made for training q | Small-fold selection can be unstable; feedback from model ranking |
| Score ties in q scan | All equal scores share one rejection boundary | Exact numeric `==` grouping after descending total order | Exact equality captures meaningful ties | NaNs form separate groups; upstream competition has already discarded tied losers incorrectly |
| q-values | Minimum attainable estimated FDP over lower cutoffs | `q_i = min_{j>=i} FDP_j`, reverse cumulative minimum, clamped to `[0,1]` | Input winner list and TDC assumptions valid | Correct scan cannot repair biased winners; internal NaN changes a finite prefix's q-values |
| CV splitting | Three spectrum-grouped outer folds | Deterministic balanced assignment by `(source,scan)`; ensemble uses scan | Group key captures all correlated candidates | Ensemble files can reuse scan numbers; global ensemble feature precedes folds |
| Preprocessing | Fit transformations on outer training rows only | Per-fold z-score matrix and optional RT residual alignment fit on training indices | Test-time transformation does not use held-out labels/statistics | Holds on canonical path; named pre-fold/global features can violate it |
| Initial direction | Choose the best signed feature within training data | Training-only rank/yield search in each outer fold | Training heuristic does not consult held-out labels | Canonical attack passes; global variants would leak and mutation test detects it |
| Semi-supervised loop | Iterate positives/decoys and fit a linear classifier | Training-only q heuristic; confident targets vs all decoys; fixed number of iterations | Decoys are suitable negatives; selected positives are sufficiently clean | Selection bias and TDC assumptions; fixed iteration convergence is not a guarantee |
| SVM | Exact L2-regularized squared-hinge optimization | `0.5||w||^2 + sum_i C_i max(0,1-y_i w'x_i)^2`; explicit active-set Hessian, Cholesky Newton step, line search | Solver reaches an adequate optimum; finite matrix | Different solver/convergence path from C++; fallback direction and tolerance can alter ranking |
| Held-out scoring/merge | Out-of-fold predictions on a comparable scale | Held-out margin standardized by **training-decoy** mean and SD | Training-decoy location/scale is stable and comparable across folds | No leakage in default, but cross-fold ordering differs from official held-out anchoring |
| Peptide statistics | One best PSM per peptide/label, then new error estimates | Best score per `(label, core peptide)`, recompute same q and PEP procedures | Peptide-level target/decoy construction is valid; PSM score is appropriate peptide score | Target and decoy forms are selected separately; all upstream competition/PEP issues propagate |
| PEP | Local posterior error probability monotone in score | For target ranks `k`: `r_k=max(0,kq_k-(k-1)q_(k-1))+0.5/N`; unit-weight PAVA, clamp; interpolate decoy outputs | q curve is accurate and differentiable enough; prior is justified; isotonic fit is calibrated | Valid FDR does not imply valid derivative PEP; prior breaks the stated cumulative identity; small/tied samples; observed optimism |
| Picked protein | Indistinguishable groups, paired target/decoy picking, protein-level error estimates | Union-find connected components; best peptide score; target wins exact tie; TDC q; best peptide PEP copied | Groups really are indistinguishable; correct pairing; fair ties; protein TDC assumptions | All four assumptions fail or are unvalidated in adversarial fixtures |
| Bayesian protein | Posterior presence inference from peptide PEPs | Configurable factor graph/noisy-OR-like likelihood; exact tree or loopy BP; cumulative mean posterior PEP gives q | Input peptide PEPs calibrated; generative parameters valid; BP converged | Input PEP fails calibration; defaults unvalidated; cyclic nonconvergence still yields output |

The reported q formula agrees structurally with official Percolator's `PosteriorEstimator::getQValues`, including the `+1` outside training and `p/(1-p)` factor. The official implementation is a compatibility reference, not a proof. The relevant source is [official `PosteriorEstimator.cpp` at `rel-3-09`](https://github.com/percolator/percolator/blob/rel-3-09/src/PosteriorEstimator.cpp).

The primary TDC+ literature supports a plus-one finite-sample correction only under exchangeability/competition assumptions; it does not license target-favoring exact ties ([TDC+ analysis](https://pmc.ncbi.nlm.nih.gov/articles/PMC6919216/)). A primary treatment explicitly resolves target/decoy score ties by a fair coin ([Granholm et al., 2011](https://pmc.ncbi.nlm.nih.gov/articles/PMC3220955/)). Percolator's CV paper describes three train/test-separated folds and held-out score normalization, but its normalization is not Rust's training-decoy map ([Spivak et al., 2009/2012 full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC3489528/); [official `Scores.cpp`](https://github.com/percolator/percolator/blob/rel-3-09/src/Scores.cpp)).

For PEP, the decisive literature point is directional: calibrated PEPs imply cumulative FDR by averaging, but calibrated FDR does not imply calibrated PEPs. Differentiating an empirical FDR curve is high variance, motivating nonparametric logistic estimation in QVALITY ([Käll et al., 2008](https://noble.gs.washington.edu/papers/kall2008nonparametric.pdf); [Käll et al., 2009](https://noble.gs.washington.edu/papers/kall2009qvality.pdf)). Rust's raw derivative and `0.5/N` resemble official `IsotonicPEP.cpp`, but official 3.09 normally uses an I-spline rather than Rust's PAVA, and source similarity is not a calibration result ([official `IsotonicPEP.cpp`](https://github.com/percolator/percolator/blob/rel-3-09/src/IsotonicPEP.cpp)).

## 4. CV leakage results

### Canonical fixed-C path

For an attacked outer fold, changing only held-out labels, feature values, target/decoy composition, order, and extreme outliers left the fold's training normalization, RT alignment, initial direction, and fitted weights bit-identical. In the explicit 200-row label-flip fixture, fixed-C remained 1:4 and 0/200 held-out scores changed. Existing internal tests also corrupt held-out data and compare fitted objects directly. Mutating either normalization or initial direction back to global fitting makes the isolation test fail.

Information flow on this path is:

| Stage | Training information | Held-out information | Assessment |
|---|---|---|---|
| Fold assignment | Spectrum/source keys, row counts, seed | Unsupervised grouping keys | Acceptable |
| Normalization and RT alignment | Outer training rows and labels needed for RT fit | Transformed only after fit | Isolated |
| Initial direction and positive selection | Outer training features/labels | None | Isolated |
| SVM fit | Outer training matrix/labels | None | Isolated |
| Held-out transform/prediction | Stored training transforms/model | Held-out features | Isolated |
| Fold merge | Training-decoy mean/SD plus predictions | No held-out labels | Isolated, method unvalidated |

Therefore the narrowly scoped claim “the canonical fixed-C outer-fold model does not train on its held-out fold” has strong evidence.

### Modes that fail the unqualified claim

| Mode | Attack/result | Conclusion |
|---|---|---|
| `--select-c` | Flipping only outer fold 0 labels changed selected weights 4:1 → 0.25:1 and all 200 fold-0 scores | Non-nested selection leakage by construction |
| `--ensemble` | Agreement feature built globally with key `(ScanNr, Label, Peptide)` before folds | Direct held-out-label-derived feature leakage |
| `--auto-model` | Source inspection shows candidate choice uses inner folds inside outer training data | Appears nested; not attacked as extensively as fixed-C |
| `--rt-features` | Both code and mutation-sensitive tests fit alignment within each fold | Current README statement that it remains outside isolation is stale |

Verdict interpretation: default fixed-C isolation is **STRONG EVIDENCE**, but the program-wide dimension is **FAILED VALIDATION** because supported named modes violate the property.

## 5. q-value/TDC adversarial tests

### Independent hand oracle

The standalone probe imports current `src/stats.rs` and compares to separately enumerated expected vectors. All finite hand cases passed.

| Case | Expected/reproduced result at default `p=0.5` |
|---|---|
| Empty | Empty |
| One target | q = 1 |
| One decoy | q = 1 |
| Target before decoy / decoy before target | Bounded q = 1 |
| Four targets, no decoys | q = 0.25 each |
| All decoys | q = 1 |
| `T,T,D,T,D,T` | `0.5,0.5,2/3,2/3,0.75,0.75` |
| Mixed target/decoy exact ties | Group-shared q = 2/3 |
| All identical scores | q = 1 |
| Extreme finite scores | Finite, bounded q = 1 in fixture |
| Imbalanced opportunity, `p=1/3` | First five q = 0.1; final decoy q = 0.2 |

These results independently verify the implemented formula, finite-sample correction, tie grouping within the q scan, counts, bounds, and monotonicity. No unjustified finite-input q=0 was reproduced.

### Attacks that break reported q-values

The structured competition fixture is decisive:

| Exact tied-pair order | Target winners | Decoy winners | Targets `q<0.01` | Minimum target q |
|---|---:|---:|---:|---:|
| Target row first | 200 | 0 | 200 | 0.005 |
| Decoy row first | 0 | 200 | 0 | — |

The score and feature multiset is identical. Only within-pair row order changes. This is not a q-scan tie-grouping bug; it is biased winner construction that invalidates the q scan's null exchangeability assumption.

An internal API boundary also remains unsafe: scores/labels `[finite T, finite D, finite T]` give q `[1,1,1]`; appending one trailing NaN target changes all four q-values to `2/3`. The PIN parser rejects NaN and Inf, so ordinary CLI input cannot trigger this directly, but `stats` itself does not fail closed.

## 6. Complete-null replication

I reran the predefined 30 current Rust null inputs without redesigning seeds or thresholds. All 30 target output hashes exactly matched the stored repaired artifacts.

| Implementation/stage | Completed runs | Runs with any discovery at each of q 0.001, 0.005, 0.01, 0.02, 0.05, 0.10 | Count distribution |
|---|---:|---:|---|
| Pre-repair Rust | 30 | 17/30 at every threshold | mean 1.37, maximum 9 |
| Post-repair Rust | 30 | 0/30 at every threshold | all zero |
| C++ | 29 | 0/29 at every threshold | all zero |

Four post-hoc design-variant arms (20 runs each) also produced zero discoveries at all thresholds. The two imbalanced arms had no power to challenge small q-values: their minimum q-values were 0.276 and 0.552.

This is strong evidence that the original complete-null pathology was repaired for these fixtures. It is not evidence for calibration at q=0.001 or q=0.01: with 0/30 run-level failures, the two-sided Clopper–Pearson 95% upper bound is approximately 0.116. It also does not cover exact tied competitions. The experiment tests exchangeable random relabeling conditional on stored data; its seeds are relabel replicates, not independent biological experiments.

## 7. Entrapment replication

All five predefined Rust seeds over the six stored datasets were rerun. Per-seed q=0.01 results exactly reproduced:

| Seed | Accepted targets | Pure entrapment | Effective entrapment fraction | Adjusted FDP |
|---:|---:|---:|---:|---:|
| 1 | 19,545 | 258 | 0.738889 | 0.0178651 |
| 2 | 19,488 | 244 | 0.723757 | 0.0172994 |
| 3 | 19,532 | 262 | 0.727778 | 0.0184313 |
| 4 | 19,514 | 262 | 0.712707 | 0.0188384 |
| 5 | 19,586 | 267 | 0.741758 | 0.0183782 |
| Aggregate | mean 19,533; median 19,532; SD 36.54; range 19,488–19,586 | mean 258.6 | — | mean 0.0181625; median 0.0183782; SD 0.0005935; range 0.0172994–0.0188384 |

The entire predefined Rust curve is anti-conservative:

| Nominal strict q threshold | Mean accepted targets | Mean adjusted FDP | FDP / nominal |
|---:|---:|---:|---:|
| 0.001 | 12,946.4 | 0.006146 | 6.15× |
| 0.005 | 17,982.4 | 0.012070 | 2.41× |
| 0.010 | 19,533.0 | 0.018162 | 1.82× |
| 0.020 | 21,139.6 | 0.027531 | 1.38× |
| 0.050 | 23,823.8 | 0.058482 | 1.17× |
| 0.100 | 26,846.2 | 0.108726 | 1.09× |

Pre-repair Rust at q=0.01 was worse: mean 0.0269832. Repaired Rust without post-processing competition gives 0.0256121. With matched competition, C++ gives mean 0.0172929 at q=0.01 and follows a similar anti-conservative curve (0.006445, 0.012775, 0.017293, 0.027170, 0.059862, 0.108434).

The adjusted estimator is `pure entrapment targets / effective foreign fraction / accepted targets`; the effective fraction is estimated from accepted nonmixed decoys. It assumes incorrect native targets fall into the foreign database at the same rate as decoys. Approximate per-seed q=0.01 intervals, even treating PSMs as independent, are roughly 0.015–0.022 and remain above nominal; that approximation is optimistic because PSMs, proteins, runs, and seeds are dependent.

Classification: empirical calibration **fails**. Attribution remains unresolved among shared TDC assumptions, entrapment adjustment assumptions, search-space effects, dependence, competition structure, and training effects. C++ similarity rules out a large Rust-only divergence in this design; it does not prove either implementation valid.

## 8. PEP calibration results

### Algorithmic audit

The repaired estimator eliminates exact zero by adding `0.5/N`, but this prior is not derived as a posterior model. For five target scores followed by one decoy, the q curve is 0.2 on targets and 0.4 on the decoy, while every target PEP becomes 0.3. Thus the asserted identity “mean PEP through rank k equals q at k” is not preserved after the prior. More fundamentally, even preserving it would not prove local posterior calibration.

Small/all-null/tied cases remain bounded and monotone due to PAVA/clamps, but bounded monotonic output is weaker than probability validity. Decoy PEPs are interpolated for display and explicitly carry no error-rate claim.

### Entrapment calibration

| Build/path | Targets | PEP=0 | Known-false PEP=0 | Weighted absolute calibration error | Signed observed−reported |
|---|---:|---:|---:|---:|---:|
| Pre-repair Rust | 435,261 | 9,155 | 87 | 0.063361 | +0.061933 |
| Repaired Rust, competed | 98,680 | 0 | 0 | 0.018514 | +0.018514 |
| Repaired Rust, no competition | 435,261 | 0 | 0 | 0.017150 | positive |
| C++, no competition | 435,261 | 0 | 0 | 0.017756 | +0.016014 |

Every populated repaired-competed bin is optimistic:

| PEP bin | N | Mean reported | Adjusted observed false fraction | Observed−reported |
|---|---:|---:|---:|---:|
| [0.0001,0.001) | 10,926 | 0.000387 | 0.003844 | +0.003457 |
| [0.001,0.005) | 2,978 | 0.003735 | 0.009850 | +0.006115 |
| [0.005,0.01) | 1,399 | 0.007166 | 0.016542 | +0.009377 |
| [0.01,0.02) | 1,703 | 0.014114 | 0.030143 | +0.016029 |
| [0.02,0.05) | 1,326 | 0.035477 | 0.060210 | +0.024733 |
| [0.05,0.10) | 1,514 | 0.072017 | 0.081345 | +0.009327 |
| [0.10,0.20) | 1,570 | 0.158629 | 0.190906 | +0.032277 |
| [0.20,0.50) | 3,796 | 0.343283 | 0.359573 | +0.016290 |
| [0.50,1.00] | 73,468 | 0.913545 | 0.934931 | +0.021386 |

The prior pathology “PEP exactly zero for known false PSMs” is fixed. The stronger claim “PEP is a calibrated posterior error probability” fails.

## 9. Mutation-test results

All mutations were performed in isolated worktrees and removed afterward. Production source remained unchanged.

| Deliberate defect | Existing suite response | Audit interpretation |
|---|---|---|
| Remove reported `+1` decoy | 7 tests fail | Good protection |
| End every score tie after one row | 1 test fails | Some protection; surprisingly narrow |
| Remove reverse cumulative minimum | 7 tests fail | Good protection |
| Reintroduce global normalization | Held-out isolation test fails | Good protection |
| Reintroduce global initial direction | Held-out isolation test fails | Good protection |
| Double PEP prior from `0.5/N` to `1/N` | All 82 tests pass | Serious probability behavior unconstrained |
| Raise PEP floor from `1e-12` to `0.01` | All 82 tests pass | Serious probability behavior unconstrained |
| Keep target-first exact competition ties | All tests pass; one test requires it | Critical invalid behavior encoded as expected |

The suite is effective for the repaired q arithmetic and canonical CV isolation. It is not a scientific validation suite for PEP, competition validity, or protein inference.

## 10. Multi-seed results

I reran all 20 current Rust runs (four datasets × five predefined seeds) and all 20 matched C++ runs. For each implementation, all four output hashes per run match the stored manifest. No seed was excluded. Protein counts were not produced by this canonical study and should not be retrofitted now, particularly given the independent protein failures.

PSM results at strict q < 0.01:

| Dataset | Implementation | Counts, seeds 1–5 | Mean | Median | SD | Min–max |
|---|---|---|---:|---:|---:|---:|
| MSFragger | Rust | 1367, 1377, 1312, 1361, 1340 | 1351.4 | 1361 | 25.85 | 1312–1377 |
| MSFragger | C++ | 1321, 1403, 1357, 1358, 1332 | 1354.2 | 1357 | 31.62 | 1321–1403 |
| Sage | Rust | 25806, 25803, 25804, 25791, 25779 | 25796.6 | 25803 | 11.46 | 25779–25806 |
| Sage | C++ | 25787, 25780, 25776, 25782, 25782 | 25781.4 | 25782 | 3.97 | 25776–25787 |
| Tide | Rust | 27640, 27665, 27616, 27624, 27624 | 27633.8 | 27624 | 19.50 | 27616–27665 |
| Tide | C++ | 27617, 27656, 27610, 27612, 27670 | 27633.0 | 27617 | 27.95 | 27610–27670 |
| Yeast | Rust | 1150, 1144, 1124, 1180, 1159 | 1151.4 | 1150 | 20.51 | 1124–1180 |
| Yeast | C++ | 1182, 1192, 1055, 1133, 1118 | 1136.0 | 1133 | 55.10 | 1055–1192 |

Peptide results at strict q < 0.01:

| Dataset | Implementation | Counts, seeds 1–5 | Mean | Median | SD | Min–max |
|---|---|---|---:|---:|---:|---:|
| MSFragger | Rust | 1060, 1067, 1030, 1065, 1056 | 1055.6 | 1060 | 14.94 | 1030–1067 |
| MSFragger | C++ | 1043, 1104, 1051, 1058, 1036 | 1058.4 | 1051 | 26.80 | 1036–1104 |
| Sage | Rust | 11245, 11253, 11247, 11246, 11253 | 11248.8 | 11247 | 3.90 | 11245–11253 |
| Sage | C++ | 11336, 11317, 11322, 11327, 11317 | 11323.8 | 11322 | 7.98 | 11317–11336 |
| Tide | Rust | 19736, 19713, 19715, 19701, 19735 | 19720.0 | 19715 | 15.13 | 19701–19736 |
| Tide | C++ | 19722, 19736, 19745, 19737, 19722 | 19732.4 | 19736 | 10.11 | 19722–19745 |
| Yeast | Rust | 904, 908, 907, 954, 935 | 921.6 | 908 | 22.01 | 904–954 |
| Yeast | C++ | 943, 943, 863, 906, 908 | 912.6 | 908 | 33.07 | 863–943 |

The canonical seed 1 is not uniformly favorable: it is neither the maximum Rust PSM yield for MSFragger, Tide, nor yeast, and the Sage spread is negligible. The data do not show favorable seed exclusion. They do show greater seed instability for MSFragger and yeast than for Sage/Tide.

## 11. Cross-dataset results

The four datasets cover different search-engine/PIN characteristics. At q=0.01, mean Rust/C++ Jaccard is 0.9235 (MSFragger), 0.9962 (Sage), 0.9930 (Tide), and 0.9201 (yeast). Sage and Tide are strong agreement cases; MSFragger and yeast are materially weaker.

Mean Jaccard over five seeds across all predefined thresholds:

| Dataset | q=.001 | q=.005 | q=.01 | q=.02 | q=.05 | q=.10 |
|---|---:|---:|---:|---:|---:|---:|
| MSFragger | 0.6599 | 0.8955 | 0.9235 | 0.9477 | 0.9319 | 0.9092 |
| Sage | 0.9941 | 0.9958 | 0.9962 | 0.9968 | 0.9949 | 0.9951 |
| Tide | 0.9879 | 0.9922 | 0.9930 | 0.9917 | 0.9900 | 0.9874 |
| Yeast | undefined/empty in at least one run | 0.9102 | 0.9201 | 0.9263 | 0.9147 | 0.9116 |

Compatibility is broad enough to accept “reads and produces plausible ranked PSM/peptide results on these four input families.” It is not a generalization or calibration experiment: there are only four curated datasets, all were available during repair, labels are target/decoy rather than ground truth, and the entrapment data already reject nominal calibration.

The forced-concatenated PXD032157 file-0001 artifact must not be presented as post-repair evidence. It is a historical pre-repair/unmatched run: 183,982 matching rows, Rust 208 vs C++ 0 discoveries at q=0.01, with score Pearson/Spearman 0.539/0.471, q Pearson/Spearman 0.607/0.464, and PEP Pearson/Spearman 0.563/0.414. It documents prior disagreement, not current corrected behavior. PXD032157 was also the development/performance dataset, so it is not untouched external validation.

No dataset-specific truth labels establish accuracy. The entrapment study pools six data sources and reveals anti-conservatism, but its current summaries do not isolate calibration separately by search engine/dataset. That remains a material gap.

## 12. Rust/C++ disagreement analysis

### Correlation and accepted sets

Mean correlations across five seeds, calculated only on matched PSMs:

| Dataset | Score Pearson / Spearman | q Pearson / Spearman | PEP Pearson / Spearman | q=.01 mean Jaccard |
|---|---:|---:|---:|---:|
| MSFragger | 0.9457 / 0.9627 | 0.9324 / 0.9570 | 0.9608 / 0.9492 | 0.9235 |
| Sage | 0.9966 / 0.9955 | 0.9953 / 0.9728 | 0.9973 / 0.9408 | 0.9962 |
| Tide | 0.9986 / 0.9988 | 0.9975 / 0.9966 | 0.9984 / 0.9825 | 0.9930 |
| Yeast | 0.9793 / 0.9498 | 0.9599 / 0.9496 | 0.9747 / 0.9325 | 0.9201 |

High correlation does not imply equal accepted sets or valid probabilities. Across five seeds at q=0.01, exclusive-set tracing gives:

| Dataset | Rust-only / C++-only | Counterpart q<.02 | Counterpart q≥.05 | Median absolute normalized-rank gap |
|---|---:|---:|---:|---:|
| MSFragger | 250 / 266 | 195 / 239 | 8 / 2 | 0.0495 / 0.0473 |
| Sage | 283 / 207 | 280 / 197 | 1 / 3 | 0.0058 / 0.0069 |
| Tide | 489 / 485 | 485 / 471 | 0 / 0 | 0.0073 / 0.0065 |
| Yeast | 242 / 175 | 196 / 134 | 3 / 5 | 0.0142 / 0.0173 |

Sage and Tide disagreements are overwhelmingly near the cutoff. MSFragger has materially larger ranking differences; yeast is intermediate. Representative exclusives were traced through matched identifiers: most have a counterpart just beyond q=0.01, while the small q≥0.05 subsets are true ranking/model disagreements rather than output-ID mismatches.

### Likely sources

The matched comparison uses equivalent post-processing competition, so its remaining differences precede the final q scan. Source inspection identifies legitimate methodological differences: Rust uses fixed 1:4 class weights by default; independent fold construction/RNG; an explicit squared-hinge Newton solver; training-decoy mean/SD fold scaling; and PAVA-derived PEPs. Official C++ uses different optimization/convergence details and held-out score anchoring. q and PEP transformations then magnify ranking differences around sparse decoy boundaries.

The audit cannot uniquely assign every exclusive PSM to one stage because the stored C++ artifacts do not expose fold-local normalized rows, initial directions, weights, or intermediate margins. The evidence supports “broad compatibility/agreement,” not bitwise, algorithmic, statistical, or probability equivalence.

## 13. Validation-design audit

### Complete null

The design is a legitimate conditional exchangeable-label complete null: labels are randomized on stored decoy-derived data and every pseudo-target is false. Under this definition, FDP is one whenever any discovery occurs. Thresholds and primary seeds were predefined. Weaknesses are low replication for small nominal rates, dependence inherited from spectra/candidates, lack of biological independence, and failure to preserve/target structured target/decoy score ties. Variant arms added after primary results are sensitivity analyses, not confirmatory evidence.

### Entrapment

Pure foreign-protein targets are credible known-false calls for the native samples; mixed native/foreign mappings are correctly excluded. The adjustment by the accepted-decoy foreign fraction is reasonable only if incorrect native targets and decoys have the same foreign-placement probability. Search-space composition, homologous peptides, protein/PSM dependence, model training, and competition can violate that condition. Only five relabel/search seeds share six biological datasets, so seed dispersion is not a biological confidence interval. The validation reports lack a predeclared uncertainty model.

Consequently:

- The observed anti-conservative curve is a **validation failure of the scientific FDR claim**.
- It is **not sufficient to diagnose a Rust implementation defect** because matched C++ shows nearly the same behavior and the adjustment's assumptions are not proven.
- Describing q=.001 FDP 0.00615 as “conservative” is incorrect; it is lower than other arms but 6.15 times nominal.

### Multi-seed/reference/cross-dataset

These studies are well manifested and highly reproducible. They test deterministic behavior, yield stability, parser compatibility, ranking similarity, and accepted-set agreement. They do not test truth, accuracy, calibration, or biological generalization. Reuse of the available datasets during development makes them confirmation on known data, not independent external validation.

### Canonical baseline provenance

The repaired canonical PXD032157 baseline is 106,795 PSMs and 35,866 peptides at strict q<0.01, seed 1, all 65 files, fixed default settings. Commit `1348b0f` records the regression gate, and the current release binary exactly matches the frozen binary. The historical 107,046 PSM / 37,469 peptide counts did not become the corrected acceptance target; documentation explicitly demotes them.

There is no evidence that seed 1, a subset of files, a weaker threshold, fewer folds, or relaxed iterations was selected to recover yield. However, PXD032157 was repeatedly used to develop and tune class weights, folds, score scaling, and performance. Therefore these counts are a reproducible development baseline, not an unbiased accuracy or sensitivity estimate.

## 14. Performance/correctness audit

The inspected optimizations mostly preserve the stated algorithm:

- Active-set SVM evaluation computes the exact squared-hinge objective/Hessian for active rows; it is not a sampling approximation. Cholesky Newton plus line search is an algorithmic solver choice, so convergence/ranking equivalence to C++ is not guaranteed.
- Packed row layout copies values exactly. Cached initial scores are used only with the weights that produced them.
- q sorting constructs an in-bounds permutation before `get_unchecked`; exact numeric tie grouping occurs after the sort. The optimized q scan passes the independent oracle, but cannot cure upstream biased competition ties.
- The fixed `dot_22` reduction order is tested bit-exact against scalar order; SIMD AXPY is elementwise and does not reassociate the dot reduction.
- Fast parsing rejects nonfinite numbers; output-format fuzzing compares 100,000 finite values with the standard formatter.
- Debug and release builds produced byte-identical four-file fixture outputs. Current reruns exactly match all frozen null and multi-seed hashes.
- `tests/regression.sh` reproduces 117 PSM targets at q<0.01 and 43 peptide targets at q<0.05, and passes its named wall/RSS gates.

The reported speedup is valid only for its named host, dataset, command, and comparator. Exact fixture/hash equality and source inspection provide moderate evidence that the hot-path optimizations did not add hidden approximations. They do not establish exhaustive floating-point equivalence, equal convergence to official C++, or correctness of the statistical method itself. The performance/correctness-equivalence verdict is therefore **MODERATE EVIDENCE**.

## 15. Claim-by-claim README and validation-document audit

Classifications use the requested claim vocabulary, distinct from the final verdict vocabulary.

| Claim or wording | Classification | Basis / required correction |
|---|---|---|
| “Reported q scan uses TDC+, ties, opportunity ratio, monotonicization” | SUPPORTED WITH CAVEATS | Standalone finite-input oracle passes; only conditional on a valid winner list |
| “FDR controlled” / “valid q-values” without scope | INCORRECT | Exact-tie competition attack and anti-conservative entrapment curve |
| “Complete-null pathology repaired” | SUPPORTED WITH CAVEATS | 17/30 → 0/30 on predefined fixtures; underpowered and misses structured ties |
| “Nothing about heldout reaches the model” on canonical fixed-C path | SUPPORTED | Direct corruption and mutation tests |
| “Leakage-free” across supported program modes | INCORRECT | `--select-c` and `--ensemble` leak labels/selection information |
| `--auto-model` is nested | SUPPORTED WITH CAVEATS | Code structure is nested; less adversarial coverage than fixed-C |
| `--rt-features` remains outside isolation | INCORRECT | Current code refits fold-local; README is stale |
| “PEP follows the Käll identity” / implies published validation | MISLEADING | Direction of implication reversed; added prior breaks exact cumulative identity; official-like code is not calibration |
| “PEP cannot be zero” | SUPPORTED | Code/oracle and 98,680 competed targets show zero exact zeros |
| “PEP calibrated” / valid posterior probability | INCORRECT | Every populated entrapment bin optimistic; ECE/signed error 0.0185 |
| “Peptide level uses the same estimators” | SUPPORTED WITH CAVEATS | It recomputes them on best peptide representatives; inherits q/PEP defects |
| “Protein level uses the same estimators” | INCORRECT | Picked output copies peptide PEP; Bayesian uses a different unvalidated posterior model |
| Protein groups are “indistinguishable” | INCORRECT | Connected-component counterexample merges distinguishable A/B proteins |
| Picked protein method is more sensitive/better calibrated | INCORRECT | Yield-only test plus target-favoring tie attack; no truth/calibration study |
| “Expected picked-FDR benefit” from the synthetic fixture | MISLEADING | Test encodes `picked >= classic`; it cannot validate scientific benefit |
| “Percolator-compatible” | SUPPORTED WITH CAVEATS | Strong PSM/peptide agreement on four datasets; material MSFragger/yeast and method differences |
| “Equivalent” to C++ statistically or algorithmically | INSUFFICIENT EVIDENCE | High correlations/overlap are not equivalence; solver/scaling/PEP differ |
| “Generalizes” | INSUFFICIENT EVIDENCE | Four known datasets, no untouched external truth set, calibration fails |
| “Improved” relative to rejected build | SUPPORTED WITH CAVEATS | Null, zero-PEP, leakage, and entrapment metrics improve; method still fails validation |
| “More accurate” | INSUFFICIENT EVIDENCE | No independent truth-based accuracy comparison |
| “More sensitive” | INSUFFICIENT EVIDENCE | Identification yield is not sensitivity; FDR differs and protein evidence is invalid |
| Repaired identification counts | SUPPORTED | Frozen baseline and regression history reproduce; development baseline only |
| 23.4× performance statement | SUPPORTED WITH CAVEATS | Named benchmark provenance; not portable and not correctness evidence |
| Parser fails closed on malformed/nonfinite input | SUPPORTED WITH CAVEATS | Nonfinite/text failures reject; non-±1 numeric labels are accepted by sign |
| Repair report: q=.001 entrapment result is “conservative” | INCORRECT | 0.00615 > 0.001 |
| Repair report: residual shared C++ failure is not `percolator-rs` | MISLEADING | Excludes a large Rust-only divergence, not a shared algorithm/method defect |
| Completion audit: scientific repair is complete | INCORRECT | Structured competition ties, PEP calibration, leakage modes, and protein defects remain |

The older `SCIENTIFIC_VALIDATION.md` and `IMPLEMENTATION_AUDIT.md` are clearly marked as rejected-build history. Their failure statements are appropriate historical records, not current claims.

## 16. Unresolved questions

1. What proportion and label composition of exact/near-exact target-decoy score ties occurs across real supported PIN inputs, including quantized features or degenerate models? The constructed attack proves invalidity even if prevalence is low.
2. Which tie rule will be used—seeded fair coin per competition, a conservative decoy win, or tie removal—and how will it be proven order-invariant and propagated to protein picking?
3. Can a proper PEP model (for example, independently validated nonparametric logistic/qvality-style estimation) calibrate across held-out entrapment data without tuning on those same data?
4. Why does the matched entrapment design remain anti-conservative, and how much comes from TDC assumptions, foreign-fraction adjustment, spectrum/protein dependence, competition, and semi-supervised training?
5. What is the scientifically intended protein grouping definition, target/decoy pairing rule, and protein-level PEP model? Current code implements none of these claims correctly enough to validate.
6. Can dataset-specific entrapment/calibration be measured for MSFragger, Sage, Tide, and yeast rather than only pooled?
7. Can an untouched external collection with predeclared thresholds, seeds, success criteria, and uncertainty estimates be reserved for generalization?
8. Can stored fold-local C++ intermediate scores/weights isolate the cause of MSFragger and yeast exclusive calls?
9. Should the public PIN contract require labels exactly `+1/-1`, and should all public statistical APIs reject nonfinite scores?
10. Can `--select-c` be removed/deprecated in favor of nested `--auto-model`, and can ensemble feature construction be made outer-fold-local?

Proposed fixes are deliberately not implemented here. The minimum scientific repair order is: unbiased/conservative competition ties; strict input/API invariants; replacement and external calibration of PEP; correct protein grouping/pairing/posteriors; elimination of nonnested/global label-derived modes; then a fresh, frozen validation on untouched data.

## 17. Final verdict table

The verdict is for the current implementation and its important documented modes, not merely the best-behaved helper function.

| Dimension | Verdict | Concrete evidence |
|---|---|---|
| IMPLEMENTATION CORRECTNESS | **FAILED VALIDATION** | 200-pair row-order competition attack; invalid protein grouping/ties/PEP; leaking modes |
| CROSS-VALIDATION ISOLATION | **FAILED VALIDATION** | Canonical fixed-C is strongly isolated, but `--select-c` changes all 200 held-out scores after label-only corruption and ensemble constructs a global label-keyed feature |
| Q-VALUE IMPLEMENTATION | **FAILED VALIDATION** | Standalone q arithmetic passes 12 hand cases, but reported pipeline q-values can be forced from zero target winners to 200 discoveries at q<.01 solely by tied-row order; internal NaN also perturbs finite-prefix q |
| FDR CALIBRATION | **FAILED VALIDATION** | Five-seed entrapment mean adjusted FDP 0.01816 at nominal 0.01 and above nominal at all six thresholds; exact-tie attack is arbitrarily anti-conservative |
| PEP IMPLEMENTATION | **FAILED VALIDATION** | Unjustified derivative+prior construction, cumulative identity broken, large PEP mutations pass all tests, protein PEP is copied peptide PEP |
| PEP CALIBRATION | **FAILED VALIDATION** | Every nonempty repaired bin optimistic; weighted absolute and signed error both +0.018514 |
| MULTI-SEED REPRODUCIBILITY | **STRONG EVIDENCE** | 20 Rust and 20 C++ reruns; four output hashes per run exactly reproduce; all predefined seeds retained |
| REFERENCE AGREEMENT | **STRONG EVIDENCE** | Mean q=.01 Jaccard 0.920–0.996 and high score/q/PEP rank correlations over four matched datasets, explicitly agreement rather than correctness |
| CROSS-DATASET GENERALIZATION | **WEAK EVIDENCE** | Four known datasets show compatibility, but MSFragger/yeast disagreement is material, no untouched truth set exists, and pooled entrapment fails calibration |
| VALIDATION-SUITE QUALITY | **WEAK EVIDENCE** | Detects q/CV mutations, but permits two major PEP mutations and explicitly requires invalid input-order competition ties; protein tests encode yield/implementation behavior |
| PERFORMANCE/CORRECTNESS EQUIVALENCE | **MODERATE EVIDENCE** | Exact debug/release fixture and frozen hashes plus source inspection show no obvious hidden approximation; solver/convergence and exhaustive floating-point equivalence remain unproven |

### Claims I would accept, qualify, or reject in a paper

**Accept:** a deterministic Rust implementation exists; its standalone finite-input q scan implements the stated TDC+ arithmetic conditional on valid competition; canonical fixed-C outer-fold model fitting is isolated; frozen results are exactly reproducible; and it has strong empirical PSM/peptide agreement with official C++ on the four tested datasets plus a reproducible named speed benchmark.

**Require qualification:** the repair substantially improves the rejected build; complete-null results are encouraging but low-powered and incomplete; canonical identification counts are development baselines, not sensitivity; “Percolator-compatible” means high empirical agreement, not equivalence; and all CV language must explicitly exclude `--select-c` and current ensemble construction.

**Reject:** valid/calibrated q-values or FDR for arbitrary supported inputs; calibrated PEPs; unqualified leakage-free behavior; valid protein grouping, q-values, FDR, or PEPs; claims of improved accuracy/sensitivity; broad biological/search-engine generalization; and statistical, algorithmic, or probability equivalence to C++.
