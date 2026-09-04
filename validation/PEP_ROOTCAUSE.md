# Root cause of optimistic PEPs in signal-present entrapment

Investigation date: 2026-08-30  
Audited commit: `e8d83d1c76e4cf651fdfcf22d98b0b499c35943a`  
Audited binary SHA-256: `be9bf670bfd69df4dc3ba3b8be6c4c164acaf56a91f4a2819d115f49604b2c45`

No production statistical methodology was changed. No PEP was tuned, shifted,
clamped, rescaled, or fitted to entrapment outcomes. The investigation added only
validation scripts and reports.

Machine-readable evidence is in `PEP_ROOTCAUSE_FINAL_RESULTS.json`. Every populated
calibration bin, pooled and separately by seed and dataset, is in
`PEP_ROOTCAUSE_CALIBRATION_TABLES.md`.

## Executive conclusion

**Primary classification: TARGET/DECOY NON-EXCHANGEABILITY — HIGH confidence.**

The reported PEP is an estimate of the local number of false targets based entirely
on the local number of decoys. In the confident tail of this search, reversed decoys
do not represent the false biological targets. Consequently the estimator has too
little null mass to assign and reports PEPs that are too small.

This is quantitatively sufficient, not merely correlated. In PEP bin `b`, let

- `N_b` be the number of target PSMs;
- `M_b = N_b * mean(PEP_b)` be the predicted false-target mass;
- `D_b` be the usable decoy count;
- `E_T,b` and `E_D,b` be pure-entrapment target and decoy counts.

Using the local entrapment fraction among decoys, `f_b = E_D,b / D_b`, the
entrapment-implied calibration ratio factorizes exactly as

```text
observed_b / predicted_b
  = [E_T,b / (f_b N_b)] / [M_b / N_b]
  = (E_T,b / E_D,b) * (D_b / M_b).
```

The second factor is approximately one because PEP mass is derived from decoy mass.
Across the nine populated Rust bins, log(`observed/predicted`) and
log(`entrapment target/decoy`) correlate at **0.9981**; their median residual is only
**5.1% on a fold scale**. The observed target/decoy imbalance therefore explains the
direction and almost all of the magnitude bin by bin.

Several secondary findings qualify that conclusion:

1. **enzN/enzC are causal contributors, not the sole cause.** Removing them reduces
   the q<0.01 internal-null ratio from 1.96 to 1.49 and pooled adjusted calibration
   error from +0.01862 to +0.00971, but does not restore exchangeability.
2. **The narrow semi-tryptic/reversal-terminus hypothesis is falsified.** A matched
   fully-tryptic re-search makes enzN/enzC essentially constant but does not improve
   either PEP calibration or the internal null; both are slightly worse in the
   confident region.
3. **Biologically structured near-homologs are a major residual channel.** At
   PEP<0.01, 16/45 distinct entrapment-target peptides have a native same-length
   substring within two I/L-aware substitutions, versus 3/27 entrapment decoys.
   Those 16 peptides account for 317/414 PSM rows and nine are actin/tubulin-like.
4. **Finite local samples matter.** Ideal-assumption simulations show modest
   finite-sample optimism for very small lists and sparse tails. This inflates
   uncertainty and can add bias, but it does not create the real 3–13-fold local
   deficit of decoys.
5. **The pooled +0.01862 is not a clean effect size.** About 87% comes from applying
   the entrapment `1/f` extrapolation to the `[0.5,1]` junk bin. Without extrapolation,
   the pooled known-false lower bound is conservative. Below PEP 0.05, however, the
   directly observed entrapment-target fraction alone exceeds predicted PEP, so the
   confident-tail failure is not an adjustment artifact.

There is **no demonstrated Rust defect affecting target PEPs**. A separate unstable
tie-order defect was found in the PEP values written for decoy rows only; it cannot
cause target calibration failure.

## 1. Exact production path and implemented mathematics

### 1.1 Scores entering the estimator

The path was traced from code rather than comments:

1. `percolator::run` creates deterministic, spectrum-grouped three-fold splits.
2. Every row is scored out of fold by a model trained without that fold.
3. Each fold's held-out margins are transformed by the increasing affine map
   `(margin - training_decoy_mean) / training_decoy_sd`. This aligns fold scales
   but cannot change within-fold order or orientation; larger remains better.
4. `main.rs` performs spectrum-level target/decoy competition on the out-of-fold
   scores using `(source, ScanNr, ExpMass)` as the spectrum key. Exact score ties
   are resolved by a seeded fair draw.
5. q-values and PEPs are recomputed on the actual reported winner list. Only then
   are rows split into target and decoy TSV files.

There is no hidden score negation, score-to-p-value conversion, logarithm, density
fit, or interpolation before the Rust target PEP estimator.

### 1.2 Ordering and ties

Let the reported winner scores be `s_i` and labels be `l_i in {+1,-1}`. Non-finite
scores are rejected. Rows are sorted by descending `f64::total_cmp`; exact numeric
equality defines a score group, so `+0.0` and `-0.0` share a group. Every target in
one exact-score group receives the same raw increment.

For a concatenated 1:1 TDC search, production uses

```text
pi0 = 1
p = P(null target wins) = 0.5
lambda = p / (1-p) = 1
initial decoys = 1
```

The initial decoy is the finite-sample `+1` safeguard. `pi0` is not estimated.

### 1.3 q-values

For score group `g`, after counting the whole tie group, let `T_g` and `D_g` be
cumulative targets and observed decoys; in the following expression the safeguard
is explicit:

```text
raw_FDP_g = min(1, pi0 * lambda * (D_g + 1) / max(T_g, 1)).
q_g       = min_{h >= g} raw_FDP_h.
```

The reverse cumulative minimum is shared by the entire tie group.

### 1.4 Target PEPs

Index only score groups that contain at least one target. Let `n_g` be the number
of targets in group `g`, and let `F_{g-1}` be the false mass already assigned at the
previous target-containing group, starting at zero. Rust computes

```text
F_g = min(
          pi0 * lambda * (D_g + 1),
          T_g,
          F_{g-1} + n_g
      )

y_g = max((F_g - F_{g-1}) / n_g, 0).
```

`y_g` is repeated once for every target in the group. A decoy-only score group
raises `D` but assigns no mass; that mass reaches the next target-containing group.
The third bound prevents any group increment from exceeding one per target.

Unit-weight, non-decreasing PAVA is then applied to the target sequence in
best-to-worst order:

```text
(pep_1, ..., pep_T) = argmin_{z_1 <= ... <= z_T} sum_i (z_i - y_i)^2.
```

Equivalently, the fitted PEPs are the slopes of the greatest convex minorant of
the cumulative raw-increment curve. PAVA preserves total mass over its complete
input, but it may redistribute mass across intermediate score thresholds.

Finally target PEPs are clamped to `[1e-12, 1]`. In finite production data the
`+1` safeguard and leading PAVA pool, not the numeric floor, determine the minimum.
Target values are not interpolated as a function of score. The TSV writer rounds
to six decimals.

### 1.5 Decoy PEPs and the one implementation defect found

After target PEPs are fixed, a row-wise pass gives a decoy the latest target PEP
encountered above it; leading decoys receive the first target PEP. These values make
the decoy output visually comparable but carry no target-error interpretation.

Because the score sort is unstable within exact ties, a decoy tied to a target can
receive either the prior or tied target's PEP depending on intra-tie order. A six-row
reproducer moves one decoy PEP from 0.25 to 1.0 while every target PEP and q-value is
unchanged. The smallest justified repair would fill decoys once per tie group. It was
not implemented here, and it is irrelevant to target calibration.

## 2. Comparison with the intended local-FDR methodology

The primary method defines local FDR/PEP as `pi0 f0(s) / f(s)` and estimates the
density ratio using decoys. The original method uses binned non-parametric logistic
regression; QVALITY pools target and null scores into up to 500 equal-sized,
tie-respecting bins, uses bin medians and decoy fractions, fits a penalized natural
cubic spline on the logit scale by IRLS, selects roughness by cross-validation, and
evaluates the fitted ratio at scores. Generic QVALITY estimates `pi0` by a Storey
bootstrap; its `-Y/--tdc-input` option fixes `pi0=1` for concatenated TDC. See the
[primary local-PEP paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC2732210/), the
[QVALITY paper](https://noble.gs.washington.edu/papers/kall2009qvality.pdf), and the
[q-value/PEP relationship paper](https://noble.gs.washington.edu/papers/kall2008posterior.html).

Current C++ Percolator 3.09 must be distinguished from standalone QVALITY. Its
default is now a score-aware monotone I-spline applied to increments of `q_i*i`,
with a `0.5/n` pseudocount; `--pava-pep` selects PAVA and `--irls-pep` selects the
historical QVALITY spline. This was confirmed from release-3.09 source and agrees
with the [3.09 release notes](https://github.com/percolator/percolator/releases).

| Stage | Rust | Original QVALITY / current C++ | Classification |
|---|---|---|---|
| Score orientation | larger is better; increasing fold standardization | orientation supplied by caller | exact in meaning |
| TDC input | competed target/decoy winners | `-Y` TDC supports the same design | exact in design |
| `pi0` | fixed 1 | QVALITY generic: Storey bootstrap; `-Y`: 1; C++ TDC: 1 | exact for TDC; intentional difference from generic mode |
| Null opportunity | explicit `lambda=p/(1-p)` | C++ TDC has same factor; historical QVALITY assumes supplied null scale | mathematically equivalent at `p=.5` |
| Local evidence | increments of cumulative decoy false mass | QVALITY: binned decoy fraction; C++ default: increments of `q_i*i` | intentional methodological difference |
| Smoother | unit-weight PAVA / GCM | QVALITY: penalized logistic spline; C++ default: monotone I-spline | intentional difference |
| Score interpolation | none for targets | both spline methods evaluate a score-aware curve | intentional difference |
| Monotonicity | PEP non-decreasing as score worsens | monotone fit / final monotone processing | equivalent constraint, different loss |
| Finite-sample boundary | TDC `+1`, leading pool, `1e-12` floor | QVALITY spline extrapolation; C++ default `0.5/n` pseudocount | intentional difference |
| Exact ties | one Rust count/increment per exact score | QVALITY bins preserve ties; C++ score fit shares score coordinate | compatible, not identical |
| Decoy output | previous target value; tie-order defect | interpolation | suspicious Rust difference, decoy-only |

Therefore Rust is **not an exact implementation of historical QVALITY**, and its
function names should not be read that way. The differences are explicit estimator
choices, not evidence of an implementation error. The empirical comparisons in
section 11 test whether any of them explains the optimism.

## 3. Synthetic calibration under valid assumptions

### 3.1 Code-to-mathematics verification

`pep_rootcause_probe.rs` links production `src/stats.rs` and compares target PEPs
with an independently written greatest-convex-minorant oracle. It passes:

- complete null, perfect separation, moderate overlap;
- exact all-score ties and repeated score blocks;
- samples of size 1 through 8;
- separate severe-imbalance fixtures with 1% targets and with 1% decoys;
- `+0`, `-0`, `MIN_POSITIVE`, and extreme finite scores;
- null-win probabilities 0.1, 0.25, 0.5, 0.75, and 0.9;
- hand-computed false-mass examples; and
- bounds, monotonicity, and mass conservation on 5,000 rows.

A separate NumPy implementation matched production target PEPs exactly on 60
randomized fixtures. All 109 Rust unit tests and all integration tests pass.

### 3.2 Analytic generative tests

Each synthetic spectrum is null with probability `pi0`. Under the null, target and
decoy scores are iid from `f0`; under signal, target is drawn from `f1` and decoy from
`f0`; the winner is retained. For a target winner at score `s`, competition cancels
from numerator and denominator, giving the analytic truth

```text
PEP_true(s) = pi0 f0(s) / [pi0 f0(s) + (1-pi0) f1(s)].
```

For quantized scores the corresponding normal-bin masses replace densities and
ties are decided by a fair coin.

| assumption-holding case | target rows | observed−predicted |
|---|---:|---:|
| complete null, n=200k ×4 | 399,255 | +0.00296 |
| all 200k scores exactly tied, complete null | 100,000 | 0.00000 |
| strong separation, n=200k ×4 | 599,164 | −0.00019 |
| moderate overlap | 568,460 | −0.05517 |
| heavy overlap | 455,759 | −0.31706 |
| target-heavy class imbalance | 748,049 | −0.01549 |
| near-balanced output | 438,910 | −0.00198 |
| repeated scores on a 0.5 grid | 567,200 | −0.05757 |
| dense ties on a 0.1 grid | 568,109 | −0.05651 |
| small n=50, 1,000 replicates | 36,964 | +0.04531 |
| sparse-tail n=2,000, 200 replicates | 299,468 | +0.01132 |
| large n=1,000,000 | 742,162 | −0.00974 |
| dense extreme tail, n=1,000,000 | 749,164 | −0.00128 |

The estimator is usually conservative because decoys that beat true targets add
null mass without a corresponding false target. It is not universally conservative:
very small lists and sparse tails show finite-sample downward bias. In the sparse-tail
case the first populated `[.001,.005)` bin is nevertheless calibrated
(predicted 0.001419, observed 0.001386); the +0.011 global error is mostly in the
poor-score bulk. Dense data recover the analytic local posterior down to roughly
`1e-5`.

The broader 72-condition grid has median signed error −0.01447. Eleven of 72 cases
are positive, mostly smaller samples. In 40 repeated n=200k runs, mean error is
−0.01072 (sd 0.00195) and none is positive.

### 3.3 Violating one assumption at a time

| scenario | weighted signed error | observed shape |
|---|---:|---|
| exchangeable control | −0.0108 | conservative |
| clustered PSMs | −0.0097 to −0.0156 | conservative |
| quantized ties | −0.0103 to −0.0114 | conservative |
| semi-supervised training with exchangeable features | −0.0103 to −0.0118 | no iteration trend |
| narrower decoy tail | −0.0034 to +0.0079 | small effect |
| global decoy shift | +0.0094 to +0.1335 | wrong shape: strongest in bulk |
| unmatched false-target/homolog component | +0.0126 to +0.1167 | large low-tail ratio decaying toward one |

Only an unmatched, high-scoring false-target component reproduces the real shape.
A synthetic partial-match homolog channel also reproduces the observed training
dose-response, whereas a uniform target/decoy provenance feature is downweighted by
training. These simulations establish possibility and mechanism; they were not fit
to the benchmark and are not proof about the biological source.

## 4. Canonical real-data calibration

The canonical result is 5 predefined seeds × 6 semi-tryptic Comet PINs, separately
trained per file, fixed SVM class weights 1:4, `maxiter=10`, and production
spectrum competition. The rerun has 838,135 reported winner rows and 493,760 target
rows. The global pure-entrapment share among non-mixed decoys is `f=0.783010`.

| PEP interval | targets | mean PEP | ent T | ent D | observed f=1 | observed adjusted | predicted−observed adjusted | observed/predicted | entT/entD |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| [1e-4,1e-3) | 54,378 | .000386 | 217 | 17 | .003991 | .005096 | −.004710 | 13.20 | 12.76 |
| [1e-3,5e-3) | 15,887 | .002707 | 113 | 32 | .007113 | .009084 | −.006377 | 3.36 | 3.53 |
| [5e-3,.01) | 5,750 | .006783 | 84 | 28 | .014609 | .018657 | −.011874 | 2.75 | 3.00 |
| [.01,.02) | 7,696 | .015852 | 193 | 85 | .025078 | .032028 | −.016175 | 2.02 | 2.27 |
| [.02,.05) | 9,091 | .033330 | 370 | 215 | .040700 | .051978 | −.018649 | 1.56 | 1.72 |
| [.05,.10) | 5,374 | .073874 | 337 | 276 | .062709 | .080088 | −.006213 | 1.08 | 1.22 |
| [.10,.20) | 9,231 | .145705 | 1,169 | 1,032 | .126639 | .161733 | −.016028 | 1.11 | 1.13 |
| [.20,.50) | 19,260 | .345587 | 5,418 | 5,030 | .281308 | .359265 | −.013679 | 1.04 | 1.08 |
| [.50,1] | 367,093 | .913066 | 268,738 | 262,592 | .732071 | .934944 | −.021878 | 1.02 | 1.02 |

Observed minus predicted, pooled and adjusted, is **+0.018624**, reproducing the
prior +0.018685 within its previously documented choice of entrapment fraction.
Every adjusted bin is optimistic. The first five bins are also optimistic under
the adjustment-free known-false lower bound; from PEP 0.05 upward they are not.

The pooled number is dominated by validation design. The `[.5,1]` bin contains
74.4% of targets and contributes +0.01627, or 87.3% of the pooled +0.01862. With
`f=1`, pooled observed minus predicted is −0.13664. Thus the defensible result is:

- **confident PEPs below 0.05 are demonstrably optimistic even without `1/f`;**
- the claim that the entire PEP range is optimistic depends on extrapolating the
  entrapment fraction and should not be interpreted as a direct ground-truth fact.

Complete per-seed and per-dataset bin tables, including Wilson intervals and all
unfavorable bins, are in `PEP_ROOTCAUSE_CALIBRATION_TABLES.md`. Seed pooled errors
range +0.01657 to +0.02018. Dataset pooled errors range +0.00979 to +0.03918; three
of six datasets do not populate the lowest PEP bin.

## 5. Extreme-tail and feature analysis

Score percentiles are computed within each seed × dataset, because fold-standardized
scores are not necessarily comparable across trained files.

| score region | ent T | ent D | ratio |
|---|---:|---:|---:|
| top 0.1% | 0 | 0 | undefined |
| top 0.5% | 2 | 0 | undefined |
| top 1% | 11 | 0 | undefined |
| top 2% | 55 | 3 | 18.33 |
| top 5% | 972 | 796 | 1.22 |

The first three regions are too sparse to estimate a ratio. The imbalance appears
abruptly where entrapment nulls first populate the high-score tail, then decays.
Bulk score distributions nearly coincide: mean score is 0.9164 for entrapment
targets and 0.9051 for entrapment decoys. At the 99th percentile they are 3.2734
and 3.1570. A global AUC near 0.5 therefore does not contradict a severe tail ratio.

Conditional feature distributions show the enzymatic channel but also its limits:

| region | enzN ent T / ent D | enzC ent T / ent D | strongest standardized difference |
|---|---:|---:|---|
| PEP<.001 | .995 / 1.000 | 1.000 / 1.000 | lnrSp, but only 17 decoys |
| PEP<.01 | .969 / .948 | 1.000 / .922 | enzC, SMD +0.733 |
| PEP<.05 | .967 / .899 | .989 / .950 | Xcorr, SMD +0.491 |
| q<.01 | .945 / .885 | .988 / .953 | Xcorr, SMD +0.482 |

Thus enzN/enzC distinguish false targets and decoys in a relevant region, but the
most extreme PEP block is dominated by repeated homolog-like PSMs and too few
decoys for stable feature attribution. Bin-wise `entT/entD` in section 4 is the
empirical local null-density ratio.

## 6. Explicit shared-cause test

The factorization in the executive conclusion was evaluated in every bin. The
ratio `(observed/predicted)/(entT/entD)` is

```text
1.034, 0.950, 0.917, 0.890, 0.906, 0.888, 0.980, 0.965, 1.001.
```

For example, in `[1e-4,1e-3)`, local usable-decoy count divided by PEP false mass
is 1.00034. Therefore an entrapment target/decoy ratio of 12.76 predicts a
calibration ratio of 12.77; observed adjusted/predicted is 13.20. In `[.5,1]`, the
same calculation predicts 1.023 and observes 1.024.

This is the quantitative bridge between q-value and PEP investigations. The PEP
estimator assumes that decoys approximate the local false-target distribution.
When false biological targets have a heavier high-score tail than reversed decoys,
the local decoy density is too low, `F_g` grows too slowly, and PAVA can only smooth
the deficient mass—it cannot invent the missing false targets.

## 7. enzN/enzC causal ablation

Only the two PIN feature columns were removed; candidate rows, labels, spectra,
seeds, training settings, competition, and PEP estimator were held fixed. This is
an investigative intervention, not a proposed production change.

| metric | canonical | drop enzN/enzC | result |
|---|---:|---:|---:|
| all reported target rows | 493,760 | 485,116 | −1.8% |
| pooled adjusted calibration error | +.018624 | +.009705 | −47.9% |
| targets at PEP<.01 | 76,015 | 56,827 | −25.2% |
| mean PEP at PEP<.01 | .001355 | .001865 | higher |
| adjusted false fraction at PEP<.01 | .006956 | .006076 | lower |
| calibration ratio at PEP<.01 | 5.13 | 3.26 | improved |
| entT/entD at PEP<.01 | 414/77 = 5.38 | 270/80 = 3.38 | improved |
| q<.01 targets | 97,682 | 73,307 | −24.9% |
| q<.01 normal targets | 93,800 | 70,980 | −24.3% |
| q<.01 entrapment targets | 1,293 | 760 | −41.2% |
| entT/entD at q<.01 | 1.956 | 1.493 | improved |
| adjusted FDP at q<.01 | .01691 | .01326 | improved |

The very lowest PEP region remains sparse and unfavorable: PEP<.001 has 109/6
entrapment targets/decoys after ablation versus 217/17 before. Its adjusted gap
improves from +.00471 to +.00326, but the internal ratio point estimate worsens
because only six decoys remain. This bin is not hidden and cannot support a precise
ablation effect by itself.

The joint improvement in cumulative FDP and local calibration supports a shared
causal mechanism: training used termini information that was associated with
sequence provenance. The incomplete improvement proves that it is only one channel.
Identification loss is reported to characterize the intervention, not optimized.

## 8. Semi-tryptic versus fully-tryptic search

All six mzML files were re-searched against the same combined FASTA with the same
Comet parameter files except `num_enzyme_termini=2`. The semi- and fully-tryptic
PINs contain 837,367 and 829,476 candidates, respectively. Both were rescored with
the same audited binary and five seeds.

| metric | semi-tryptic | fully-tryptic |
|---|---:|---:|
| pooled adjusted calibration error | +.01862 | +.04869 |
| PEP<.01 targets | 76,015 | 68,160 |
| PEP<.01 mean predicted | .001355 | .001203 |
| PEP<.01 adjusted false fraction | .006956 | .007643 |
| PEP<.01 calibration ratio | 5.13 | 6.35 |
| PEP<.01 entT/entD | 5.38 | 6.17 |
| q<.01 targets | 97,682 | 90,316 |
| q<.01 normal targets | 93,800 | 87,003 |
| q<.01 entrapment targets | 1,293 | 1,325 |
| q<.01 entT/entD | 1.956 | 2.127 |
| q<.01 adjusted FDP | .01691 | .01906 |

In the fully-tryptic confident sets, enzN and enzC are essentially one for both
classes. Calibration nevertheless worsens. This rejects the prediction that merely
removing the semi-tryptic/reversal enzymatic asymmetry restores the null model.
Other features (`lnExpect`, `deltCn`, internal-cleavage count) and biological
near-homology still distinguish the selected populations. The fully-tryptic arm is
a new search, not a row-wise ablation, so its much larger pooled error also reflects
changed target composition and entrapment extrapolation; the low-PEP internal ratio
is the cleaner comparison.

## 9. Training-iteration dose response

Seeds 1–3 were rerun at predefined `maxiter` 0, 1, 2, 3, and 10. `maxiter=0` is
Percolator's selected initial direction, not an untouched raw score.

| design | iter | pooled error | PEP<.01 observed / predicted | PEP<.01 entT/entD | q<.01 entT/entD | q<.01 adjusted FDP |
|---|---:|---:|---:|---:|---:|---:|
| semi | 0 | +.00623 | .004965 / .001326 | 6.00 | 1.290 | .01157 |
| semi | 1 | +.01231 | .006468 / .001439 | 4.60 | 1.514 | .01418 |
| semi | 2 | +.01497 | .006903 / .001278 | 5.43 | 1.826 | .01633 |
| semi | 3 | +.01677 | .006898 / .001246 | 5.67 | 1.915 | .01677 |
| semi | 10 | +.01842 | .006977 / .001398 | 5.10 | 1.924 | .01666 |
| fully | 0 | +.02336 | .008082 / .001546 | 6.15 | 1.841 | .01615 |
| fully | 1 | +.04802 | .008908 / .001496 | 6.04 | 1.992 | .01842 |
| fully | 2 | +.04980 | .008796 / .001435 | 5.98 | 2.187 | .01945 |
| fully | 3 | +.05032 | .008307 / .001262 | 6.23 | 2.159 | .01937 |
| fully | 10 | +.04946 | .008125 / .001251 | 6.50 | 2.165 | .01920 |

In both designs, the cumulative internal-null ratio and pooled PEP error rise with
training and saturate around iteration 2–3. The very local PEP ratio is noisier
because it has tens of decoys, but remains several-fold throughout. Exchangeable-
feature synthetic training has no such trend. This supports amplification of a
pre-existing structured false-target channel rather than an error in later PEP code.
All iteration-specific bin tables are in the calibration appendix.

## 10. Raw Comet control

The same spectrum key and competition were applied to raw Comet XCorr and
`-lnExpect`, followed by the same analysis estimator. This does not turn a search
score into a PEP; it tests the ranking before semi-supervised rescoring.

At the matched depth of 133 pure-entrapment decoys:

| ranking | depth | ent T | ent D | ratio |
|---|---:|---:|---:|---:|
| raw XCorr | 6,117 | 134 | 133 | **1.008** |
| raw `-lnExpect` | 12,557 | 147 | 133 | **1.105** |
| Percolator maxiter 10, canonical prior audit | — | 258 | 133 | **1.940** |

Thus the raw XCorr ranking is approximately exchangeable at the depth where the
rescored ranking is strongly imbalanced. Raw `-lnExpect` does expose a much sparser
PEP<.001 block (22/2), showing that isolated raw-score extremes can already contain
structured matches. The matched-depth result is the stable control: most of the
broad extreme-tail amplification appears after Percolator training.

## 11. Rust, standalone QVALITY, and C++ Percolator

These are comparison controls, not correctness oracles.

### 11.1 Identical Rust scores and labels: estimator only

Seed-1 production winners were passed unchanged to standalone QVALITY 3.09. Both
generic pi0 estimation and TDC `-Y` were run.

| metric | Rust | QVALITY pi0 estimated | QVALITY `-Y`, pi0=1 |
|---|---:|---:|---:|
| target rows | 98,692 | 98,692 | 98,692 |
| target PEP sum | 68,867.996 | 68,547.218 | 68,530.748 |
| pooled adjusted error | +.018505 | +.021756 | +.021923 |
| target PEP<1e-4 | 0 | 6,739 | 6,748 |
| rank correlation with Rust | 1 | .8973 | .8970 |
| median absolute PEP difference | 0 | .02357 | .02348 |

The near-identity of QVALITY's two pi0 arms excludes pi0 as the explanation here.
Changing the local smoother does not fix calibration and its boundary extrapolation
creates thousands of much smaller PEPs.

### 11.2 Identical PINs, separate rescorers

C++ Percolator 3.09 was run with `--post-processing-tdc`; its default PEP is the
current I-spline, not QVALITY IRLS. Across five seeds:

| metric | Rust | C++ 3.09 |
|---|---:|---:|
| targets | 493,760 | 493,801 |
| pooled adjusted error | +.018624 | +.018258 |
| PEP<.01 entT/entD | 5.38 | 6.97 |
| q<.01 entT/entD | 1.956 | 1.868 |
| targets below PEP 1e-4 | 0 | 4,647 |

On the 88,825 seed-1 target PSMs selected by both programs, PEP rank correlation is
0.9418; median, 90th, and 99th percentile absolute differences are 0.01265, 0.08051,
and 0.22403. Of 47,778 common entrapment targets, Rust and C++ assign PEP<.01 to 81
and 77, with 64 in common. Differences are expected because the training procedures
and winner sets differ. The key result is that two different rescoring and local-fit
implementations have nearly identical pooled optimism and the same extreme-tail
internal-null failure. That argues strongly against a Rust-specific target-PEP bug.

## 12. Local-FDR assumption audit

| Assumption | Status | Evidence |
|---|---|---|
| concatenated equal-size target/decoy design | SUPPORTED | Comet decoys generated for the same combined database |
| correct score orientation | SUPPORTED | code trace and monotone affine fold transform |
| one spectrum-level winner | SUPPORTED | production competition rerun on `(source,scan,mass)` |
| `pi0=1` valid for direct TDC | LIKELY SUPPORTED | conservative setting; QVALITY pi0 arms agree |
| incorrect target and decoy exchangeability | **VIOLATED** | local entT/entD 12.76, 5.38, 2.59 below PEP .001, .01, .05 |
| decoys represent structured false biological targets | **VIOLATED** | near-homolog enrichment and quantitative bin factorization |
| semi-supervised score does not exploit provenance | **VIOLATED** | iteration dose response and enz ablation |
| score monotonicity of true PEP | UNCERTAIN | biologically distinct null components can cross; PAVA enforces one curve |
| local density estimate has enough tail data | LIKELY VIOLATED IN EXTREME TAIL | 17 entrapment decoys in the first pooled bin, 2–4 per seed |
| PAVA implementation is correct | SUPPORTED | exact independent GCM agreement |
| PAVA is unbiased in finite samples | LIKELY VIOLATED FOR SMALL/SPARSE LISTS | synthetic +.045 at n=50 and +.011 at n=2,000 |
| exact ties are handled for target PEPs | SUPPORTED | group-level increments and oracle tests |
| decoy PEP tie interpolation is deterministic | **VIOLATED, DECOY OUTPUT ONLY** | minimal unstable-sort reproducer |
| PSM independence | VIOLATED | repeated peptides and seeds; affects uncertainty, not observed direction in controls |
| class balance needs empirical count rescaling | SUPPORTED AS NOT REQUIRED | TDC opportunity factor is design-based; synthetic imbalance tests conservative |
| entrapment `1/f` gives total FDP | **LIKELY VIOLATED IN BULK** | foreign near-homologs are not representative of all false targets |

## 13. q-values and PEPs are related but not interchangeable

A q-value estimates cumulative false-discovery proportion for an accepted prefix.
A PEP estimates a local posterior for one score location. Correct q-value arithmetic
does not prove calibrated PEPs: local smoothing can place correct total mass at the
wrong scores, and sparse differentiation is noisier than cumulative counting.

In Rust they share the same TDC counts but are not literally identical algorithms.
q-values are reverse minima of cumulative `(D+1)/T`; PEPs are PAVA slopes of a
cumulative false-mass curve. PAVA preserves total mass over the full target list,
not every prefix. Therefore cumulative and local behavior can differ in degree.

Here, however, the evidence supports alternative **A**: both q-values and PEPs are
affected by the same null-model failure. The q<.01 internal ratio rises from 1.29 to
1.92 with training, while local PEP bins show the corresponding stronger undiluted
ratios. The synthetic finite-sample PEP limitation is secondary; there is no extra
Rust implementation defect and changing to QVALITY does not remove the failure.

## 14. Near-homolog analysis

All distinct pure-entrapment peptides in confident regions were compared with every
native FASTA protein using an exhaustive same-length substring search at Hamming
distance <=2 after I/L canonicalization. Three disjoint exact anchors guarantee that
every match within two substitutions is found.

| class and region | rows | distinct peptides | native distance <=2 | distances 0/1/2 |
|---|---:|---:|---:|---:|
| entrapment target, PEP<.01 | 414 | 45 | 16 (35.6%) | 0 / 9 / 7 |
| entrapment decoy, PEP<.01 | 77 | 27 | 3 (11.1%) | 0 / 1 / 2 |
| entrapment target, PEP<.05 | 977 | 156 | 37 (23.7%) | 1 / 18 / 18 |
| entrapment decoy, PEP<.05 | 377 | 111 | 18 (16.2%) | 0 / 7 / 11 |

At PEP<.01, near-homolog entrapment-target rows have mean PEP 0.00170 and mean score
5.480, versus 0.00474 and 4.482 for other entrapment-target rows. Nine of the 16
near-homolog target peptides match actin/tubulin headers. Examples include one- or
two-residue variants of `VAPEEHPVLLTEAPLNPK`, actin's
`TTGLVLDSGDGVSHTVPLYEGYALPHALLR`, and conserved tubulin segments.

Rows are repeated over seeds and spectra, so 317/414 is not an independent-sample
proportion; the distinct-peptide comparison is the relevant unit. The result shows
why a reversed sequence is a poor null for this population: a foreign biological
peptide can preserve almost all fragment evidence and ordinary enzymatic structure,
whereas its reversed decoy cannot. The one I/L-equivalent exact native match at
PEP<.05 also illustrates an entrapment-validation ambiguity—some “known false”
sequence assignments can be spectrally indistinguishable from native peptides.

## 15. Calibration uncertainty

The pooled adjusted error has a dataset-cluster bootstrap 95% interval of
**[+0.01436, +0.02586]** over six LC-MS/MS runs. A seed-resampling interval is
**[+0.01754, +0.01957]**, but this measures algorithmic reproducibility only because
all seeds reuse the same spectra.

Dataset-cluster intervals for observed-minus-predicted by bin are:

| bin | point | 95% interval |
|---|---:|---:|
| [1e-4,1e-3) | +.00471 | [.00208, .00688] |
| [1e-3,5e-3) | +.00638 | [.00204, .00898] |
| [5e-3,.01) | +.01187 | [.00952, .01845] |
| [.01,.02) | +.01618 | [.01161, .02519] |
| [.02,.05) | +.01865 | [.00615, .02459] |
| [.05,.10) | +.00621 | [−.00905, .02845] |
| [.10,.20) | +.01603 | [.00883, .03866] |
| [.20,.50) | +.01368 | [.00602, .03325] |
| [.50,1] | +.02188 | [.01703, .02881] |

Wilson intervals in the calibration appendix describe only the raw entrapment-target
fraction under an iid PSM approximation. They do not account for peptide clustering,
repeated seeds, or uncertainty in `f`. With only one realization per dataset, no
fully defensible dataset-specific sampling interval exists; dataset points and their
seed reproducibility are reported without pretending otherwise.

At PEP<.01 the internal ratio is 414/77=5.38; prior run-cluster bootstrap gives
[3.44,7.85], and peptide resampling [2.06,11.54]. At PEP<.001 it is 217/17, so its
12.76 point estimate is intrinsically imprecise. Tail sparsity limits the numerical
magnitude, not the conclusion that the ratio is above one.

## 16. Pre-registered alternatives: evidence for and against

| Hypothesis | Evidence for | Evidence against | Assessment |
|---|---|---|---|
| A. Rust implementation error | decoy-row tie-fill defect exists | target PEP exact oracle; ideal large samples; C++/QVALITY same failure | rejected for target calibration |
| B. QVALITY implementation difference | Rust is intentionally PAVA/GCM, not historical QVALITY | identical-score QVALITY is slightly more optimistic | not causal |
| C. pi0 error | generic local FDR depends on pi0 | TDC pi0=1 is conservative; QVALITY `-Y` and fitted arms nearly identical | rejected |
| D. PAVA/isotonic error | local differentiation can be high variance | production equals independent GCM exactly | no implementation error; finite-sample limitation remains |
| E. interpolation/boundary error | boundary behavior differs across estimators | Rust has no target interpolation; `+1` is conservative; QVALITY tail is worse | rejected as primary |
| F. target/decoy non-exchangeability | bin ratios, dose response, raw control, quantitative identity | none that explains the data comparably | **primary, strongly supported** |
| G. semi-tryptic/reversal provenance mismatch | enz ablation improves both q and PEP | fully tryptic does not improve; uniform-shift simulation has wrong trend | contributor, narrow form rejected |
| H. extreme-tail sparsity | 17 decoys in first bin; sparse synthetic bias | mid-PEP bins fail; dense shared-cause ratio; consistent dose response | secondary contributor |
| I. near-homolog false targets | direct native matches, lower PEP/higher scores, actin/tubulin enrichment | not all false targets are near homologs | major biological channel |
| J. entrapment validation bias | `1/f` creates 87% of pooled error; one I/L-equivalent exact match | raw known-false fraction alone exceeds PEP below .05; internal ratio needs no `1/f` | distorts pooled magnitude, not tail diagnosis |
| K. multiple contributing causes | F + I + termini + sparsity + validation extrapolation all measurable | one central identity explains local direction | true as secondary structure; F is primary |

## 17. Final classification and production-code decision

**Primary cause: TARGET/DECOY NON-EXCHANGEABILITY. Confidence: HIGH.**

The operational estimator limitation is that any decoy-based local-FDR method can
only be calibrated if decoys represent false targets locally. That assumption is
violated by a selected mixture of near-homolog biological targets and sequence-
provenance features. Semi-supervised training amplifies the mixture. Extreme-tail
sparsity adds uncertainty and modest estimator bias. The entrapment adjustment makes
the pooled all-range number look more uniformly optimistic than direct evidence
supports.

**Production statistical methodology should not change on this evidence.** A method
change that forces this benchmark to calibrate would be benchmark-specific
recalibration, not a principled repair. Documentation should state the TDC
exchangeability requirement, warn about homolog-rich foreign search spaces and
semi-tryptic reversed decoys, report local internal-null ratios, and separate direct
known-false lower bounds from `1/f` extrapolation.

The decoy-only tie-fill defect deserves an independent, minimal reproducibility fix,
but it must not be presented as a PEP-calibration repair and was not changed here.

## 18. Remaining uncertainty and the single next experiment

Magnitude is based on one dataset family, one search engine, one entrapment design,
and six runs. Fully tryptic results show that enzymatic termini do not exhaust the
provenance problem, but they do not identify every residual feature. Near-homology
is strongly enriched yet not universal. Local tail ratios have small decoy
denominators. The entrapment `1/f` conversion is not validated under this structured
foreign-target population.

The single most informative next experiment is a **pre-registered homology-depleted
entrapment re-search**: before searching, remove every entrapment peptide with a
same-length native peptide within two I/L-aware substitutions, regenerate decoys
from that filtered database, and rerun the identical semi-tryptic pipeline and seeds.
Primary endpoints should be PEP-bin `entT/entD`, direct f=1 calibration below .05,
and the iteration dose response. If those collapse toward one/calibration while raw
native yield behavior remains otherwise comparable, the near-homolog channel is
causally isolated. If not, the residual sequence-provenance features become the
next target of investigation.

## Direct answer

`percolator-rs` PEPs are systematically optimistic in the confident entrapment
region because the score-local decoy count substantially underestimates the number
of high-scoring false biological targets. The strongest false targets are enriched
for native near-homologs and ordinary sequence/enzymatic structure that reversed
decoys do not reproduce; semi-supervised training exploits this difference. The
implemented Rust target-PEP arithmetic is correct, and historical QVALITY and
current C++ Percolator remain optimistic on the same evidence. Sparse-tail local
estimation contributes bias and uncertainty, and the entrapment `1/f` adjustment
dominates the pooled +0.01862, but neither explains the adjustment-free low-PEP
failure. The root cause is therefore a local null-model failure—primarily
target/decoy non-exchangeability—with secondary extreme-tail and validation-design
limitations, not a production target-PEP implementation bug.
