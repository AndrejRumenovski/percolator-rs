# Second repair: competition ties, PEP semantics, protein inference, and the leaking CV modes

Repair date: 2026-08-26
Responds to: [`INDEPENDENT_AUDIT.md`](INDEPENDENT_AUDIT.md), which rejected commit `b38c0db`
Build under repair: `b38c0db` (statistically identical to the frozen `1348b0f`)
Reference implementation: Percolator 3.09.0, `~/opt/percolator-root/usr/bin/percolator`
Host: 12-core Ryzen 5 5600G, Linux 7.0.0-30

Every number below was produced by rerunning a predeclared experiment or by an attack script
committed alongside it. Nothing was tuned against a validation dataset, no seed or dataset was
dropped, and no adversarial fixture was removed or weakened. Where the repair did not help, that is
said plainly.

---

## 1. What changed

| # | Correction | Commit |
|---|---|---|
| 1 | Order-invariant tie-break primitive (seeded fair coin over competition-unit identity) | `d0fecb3` |
| 2 | PEPs estimated from the competition scan instead of differencing monotonized q-values; estimators fail closed on non-finite scores | `cb4b11f` |
| 3 | Exact target/decoy PSM score ties drawn with a fair coin; `--null-target-win-prob` contract closed at both ends | `651960e` |
| 4 | Protein grouping by identical peptide evidence; fair picked-protein ties; picked protein PEP removed (`NA`) | `a2bed3e` |
| 5 | `--select-c` made nested; ensemble agreement features de-labelled | `8e49e28` |
| 6 | Tie-grouping requirement tests; twelve-mutation test of the suite; entrapment root-cause probe | `7db2ee3` |
| 7 | Protein gate no longer asserts picked > classic; requires the `NA` column | `fd5524d` |

The test suite goes from 82 to **116** tests.

---

## 2. Target/decoy tie handling

### 2.1 Root cause

`competition_winners` kept the earlier row whenever two candidates of one precursor scored exactly
the same (`score[previous] >= score[i]`). The surviving label was therefore decided by the input
file's row order, which carries no information about the measurement. A unit test asserted this
behaviour rather than catching it.

The defect is upstream of the q-value scan, which is why the scan's own arithmetic passes every hand
oracle while the reported result is invalid: target-decoy competition converts one observed decoy win
into one expected incorrect target, and that conversion needs an incorrect target to beat its decoy
with the declared probability `p`. A rule that hands every tie to whichever label the file happens to
list first sets that probability to 0 or 1.

The same defect existed at the protein level (`target >= decoy`, section 5).

### 2.2 Correction

Ties are drawn **uniformly among the tied candidates**, with the draw taken from a SplitMix64 hash of
the competition unit's own identity — file, scan, experimental mass — mixed with the run seed. This
is the fair coin the primary target-decoy literature prescribes (Granholm, Noble & Käll 2011). The
hash sees no labels and no row indices, so:

* permuting the input cannot move any winner;
* a `k`-way tie holding `t` targets is won by a target with probability `t/k`, which is 1/2 in
  expectation under a 1:1 database for any `k`;
* the same seed on the same content always gives the same result, and a different seed re-flips every
  coin, so tie sensitivity is measurable across seeds instead of hidden in a fixed ordering.

Order-sensitivity was audited at every point where score order decides an outcome: PSM competition
(fixed), the q-value scan (already grouped by exact score), peptide representative selection (ties are
between same-label rows, so the statistics cannot change; the representative is now chosen by content
anyway), protein picking (fixed), protein group member order and output row order (now content-keyed),
and the fold assignment (already keyed on sorted spectrum identity).

### 2.3 Result — `validation/adversarial_competition.py`

200 spectra, one target and one decoy candidate each, identical constant features. Nine permutations
of the same rows:

| arm | pre-repair targets/decoys | pre-repair q<0.01 | **repaired** targets/decoys | **repaired** q<0.01 |
|---|---:|---:|---:|---:|
| target row first | 200 / 0 | **200** | 100 / 100 | **0** |
| decoy row first | 0 / 200 | 0 | 100 / 100 | 0 |
| file reversed | 0 / 200 | 0 | 100 / 100 | 0 |
| targets grouped first | 200 / 0 | **200** | 100 / 100 | 0 |
| decoys grouped first | 0 / 200 | 0 | 100 / 100 | 0 |
| shuffles 1–4 | 99–104 / 96–101 | 0 | 100 / 100 | 0 |

**Distinct winner sets across the nine permutations: 6 before, 1 after.** Minimum reported target
q-value goes from 0.005 to 1.0 — the correct answer for a file that contains no information.

A four-way tie fixture (three tied targets, one tied decoy per spectrum, 1000 spectra) gave 1000
target winners and 1000 discoveries at q<0.01 before; it now gives 764/236, against an expected
750/250, and no discovery.

Reseeding to seed 2 changes the split to 114/86 and still reports nothing, which is what a fair coin
should do.

### 2.4 Prevalence on real data

The attack proves invalidity; it does not show that real data triggers it. It largely does not.
Across the four reference datasets and five seeds, **every PSM and peptide count at q<0.01 is
bit-identical before and after the tie repair** (section 8). On the six-dataset entrapment study, three
of five seeds reproduce exactly and two move by 13 and 4 accepted PSMs. Continuous SVM scores rarely
tie exactly.

That is the honest scope: this was a **validity** defect with a near-zero measured effect on the
datasets at hand, and an unbounded effect on inputs that produce ties.

---

## 3. q-value revalidation after the tie repair

### 3.1 Hand oracles

`validation/independent_stats_probe.rs` compares `src/stats.rs` against vectors enumerated from the
formula in the script itself. All twelve cases pass unchanged: empty, single target, single decoy,
target-then-decoy, decoy-then-target, four targets, all decoys, `T,T,D,T,D,T`, mixed exact ties, all
identical scores, extreme finite scores, and an imbalanced opportunity ratio at `p = 1/3`. The
q-value arithmetic was never the defect and is unchanged.

One behaviour did change: a non-finite score is now **refused** rather than admitted. Previously a
trailing NaN target changed the q-values of every finite row above it through the reverse cumulative
minimum.

### 3.2 Tie grouping is now a requirement, not a side effect

Three targets and one decoy at score 5, one decoy at score 1: the tie group is evaluated once, after
all four rows have been counted, giving `1 * (1 + 1) / 3 = 2/3` for all four. Ending the group after
each row instead scores the leading targets before their tied decoy has been counted, and the reverse
cumulative minimum hands them 1/3 — half the estimate, from identical data. Both the q-value and the
PEP form of this requirement are now tested, under every rotation of the tie group.

### 3.3 Complete null — `validation/run_null.py`

Ten exact-balance relabelings of each of three PXD032157 PINs, 30 runs, thresholds predeclared.

| implementation | runs | runs with any false discovery at q<0.001 … q<0.10 |
|---|---:|---:|
| pre-repair (`d83a7ba`) | 30 | 17/30 at every threshold |
| audited build (`b38c0db`) | 30 | 0/30 |
| **this repair** | 30 | **0/30 at every threshold** |
| C++ 3.09 | 29 (1 crash) | 0/29 |

**No regression.** Zero rejections in 30 replicates still bounds the complete-null rejection
probability above by about 0.12 at 95% confidence; that is consistent with control, not a
demonstration of it at q = 0.001.

### 3.4 Identification counts across thresholds

Recorded, not targeted. See sections 8 and 10.

---

## 4. Posterior error probabilities

### 4.1 Root cause — implementation error, then estimator limitation

The previous estimator computed `raw PEP_k = k*q_k - (k-1)*q_{k-1}` from the **monotonized** q-value
curve and then added `0.5/N` of prior mass to every target. Three things are wrong with that, in
increasing order of importance.

1. `q_k` is a reverse cumulative minimum, so `k*q_k` is not the estimated number of false discoveries
   among the top `k` and its differences are not local error probabilities.
2. The added constant broke the identity the code cited as its justification. For five targets above
   one decoy the q-value is 0.2 and every reported PEP was 0.3.
3. Neither point is why the values are optimistic on entrapment data. That is inherited (§4.4).

### 4.2 Correction

The same scan that produces `FDP` also produces the estimated number of incorrect targets at or above
a score,

```text
F(s) = min(T(s), pi0 * lambda * (D(s) + 1))
```

which is exactly `T(s) * FDP(s)`. `F` is a non-decreasing step function, and the PEP of a target is
the increment of `F` it is responsible for: each tie group's increment is shared equally among the
targets it contains, decoy-only groups carry their increment forward, and a group of `g` targets can
absorb at most `g` further false discoveries because a probability cannot exceed one. PAVA then
enforces monotonicity in score, redistributing mass only inside a pooled block.

Consequences, none of them tuned:

* **The reported PEPs sum exactly to the reported estimated false-discovery count.** Tested.
* **No finite input gives PEP = 0.** The finite-sample safeguard decoy supplies one estimated false
  discovery spread across the targets that outrank every decoy — for five targets above one decoy,
  0.2 each, which is also their q-value. No constant is added anywhere.
* **The declared opportunity ratio scales PEPs exactly as it scales q-values.** Tested at `p = 1/3`.

### 4.3 Frozen oracles

| case | expected | pre-repair | repaired |
|---|---|---|---|
| 5 targets, 1 decoy | every PEP 0.2, sum 1.0 | 0.3 each | **0.2 each** |
| 10 targets, 1 decoy, 10 targets | every PEP 0.1, sum 2.0 | sum 2.5 | **sum 2.0** |
| 20 targets above every decoy | 0.05 each | 0.075 | **0.05** |
| `p = 1/3` vs `p = 1/2` | exactly half | not proportional | **exactly half** |
| summation identity, 400-row case | sum = estimate | 48.46 vs 50 | **exact** |

All five failed before the change.

### 4.4 Calibration — `validation/pep_entrapment.py`

Six-dataset entrapment, all five predeclared seeds, ~98,700 targets per seed:

| seed | targets | PEP = 0 | known-false at PEP = 0 | weighted absolute error | signed (observed − reported) |
|---:|---:|---:|---:|---:|---:|
| 1 | 98,692 | 0 | 0 | 0.018825 | +0.018825 |
| 2 | 98,810 | 0 | 0 | 0.020400 | +0.020400 |
| 3 | 98,656 | 0 | 0 | 0.016540 | +0.016540 |
| 4 | 98,825 | 0 | 0 | 0.019556 | +0.019556 |
| 5 | 98,777 | 0 | 0 | 0.018467 | +0.018467 |

Audited build, seed 1: 0.018514. **This repair: 0.018825. Calibration did not improve.**

Per bin, seed 1:

| PEP bin | targets | mean reported | adjusted observed | observed / reported |
|---|---:|---:|---:|---:|
| [1e-4, 1e-3) | 10,926 | 0.000366 | 0.007688 | 21.0 |
| [1e-3, 5e-3) | 2,978 | 0.003694 | 0.008126 | 2.20 |
| [5e-3, 0.01) | 1,399 | 0.007148 | 0.016542 | 2.31 |
| [0.01, 0.02) | 1,703 | 0.014092 | 0.030143 | 2.14 |
| [0.02, 0.05) | 1,326 | 0.035445 | 0.060210 | 1.70 |
| [0.05, 0.10) | 1,514 | 0.071995 | 0.081345 | 1.13 |
| [0.10, 0.20) | 1,571 | 0.157861 | 0.189526 | 1.20 |
| [0.20, 0.50) | 3,796 | 0.343256 | 0.359661 | 1.05 |
| [0.50, 1.00] | 73,479 | 0.913349 | 0.934651 | 1.02 |

Every populated bin is optimistic. Per dataset (seed 1) the weighted absolute error ranges from
0.0131 to 0.0307, and the sign is positive on all six.

### 4.5 Classification

The differentiation defect and the arbitrary prior were **implementation errors** and are fixed. The
remaining optimism is an **inherited property of the cumulative estimator**, not a further defect in
the PEP step: the observed/reported ratios above (2.2, 2.3, 2.1, 1.7, 1.1) track the cumulative
q-value ratios on the same data (2.41 at nominal 0.005, 1.81 at 0.01, 1.37 at 0.02, 1.17 at 0.05,
1.09 at 0.10) in both shape and magnitude. Differencing an anti-conservative cumulative curve produces
an anti-conservative local one, and no property of the isotonic fit repairs that.

No correction factor, clamp or smoothing parameter was introduced. **These values are not calibrated
posterior error probabilities and must not be described as such.**

---

## 5. Protein inference

### 5.1 Grouping — root cause and correction

Proteins were collapsed by connected components of the peptide-sharing graph. Sharing a peptide is
not indistinguishability: if peptide `AB` maps to `A` and `B` while peptide `A1` maps only to `A`,
the data separates them, and a chain of shared peptides collapses an entire component.

Two proteins are now grouped only when their **observed peptide sets are identical**, the standard
definition. Subset proteins stay separate because they are distinguishable. A shared peptide counts
as evidence for every group it maps to; no parsimony step is applied, and that is now documented as a
limitation rather than left implicit.

Graph-level unit tests cover: identical evidence groups; a unique peptide keeps two proteins apart; a
subset is not absorbed by its superset; a sharing chain `A–B–C–D` yields four groups, not one; shared
peptides support every group they map to; and grouping is invariant to entry order.

On real single-organism data (`data/F_3.pin`, 105,560 PSMs, seed 1):

| | pre-repair | repaired |
|---|---:|---:|
| protein groups (target + decoy) | 4,182 | 4,282 |
| reported target groups | 1,793 | 1,844 |
| largest reported group, members | **13** | **2** |
| groups with >1 member | 31 | 1 |
| target proteins q<0.01, picked | 1,409 | 1,456 |
| target proteins q<0.01, classic | 1,368 | 1,417 |

Connected-component grouping was merging up to thirteen distinguishable proteins into one
"indistinguishable" group on this dataset. Exactly one genuinely indistinguishable pair remains.

### 5.2 Competition — root cause and correction

`groups[ti].score >= groups[di].score` gave every exact target/decoy protein tie to the target. 200
tied pairs produced 200 target wins, all at q<0.01, minimum q 0.005.

Picking now uses the same fair coin, keyed on the pairing key — the sorted, decoy-stripped member
names, which carry no label. Within one slot, ties between candidate groups are broken by member
names, which are unique per group.

`validation/adversarial_protein.rs`, 200 exactly tied pairs:

| arm | pre-repair | repaired |
|---|---|---|
| target rows first | 200 target / 0 decoy, 200 at q<0.01 | 107 / 93, **0** at q<0.01 |
| decoy rows first | (not run; rule is order-independent by construction) | 107 / 93, identical pick signature |
| seed 2 | 200 / 0 | 115 / 85, different signature |

Minimum target q goes from 0.005 to 0.879. A strictly better group still always wins its pair, at
every seed.

### 5.3 Protein PEP — resolution

The picked path copied the best peptide's PEP into a column named `posterior_error_prob`. A fixture
with peptide PEP 0.123456 produced protein PEP 0.123456.

Picked-protein FDR estimates a **cumulative** error rate over protein groups and no posterior. There
is no protein-level posterior to report on that path, so `ProtGroup::pep` is now `Option<f64>`, it is
`None` for picked inference, and the column reads **`NA`**. On `data/F_3.pin` the top group's reported
protein PEP goes from `0.000078` (the peptide value) to `NA`.

The Bayesian path (`--protein-inference bayesian`) does estimate a protein-level posterior from a
noisy-OR factor graph and still fills the column in. That posterior is **unvalidated**: its inputs are
the peptide PEPs of section 4, which fail calibration, and its α/β/γ parameters were selected on a
PrEST dataset that did not demonstrate nominal 1% protein FDR. It remains available and is marked
experimental.

### 5.4 What is still not established

No truth-based protein-level validation exists. The grouping is now correct by a stated definition
and the competition is now fair, both under unit and adversarial tests, but **no protein q-value in
this repository has been checked against protein-level ground truth.** Protein output should be read
as a ranked list with a target-decoy error estimate whose calibration is untested.

---

## 6. Cross-validation

### 6.1 `--select-c`

**Root cause.** The Cpos/Cneg grid was scored on the same out-of-fold predictions that were later
reported. Every held-out row helped choose the hyperparameters of the model that scored it.

**Correction.** Selection is nested. For each outer fold: split its own training partition again,
train each grid candidate on the inner training rows, score the inner validation rows, take the
highest pooled inner-validation yield, then train the final model on the whole outer training
partition and score the untouched held-out fold.

### 6.2 `--ensemble`

**Root cause.** The per-candidate cross-engine agreement feature was built from a map keyed on
`(ScanNr, Label, Peptide)` over every row, before folds existed.

**Correction.** The key is `(ScanNr, Peptide)`. Whether two engines reported the same candidate is a
property of the searches; the peptide sequence already determines which database it came from.

### 6.3 Results — `validation/adversarial_cv.py`

Two independent attacks per mode. `fold` flips every label in outer fold 0 and requires that fold's
scores and its own selection to be unchanged. `row` flips one row's label and requires that row's own
score to be unchanged; it needs no knowledge of the fold assignment.

| mode | attack | pre-repair | **repaired** |
|---|---|---|---|
| fixed-C | fold | 0/200 scores changed, selection stable | 0/200, stable |
| fixed-C | row | 0/6 own scores changed | 0/6 |
| `--select-c` | fold | **200/200 changed; fold-0 weights 4:1 → 0.25:1** | **0/200; fold 0 keeps 0.25:1** |
| `--select-c` | row | **1/6 own scores changed** | **0/6** |
| `--ensemble` | fold | 0/195 changed | 0/195 |
| `--ensemble` | row | 0/6 changed | 0/6 |

An end-to-end Cargo test (`tests/cv_leakage.rs`) runs the per-row attack for all three modes on every
build.

### 6.4 A correction to the audit's ensemble finding

The audit reported `--ensemble` as leaking. The source-level dependence was real — a global map keyed
on labels — and is now gone. But **no end-to-end attack reproduces it on well-formed input, on either
build.** The audit's earlier ensemble arm used an equal-size fold map, while ensemble input groups by
`ScanNr` alone and deals folds by candidate count; with the correct map the attacked rows are not the
held-out ones. With the right map, both builds show 0 changed scores.

The reason is that on a concatenated search the peptide sequence already determines the label, so the
`Label` component of that key was redundant. A unit test that sweeps all 64 relabelings of a small
ensemble fixture does fail on the previous code and pass on this one, so the dependence was genuine;
its measured exploitability was nil.

Both facts are reported. The property "no ensemble feature depends on any label" now holds
unconditionally and is checkable, which is what the claim needs.

---

## 7. The entrapment FDP — root-cause investigation

### 7.1 The repair did not move it

Six mzML runs re-searched against the native database plus an equally sized foreign plant proteome;
pure foreign assignments are known errors. Five predeclared seeds, at reported q<0.01:

| seed | audited build accepted / FDP | **this repair** accepted / FDP |
|---:|---|---|
| 1 | 19,545 / 0.017865 | 19,545 / 0.017865 |
| 2 | 19,488 / 0.017299 | 19,501 / 0.017288 |
| 3 | 19,532 / 0.018431 | 19,536 / 0.018150 |
| 4 | 19,514 / 0.018838 | 19,514 / 0.018838 |
| 5 | 19,586 / 0.018378 | 19,586 / 0.018378 |
| **mean** | **19,533.0 / 0.018162** | **19,536.4 / 0.018104** |

Full curve, this repair, five-seed means:

| nominal q | accepted | adjusted FDP | FDP / nominal |
|---:|---:|---:|---:|
| 0.001 | 12,946.4 | 0.006146 | 6.15 |
| 0.005 | 17,985.8 | 0.012052 | 2.41 |
| 0.010 | 19,536.4 | 0.018104 | 1.81 |
| 0.020 | 21,142.0 | 0.027455 | 1.37 |
| 0.050 | 23,828.4 | 0.058509 | 1.17 |
| 0.100 | 26,848.4 | 0.108723 | 1.09 |

**Unchanged, and still anti-conservative at every threshold.** That is the expected outcome once
sections 2 and 6 showed the tie and leakage defects have no measurable effect on this data, and
section 4 showed PEPs do not feed back into q-values.

### 7.2 Attribution — `validation/entrapment_rootcause.py`

The probe applies the identical entrapment accounting to scores produced with **no semi-supervised
rescoring at all**: target-decoy competition and q-values computed directly on the best single PIN
feature (`-lnExpect` on all six datasets), using an estimator written out in the script rather than
called from the crate. Seed 1, pooled over the same six datasets:

| nominal q | **raw score** accepted / FDP / ratio | maxiter 1 | maxiter 3 | **maxiter 10 (default)** |
|---:|---|---|---|---|
| 0.001 | 7,203 / 0.0117 / 11.7 | 10,627 / 0.0083 / 8.3 | 12,729 / 0.0047 / 4.7 | 12,454 / 0.0058 / 5.8 |
| 0.010 | 11,728 / 0.0121 / **1.21** | 16,947 / 0.0140 / 1.40 | 18,967 / 0.0184 / 1.84 | 19,545 / 0.0179 / **1.79** |
| 0.050 | 15,876 / 0.0503 / **1.01** | 21,153 / 0.0547 / 1.09 | 23,211 / 0.0571 / 1.14 | 23,813 / 0.0589 / **1.18** |
| 0.100 | 18,656 / 0.0992 / **0.99** | 24,096 / 0.1059 / 1.06 | 26,256 / 0.1084 / 1.08 | 26,796 / 0.1098 / **1.10** |

Three findings.

**A baseline anti-conservatism exists with no rescoring.** At nominal 1% the raw search score is
already 1.21× over. At nominal 5% and 10% it is essentially calibrated (1.01×, 0.99×). So the excess
at loose thresholds is entirely attributable to rescoring, and roughly a third of the excess at 1%
is present before percolator-rs does anything.

**Semi-supervised training accounts for the rest, and scales with the training budget.** 1.21× → 1.40×
→ 1.84× → 1.79× as iterations go 0 → 1 → 3 → 10. The mechanism is not fold leakage — the default path
passes every leakage attack — but the training objective itself: the model is fitted to separate the
target database from the decoy database, and any feature that does so systematically (composition,
length, enzymatic termini after reversal) promotes *all* target-database matches, foreign entrapment
targets included, while promoting no decoys. That breaks the exchangeability target-decoy competition
assumes. It is a property of the Percolator method, not of this implementation; the reference shows
the same behaviour on the same data.

**The entrapment adjustment's own assumption is worth as much as the effect being measured.** The
adjusted FDP divides pure-entrapment target counts by an estimate of the probability that an
incorrect target lands in the foreign database. Estimating it from accepted non-mixed decoys gives
≈0.75; the study's declared database share is 0.50. At nominal 1% the raw arm reads 0.0121 under the
first and 0.0181 under the second:

| nominal q | raw, decoy-estimated *f* | raw, whole-file decoy *f* | raw, declared *f* = 0.5 |
|---:|---:|---:|---:|
| 0.001 | 0.0117 | 0.0050 | 0.0078 |
| 0.010 | 0.0121 | 0.0116 | 0.0181 |
| 0.050 | 0.0503 | 0.0508 | 0.0789 |
| 0.100 | 0.0992 | 0.1029 | 0.1598 |

The decoy-estimated value is the defensible one — decoys sample the same random-match opportunity —
and the study already uses it. But that estimate assumes incorrect **targets** place into the foreign
database at the same rate as decoys. Incorrect target PSMs can match homologous native peptides that
decoys cannot, which would push the true rate *below* 0.75 and make the reported FDP an
**underestimate**. The direction of the residual bias is therefore not established, and the q<0.001
row is unusable in both arms (3 and 8 accepted decoys support the fraction estimate).

### 7.3 Per dataset

Pooling hid real spread. Seed 1, rescored, per dataset at nominal 1%:

| dataset | accepted | adjusted FDP |
|---|---:|---:|
| 09Dec2015-…-5-atrium-P-12hpm-3rd-01 | 814 | 0.0189 |
| 22Oct2014-…-8-MAGs-S-01 | 7,410 | 0.0179 |
| 28May2015-…-22-atrium-S-24H-2nd-01 | 345 | no accepted decoy |
| 28May2015-…-23-atrium-P-24H-2nd-01 | 307 | 0.0000 |
| 28May2015-…-38-MAGs-P-3rd-02 | 4,988 | 0.0151 |
| 9March2015-…-14N-male-02 | 5,681 | 0.0202 |

### 7.4 Conclusion

**Nominal 1% FDR control is not demonstrated on this experiment and the documentation says so.** The
q-value estimator's mathematical guarantee is conditional on target/decoy exchangeability among
competition winners; the entrapment experiment measures a foreign-proteome false-discovery proportion
under its own plug-in assumption. The gap between them is now attributed rather than merely reported:
roughly a third of the excess at nominal 1% is present without any rescoring, the rest tracks the
semi-supervised training budget, and the adjustment's foreign-fraction assumption can move the number
by a factor of 1.5 in either direction.

---

## 8. Multi-seed and cross-dataset results

Four datasets, five seeds, matched Rust and C++ runs. Target PSMs at q<0.01:

| dataset | Rust, seeds 1–5 | mean | C++, seeds 1–5 | mean |
|---|---|---:|---|---:|
| MSFragger | 1367, 1377, 1312, 1361, 1340 | 1351.4 | 1388, 1477, 1399, 1409, 1383 | 1411.2 |
| Sage | 25806, 25803, 25804, 25791, 25779 | 25796.6 | 25795, 25788, 25784, 25790, 25790 | 25789.4 |
| Tide | 27640, 27665, 27616, 27624, 27624 | 27633.8 | 27617, 27656, 27610, 27612, 27670 | 27633.0 |
| yeast | 1150, 1144, 1124, 1180, 1159 | 1151.4 | 1147, 1137, 1059, 1105, 1111 | 1111.8 |

Target peptides at q<0.01:

| dataset | Rust mean | C++ mean |
|---|---:|---:|
| MSFragger | 1055.6 | 1080.6 |
| Sage | 11248.8 | 11323.8 |
| Tide | 19720.0 | 19732.4 |
| yeast | 921.6 | 897.0 |

**Every one of these 40 Rust counts is identical to the audited build's**, and all 80 C++ artifacts
reproduce byte-for-byte. All 80 Rust artifacts differ by exactly one column: `posterior_error_prob`.
On Tide seed 1, 40,487 of 42,330 PEP cells changed while `PSMId`, `score`, `q-value` and `peptide`
changed in zero cells.

### 8.1 A post-processing mismatch in the reference comparison

Percolator 3.09's `--post-processing-tdc` does **not** perform spectrum-level competition on
concatenated input; it changes how q-values are assigned when the input came from separate searches.
On concatenated multi-candidate PINs the reference therefore reports every candidate while
percolator-rs reports one winner per precursor. Candidates per precursor:

| dataset | candidates/precursor | comparison is like-for-like |
|---|---|---|
| Tide | 1 | yes |
| Sage | 1 (9 precursors with 2) | yes |
| MSFragger | mostly 3 | **no** |
| yeast | mostly 2 | **no** |
| PXD032157 | 5 | **no** |

This is exactly the split between the strong-agreement and weak-agreement datasets the audit
observed. It is not, however, the whole explanation. Running percolator-rs with
`--no-psm-competition` to match the reference's list (seed 1, q<0.01 Jaccard):

| dataset | competed (default) | uncompeted (matched to C++) |
|---|---:|---:|
| MSFragger | 0.877 | 0.914 |
| yeast | 0.934 | 0.909 |

Aligning post-processing helps MSFragger and hurts yeast. The residual 0.88–0.93 agreement on these
two datasets is a genuine model/solver difference, not a bookkeeping artefact.

### 8.2 PXD032157, the pathological case

The audit required this dataset to be reported rather than excluded. It is a **concatenated** search
reporting five candidates per precursor — verified directly: every precursor has exactly five
candidates split between labels, 93.4% of precursors carry both. percolator-rs's input contract holds;
the reference's auto-detection of "separate" input on this file is a misdetection on its side.

Per-file yield is extremely skewed: across the 65 files, target PSMs at q<0.01 range from 142 to
9,938. File `03May2016-…-12-WB-virgin-neg-4-18-16-01` — the file the earlier study used — is one of
the five worst.

| file | Rust q<0.01 | C++ q<0.01 | matching PSMs | score Spearman | Jaccard |
|---|---:|---:|---:|---:|---:|
| `03May2016-…-01` (low yield) | 152 | 0 | 36,802 | 0.413 | 0.000 |
| `22Oct2014-…-9-MAGs-P-01` (high yield) | 9,938 | 9,762 | 52,804 | 0.991 | 0.962 |

On the low-yield file the reference reported 92,989 target rows to percolator-rs's 19,082 — it did not
compete candidates — so the two are scoring different lists and the disagreement is not a modelling
result. On a high-yield file from the same dataset, agreement is 0.96.

**PXD032157 is a development dataset.** It was used to tune class weights, folds, score scaling and
performance. Its counts are a reproducible baseline, not an accuracy estimate.

---

## 9. Mutation test of the validation suite — `validation/mutation_test.py`

Twelve known scientific defects reintroduced in a throwaway worktree. **All twelve are caught.**

| mutation | caught by |
|---|---|
| input-order-dependent PSM ties | 4 tests |
| missing finite-sample `+1` on the reported path | 15 tests |
| broken tie grouping (every row its own group) | 3 tests |
| PEP prior added back (`+0.5/N`) | 6 tests |
| PEP floor raised to 0.01 | 1 test |
| normalization fitted on all rows | 2 tests |
| initial direction chosen with held-out labels | 1 test |
| C selection given the held-out fold | 1 test (`tests/cv_leakage.rs`) |
| ensemble agreement feature keyed on labels | rejected at compile time; the 64-relabeling unit test fails when it does compile |
| target-favouring protein ties | 2 tests |
| connected-component protein grouping | 5 tests |
| peptide PEP emitted as protein PEP | 2 tests |

The two PEP mutations that passed all 82 tests in the audit are now caught, and tie grouping went from
1 detecting test to 3 after a direct requirement was added.

---

## 10. New canonical baseline and performance

PXD032157, all 65 files, seed 1, `--canonical` defaults, outputs written to local ext4.

| metric | audited build | **this repair** | change |
|---|---:|---:|---|
| target PSMs q<0.01 | 106,795 | **106,823** | +28 (+0.03%) |
| target peptides q<0.01 | 35,866 | **35,886** | +20 (+0.06%) |
| wall, N=4 (3-run median) | 16.1 s | **15.35 s** | −4.7% |
| wall, sequential (3-run median) | 51.1 s | **49.36 s** | −3.4% |
| peak RSS, N=4 | 0.76 GiB | **0.75 GiB** | −1.3% |
| largest file, 1 thread (3-run median) | 2.05 s | **1.56 s** | — |
| largest file, 3 threads (3-run median) | 1.39 s | **0.98 s** | — |

The count change comes from the fair coin resolving a small number of exact ties differently. The
timing change was not sought; the largest single contributor is the hash map that replaced an ordered
map in competition. The single-file numbers are not directly comparable to the previously published
ones — the input was warm in page cache after repeated benchmark runs — so they are reported as
"no regression", not as a speedup.

Regression gates: `tests/regression.sh` 117 PSM / 43 peptide (unchanged), `model_regression.sh`,
`selection_regression.sh`, `ensemble_regression.sh` and `protein_regression.sh` all pass.

---

## 11. What still fails, and what remains untested

**Fails.**

* Nominal FDR control on the signal-present entrapment experiment, at every threshold. Attributed in
  §7 but not fixed, and not fixable without changing the method.
* PEP calibration on the same data. Every populated bin is optimistic; weighted absolute error 0.017–0.020
  across five seeds.

**Untested.**

* Protein-level q-values against protein ground truth. The grouping and competition are now correct
  under their stated definitions, but no truth-based protein study exists.
* The Bayesian protein posterior. Its α/β/γ come from a PrEST study that did not establish nominal
  protein FDR.
* `--rescore-model mlp`, `--join`, `--rt-features` beyond fold-isolation tests, and `--auto-model`
  beyond source inspection and the fold-attack it passes.
* Generalization. Four reference datasets plus six entrapment runs, all available during development.
  No untouched external collection with predeclared success criteria is reserved.

**Known limitations of the corrected methods.**

* Protein inference applies no parsimony. Shared peptides are evidence for every group they map to,
  and subset proteins are reported separately.
* Ties are resolved by a seeded coin, so a dataset with many exact ties will show seed-to-seed
  variation in its results. That variation is real and was previously hidden by row order.
* `pi0` is fixed at 1. Separate target/decoy searches (mix-max) are not supported and produce an
  unvalidated number rather than an error.
