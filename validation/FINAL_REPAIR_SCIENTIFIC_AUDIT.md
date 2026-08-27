# Fresh adversarial scientific audit after the joined-input and protein-grouping repairs

Audit date: 2026-08-27<br>
Audited commit: `b2501f03ed727a74a633b176d3ef9efc4d41d88c`<br>
Latest repair commits: `2dcba362ca66ff47c22f05d31db7ebe71a01d921`, `b2501f03ed727a74a633b176d3ef9efc4d41d88c`<br>
Frozen binary SHA-256: `be9bf670bfd69df4dc3ba3b8be6c4c164acaf56a91f4a2819d115f49604b2c45`

## Executive verdict

**No. `percolator-rs` is still not scientifically defensible as an unqualified source of calibrated
PSM, PEP, or protein-level confidence.**

The two latest repairs are real and their narrow claims survive independent attacks:

- fixed file names plus file-order, row-order, target/decoy-order, reverse, deterministic-shuffle,
  exact-tie, near-tie, and four-file joined permutations now give one identical result;
- the repaired protein grouping relation correctly handles new indistinguishable, distinguishable,
  subset, shared-peptide, target/decoy, and insertion-order graphs;
- mutations that remove joined row canonicalization, joined file canonicalization, complete reported
  peptide mapping, target/decoy class separation, or label-free ensemble features are caught by
  behavioral tests.

Those facts do not rescue scientific validity. Fresh counterexamples establish all of the following:

1. **PSM competition does not preserve the declared null probability under duplicate candidates.**
   Exact duplication of an otherwise identical target row changes a 56 target / 45 decoy null result
   with zero q<0.01 discoveries into 101 target / 0 decoy winners and 101 false q<0.01 discoveries.
2. **Joined results remain dependent on the lexical source name.** The same file bytes, accessed via
   differently named symlinks, can flip an exact tie and move strict q<0.01 discoveries from 0 to 101.
3. **Protein evidence can still be lost before the repaired union runs.** The union uses only PSM
   competition winners; two identical peptide candidates with complementary protein mappings retain
   only the mapping selected by the tie seed.
4. **Picked-protein group pairing has a serialization collision.** The distinct groups
   `{LEFT, RIGHT}` and `{LEFT|RIGHT}` share the same unescaped key, so one target group is silently
   removed from competition.
5. **Current protein truth validation fails severely.** In the held-out PrEST A and B samples,
   picked protein q≤0.01 has raw known-absent FDP 45.92% and 48.08%; the predefined count-adjusted
   values are 53.32% and 55.78%.
6. **Signal-present PSM FDR remains anti-conservative without methodology tuning.** Mean adjusted FDP
   is 1.8104% at nominal q<1%, and exceeds nominal at all six tested thresholds.
7. **PSM PEPs remain optimistic in every populated entrapment bin.** Weighted signed and absolute
   calibration error are both +0.018685; 217 known-false PSMs have PEP<0.001.
8. **The PEP implementation still violates its exact mass claim at an accepted parameter value.** At
   `p=1e-15`, target PEP mass is 5,000 times the implemented false-count estimate because of the
   numerical floor.
9. **The committed protein evaluator cannot read the current picked-protein `NA` schema.** Therefore
   current protein output can pass the ordinary suite while its truth benchmark is not runnable.

No C++ Percolator result was used as a correctness oracle in this audit.

## Audit discipline and provenance

The primary checkout was already dirty when the audit began. Existing modifications and prior audit
files were preserved. No production source under `src/`, no test under `tests/`, and no build
configuration was modified. New audit drivers and probes live under `validation/`; all production
mutations were confined to a detached throwaway Git worktree and restored after each experiment.

The audit reviewed:

- `validation/POST_REPAIR_AUDIT.md` and the earlier independent/repair records;
- every prior adversarial probe named by those records;
- the complete diffs of repair commits `2dcba36` and `b2501f0`;
- the current competition, TDC/PEP, CV, joined parsing, grouping, picked-protein, Bayesian-protein,
  output, and benchmark-evaluator paths.

The audited binary was frozen before adversarial execution. `cargo test --release --all-targets`
passed, and the current test graph lists 126 tests. A passing baseline is therefore not being
confused with absence of defects.

## Failure classification

| ID | classification | minimized failure |
|---|---|---|
| I1 | **IMPLEMENTATION DEFECT** | PSM exact-tie selection samples candidate rows, so duplicating one label changes its win probability from 1/2 to 94/95 and creates 101/101 false discoveries |
| I2 | **IMPLEMENTATION DEFECT** | joined numeric source identity is assigned after sorting by lexical file name; renaming/symlinking identical bytes changes folds and tie coins, including a 0 → 101 cutoff change |
| I3 | **IMPLEMENTATION DEFECT** | protein mappings are unioned over `reported_indices` after PSM competition; a losing equivalent row's mapping disappears |
| I4 | **IMPLEMENTATION DEFECT** | picked-protein member sets are serialized with unescaped `|`; distinct sets collide and one group is not picked |
| I5 | **IMPLEMENTATION DEFECT** | the public PEP conservation/scaling claim is false at accepted `p=1e-15` because a fixed per-target floor dominates the estimate |
| I6 | **IMPLEMENTATION DEFECT** | decoy PEP display values change under single-file row permutations even when winners, scores, q-values, and target PEPs do not; no decoy posterior claim is made, so impact is presentation-only |
| M1 | **ESTIMATOR/METHODOLOGY LIMITATION** | direct TDC arithmetic is correct but the untuned signal-present design has adjusted FDP above nominal at every threshold |
| M2 | **ESTIMATOR/METHODOLOGY LIMITATION** | PSM PEP construction has ordinary algebraic consistency but is not a calibrated posterior estimator on entrapment truth |
| M3 | **ESTIMATOR/METHODOLOGY LIMITATION** | grouping can produce unmatched target and decoy group topologies; the picked method retains unpaired groups yet applies a 1:1 TDC estimator, and current PrEST FDR fails badly |
| M4 | **ESTIMATOR/METHODOLOGY LIMITATION** | default Bayesian protein probabilities are severely miscalibrated on held-out PrEST truth; frozen selected parameters transfer much better but are benchmark-specific |
| V1 | **VALIDATION-DESIGN LIMITATION** | the predefined complete-null design has only ten relabels per input and cannot establish calibration at small nominal rates despite observing 0/30 rejections |
| V2 | **VALIDATION-DESIGN LIMITATION** | the green ordinary suite has no gates for I1-I6, empirical FDR, PSM PEP calibration, or current protein truth calibration |
| V3 | **VALIDATION-DESIGN LIMITATION** | `bench/protein_calibration/report.py` calls `float("NA")`; the committed truth workflow is stale relative to the repaired output schema |
| V4 | **VALIDATION-DESIGN LIMITATION** | the historical mutation harness has drifted: one protein mutation is no longer applicable and its old ensemble mutation only fails to compile |
| V5 | **VALIDATION-DESIGN LIMITATION** | all four compact compatibility datasets and both truth designs were available during development; there is no untouched biological generalization set |
| U1 | **UNRESOLVED** | the relative causal contributions of adaptive training, TDC exchangeability violations, candidate reporting, and entrapment adjustment to the real-data FDR excess are not identified |
| U2 | **UNRESOLVED** | the frequency and practical impact of exact candidate multiplicity and source-alias sensitivity in external PIN producers are not quantified |
| U3 | **UNRESOLVED** | the PrEST picked-protein failure is consistent with asymmetric/unpaired group topology, but the benchmark failure has not been decomposed into one unique causal mechanism |

## 1. Single-file PSM competition and tie invariance

A new 233-spectrum fixture contained one exactly tied target/decoy pair per spectrum. Six independent
layouts were used: original, reversed, target-first, decoy-first, and deterministic shuffles 1701 and
9907.

All arms produced exactly 119 target and 114 decoy winners. Winner identities, printed scores,
q-values, target PEPs, and q<0.01 discoveries were invariant. This independently confirms the narrow
single-file permutation repair.

The complete output was not invariant because 2–4 **decoy** PEP fields changed in five arms. The
largest observed change was approximately 0.276. The implementation makes no posterior-error claim
for decoy rows, so this is I6 rather than evidence that target confidence changed.

The broader PSM competition contract fails I1. A second new fixture had 101 spectra:

| fixture | unique scientific candidates per spectrum | emitted rows per spectrum | target / decoy winners | minimum target q | targets q<0.01 |
|---|---:|---:|---:|---:|---:|
| balanced | 1 target + 1 decoy | 2 | 56 / 45 | 0.821429 | 0 |
| target duplicated | same 1 target + 1 decoy | 94 identical target + 1 decoy | 101 / 0 | 0.009901 | 101 |

The tie code draws uniformly over `k` rows. It therefore gives a target probability `t/k`, not the
declared `p=0.5`, when a candidate is duplicated. The downstream q scan correctly computes q-values
for a winner list that the competition implementation has made non-null.

## 2. Joined-input permutation invariance

A new four-file trained fixture included ordinary rows, exact target/decoy ties, and one-ULP near
ties. The following were permuted jointly and independently: file order, row order, target/decoy
layout, full reversal, and deterministic shuffles 27182 and 31415.

All six arms produced one result: 218 target and 114 decoy winners, byte-identical target/decoy
statistics, and zero q<0.01 discoveries. The historical joined-input counterexample also passes on
the current build. Mutations removing row canonicalization or file canonicalization are killed by
behavioral joined tests.

The repair is nevertheless path-sensitive. The same four file bytes, reached through symlinks whose
lexical names change the sort order, changed all 332 printed PSM records and three winner identities.
A minimized no-training cutoff fixture isolates the consequence:

| access path | boundary winner | targets q<0.01 |
|---|---|---:|
| original names | decoy | 0 |
| renamed symlinks to the same bytes | target | 101 |

Thus argument and row permutations of a fixed named multiset are repaired; content identity under a
scientifically irrelevant source rename is not. This is why the verdict is moderate rather than
strong.

## 3. Q-value/TDC arithmetic

A standalone probe, outside the Cargo test graph, compared current output to a separately written
TDC+ oracle over:

- every target/decoy label pattern for lengths 0 through 10;
- five score-partition families, including all ties, mixed ties, strict order, and one-ULP near ties;
- `p` values 0.2, 1/3, 0.5, and 0.8;
- 40,928 full q-value cases;
- 327,424 optimized count/mask comparisons against materialized q-values;
- exact strict boundaries at 0.001, 0.005, 0.01, 0.02, 0.05, and 0.10.

Every q-value, fast count, fast mask, tie boundary, finite-sample `+1`, opportunity ratio, and reverse
cumulative minimum matched. Q-value arithmetic is the strongest component, conditional on receiving
a scientifically valid winner list and a valid null-target probability.

## 4. Complete-null behavior

The predefined design was rerun unchanged: three specified PXD032157 inputs, relabel seeds
1001–1010, model seed 1, and thresholds 0.001, 0.005, 0.01, 0.02, 0.05, and 0.10.

- 30/30 executions completed;
- 0/30 had any false target at any threshold;
- every false-target count was zero.

This is a successful conservative complete-null result, not evidence of calibrated small FDR. Each
input has only ten relabels, whose two-sided Wilson upper bound after 0/10 events is 27.75%. Even an
unqualified 0/30 pool has an upper bound around 11.35%, before accounting for within-dataset
dependence. It also does not exercise the I1 duplicate-candidate construction.

## 5. Signal-present entrapment and FDR calibration

The six predefined Comet PINs, seeds 1–5, fixed canonical parameters, fixed-C mode, thresholds, and
the predeclared entrapment-fraction calculation were rerun without tuning or excluding a seed.

| nominal q | mean accepted targets | mean pure entrapment | mean adjusted FDP | FDP / nominal |
|---:|---:|---:|---:|---:|
| 0.001 | 12,946.4 | 60.6 | 0.6146% | 6.15× |
| 0.005 | 17,985.8 | 170.8 | 1.2052% | 2.41× |
| 0.010 | 19,536.4 | 258.6 | 1.8104% | 1.81× |
| 0.020 | 21,142.0 | 430.4 | 2.7455% | 1.37× |
| 0.050 | 23,828.4 | 1,046.6 | 5.8509% | 1.17× |
| 0.100 | 26,848.4 | 2,215.2 | 10.8723% | 1.09× |

At q<0.01, the five seed-specific adjusted FDPs were 1.7865%, 1.7288%, 1.8150%, 1.8838%, and
1.8378%. This is a stable failure of the available signal-present calibration design. It is M1, not
a falsification of the q-scan arithmetic proved in section 3. U1 remains because the design cannot
uniquely attribute the excess.

## 6. PSM PEP implementation and calibration by bin

Ordinary PEP arithmetic passed finite/range/monotonicity and mass checks. In the standalone mixed
case, target PEP mass was 2.99999999999999956 for an estimated false count of three.

The public edge case still fails. Five all-target scores at `p=1e-15` give estimated false count
`1.00000000000000106e-15`, but the fixed `1e-12` floor gives PEP sum
`4.99999999999999970e-12`, a ratio of 5,000.

Pooled signal-present calibration uses 493,760 target PSMs from all five predefined seeds. Empty bins
are reported explicitly.

| PEP bin | targets | mean reported PEP | adjusted observed error | observed − reported |
|---|---:|---:|---:|---:|
| [0, 1e-12) | 0 | — | — | — |
| [1e-12, 1e-6) | 0 | — | — | — |
| [1e-6, 1e-4) | 0 | — | — | — |
| [1e-4, 0.001) | 54,378 | 0.000386 | 0.004930 | +0.004543 |
| [0.001, 0.005) | 15,887 | 0.002707 | 0.009335 | +0.006629 |
| [0.005, 0.01) | 5,750 | 0.006783 | 0.017739 | +0.010956 |
| [0.01, 0.02) | 7,696 | 0.015852 | 0.031274 | +0.015421 |
| [0.02, 0.05) | 9,091 | 0.033330 | 0.056601 | +0.023271 |
| [0.05, 0.10) | 5,374 | 0.073874 | 0.090656 | +0.016782 |
| [0.10, 0.20) | 9,231 | 0.145705 | 0.166029 | +0.020324 |
| [0.20, 0.50) | 19,260 | 0.345587 | 0.372132 | +0.026546 |
| [0.50, 1.0000001) | 367,093 | 0.913066 | 0.934018 | +0.020952 |

Every populated bin is optimistic. Weighted signed error and weighted absolute error are both
0.01868545. There are 217 known-false PSMs below PEP 0.001; the smallest known-false PEP is 0.000262.
The values are therefore not defensible as calibrated posterior error probabilities.

## 7. Cross-validation isolation

A fresh 411-spectrum fixture attacked one complete held-out fold (137 spectra). All labels in the
fold were flipped; all but 11 held-out feature vectors were replaced by ±1e12 outliers; reverse and
deterministic-shuffle arms were also run. In ensemble mode, three overlapping engine inputs were
used. Tests compared only the held-out rows whose model must not have seen the attack.

| mode | held-out score changes after label attack | sentinel changes after other held-out outliers | reverse/shuffle changes | fold-1 selection changed |
|---|---:|---:|---:|---:|
| fixed-C | 0 | 0/11 | 0 / 0 | n/a |
| `--select-c` | 0 | 0/11 | 0 / 0 | no |
| `--ensemble` | 0 | 0/11 | 0 / 0 | n/a |

Source inspection confirms normalization/direction are fitted on outer-training rows, C selection is
nested inside each outer fold, ensemble agreement features are label-free, and overlapping engine
copies are assigned together. Fixed-C, select-C, and ensemble CV isolation receive strong evidence.

## 8. Protein grouping

A new hand graph exercised all requested structures:

- indistinguishable proteins `IX/IY` with identical two-peptide evidence;
- a distinguishable `A/B/C` shared-peptide component;
- strict `SUB/SUPER` evidence;
- paired target/decoy evidence;
- an independent target/decoy exact score tie;
- an unconnected protein;
- original, reverse, and two deterministic shuffle orders.

All four layouts produced the same 11 groups and complete signatures. Mixed target/decoy proteins
did not co-group. The core repaired grouping relation is correct on this new graph.

End-to-end grouping still fails I3. With `--no-psm-competition`, two rows for one peptide mapping to
`LOSS_A_PROT` and `LOSS_B_PROT` produce their complete union. Under default competition the rows are
equivalent candidates for one spectrum, so only one reaches `reported_indices`: seed 1 retains A and
loses B; seed 3 loses A and retains B. Reversing insertion order at a fixed seed does not change the
result, but changing a legitimate tie seed changes the biological protein graph. The latest
integration repair test switches competition off and therefore cannot detect this stage-order gap.

## 9. Protein target/decoy competition

Fairness of the exact tie coin passes. A new population of 509 paired groups gave 232 target wins at
seed 37 and 25,415 target wins out of 50,900 competitions over 100 seeds (49.93%). Reversing entries
did not change a fixed-seed signature; changing the seed did.

The picked grouping key fails I4. `{LEFT, RIGHT}`, `{LEFT|RIGHT}`, and
`{DECOY_LEFT, DECOY_RIGHT}` create three valid groups, but the two target sets serialize to the same
`LEFT|RIGHT` key. Only one target group is marked picked.

More importantly, current PrEST truth validation fails protein FDR. The repaired current binary was
run on all 12 predefined samples, with replicate 1 calibration, replicate 2 validation, replicate 3
test, and no parameter retuning. At q≤0.01 on final test:

| vial | method | accepted | known absent | raw known-absent FDP | predefined adjusted FDP |
|---|---|---:|---:|---:|---:|
| A | picked | 331 | 152 | 45.92% | 53.32% |
| B | picked | 364 | 175 | 48.08% | 55.78% |
| A+B | picked | 389 | 14 | 3.60% | 4.98% |
| blank | picked | 0 | 0 | — | — |
| A | Bayesian fixed | 197 | 22 | 11.17% | 12.97% |
| B | Bayesian fixed | 203 | 27 | 13.30% | 15.43% |
| A+B | Bayesian fixed | 350 | 5 | 1.43% | 1.98% |
| A | Bayesian selected/frozen | 173 | 3 | 1.73% | 2.01% |
| B | Bayesian selected/frozen | 179 | 1 | 0.56% | 0.65% |
| A+B | Bayesian selected/frozen | 302 | 0 | 0 | 0 |

The picked A run, for example, forms 1,662 groups (900 target, 762 decoy) but only 1,092 picked
buckets. Whole-group pairing retains unmatched groups and then uses a fixed 1:1 TDC estimator. The
truth failure is definitive; U3 records that this audit has not assigned all of it to a single
mechanism.

## 10. Protein PEP semantics and calibration

Picked protein output consistently reports `NA`; it no longer mislabels a best-peptide PEP as a
protein posterior. Bayesian output is numeric. Those schema semantics are correct.

Scientific probability calibration is mixed and configuration-dependent. On final-test PrEST truth:

| Bayesian method | groups | 10-bin ECE | Brier | [0,0.1) mean PEP / observed error |
|---|---:|---:|---:|---:|
| fixed defaults | 3,234 | 0.39445 | 0.22444 | 0.01209 / 0.09935 |
| frozen selected parameters | 3,234 | 0.02264 | 0.02422 | 0.00154 / 0.00315 |

The selected configuration transfers substantially better to the held-out replicate, but middle
bins have only 7–31 groups and the parameters were selected on this same benchmark family. Default
Bayesian PEP is not calibrated. Picked inference supplies no PEP at all. A general protein-PEP claim
therefore fails even though the current column meanings are honest.

The committed evaluator is currently unusable: it attempts `float(row["posterior_error_prob"])` and
terminates on the first picked `NA`. The audit's truth accounting was separately implemented from the
manifest and ground-truth table.

## 11. Multi-seed reproducibility and cross-dataset behavior

Four compact datasets (Tide, MSFragger, Sage, and yeast) were run at seeds 1–5, twice per seed.
Every same-seed target and decoy output was byte-identical. Debug/release and thread-count checks are
separate in section 12.

Different seeds produce real cutoff variability:

| dataset | q<0.01 target range over five seeds | SD | q<0.001 range |
|---|---:|---:|---:|
| Tide | 27,616–27,665 | 19.50 | 24,601–24,849 |
| MSFragger | 1,312–1,377 | 25.85 | 0–1,078 |
| Sage | 25,779–25,806 | 11.46 | 24,747–24,829 |
| yeast | 1,124–1,180 | 20.51 | 0–0 |

The implementation is deterministic conditional on bytes and seed. It is not seed-insensitive,
especially at a sparse threshold (MSFragger q<0.001). The four parsers/datasets execute successfully,
but without untouched truth they establish compatibility and numerical stability, not biological
generalization.

## 12. Performance/correctness equivalence

Two independent checks addressed optimized paths without using another implementation as an oracle:

1. Section 3 compared optimized q count/mask/reused-buffer paths to independently materialized
   q-values in 327,424 cases.
2. A 14.6 MB Tide PIN was run with the frozen release binary serially and with three threads, and
   with a current debug build serially. All target and decoy files were byte-identical. Observed wall
   times were 0.364 s release/serial, 0.219 s release/three-thread, and 4.828 s debug/serial.

This is moderate evidence that current optimization and parallel paths preserve this computation.
It is not a proof on other CPUs, modes, or datasets, and byte equivalence to a debug build cannot make
an invalid estimator scientifically correct.

## 13. Mutation testing and validation-suite quality

Five compilable latest-repair mutants were applied one at a time in a detached worktree. The final
result was evaluated only by behavioral test failure, not by compilation failure:

| mutation | behavioral result |
|---|---|
| remove joined row canonicalization | caught by trained joined permutation test |
| remove joined file canonicalization | caught by all three joined permutation tests |
| revert to representative-only protein mapping | caught by end-to-end protein grouping test |
| collapse target/decoy proteins with identical evidence | caught by two protein graph tests |
| put held-out labels back in ensemble agreement features | caught by independent-label ensemble feature property |

The legacy harness reports 11/12 “caught,” but one result is only a compiler rejection and the
connected-component mutation is no longer applicable after the source type changed. Fixed-C and
select-C leakage mutants are killed behaviorally. The fresh five-mutant audit is the defensible
latest-repair result.

Suite quality still fails scientific validation. The unmutated 126-test suite passes with I1-I6,
M1-M4, the broken protein evaluator, no empirical calibration gate, and no untouched truth set. The
suite is good at protecting known local repairs but not sufficient to certify scientific validity.

## 14. Final verdicts

| requested area | verdict | reason |
|---|---|---|
| PSM COMPETITION | **FAILED VALIDATION** | narrow two-way permutation repair passes, but exact duplicate target candidates change null wins to 101/0 and create 101 false q<0.01 discoveries |
| JOINED-INPUT INVARIANCE | **MODERATE EVIDENCE** | every requested fixed-name permutation passes; same bytes under renamed paths can still flip 0 → 101 discoveries |
| Q-VALUE IMPLEMENTATION | **STRONG EVIDENCE** | 40,928 independent-oracle cases, 327,424 fast-path checks, all tie/boundary/p cases pass; conditional on valid winners and `p` |
| COMPLETE-NULL BEHAVIOR | **MODERATE EVIDENCE** | 0/30 at all thresholds, but narrow and underpowered; it misses the multiplicity null |
| FDR CALIBRATION | **FAILED VALIDATION** | untuned signal-present FDP is above nominal at all thresholds, including 1.8104% at nominal 1% |
| PEP IMPLEMENTATION | **FAILED VALIDATION** | ordinary algebra passes, but accepted extreme `p` breaks exact mass by 5,000× and decoy display PEPs are permutation-unstable |
| PEP CALIBRATION | **FAILED VALIDATION** | every populated PSM bin is optimistic; signed/absolute error +0.018685; 217 known false below 0.001 |
| FIXED-C CV | **STRONG EVIDENCE** | complete-held-out-fold labels, outliers, permutations, source inspection, and mutation all support isolation |
| SELECT-C CV | **STRONG EVIDENCE** | nested selection and held-out scores/selections survive the new attacks and leakage mutation |
| ENSEMBLE CV | **STRONG EVIDENCE** | label-free features, grouped copies, new three-engine held-out attacks, and a compilable label-leak mutant pass |
| PROTEIN GROUPING | **FAILED VALIDATION** | core graph relation passes, but supported default competition loses peptide-to-protein mappings before the union |
| PROTEIN COMPETITION | **FAILED VALIDATION** | fair exact ties pass; unescaped keys lose groups and current picked q≤0.01 fails held-out PrEST truth by tens of percentage points |
| PROTEIN PEP | **FAILED VALIDATION** | picked `NA` semantics are correct, but default Bayesian PEP is severely miscalibrated and selected transfer is benchmark-specific |
| MULTI-SEED REPRODUCIBILITY | **MODERATE EVIDENCE** | exact conditional reproducibility passes; cross-seed cutoff sensitivity is material in MSFragger/yeast |
| CROSS-DATASET GENERALIZATION | **WEAK EVIDENCE** | four known formats run, but there is no untouched truth-bearing dataset and both truth studies show limitations/failures |
| PERFORMANCE/CORRECTNESS EQUIVALENCE | **MODERATE EVIDENCE** | debug/release/thread output identity and independent fast-path equivalence pass, but coverage is finite and estimator validity is separate |
| VALIDATION-SUITE QUALITY | **FAILED VALIDATION** | local mutation protection is good, yet the green suite misses current discrete defects and failed calibration; protein truth evaluator is broken |

## What still prevents scientific defensibility?

Four independent barriers are sufficient on their own:

1. The winner list is not guaranteed to satisfy the null probability assumed by TDC. I1 supplies a
   complete-null construction with 101 false discoveries at nominal q<0.01.
2. On the fixed real signal-present design, reported PSM q-values and PEPs are empirically
   anti-conservative. Correct q arithmetic does not imply a valid estimator when its upstream and
   exchangeability assumptions fail.
3. Protein output is neither structurally valid for all supported inputs nor calibrated on current
   truth. The latest mapping union happens too late, the picked key is non-injective, and current
   picked protein q-values fail the held-out PrEST standard catastrophically.
4. Validation is not an enforceable scientific gate. It lacks the discovered properties and
   calibration criteria, has no untouched generalization set, and its committed protein evaluator
   cannot consume current output.

The defensible narrow statement is:

> `percolator-rs` has a strongly validated standalone direct-TDC q scan, robustly isolated tested CV
> modes, and substantive fixes for fixed-name joined permutations and the core protein grouping
> equivalence relation. Its general PSM competition contract, calibrated FDR/PEP interpretation, and
> end-to-end protein inference remain scientifically invalid or unvalidated.

## Preserved evidence

Raw evidence root:
`/run/media/andrej-rumenovski/New Volume/percolator_rs_final_audit_20260827`

| evidence | SHA-256 |
|---|---|
| `fresh-adversarial-v3.json` | `3a32a2e96d7707b86cfda31c89be76e8a36eb704ab2ab399532028a32237e4b5` |
| `final-repair-stats-probe.log` | `1f1dc526c847d25d5f2dd6c35221b714f768f4d166700c57301aef18d8526206` |
| `final-repair-protein-probe-v2.log` | `6fce227b8fdd0b841308f66d53c80c1da494932b1b59d0997bbfea780f2c6780` |
| `null-predefined-current/manifest.json` | `0f56e9e6cfe2ef16c3150a28b859a2363a079fde1a228b7ef97d1b2b68fd57b0` |
| `empirical-current/manifest.json` | `7fce5d2e51e96b1ab37515aedeb19b046a6dbf9d8591721e841c9c9aa672d3d4` |
| `performance-equivalence/manifest.json` | `1259d457c303c08060faa7b055d6eff12284de56ba7a4c436e85b90c726a24d9` |
| `final-repair-mutations-v4.json` | `4020fc4455c7275d1c5068963c74d937d49f47bb70039b5f8dd9650bb966c295` |
| `historical-post-repair-replay.json` | `55218af58b183f1eb70728099c3e838e496a354552ff11298a3e5c967e774081` |
| `historical-mutations-current.json` | `21aeecfa5500745d07af913e5cd4fc213a3f163b090b416d2ae0b1a3b5c9a568` |

Audit probe sources are `final_repair_adversarial.py`, `final_repair_stats_probe.rs`,
`final_repair_protein_probe.rs`, `final_repair_empirical.py`, `final_repair_performance.py`, and
`final_repair_mutation.py`.
