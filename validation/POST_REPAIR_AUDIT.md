# Independent post-repair scientific audit

Audit date: 2026-08-26  
Candidate: `63d4783748356d65a0be2e065859a13bd987027b` (`main`)  
Immediate pre-repair baseline: `b38c0db`  
Frozen candidate binary SHA-256: `5c1267e85e09b3b0ef493d33dd401a839ab37ee14de9b082471ecb7ff902486b`  
Reference used only for agreement, never as an oracle: C++ Percolator 3.09.0  
Production changes during this audit: **none**

This audit was performed after, and independently of, the repair report in
[`SECOND_REPAIR.md`](SECOND_REPAIR.md). Its claims and fixtures were treated as hypotheses to attack,
not as evidence. All long-running experiments used the frozen binary above. Temporary diagnostic and
mutation builds lived in isolated Git worktrees and were removed after use.

### Evidence and reproduction

Machine-readable evidence is preserved outside the repository under
`/run/media/andrej-rumenovski/New Volume/percolator_rs_postrepair_audit_20260826`. The principal
artifacts and SHA-256 digests are:

| evidence | relative path | SHA-256 |
|---|---|---|
| frozen binary | `frozen/percolator-rs` | `5c1267e85e09b3b0ef493d33dd401a839ab37ee14de9b082471ecb7ff902486b` |
| original tie rerun | `original-200-tie-rerun.json` | `a1beb3e22c0097c1f57361cb3a6cc0e222cd2cc5834f880b396428b2213e2e67` |
| new end-to-end attacks | `post-repair-adversarial-v6.json` | `97b2505cc8d861081089d64967cbd68be1a6cc2370e75d46fcbda0cece65c4bd` |
| complete-null rerun | `null-predefined/manifest.json` | `252bd59a9db4943084ae484d90c80af15b6aa5bae65e4549e0c2597fa81cb703` |
| entrapment rerun | `entrapment-predefined/manifest.json` | `1e74a1da53d11b35e3961a639831aa45b8a79c689e561ed12a902e7b4987a066` |
| pooled PEP calibration | `pep-rust-pooled.json` | `ff8d223feb55107a90b2bb1630f617a515888a5176f04ea4844e4507b4511c95` |
| multi-dataset rerun | `multidataset-predefined/manifest.json` | `c1134c78996d41dfcf2435d0dabb00ac3697067500a392aec0936538cda196aa` |
| mutation rerun | `mutation-rerun.json` | `29228ce4217ecce600182bf570b7b512d11b918d43dea82ff5a4c595bddbc53e` |
| 65-file PXD agreement | `pxd032157-agreement.json` | `fafb0e3db6ade7dd8ef6f2874acc956c722946dea898fb85f66f8db44257c33b` |

The three new probes are in this directory and accept either the frozen binary or compile directly
against the audited source. The full end-to-end rerun command was:

```bash
python3 validation/post_repair_adversarial.py \
  --binary "/run/media/andrej-rumenovski/New Volume/percolator_rs_postrepair_audit_20260826/frozen/percolator-rs" \
  --json "/run/media/andrej-rumenovski/New Volume/percolator_rs_postrepair_audit_20260826/post-repair-adversarial-v6.json"
```

## Executive answer

**No. The repairs fixed several of the exact defects previously demonstrated, but `percolator-rs` is
not yet scientifically defensible as an unqualified source of calibrated PSM, PEP, or protein-level
confidence.**

The strongest positive findings are:

1. The standalone TDC+ q-value scan matches an independent oracle, groups exact score ties, applies
   the finite-sample correction, obeys strict threshold boundaries, and is deterministic.
2. The original single-file 200-spectrum target/decoy tie attack is fixed: all nine row permutations
   give the same 100 target / 100 decoy winners and no q<0.01 discovery.
3. Fixed-C, nested `--select-c`, and `--ensemble` pass new whole-held-out-fold label, feature-outlier,
   and reorder attacks. The repaired `--select-c` is genuinely nested in the inspected code.
4. Ordinary protein evidence-set grouping and fair picked-protein competition pass new hand graphs.
   Picked protein PEP is no longer fabricated: it is consistently `NA`.
5. The complete-null experiment independently reproduces 0/30 runs with a rejection at every tested
   threshold. All 120 entrapment PSM artifacts and all 160 multi-dataset artifacts reproduce earlier
   hashes exactly.

The repair is nevertheless rejected because independent attacks found:

1. **Joined-input tie handling is still input-order dependent.** Swapping two `--join` arguments
   changes 294 winner identities and 135/151 target/decoy winners to 156/130. A cutoff fixture changes
   strict q<0.01 discoveries from **0 to 101**. The tie coin uses numeric `source` position rather than
   stable file identity.
2. **PEPs remain systematically optimistic.** Pooled five-seed weighted absolute and signed
   calibration error are both **+0.018685**; every populated bin is optimistic. There are 217 known-
   false PSMs at PEP<0.001, with a minimum of 0.000262.
3. **PEPs are permutation unstable.** Row-order-induced score changes of only 4.4e-16 to 1.33e-15
   split mathematically equivalent ties and change 2–3 PEPs by as much as **0.364865**, while q-values
   and winner identities remain unchanged.
4. **The public PEP conservation claim has an allowed-parameter counterexample.** At
   `p=1e-15`, five targets imply an estimated false count of approximately 1e-15, but the PEP floor
   makes the PEP sum approximately 5e-12.
5. **Protein output remains order dependent upstream of the repaired grouping.** Two equal-scoring
   rows for the same peptide but different protein mappings select the earlier mapping; reversing them
   changes the peptide representative from `PROT_A` to `PROT_B` and changes 242 protein groups to 243.
6. **Mixed target/decoy proteins are improperly collapsed.** A peptide mapping to `MIXED` and
   `DECOY_MIXED` produces one group containing both and classifies the whole group as decoy.
7. **Signal-present FDR control still fails.** At nominal q<0.01, current Rust has mean adjusted FDP
   **1.8104%** across the five predefined seeds. It is above nominal at every tested threshold.

These are reproducible implementation failures (items 1, 3–6) plus a calibration/estimator failure
(items 2 and 7). Entrapment-design limitations prevent uniquely attributing item 7 to a Rust coding
bug, but they do not turn the failed scientific calibration claim into a pass.

## 1. Exact repairs identified

The production diff from `b38c0db` through the last implementation commit `8e49e28` was inspected
directly. Later commits add tests, reports, expectations, and stored results.

| repair | code and mathematical/statistical change | tests added by repair | intended restored claim | independent result |
|---|---|---|---|---|
| Order-invariant coin, `d0fecb3` | Adds `src/tiebreak.rs`: SplitMix64 identity hash and uniform draw over a tie group | six primitive distribution/reproducibility tests | A tie can be randomized reproducibly without row-index dependence | Primitive works, but its PSM key later uses unstable numeric source order under `--join` |
| PEP reconstruction, `cb4b11f` | Replaces differences of monotonized q-values plus prior with per-score increments of `min(T, lambda(D+1))`, then PAVA; rejects nonfinite scores | closed forms, mass sums, monotonicity, opportunity-ratio tests; an independent repair probe | Reported PEP has cumulative-count semantics, no arbitrary prior, and no exact zero | Ordinary closed forms pass; posterior calibration still fails, permutation instability is severe, and extreme allowed `p` breaks exact conservation |
| PSM tie lottery, `651960e` | `competition_winners` accumulates all equal best candidates, canonically sorts content, draws using `(source, scan, ExpMass, seed)`; CLI restricts `p` to `(0,1)` | original 200-pair permutations, wide ties, seed sensitivity, strict winner | Scientifically equivalent row permutations cannot change PSM competition or q conclusions | Fixed for row permutations within one input; failed for permutations of joined input files (0 versus 101 q<.01) |
| Protein repair, `a2bed3e` | Groups proteins by identical observed peptide-index sets rather than connected components; fair target/decoy ties; `pep: Option<f64>`; picked output writes `NA` | evidence graphs, tied protein pairs, output PEP, CLI fixture | Distinguishable proteins stay separate, protein ties are fair, peptide PEP is not mislabeled | Core graph and tie logic pass; upstream peptide representative and mixed target/decoy grouping still invalidate some protein output |
| CV repair, `8e49e28` | Selects class weights separately inside each outer training fold with inner splits; removes `Label` from ensemble agreement keys | binary held-out-label tests plus repair attack script | `--select-c` is nested and ensemble has no label-derived global feature | Passed new label, outlier, composition, and reorder attacks in all three modes |
| Tie/mutation gates, `7db2ee3` | Adds explicit shared rejection-boundary tests and a 12-mutation runner; no new production method | tie rotations, q/PEP boundary tests, mutation harness | Repairs cannot silently regress | All 12 historical mutations caught, but five newly observed failures pass the unmutated suite |
| Protein shell gate, `fd5524d` | Removes `picked >= classic` yield assertion; requires `NA` protein PEP | shell regression | Gate no longer mistakes yield for validity | Correct change; still not a truth-based protein validation |

Documentation and result commits were `eb40377`, `ec31b74`, and `63d4783`. They do not add another
scientific repair.

## 2. New adversarial tests

The independent additions are:

- [`post_repair_stats_probe.rs`](post_repair_stats_probe.rs): separate q-value oracle, 13 score/label
  patterns, all 24 permutations of a four-row mixed tie, six exact threshold boundaries plus the
  just-below cases, closed-form PEPs, and an extreme-`p` conservation counterexample.
- [`post_repair_adversarial.py`](post_repair_adversarial.py): a new 311-spectrum nonconstant partial-
  tie fixture under seven permutations, joined-input ordering and cutoff attacks, new whole-fold CV
  attacks, and an equal-score peptide/protein mapping attack.
- [`post_repair_protein_probe.rs`](post_repair_protein_probe.rs): a two-peptide indistinguishable pair,
  a six-node distinguishable sharing cycle, subset/superset and independent groups, 317 new tied
  target/decoy pairs, and a mixed target/decoy evidence case.

The original 200-spectrum fixture was also rerun to answer the exact historical failure, but it was
not counted as a new attack.

## 3. Tie-handling results

### Original 200-spectrum fixture

One target and one decoy candidate per spectrum, identical features, candidate seed 1:

| input arrangement | target winners | decoy winners | targets q<0.01 | minimum target q |
|---|---:|---:|---:|---:|
| target first | 100 | 100 | 0 | 1.0 |
| decoy first | 100 | 100 | 0 | 1.0 |
| whole file reversed | 100 | 100 | 0 | 1.0 |
| all targets first | 100 | 100 | 0 | 1.0 |
| all decoys first | 100 | 100 | 0 | 1.0 |
| deterministic shuffles 1–4 | 100 each | 100 each | 0 each | 1.0 |

The winner set is identical in all nine permutations. Seed 2 changes the split to 114/86, as a
seeded lottery should. A 1,000-spectrum four-way tie (three target, one decoy per spectrum) gives
764/236 and no q<0.01 discovery, consistent with uniform candidate selection.

### New nonconstant partial/wide ties

The new fixture has 311 spectra, variable score levels, two- to five-way best-score ties, and one
strictly lower candidate per spectrum. Canonical, reversed, reverse-spectrum-block, and five
deterministic shuffle arrangements all produce:

- the same 187 target and 124 decoy winners;
- the same winner identities;
- the same reported q-values and discovery counts;
- zero targets below every requested q threshold.

PEPs are not invariant. Depending on the permutation, 2–3 winner PEPs move. A precision-only
diagnostic build (serialization change only) found that 231–266 scores moved by at most 1.33e-15.
Examples:

| PSM | canonical score | reversed score | canonical PEP | reversed PEP |
|---|---:|---:|---:|---:|
| `V219_1` | 0.262024109314653364 | 0.262024109314653308 | 0.611111 | 0.635135 |
| `V243_1` | -0.641229594081162202 | -0.641229594081162313 | 0.635135 | 1.000000 |
| `V272_1` | -0.641229594081162202 | -0.641229594081162313 | 0.635135 | 1.000000 |

Thus the competition and q repairs handle the new ties, but exact floating-point grouping makes local
PEPs discontinuous under scientifically irrelevant one-ULP changes.

### Joined input order: remaining PSM failure

Two constant-score PINs containing 137 and 149 tied target/decoy spectra were passed to `--join` in
both argument orders.

| order | target winners | decoy winners | changed winner identities | q<0.01 |
|---|---:|---:|---:|---:|
| alpha, beta | 135 | 151 | — | 0 |
| beta, alpha | 156 | 130 | 294 symmetric-difference entries | 0 |

The reason is direct in the implementation: the coin key contains `Dataset::source`, an integer
assigned from argument order. It does not contain stable source filename/content identity.

The cutoff version contains a leading target-rich block, one exact target/decoy tie, and low decoys:

| order | tied winner | target/decoy winners | targets q<0.01 | tied-winner q |
|---|---|---:|---:|---:|
| alpha, beta | decoy | 100 / 101 | **0** | 0.020000 |
| beta, alpha | target | 101 / 100 | **101** | 0.009901 |

This is the same scientific input and run seed. Only the order of joined files changed. The broad
README claim that input permutation cannot change a single winner is false.

### Peptide tie propagation

Two target PSM rows have the same scan, peptide, features, and reported score but map to different
proteins. Peptide q=0.008130 and PEP=0.008130 are invariant, but the representative is not:

| order | peptide representative | protein mapping | protein groups |
|---|---|---|---:|
| `AMB_A` first | `AMB_A` | `PROT_A` | 242 |
| `AMB_B` first | `AMB_B` | `PROT_B` | 243 |

The peptide-level confidence statistics survive; peptide identity/protein evidence and downstream
protein conclusions do not.

## 4. Q-value implementation

The independent oracle implements, separately from the crate:

`raw FDP = min(1, lambda * (D + 1) / max(1,T))`, once after each exact score group, followed by a
reverse cumulative minimum. The implementation matches it in every tested case.

| hand case | expected/observed result |
|---|---|
| empty | empty |
| one target | q=1 |
| one decoy | q=1 |
| one target/one decoy in either order | both q=1 |
| six targets, no decoy | every q=1/6 |
| six decoys, no target | every q=1 |
| alternating labels | exact oracle match, monotone by worsening score |
| mixed three-target/one-decoy tie | whole tie evaluated together; no within-tie boundary |
| all identical scores | one q-value for the complete group |
| target-heavy and decoy-heavy lists | exact oracle match, bounded |
| `f64::MAX`, `1e300`, `0`, `-0`, `-1e300`, `-f64::MAX` | finite, bounded, deterministic; `-0` and `0` share a group |
| nonfinite score | estimator fails closed |

All 24 permutations of a four-row mixed tie restore exactly the same q-values by row identity.

Strict thresholds were tested with one observed decoy placed above exactly enough targets to make
`(D+1)/T` equal the threshold. The exact boundary is rejected; adding one target crosses it:

| threshold | targets at exact boundary | accepted at exact q | accepted after one extra target |
|---:|---:|---:|---:|
| 0.001 | 2,000 | 0 | 2,001 |
| 0.005 | 400 | 0 | 401 |
| 0.01 | 200 | 0 | 201 |
| 0.02 | 100 | 0 | 101 |
| 0.05 | 40 | 0 | 41 |
| 0.10 | 20 | 0 | 21 |

Conclusion: the standalone q-value implementation has strong evidence under its direct-TDC input
contract. It cannot repair an invalid/order-dependent upstream winner list, as the joined-input
attack demonstrates.

## 5. Complete-null rerun

The predefined design was not changed: ten exact-balance relabelings of original decoy rows in each
of three PXD032157 PINs, model seed 1, strict thresholds.

Every cell below is the false-target count at q<0.001 / 0.005 / 0.01 / 0.02 / 0.05 / 0.10:

| input | relabel seed | six false-target counts |
|---|---:|---|
| PXD032157-01 | 1001–1010 (each) | 0 / 0 / 0 / 0 / 0 / 0 |
| PXD032157-02 | 1001–1010 (each) | 0 / 0 / 0 / 0 / 0 / 0 |
| PXD032157-03 | 1001–1010 (each) | 0 / 0 / 0 / 0 / 0 / 0 |

Thus every one of the 30 individual runs has zero discovery at all six thresholds: **0/30 failing
runs**, unchanged from the immediate pre-second-repair result.

This is positive but low-resolution evidence. With zero events in 30 Bernoulli replicates, the exact
one-sided 95% upper bound is about 9.5% (the stored two-sided Wilson upper bound is about 11.4%). It
cannot demonstrate calibration at 0.1%, 0.5%, or 1%.

## 6. Signal-present entrapment

### Every current Rust run

| seed | q | accepted PSMs | pure entrapment PSMs | adjusted FDP |
|---:|---:|---:|---:|---:|
| 1 | .001 | 12,454 | 54 | 0.5781% |
| 1 | .005 | 17,863 | 169 | 1.1632% |
| 1 | .01 | 19,545 | 258 | 1.7865% |
| 1 | .02 | 21,064 | 435 | 2.7671% |
| 1 | .05 | 23,813 | 1,045 | 5.8445% |
| 1 | .10 | 26,796 | 2,219 | 10.8628% |
| 2 | .001 | 13,019 | 63 | 0.6452% |
| 2 | .005 | 17,813 | 159 | 1.1588% |
| 2 | .01 | 19,501 | 244 | 1.7288% |
| 2 | .02 | 21,277 | 447 | 2.8358% |
| 2 | .05 | 23,794 | 1,048 | 5.8246% |
| 2 | .10 | 26,914 | 2,218 | 10.9226% |
| 3 | .001 | 13,101 | 64 | 0.6281% |
| 3 | .005 | 18,063 | 183 | 1.2664% |
| 3 | .01 | 19,536 | 262 | 1.8150% |
| 3 | .02 | 21,141 | 423 | 2.7212% |
| 3 | .05 | 23,775 | 1,032 | 5.8323% |
| 3 | .10 | 26,853 | 2,223 | 11.0048% |
| 4 | .001 | 13,174 | 61 | 0.5953% |
| 4 | .005 | 17,950 | 157 | 1.1461% |
| 4 | .01 | 19,514 | 262 | 1.8838% |
| 4 | .02 | 21,108 | 420 | 2.6444% |
| 4 | .05 | 23,854 | 1,028 | 5.6865% |
| 4 | .10 | 26,833 | 2,194 | 10.6800% |
| 5 | .001 | 12,984 | 61 | 0.6264% |
| 5 | .005 | 18,240 | 186 | 1.2917% |
| 5 | .01 | 19,586 | 267 | 1.8378% |
| 5 | .02 | 21,120 | 427 | 2.7588% |
| 5 | .05 | 23,906 | 1,080 | 6.0665% |
| 5 | .10 | 26,846 | 2,222 | 10.8912% |

### Every raw C++ run

The raw reference does not apply the same explicit spectrum competition on these concatenated
multi-candidate PINs, so this is an agreement/context arm, not an oracle.

| seed | q | accepted PSMs | pure entrapment PSMs | adjusted FDP |
|---:|---:|---:|---:|---:|
| 1 | .001/.005/.01/.02/.05/.10 | 13,076 / 17,600 / 19,122 / 20,619 / 23,255 / 26,187 | 143 / 277 / 368 / 537 / 1,132 / 2,236 | 1.4061 / 2.0007 / 2.6149 / 3.5288 / 6.6590 / 11.5360% |
| 2 | .001/.005/.01/.02/.05/.10 | 13,053 / 17,252 / 19,056 / 20,511 / 23,131 / 26,022 | 137 / 264 / 369 / 537 / 1,136 / 2,251 | 1.1808 / 1.8873 / 2.5625 / 3.6414 / 6.7804 / 11.7404% |
| 3 | .001/.005/.01/.02/.05/.10 | 12,852 / 17,506 / 19,218 / 20,744 / 23,311 / 26,304 | 130 / 257 / 365 / 546 / 1,143 / 2,275 | 1.3005 / 1.8050 / 2.4718 / 3.6465 / 6.6615 / 11.7468% |
| 4 | .001/.005/.01/.02/.05/.10 | 12,670 / 17,814 / 19,262 / 20,757 / 23,302 / 26,256 | 132 / 278 / 369 / 544 / 1,136 / 2,261 | 1.3395 / 1.8044 / 2.5259 / 3.6527 / 6.6253 / 11.6218% |
| 5 | .001/.005/.01/.02/.05/.10 | 12,808 / 17,822 / 19,232 / 20,732 / 23,366 / 26,282 | 130 / 274 / 362 / 542 / 1,148 / 2,268 | 1.1419 / 1.9804 / 2.4319 / 3.6094 / 6.7482 / 11.7143% |

### Five-seed summaries

Each cell is mean / median / sample SD / min–max.

| method | q | accepted | entrapment | adjusted FDP |
|---|---:|---:|---:|---:|
| Rust | .001 | 12,946.4 / 13,019 / 285.0 / 12,454–13,174 | 60.6 / 61 / 3.9 / 54–64 | 0.6146 / 0.6264 / 0.0272 / 0.5781–0.6452% |
| Rust | .005 | 17,985.8 / 17,950 / 170.9 / 17,813–18,240 | 170.8 / 169 / 13.3 / 157–186 | 1.2052 / 1.1632 / 0.0682 / 1.1461–1.2917% |
| Rust | .01 | 19,536.4 / 19,536 / 32.7 / 19,501–19,586 | 258.6 / 262 / 8.8 / 244–267 | **1.8104 / 1.8150 / 0.0579 / 1.7288–1.8838%** |
| Rust | .02 | 21,142.0 / 21,120 / 80.5 / 21,064–21,277 | 430.4 / 427 / 10.9 / 420–447 | 2.7455 / 2.7588 / 0.0700 / 2.6444–2.8358% |
| Rust | .05 | 23,828.4 / 23,813 / 52.3 / 23,775–23,906 | 1,046.6 / 1,045 / 20.5 / 1,028–1,080 | 5.8509 / 5.8323 / 0.1365 / 5.6865–6.0665% |
| Rust | .10 | 26,848.4 / 26,846 / 42.8 / 26,796–26,914 | 2,215.2 / 2,219 / 12.0 / 2,194–2,223 | 10.8723 / 10.8912 / 0.1199 / 10.6800–11.0048% |
| C++ raw | .001 | 12,891.8 / 12,852 / 171.5 / 12,670–13,076 | 134.4 / 132 / 5.6 / 130–143 | 1.2737 / 1.3005 / 0.1102 / 1.1419–1.4061% |
| C++ raw | .005 | 17,598.8 / 17,600 / 237.2 / 17,252–17,822 | 270.0 / 274 / 9.1 / 257–278 | 1.8956 / 1.8873 / 0.0933 / 1.8044–2.0007% |
| C++ raw | .01 | 19,178.0 / 19,218 / 86.0 / 19,056–19,262 | 366.6 / 368 / 3.0 / 362–369 | 2.5214 / 2.5259 / 0.0723 / 2.4319–2.6149% |
| C++ raw | .02 | 20,672.6 / 20,732 / 105.8 / 20,511–20,757 | 541.2 / 542 / 4.1 / 537–546 | 3.6158 / 3.6414 / 0.0514 / 3.5288–3.6527% |
| C++ raw | .05 | 23,273.0 / 23,302 / 88.6 / 23,131–23,366 | 1,139.0 / 1,136 / 6.4 / 1,132–1,148 | 6.6949 / 6.6615 / 0.0659 / 6.6253–6.7804% |
| C++ raw | .10 | 26,210.2 / 26,256 / 114.0 / 26,022–26,304 | 2,258.2 / 2,261 / 15.3 / 2,236–2,275 | 11.6719 / 11.7143 / 0.0910 / 11.5360–11.7468% |

### Before/current/reference comparison

Mean adjusted FDP:

| arm | q<.001 | q<.005 | q<.01 | q<.02 | q<.05 | q<.10 |
|---|---:|---:|---:|---:|---:|---:|
| original rejected Rust `d83a7ba` | 1.174% | 1.991% | 2.698% | 3.909% | 7.209% | 12.604% |
| immediate pre-second Rust `b38c0db` | 0.615% | 1.207% | **1.816%** | 2.753% | 5.848% | 10.873% |
| current Rust | 0.615% | 1.205% | **1.810%** | 2.745% | 5.851% | 10.872% |
| C++ raw list | 1.274% | 1.896% | 2.521% | 3.616% | 6.695% | 11.672% |
| C++ with matched external competition | 0.644% | 1.277% | 1.729% | 2.717% | 5.986% | 10.843% |

The second repair did not materially improve the previous approximately 1.82% result. The curve is
anti-conservative at every threshold. Similar matched-reference behavior argues against a uniquely
Rust-specific q-scan bug on these six files; it does not establish nominal control.

## 7. PEP implementation and calibration

### Reconstructed estimator and hand behavior

For ordinary `p=0.5` inputs, the repaired estimator behaves as documented:

- five targets above one decoy: q=PEP=0.2 for every target, target PEP sum=1;
- one target/one decoy: target PEP=1;
- decoy-only: placeholder PEP=1 on every decoy;
- two estimated false discoveries over five targets: target PEP sum=2;
- all target PEPs are finite, in [0,1], monotone after PAVA, and nonzero;
- exact score ties receive one group increment in the direct estimator probe.

These are algebraic properties of the implemented heuristic. They do not establish that the values
are posterior probabilities.

Across the 30 complete-null target outputs, all **401,711** reported target winners are known false.
Their mean PEP is 0.991691 (observed false fraction 1.0), the minimum is 0.2, 187,833 have PEP=1,
and none has PEP<0.01, <0.001, <0.0001, or exactly zero. Thus this null has no misleading near-zero
PEP tail, although the aggregate observed-minus-predicted gap is still +0.008309. This is consistent
with the zero q-threshold discoveries and does not rescue the signal-present calibration below.

At the accepted boundary `p=1e-15`, five all-target scores produce:

- `5*q = 1.000000000000001e-15` estimated false discoveries;
- PEP sum `4.9999999999999997e-12` after the `1e-12` per-target floor.

The README's unrestricted “sum exactly” and “opportunity ratio scales exactly” claims therefore have
a public parameter counterexample.

### Entrapment calibration

Pooled over all five predefined seeds (493,760 target outputs):

- exact PEP=0: **0**;
- PEP<1e-4: **0**;
- PEP<1e-3: **54,378**;
- exact PEP=1: **4,909**;
- minimum PEP: **0.000262**;
- pure known-false targets: **278,413**;
- pure known-false with PEP<0.001: **217**;
- pure known-false with PEP<0.01: **414**;
- weighted absolute calibration error: **0.01868545**;
- weighted signed observed-minus-reported error: **+0.01868545**.

Every populated bin is optimistic:

| PEP bin | PSMs | mean predicted | adjusted observed false fraction | observed − predicted |
|---|---:|---:|---:|---:|
| [1e-4,1e-3) | 54,378 | 0.000386 | 0.004930 | +0.004543 |
| [1e-3,.005) | 15,887 | 0.002707 | 0.009335 | +0.006629 |
| [.005,.01) | 5,750 | 0.006783 | 0.017739 | +0.010956 |
| [.01,.02) | 7,696 | 0.015852 | 0.031274 | +0.015421 |
| [.02,.05) | 9,091 | 0.033330 | 0.056601 | +0.023271 |
| [.05,.10) | 5,374 | 0.073874 | 0.090656 | +0.016782 |
| [.10,.20) | 9,231 | 0.145705 | 0.166029 | +0.020324 |
| [.20,.50) | 19,260 | 0.345587 | 0.372132 | +0.026546 |
| [.50,1.00] | 367,093 | 0.913066 | 0.934018 | +0.020952 |

Per-seed absolute/signed errors are 0.018825, 0.020400, 0.016540, 0.019556, and 0.018467; the sign is
positive for every seed. The previous systematic optimism did **not** disappear.

Scientific classification:

- replacing differences of monotonized q-values and removing an arbitrary prior fixed a real
  implementation defect;
- interpreting increments of a cumulative target/decoy estimate as individual posterior
  probabilities is an estimator limitation without an independently calibrated local model;
- the observed optimism is a failed calibration result, not repaired by boundedness, conservation,
  or isotonic monotonicity;
- the one-ULP permutation amplification and extreme-`p` floor are additional implementation issues.

## 8. Cross-validation leakage

The new fixture contains 360 spectra, with 120 outer-fold-0 spectra. The attack changes only outer
fold 0 while checking models/scores for that fold:

1. flip every held-out label (target/decoy composition attack);
2. replace most held-out features by +/-1e12 while leaving 12 held-out spectra untouched as
   sentinels;
3. reorder held-out examples while preserving content;
4. repeat under fixed-C, `--select-c`, and a two-engine ensemble.

| mode | held-out rows checked | label-changed scores | untouched sentinel scores changed by outliers | scores changed by reorder | fold-0 selection changed |
|---|---:|---:|---:|---:|---:|
| fixed-C | 240 | 0 | 0/24 | 0 at reported precision | not applicable |
| `--select-c` | 240 | 0 | 0/24 | 0 at reported precision | no |
| `--ensemble` | 240 | 0 | 0/24 | 0 at reported precision | not applicable |

For `--select-c`, fold 0 independently retains
`C=1, positive:negative=0.2:1.0, 3 features, tolerance=1e-5, inner yield=0` under every attack.

Source inspection confirms true nesting:

`outer training rows -> two inner folds -> candidate training/validation -> per-outer-fold class
weights -> train on complete outer training partition -> score untouched outer fold`.

The ensemble has no member-selection or weighting layer: engine rows and one-hot feature blocks are
merged into one fold-local learner. Agreement keys are now `(ScanNr, peptide)`, not label keyed, and
all rows sharing an ensemble scan are assigned to the same outer fold. No held-out information was
found to influence training, combination, or selection.

A precision serialization diagnostic observed sub-1e-6 row-order numerical differences in some
scores, but label/outlier attacks did not move the model at any reported precision and no selection
changed. This is numerical reproducibility noise, not evidence of held-out leakage.

## 9. Protein inference

### Group definition

New hand graphs give the expected result:

- two proteins supported by the same two peptides -> one indistinguishable group;
- six proteins in a sharing cycle with alternating unique evidence -> six groups, not one connected
  component;
- strict subset/superset evidence -> separate groups;
- a shared peptide supports every mapped group; no hidden razor/parsimony assignment;
- multiple independent groups coexist;
- entry permutation leaves these core groups invariant.

The evidence-set repair itself is correct for its supplied `entries`. The end-to-end pipeline does
not supply an order-invariant protein mapping when equal-score peptide representatives disagree, as
shown in section 3.

### Protein competition

On 317 new exact target/decoy protein pairs at seed 31:

- exactly one group is picked per pair;
- 166/317 picks are targets, compatible with a fair split;
- reversing every entry preserves the full pick/q signature;
- seed 32 changes the signature;
- a strict score winner is not displaced by the tie rule.

This is good conditional evidence for the repaired picking step. It is not protein-FDR calibration.

### Remaining grouping defects

1. Equal-score peptide mapping: reversing `AMB_A`/`AMB_B` makes `PROT_B` appear as an additional
   group (242 -> 243) even though confidence statistics are unchanged.
2. Mixed target/decoy evidence: one peptide mapped to `MIXED DECOY_MIXED` is grouped into
   `[DECOY_MIXED,MIXED]`. Because `is_decoy` is `any(member is decoy)`, the target is emitted as part
   of a decoy group and disappears from target output.
3. No parsimony is implemented. This is documented, but it limits biological interpretation.
4. No protein-level truth study validates picked q-values after this repair.

### Protein PEP

Every picked group in both CLI orders and every direct graph has `pep=None`; output writes `NA`.
README and current repair documentation explicitly say picked FDR estimates a cumulative rate and no
protein posterior. The previous peptide-as-protein PEP defect is fixed.

The separate Bayesian path still reports a model-derived protein posterior. It was not revalidated
here, consumes systematically optimistic peptide PEPs, and its prior PrEST study did not establish
nominal protein FDR. No scientific validity is inferred for that posterior.

## 10. Mutation testing

An isolated worktree at `HEAD` first passed all tests. Each historical defect was then introduced
alone, the release suite was run, and the original source restored before the next mutation. The
worktree was removed afterward.

| mutation | result | detecting evidence |
|---|---|---|
| input-order-dependent PSM ties | caught | 4 competition tests |
| missing finite-sample decoy | caught | 15 q/PEP/protein tests |
| missing tie grouping | caught | 3 q/PEP/rank tests |
| arbitrary PEP prior restored | caught | 6 PEP identity/closed-form tests |
| PEP floor raised to .01 | caught | PEP sum test |
| globally fitted normalization | caught | 2 fold-local tests |
| globally chosen initial direction | caught | held-out model test |
| leaked `--select-c` | caught | binary selected-C leakage test |
| label-keyed ensemble feature | caught | compile/type failure; exhaustive relabel unit test protects compilable equivalent |
| target-favoring protein ties | caught | 2 protein tie tests |
| connected-component protein grouping | caught | 5 grouping tests |
| peptide PEP as protein PEP | caught | 2 protein PEP tests |

Result: **12/12 historical mutations caught**.

The current suite nevertheless passes while joined-file order, equal-score peptide mapping, mixed
target/decoy groups, one-ULP PEP instability, and extreme-`p` PEP conservation fail this audit.
Those five properties have no production regression gate. Mutation strength is therefore good for
known repairs, incomplete for the actual advertised scientific surface.

## 11. Multi-seed results

### Exact reproducibility

- 120 current entrapment PSM TSVs versus the repair run: **120/120 identical SHA-256**.
- 160 current compact-dataset Rust/C++ PSM/peptide TSVs versus the repair run: **160/160 identical
  SHA-256**.
- No seed or unfavorable dataset was excluded.

All five entrapment seeds are reported individually in section 6, with mean, median, sample standard
deviation, and range at every threshold. All five seeds for each compact dataset are reported in the
next section. Seed variability is small for Tide/Sage rankings, larger for MSFragger/yeast accepted
sets, and does not change the failed entrapment-calibration conclusion: every seed is optimistic in
aggregate PEP calibration and every seed exceeds 1% adjusted FDP at q<0.01.

## 12. Multi-dataset results

### Per-seed q<0.01 results

Columns are Rust/C++ PSMs, Rust/C++ peptides, Rust-only/C++-only accepted PSMs, and Jaccard.

| dataset | seed | PSMs R/C | peptides R/C | exclusive R/C | Jaccard |
|---|---:|---:|---:|---:|---:|
| Tide | 1 | 27,640 / 27,617 | 19,736 / 19,722 | 106 / 83 | .9932 |
| Tide | 2 | 27,665 / 27,656 | 19,713 / 19,736 | 100 / 91 | .9931 |
| Tide | 3 | 27,616 / 27,610 | 19,715 / 19,745 | 115 / 109 | .9919 |
| Tide | 4 | 27,624 / 27,612 | 19,701 / 19,737 | 97 / 85 | .9934 |
| Tide | 5 | 27,624 / 27,670 | 19,735 / 19,722 | 71 / 117 | .9932 |
| MSFragger | 1 | 1,367 / 1,388 | 1,060 / 1,062 | 80 / 101 | .8767 |
| MSFragger | 2 | 1,377 / 1,477 | 1,067 / 1,128 | 45 / 145 | .8752 |
| MSFragger | 3 | 1,312 / 1,399 | 1,030 / 1,075 | 37 / 124 | .8879 |
| MSFragger | 4 | 1,361 / 1,409 | 1,065 / 1,085 | 57 / 105 | .8895 |
| MSFragger | 5 | 1,340 / 1,383 | 1,056 / 1,053 | 63 / 106 | .8831 |
| Sage | 1 | 25,806 / 25,795 | 11,245 / 11,336 | 59 / 48 | .9959 |
| Sage | 2 | 25,803 / 25,788 | 11,253 / 11,317 | 66 / 51 | .9955 |
| Sage | 3 | 25,804 / 25,784 | 11,247 / 11,322 | 59 / 39 | .9962 |
| Sage | 4 | 25,791 / 25,790 | 11,246 / 11,327 | 54 / 53 | .9959 |
| Sage | 5 | 25,779 / 25,790 | 11,253 / 11,317 | 45 / 56 | .9961 |
| yeast | 1 | 1,150 / 1,147 | 904 / 928 | 41 / 38 | .9335 |
| yeast | 2 | 1,144 / 1,137 | 908 / 890 | 45 / 38 | .9298 |
| yeast | 3 | 1,124 / 1,059 | 907 / 863 | 91 / 26 | .8983 |
| yeast | 4 | 1,180 / 1,105 | 954 / 896 | 98 / 23 | .8994 |
| yeast | 5 | 1,159 / 1,111 | 935 / 908 | 70 / 22 | .9221 |

Five-seed q<0.01 summaries (mean / median / SD / range):

| dataset | method | PSMs | peptides |
|---|---|---:|---:|
| Tide | Rust | 27,633.8 / 27,624 / 19.5 / 27,616–27,665 | 19,720.0 / 19,715 / 15.1 / 19,701–19,736 |
| Tide | C++ | 27,633.0 / 27,617 / 27.9 / 27,610–27,670 | 19,732.4 / 19,736 / 10.1 / 19,722–19,745 |
| MSFragger | Rust | 1,351.4 / 1,361 / 25.9 / 1,312–1,377 | 1,055.6 / 1,060 / 14.9 / 1,030–1,067 |
| MSFragger | C++ | 1,411.2 / 1,399 / 38.1 / 1,383–1,477 | 1,080.6 / 1,075 / 29.2 / 1,053–1,128 |
| Sage | Rust | 25,796.6 / 25,803 / 11.5 / 25,779–25,806 | 11,248.8 / 11,247 / 3.9 / 11,245–11,253 |
| Sage | C++ | 25,789.4 / 25,790 / 4.0 / 25,784–25,795 | 11,323.8 / 11,322 / 8.0 / 11,317–11,336 |
| yeast | Rust | 1,151.4 / 1,150 / 20.5 / 1,124–1,180 | 921.6 / 908 / 22.0 / 904–954 |
| yeast | C++ | 1,111.8 / 1,111 / 34.3 / 1,059–1,147 | 897.0 / 896 / 23.9 / 863–928 |

### Threshold behavior and ranking agreement

Mean target PSM counts and mean Jaccard over the five seeds:

| dataset | metric | .001 | .005 | .01 | .02 | .05 | .10 |
|---|---|---:|---:|---:|---:|---:|---:|
| Tide | Rust / C++ | 24,749.6/24,742.8 | 26,909.8/26,920.0 | 27,633.8/27,633.0 | 28,534.8/28,526.8 | 30,102.6/30,099.0 | 32,194.8/32,209.0 |
| Tide | Jaccard | .9879 | .9922 | .9930 | .9917 | .9900 | .9874 |
| Sage | Rust / C++ | 24,782.0/24,822.4 | 25,478.4/25,483.8 | 25,796.6/25,789.4 | 26,112.6/26,115.4 | 26,638.0/26,644.8 | 27,390.8/27,411.8 |
| Sage | Jaccard | .9938 | .9955 | .9959 | .9965 | .9947 | .9948 |
| MSFragger | Rust / C++ | 641.6/1,100.0 | 1,218.8/1,345.0 | 1,351.4/1,411.2 | 1,473.8/1,545.2 | 1,605.0/1,685.0 | 1,760.2/1,861.6 |
| MSFragger | Jaccard | .5072 | .8583 | .8825 | .8969 | .8745 | .8502 |
| yeast | Rust / C++ | 0/0 | 1,021.8/1,017.0 | 1,151.4/1,111.8 | 1,253.2/1,178.2 | 1,473.0/1,375.8 | 1,706.0/1,551.2 |
| yeast | Jaccard | undefined | .9072 | .9166 | .9080 | .8952 | .8757 |

Mean score/q/PEP Spearman on matched PSMs is:

| dataset | score | q | PEP | q<.01 Rust-only / C++-only over five seeds |
|---|---:|---:|---:|---:|
| Tide | .9988 | .9966 | .9825 | 489 / 485 |
| Sage | .9955 | .9728 | .9566 | 283 / 247 |
| MSFragger | .9602 | .9548 | .9458 | 282 / 581 |
| yeast | .9424 | .9423 | .6991 | 345 / 147 |

Tide and Sage are strong agreement cases. MSFragger and yeast are unfavorable and retained; they use
multiple candidates per precursor, while raw C++ reports all candidates and Rust competes them, and
they also retain genuine solver/model differences. These studies have target/decoy labels, not truth,
so they do not measure accuracy or calibration. PEP calibration is unavailable on these compact
datasets because they have no known-false target stratum.

Protein counts are not reported as valid cross-dataset discoveries: current picked protein q-values
have no truth validation, and this audit found remaining grouping/input-order defects.

### PXD032157 and the forced-concatenated pathology

Independent 65-file seed-1 Rust rerun:

- 65/65 valid files;
- 106,823 target PSMs and 35,886 target peptides at q<0.01;
- four-worker wall 17.70 s; peak RSS 768,928 KiB.

Direct comparison with the stored raw C++ reference:

- Rust reports 1,729,858 competed rows; C++ reports 8,639,746 uncompeted candidate rows;
- 1,729,858 unambiguous matching PSMs;
- score/q/PEP Spearman .7785/.9240/.7368;
- q<.01 Rust/C++ 106,823/103,038;
- intersection 99,926, Rust-only 6,897, C++-only 3,112, Jaccard .9090.

This is not like-for-like post-processing and is not an accuracy comparison.

The previously pathological low-yield file was rerun and not excluded:

| case | Rust/C++ q<.01 | intersection | Rust-only/C++-only | score Spearman | Jaccard |
|---|---:|---:|---:|---:|---:|
| low-yield forced-concatenated file | 152 / 0 | 0 | 152 / 0 | .4125 | 0 |
| high-yield file from same dataset | 9,938 / 9,764 | 9,662 | 276 / 102 | .9914 | .9624 |

The low-yield divergence remains reproducible. C++ keeps roughly five candidates per precursor while
Rust keeps one; it is evidence of incompatible post-processing/auto-detection on this dataset, not
proof that either result is correct. PXD032157 was development data and is not an untouched
generalization set.

## 13. Validation-design limitations

### Complete null

What is valid:

- only original decoy-search matches are retained, so every pseudo-target is false;
- exact target/decoy balance and random relabeling make labels exchangeable conditional on stored
  rows;
- under the complete null, any discovery gives FDP=1, so the event rate is empirical FDR;
- thresholds and seeds were fixed before this rerun.

Limitations:

- 30 dependent relabelings have almost no resolution for nominal rates <=1%;
- rows inherit spectrum/candidate dependence from three source files;
- relabeling stored decoy matches is not a fresh null search and may omit search-stage interactions;
- the original design did not target joined-file identity or structured partial score ties.

Classification: useful null pathology test, not a proof of small-q control.

### Entrapment

What is valid:

- pure plant-proteome assignments in native Anopheles samples are credible known-false targets;
- mixed native/foreign mappings are excluded from the known-false numerator;
- the identical, predefined accounting is applied to every seed and method;
- the observed anti-conservative curve is a legitimate failure of the empirical calibration claim
  under this design.

Limitations and assumptions:

- adjusted false count is `pure entrapment targets / foreign fraction among accepted nonmixed
  decoys`; it assumes incorrect native targets and decoys have the same foreign-placement
  probability after scoring and competition;
- the plug-in fraction is itself threshold- and seed-dependent and unstable when accepted decoys are
  sparse;
- foreign and native peptides differ in homology, detectability, database multiplicity, and score
  distribution;
- pooling six runs creates PSM/protein/spectrum dependence; five model seeds are not five biological
  replicates or a confidence interval;
- semi-supervised training explicitly uses target/decoy labels, so exchangeability after adaptive
  fitting needs assumptions beyond correct q arithmetic;
- C++ raw and externally competed lists answer different post-processing questions.

Classification of the residual 1.81%:

- **not demonstrated to be a q-scan implementation bug** on these six ordinary single-file runs;
- **a failed estimator/method calibration result** under the predefined search design;
- potentially influenced by violated TDC assumptions, dataset/search-space effects, and the plug-in
  entrapment adjustment;
- not exonerated by similar C++ behavior.

The repair report's no-rescoring and iteration-budget patterns are useful causal clues. They do not
uniquely prove that semi-supervised training is the complete cause, because search-space and
adjustment limitations remain.

### Agreement/generalization studies

These establish parser compatibility, deterministic reproducibility, ranking similarity, and
accepted-set overlap. They do not establish truth, accuracy, sensitivity at equal true FDR,
probability calibration, protein validity, DIA behavior, or biological generalization. All datasets
were available during development.

### Performance/correctness check

Source inspection and tests support these points:

- q sorting builds an in-bounds index permutation; `get_unchecked` is confined to those indices;
- `total_cmp` supplies deterministic finite ordering and exact numeric equality defines tie groups;
- target-count/mask fast paths evaluate only after complete tie groups and match materialized q-values;
- the SVM active set is the exact set with positive squared-hinge residual, not an approximation;
- gradient and Hessian sum every active row; feature masks exclude weights exactly;
- the 22-element dot product preserves sequential addition order; SIMD `axpy` is elementwise;
- fold preprocessing and score standardization are trained only on training rows;
- debug/release and thread-count regression tests remain deterministic for identical input bytes;
- the 65-file performance gate reproduces expected counts and stays inside runtime/RSS budgets.

No performance optimization was found to deliberately omit tie grouping, approximate the active
set, or bypass the repaired q scan.

Correctness is not equivalent to C++ behavior: the solvers, folds, score anchoring, candidate lists,
and PEP methods differ. The row-order precision experiment also shows that ordinary floating-point
accumulation can perturb scores by one ULP and the local PEP estimator can amplify that into a 0.365
change. Therefore performance-path equivalence receives only weak evidence, not a proof.

## 14. Remaining defects

1. `--join` tie coins use numeric input position; reorder can change discoveries from 0 to 101.
2. Equal-score peptide representatives use first-row wins, changing protein evidence/output.
3. Mixed target/decoy proteins can be grouped and classified wholly decoy.
4. PEP exact-tie behavior is discontinuous under one-ULP row-order perturbations.
5. Extreme valid `--null-target-win-prob` breaks claimed PEP mass conservation due to the floor.
6. PEPs fail calibration in every populated entrapment bin.
7. Signal-present q/FDR calibration fails at every requested threshold.
8. Picked protein q-values and Bayesian protein posteriors lack successful truth-based validation.
9. No untouched dataset establishes cross-dataset biological/generalization claims.
10. The validation suite has no gates for defects 1–5.

## 15. Claim audit

Classifications use only `SUPPORTED`, `SUPPORTED WITH CAVEATS`, `INSUFFICIENT EVIDENCE`,
`MISLEADING`, and `INCORRECT` as requested.

| current major claim | classification | evidence |
|---|---|---|
| “The independent audit findings have been repaired” | **MISLEADING** | Several exact repairs work, but joined tie order, PEP calibration/stability, and protein pipeline defects remain |
| Single-file 200-spectrum row-order tie failure is fixed | **SUPPORTED** | Nine permutations have identical 100/100 winners and zero q<.01 |
| Exact tie coin is independent of input permutation and cannot change a single winner | **INCORRECT** | Swapping `--join` inputs changes 294 identities and q<.01 0 -> 101 |
| TDC+ q formula, finite correction, tie grouping, monotonicization | **SUPPORTED WITH CAVEATS** | Independent oracle/boundaries pass; valid only with a scientifically valid upstream winner list and direct-TDC contract |
| q-values/FDR are calibrated | **INCORRECT** | README correctly denies this elsewhere; entrapment mean FDP 1.810% at nominal 1% and above nominal at every threshold |
| Complete-null pathology is repaired | **SUPPORTED WITH CAVEATS** | 0/30 independently reproduced; too underpowered to demonstrate small nominal rates |
| Fixed-C has no held-out leakage | **SUPPORTED** | label, feature-outlier, reorder, source inspection, and mutation attacks pass |
| `--select-c` is nested/leakage-free | **SUPPORTED** | per-outer-fold inner selection; selection and held-out scores invariant under attacks |
| Ensemble is leakage-free | **SUPPORTED WITH CAVEATS** | label-free agreement and grouped scans pass attacks; no truth-based accuracy study |
| PEPs are increments of the competition false-count estimate | **SUPPORTED WITH CAVEATS** | ordinary hand cases pass; interpretation as posterior probability is not established |
| PEPs sum exactly to the reported estimate for the public parameter range | **INCORRECT** | `p=1e-15` floor counterexample: about 5e-12 versus 1e-15 |
| PEP=0 is unreachable | **SUPPORTED** | no exact zeros in hand cases or 493,760 entrapment targets; explicit floor also ensures it |
| PEPs are calibrated posterior probabilities | **INCORRECT** | every pooled bin optimistic; signed error +.018685; 217 known-false PEP<.001 |
| Protein groups are proteins with identical observed peptide evidence | **SUPPORTED WITH CAVEATS** | core `infer` hand graphs pass; upstream equal-score mapping and mixed target/decoy grouping invalidate some end-to-end groups |
| Picked protein competition treats exact ties fairly | **SUPPORTED WITH CAVEATS** | 166/317 target wins, reversal invariant; conditional on correct group construction/pairing |
| Picked protein PEP is unavailable and reported `NA` | **SUPPORTED** | direct objects, CLI output, shell mutation gate, README all agree |
| Bayesian protein PEP is a validated protein posterior | **INSUFFICIENT EVIDENCE** | optimistic peptide inputs and no successful post-repair truth calibration |
| Protein output has not been validated against protein ground truth | **SUPPORTED** | candid and accurate limitation |
| Four-dataset Rust/C++ agreement is broad | **SUPPORTED WITH CAVEATS** | Tide/Sage strong; MSFragger/yeast weaker; candidate-list mismatch on multi-candidate PINs |
| Rust and C++ are statistically/algorithmically equivalent | **INSUFFICIENT EVIDENCE** | not claimed as headline; PXD score Spearman .779 and material exclusive sets/method differences |
| Identification yield establishes accuracy or sensitivity | **INCORRECT** | README appropriately says it does not; no equal-true-FDR truth comparison |
| Cross-dataset generalization | **INSUFFICIENT EVIDENCE** | four known compact datasets, one development dataset, one pooled entrapment design; no untouched set |
| 106,823 / 35,886 is the current PXD development baseline | **SUPPORTED** | independently reproduced on 65/65 files |
| Named large-case performance claim | **SUPPORTED WITH CAVEATS** | independent Rust run 17.70 s/769 MB; historical 15.3 s is a named three-run median, not correctness or portability evidence |
| Results are bit-deterministic for identical input bytes and seed | **SUPPORTED WITH CAVEATS** | exact rerun hashes match; scientifically equivalent row/file permutations need not match |
| Residual entrapment excess is fully attributed to semi-supervised training | **MISLEADING** | iteration/no-rescoring evidence is suggestive, but TDC and entrapment-design assumptions prevent unique attribution |
| “Validated configuration” for canonical profile | **MISLEADING** | implementation/regression validation exists, but signal-present FDR and PEP calibration fail |

## 16. Final verdicts

| dimension | verdict | decisive evidence |
|---|---|---|
| PSM TIE HANDLING | **FAILED VALIDATION** | joined input reorder changes q<.01 discoveries 0 -> 101 |
| Q-VALUE IMPLEMENTATION | **STRONG EVIDENCE** | independent formula oracle, 24 tie permutations, six exact boundaries, finite correction and fast-path equivalence; conditional on valid winners |
| COMPLETE-NULL BEHAVIOR | **MODERATE EVIDENCE** | 0/30 at all six thresholds, but low power and narrow construction |
| FDR CALIBRATION | **FAILED VALIDATION** | mean adjusted FDP 1.810% at nominal 1%; above nominal at every tested threshold |
| PEP IMPLEMENTATION | **FAILED VALIDATION** | one-ULP permutation amplification, extreme-`p` conservation failure, and no valid posterior derivation |
| PEP CALIBRATION | **FAILED VALIDATION** | every populated bin optimistic; pooled signed/absolute error .018685 |
| FIXED-C CV | **STRONG EVIDENCE** | held-out label, composition, outlier, reorder, source, and mutation attacks pass |
| SELECT-C CV | **STRONG EVIDENCE** | true nested source path; fold selection and scores invariant to all attacks |
| ENSEMBLE CV | **STRONG EVIDENCE** | label-free features, spectrum-grouped folds, all new leakage attacks pass |
| PROTEIN GROUPING | **FAILED VALIDATION** | equal-score peptide order changes group count; mixed target/decoy group misclassification |
| PROTEIN COMPETITION | **MODERATE EVIDENCE** | fair/invariant on 317 hand pairs, but conditional on flawed upstream groups and no protein truth |
| PROTEIN PEP | **MODERATE EVIDENCE** | picked-path mislabeling is fixed and `NA` is unambiguous; Bayesian posterior remains unvalidated |
| MULTI-SEED REPRODUCIBILITY | **STRONG EVIDENCE** | 120/120 entrapment and 160/160 compact artifacts hash-identical; no seed excluded |
| CROSS-DATASET GENERALIZATION | **WEAK EVIDENCE** | compatibility on known datasets, weaker MSFragger/yeast/PXD behavior, no untouched truth set |
| PERFORMANCE/CORRECTNESS EQUIVALENCE | **WEAK EVIDENCE** | hot paths appear exact, but PEP amplifies one-ULP ordering effects and no C++ method equivalence exists |
| VALIDATION-SUITE QUALITY | **MODERATE EVIDENCE** | catches 12/12 historical mutations, yet passes with five new implementation defects and failed calibration |

## Final scientific judgment

The second repair is substantive, not cosmetic. It fixes the original single-file tie attack, the
specific select-C and ensemble leakage paths, connected-component protein grouping on ordinary hand
graphs, target-favoring protein ties, and peptide PEP mislabeled as picked protein PEP. The q-value
scan itself is now the strongest validated component.

It does **not** fix the previously demonstrated PEP optimism or signal-present FDR failure. It also
leaves new reproducible failures in joined-input competition, tie-sensitive PEPs, and end-to-end
protein construction. Therefore the scientifically defensible statement is narrow:

> `percolator-rs` is a deterministic, fast experimental reimplementation with a well-supported
> standalone direct-TDC q scan and leakage-resistant tested CV modes. Its reported PSM q-values fail
> the available signal-present calibration experiment, its PEPs are not calibrated posterior
> probabilities, and its protein output remains invalid for some supported inputs and unvalidated
> against protein truth.

Any broader claim of calibrated FDR, calibrated PEP, generally valid protein inference, or overall
scientific defensibility is rejected by this audit.
