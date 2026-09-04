# Why the signal-present entrapment experiment reports ~1.81% adjusted FDP at nominal q<0.01

Date: 2026-08-27. Audited binary `be9bf670…04b2c45` (commit `b2501f0`), unmodified.
No production code was changed during this investigation.
Machine-readable results: `ENTRAPMENT_ROOTCAUSE_RESULTS.json`.

## Short answer

The 1.81% is **not** a percolator-rs implementation defect, and it is not primarily an
artefact of the adjustment formula. It is a **real target-decoy exchangeability failure
created by the semi-supervised rescoring step**, which the reference C++ implementation
reproduces at the same magnitude, and which the entrapment design then **amplifies by a
factor of ~1.35** through an extrapolation the data does not support.

The whole phenomenon reduces to a single, adjustment-free fact:

> Among PSMs accepted at q<0.01, entrapment **targets** outnumber entrapment **decoys**
> 258 to 133 (z = +6.3). Both populations are certainly false. An exchangeable estimator
> would put them at 1:1.

## 1. Exact experimental reconstruction

| Element | What it actually is |
|---|---|
| Sample | 6 LC-MS/MS runs from PRIDE **PXD032157** — *Anopheles gambiae* atrium / male accessory glands / whole body, Q-Exactive HF |
| True target DB | `AnogambiaeVB54_AnocoluzziiMaliVB54_SacCereUP2021_RNAseq2016-3FT_contam.fasta` — 139,191 proteins, 145,351,640 aa. Two sibling *Anopheles* proteomes + *S. cerevisiae* + 3-frame RNAseq translations + contaminants |
| Entrapment DB | 9 plant proteomes (arabidopsis, rice, maize, soybean, wheat, tomato, potato, grape, apple), `ENT_` prefix, appended until entrapment aa ≥ native aa → 389,504 proteins, 145,351,799 aa, `entrapment_fraction = 0.50000027` |
| Decoys | Comet internal, `decoy_search=1`, `decoy_prefix=DECOY_`, peptide reversal; decoys generated for **both** halves |
| Search | Crux Comet 4.0, **`num_enzyme_termini=1` (semi-tryptic)**, 2 missed cleavages, 10 ppm precursor, `isotope_error=2`, no variable mods, `num_output_lines=5` |
| Search output | PIN with **5 ranked candidates per spectrum** (167,627 spectra, 837,367 rows over the 6 files) |
| Competition | percolator-rs performs `--post-processing-tdc`-equivalent spectrum-level competition **on the rescored values**, emitting exactly one PSM per spectrum |
| Entrapment counting | a target PSM is "known false" iff **every** reported protein starts with `ENT_`; PSMs mixing `ENT_` and native proteins are excluded from the numerator but retained in the denominator |
| Adjustment | `adjusted_fdp = (R_ent / f) / R`, with `f` = fraction of accepted non-mixed **decoys** that are entrapment decoys |
| Run config | `--canonical --no-select-c --seed {1..5} --num-threads 1`, each PIN run separately, results pooled |

## 2. Independent derivation of the adjusted FDP

With `f = D_ent/(D_ent + D_nat)` estimated from accepted decoys, the audit's formula is
algebraically identical to

```
adjusted_fdp = (D_ent + D_nat)/R  ×  R_ent/D_ent
             ≈        q           ×  R_ent/D_ent
```

Verified exactly on the audit outputs (seed 1, q<0.01): R=19,545, R_ent=258, D_ent=133,
D_nat=47, f=0.7389 → 0.01787 from either expression.

**Consequence.** The reported anti-conservatism is *entirely* the ratio of accepted
entrapment targets to accepted entrapment decoys. Everything about database sizes,
opportunity ratios and the 1/f extrapolation cancels out of that ratio. That makes
`R_ent/D_ent` an assumption-free internal null, and it is the statistic this
investigation drives.

**Does the adjustment assume equal opportunity, and does it hold?** The formula does not
assume 1:1 — it estimates `f` empirically. But the *design* claims equal opportunity and
does not achieve it:

| | native half | entrapment half |
|---|---|---|
| amino acids | 145,351,640 | 145,351,799 |
| distinct 7-mers (I/L collapsed) | 40,633,136 | 61,699,335 |
| distinct 7-mers per aa | 0.2796 | **0.4245** |
| realised share of decoy wins | 21.8% | **78.2%** |

The native half is heavily redundant (two sibling *Anopheles* proteomes plus 3-frame
translations of the same transcripts), so balancing amino acids produced a ~3.5:1
effective opportunity ratio, not 1:1. The decoy-derived `f` absorbs most of this, but the
design does not do what its documentation says.

## 3. TDC assumptions

1. Concatenated 1:1 target/decoy search — **holds** (decoys built for both halves).
2. One PSM per spectrum, the winner of a target-decoy competition — **holds** (competition on rescored values, 167,627 output PSMs = 167,627 spectra).
3. Under the null, an incorrect target beats its paired decoy with probability 0.5 — **fails** (see §5).
4. Decoy scores are exchangeable with incorrect-target scores — **fails** for the rescored score, **holds** for the raw search score.
5. Independence / no shared information across spectra — **partially fails**: 258 entrapment target PSMs come from only 150 distinct peptides, one peptide contributing 26.
6. Ties resolved without a label preference — **holds** (253 target/decoy ties at the top of 167,627 spectra, drawn by a seeded coin).
7. Entrapment hits are true negatives — **holds** (see §4).

## 4. Assumption checks

**Are the entrapment hits legitimately false?** Yes, established four ways:

- None of the top entrapment peptides occurs as an exact substring anywhere in the native FASTA (0/11 checked, and `pure` requires no native protein annotation).
- The 258 accepted PSMs touch 215 distinct entrapment proteins; **211 have exactly one distinct peptide and none has three or more.** Real protein presence produces protein-level coherence; this is its opposite.
- **No abundant plant-specific marker** (RuBisCO large/small, chlorophyll a/b binding, photosystem I/II) appears at all.
- The reviewed (SwissProt) accepted hits are dominated by **alpha-tubulin (9) and actin (6)** — universally conserved proteins.

So the labels are sound. But the hits are of a specific kind: **homologs of abundant true
peptides.** `IAPEEHPVLLTEAPLNPK` (26 PSMs) shares 17 of its 18 residues with a native
actin; `VVPEEHPVLLTEAPLNPK` (7 PSMs) shares 16 of 18; `AVFVDLEPTVIDEVR` (22 PSMs, alpha-tubulin)
shares 10. Across all 150 accepted peptides, 8.7% have a ≥12 aa substring shared with the
native proteome versus 4.7% of unaccepted entrapment peptides and 4.3% of accepted
entrapment decoys.

## 5. Hypotheses and outcomes

| # | Hypothesis | Prediction | Result | Verdict |
|---|---|---|---|---|
| A | percolator-rs implementation defect | C++ on the same PIN would be calibrated | C++ `--post-processing-tdc` gives R_ent/D_ent = 1.83, adjFDP 1.69% (rs: 1.94 / 1.79%) | **Falsified** |
| A2 | q-value arithmetic bug | the audit's own 40,928-case q-value oracle would fail; complete-null would fail | both pass | **Falsified** |
| C | Entrapment-adjustment bias (the 1/f extrapolation) | headline moves a lot with f | 1.32% (f=1) → 1.79% (f=0.739) → 2.64% (f=0.5) | **Supported as an amplifier**, not the origin |
| D | Search-engine / decoy-construction effect | raw Comet TDC would show the same asymmetry | Comet XCorr: R_ent/D_ent = 1.33, z = +1.25 (n.s.); at matched depth 1.01 | **Falsified as the origin** |
| E | Semi-supervised learning effect | asymmetry grows with training | monotone 1.29 → 1.92 over 0→3 iterations | **Supported — primary cause** |
| F | Competition mismatch (5 candidates/spectrum) | rank-1-only input would fix it | rank-1 only: 1.85 / 1.72% vs top-5 1.94 / 1.79% | **Minor contributor (~10%)** |
| G | Database construction / unequal opportunity | `f` would be misestimated | opportunity is 3.5:1 not 1:1, but `f` measures it | **Real design flaw, largely absorbed** |
| H | Exchangeability violation between real and reversed sequences | features would carry provenance information | leave-one-file-out AUC 0.521 vs permuted 0.502, driven by `enzN`/`enzC` | **Supported — the mechanism** |

## 6. Controlled experiments

All on the identical six PINs, all with the frozen audited binary.

**6.1 Search engine only (no Percolator).** Top-1 per spectrum, own TDC q-values:

| score | R at q<0.01 | R_ent/D_ent | binomial z | adjusted FDP |
|---|---|---|---|---|
| Comet XCorr | 4,422 | 1.333 | +1.25 | 1.15% |
| Comet ln(E-value) | 11,728 | 1.309 | +1.83 | 1.21% |
| **percolator-rs** | **19,545** | **1.940** | **+6.32** | **1.79%** |
| C++ percolator `--post-processing-tdc` | 19,404 | 1.831 | +5.76 | 1.69% |
| percolator-rs, rank-1-only PIN | 19,632 | 1.852 | +5.86 | 1.72% |

**6.2 Internal null at matched depth.** Walking each list to the same number of accepted
entrapment decoys removes the sensitivity difference:

| entrapment decoys accepted | Comet XCorr | Comet E-value | perc-rs | perc-rs rank-1 | C++ perc |
|---|---|---|---|---|---|
| 25 | 1.200 | 2.080 | 4.280 | 4.080 | 4.440 |
| 50 | 1.120 | 1.520 | 3.100 | 3.060 | 3.120 |
| 133 | **1.015** | 1.128 | **1.932** | 1.857 | 1.857 |
| 400 | 0.968 | 0.993 | 1.350 | 1.327 | 1.373 |
| 1500 | 1.001 | 0.994 | 1.135 | 1.113 | 1.114 |

The search score is exchangeable on this null at every depth. Both Percolator
implementations are not, and by the same amount.

**6.3 Training dose-response** (`--maxiter`, mean of seeds 1–3, q<0.01):

| maxiter | accepted targets | R_ent/D_ent | adjusted FDP | binomial z |
|---|---|---|---|---|
| 0 (no training) | 11,717 | 1.290 | 1.196% | +1.73 (n.s.) |
| 1 | 16,769 | 1.514 | 1.429% | +3.59 |
| 2 | 18,410 | 1.828 | 1.717% | +5.58 |
| 3 | 18,998 | 1.918 | 1.783% | +6.11 |
| 5 | 19,374 | 1.858 | 1.720% | +5.88 |
| 10 (default) | 19,527 | 1.924 | 1.777% | +6.22 |

With training disabled the estimator is calibrated on its own internal null. Each
iteration increases the violation until it saturates by iteration 3.

**6.4 Provenance leakage in the feature space.** A logistic model on the 21 PIN features,
trained to separate entrapment **targets** from entrapment **decoys** — two populations
that are both certainly false — reaches leave-one-file-out AUC **0.5209** (sd 0.0051)
against a label-permuted control of 0.5024. The largest standardised weights are
`enzN` (+0.139) and `enzC` (+0.104): Comet's peptide-reversal decoys, in a **semi-tryptic**
search, have systematically fewer tryptic termini than the real sequences they mirror.
That difference has nothing to do with whether a match is correct, and a discriminative
rescorer is free to use it.

Note that percolator's own score has bulk AUC 0.5012 on the same contrast. **The violation
is not a uniform distributional shift; it lives in the extreme upper tail**, where
homology-driven matches that look like well-behaved real identifications sit and where
reversed decoys structurally cannot follow.

**6.5 Candidate-pool multiplicity.** The PIN offers 2.66 target and 2.34 decoy candidates
per spectrum on average (measured on `22Oct2014-…-8-MAGs-S-01`). Competing on rescored values therefore gives the target side more
draws. Stratifying the accepted set: spectra with n_T > n_D give R_ent/D_ent = 4.66 and
adjFDP 2.00%; spectra with n_T < n_D give 0.60 and 1.21%. But feeding rank-1-only PINs
moves the pooled result only from 1.79% to 1.72% (seeds 1–3: 1.72%, 1.83%, 1.58%). Real,
small, within seed noise.

## 7. Threshold curves (mean of 5 seeds)

| nominal q | adjusted FDP | ratio to nominal | unadjusted entrapment FDP (f=1, lower bound) | ratio | R_ent/D_ent | accepted targets |
|---|---|---|---|---|---|---|
| 0.001 | 0.615% | 6.15 | 0.468% | 4.68 | 9.51 | 12,946 |
| 0.005 | 1.205% | 2.41 | 0.949% | 1.90 | 2.88 | 17,986 |
| **0.01** | **1.810%** | **1.81** | 1.324% | 1.32 | 1.96 | 19,536 |
| 0.02 | 2.745% | 1.37 | 2.036% | 1.02 | 1.42 | 21,142 |
| 0.05 | 5.851% | 1.17 | 4.392% | 0.88 | 1.19 | 23,828 |
| 0.10 | 10.872% | 1.09 | 8.251% | 0.83 | 1.10 | 26,848 |

**Not multiplicative. Not threshold-uniform. Concentrated in the high-confidence tail**,
decaying smoothly to ~1.09 at q<0.10. Under the strict lower-bound reading (no
extrapolation) the method is anti-conservative only at q ≤ 0.02 and conservative at
q ≥ 0.05. The q<0.001 row rests on 6–9 accepted decoys in total and carries almost no
information.

## 8. Seed and dataset effects

Across 5 seeds the pooled estimate is stable (1.729%–1.884%, sd 0.058%) — but the seeds
share the same six files, so this measures algorithmic reproducibility, **not** sampling
error. Across the six LC-MS/MS runs, which are the real experimental units:

| run | accepted targets (seed 1) | mean adjusted FDP over seeds |
|---|---|---|
| 09Dec2015-…-5-atrium-P-12hpm | 814 | 2.09% |
| 22Oct2014-…-8-MAGs-S | 7,410 | 1.84% |
| 28May2015-…-22-atrium-S-24H | 345 | 1.20% (defined in 2/5 seeds) |
| 28May2015-…-23-atrium-P-24H | 307 | 0.58% |
| 28May2015-…-38-MAGs-P | 4,988 | 1.54% |
| 9March2015-…-29-MAGs-pellet | 5,681 | 2.05% |

The pooled 1.81% is PSM-weighted and dominated by three files. Only one dataset family,
one search engine and one entrapment construction were tested, so **dataset specificity
cannot be excluded** — though the mechanism (§6.4) is generic.

## 9. Statistical uncertainty

| basis | estimate | 95% interval |
|---|---|---|
| audit's 5-seed spread | 1.810% | ±0.058% (reproducibility only — not a CI) |
| cluster bootstrap over the 6 runs | 1.799% | [1.54%, 2.07%] |
| cluster bootstrap over peptide sequences | 1.794% | [1.25%, 2.47%] |
| unweighted mean over 6 runs, t-interval | 1.549% | [0.89%, 2.21%] |

Sensitivity to the entrapment fraction, seed 1, q<0.01: f=1 → 1.32%; f=0.7389 (as
audited) → 1.79%; f=0.7823 (all decoys, less tail-selection noise) → 1.69%; f=0.5 (the
design's stated balance) → 2.64%.

**1.8104% should be read as roughly 1.3%–2.5%.** The internal-null violation itself
(R_ent > D_ent) is far more robust than the calibrated number: PSM-level z = +6.3;
collapsed to distinct peptides, ratio 1.35–1.54 with z = +2.5 to +3.7 across the five
seeds. Significant, but weaker once PSM-level clustering is removed.

## 10. Most likely causal explanation

**Multiple contributing causes, dominated by an estimator limitation.**

1. **Primary (~0.5 pp of the 0.8 pp excess).** Percolator's semi-supervised training breaks target-decoy exchangeability in the high-confidence tail. Reversed decoys are not exchangeable with real-but-incorrect peptides in the PIN feature space (AUC 0.521, driven by `enzN`/`enzC` under a semi-tryptic search), and homology-driven false matches to conserved proteins look like well-behaved real identifications on exactly the features the learner rewards. The learner converts a signal that is invisible in the bulk into a tail excess. Causally demonstrated by the monotone `--maxiter` dose-response and by the reference implementation reproducing it.
2. **Secondary (~0.35× multiplier).** The entrapment adjustment divides by f ≈ 0.74, i.e. it posits 0.35 unobserved native false discoveries for every observed entrapment one. The homology channel produces entrapment false discoveries with no native counterpart (the native counterpart of a homolog hit is the *correct* answer, not another false discovery), so this extrapolation over-counts. Without it the excess is 1.32% rather than 1.81%.
3. **Minor.** Feeding 5 candidates per spectrum and re-competing on rescored values contributes ~0.07 pp. PSM-level clustering (150 distinct peptides behind 258 PSMs) inflates apparent significance. The database's unequal opportunity ratio (3.5:1, not the designed 1:1) is a genuine design flaw but is largely absorbed by the decoy-derived f.

**Why the complete-null test passed and this one failed.** `run_null.py` builds its null
by taking only the **decoy** rows of a real PIN and randomly relabelling half of them as
targets. Both classes are then reversed sequences, exactly exchangeable by construction.
That design is structurally blind to the failure mode found here. The two results are
fully consistent, and 0/30 carries no information about this question.

## 11. Alternatives not ruled out

- **Dataset / entrapment-construction specificity.** One organism, one search engine, one plant-proteome entrapment, one decoy scheme. The mechanism generalises in principle; the magnitude may not.
- **Semi-tryptic search as an aggravating factor.** `enzN`/`enzC` dominate the provenance signal. A fully-tryptic search might show much less. Not tested.
- **Comet-specific decoy generation.** Peptide reversal keeping the C-terminal residue is what creates the terminus asymmetry. Shuffled or Markov decoys might behave differently. Not tested.
- **Magnitude in an ordinary (non-entrapment) search.** Doubling the database with a foreign proteome makes 56% of target PSMs foreign, which changes the learning problem. The direction should carry over; the size may not.
- **Cross-fold peptide leakage.** Percolator splits folds by spectrum, not by peptide; one peptide contributed 26 PSMs across folds. Not isolated here.

## 12. Should production code change?

**No — not on this evidence.**

- No implementation defect was found. The q-value estimator passed its 40,928-case oracle; the competition matches the reference's `weedOutRedundantTDC`; the reference C++ implementation reproduces the effect at the same magnitude on the same inputs (1.69% vs 1.79%).
- The residual anti-conservatism is a property of semi-supervised target-decoy FDR estimation with reversed decoys, not of this implementation of it. Changing the estimator to hit 1% on this benchmark would be tuning against the benchmark.
- What **should** change is documentation and validation, not methodology:
  - state that reported PSM q-values are anti-conservative by roughly 1.3–2.0× at q ≤ 0.01 in the presence of a homologous foreign search space, and near-calibrated at q ≥ 0.05;
  - record that the complete-null test cannot detect decoy non-exchangeability, and add an internal-null check (accepted entrapment targets vs accepted entrapment decoys — no adjustment needed) as a standing regression;
  - correct `bench/ENTRAPMENT.md`: the database balances amino acids, not searchable sequence, and the realised opportunity ratio is ~3.5:1;
  - report the entrapment result with run-level or peptide-level uncertainty, not seed spread.

If a methodological change is ever pursued, the smallest scientifically justified one is
to **exclude the enzymatic-terminus features (`enzN`, `enzC`, `enzInt`) from the model
when the search is semi-tryptic and decoys are generated by peptide reversal**, since
those features carry sequence provenance rather than match correctness. That should be
validated on an independent entrapment construction first, and is **not implemented here**.

## 13. The experiment that would most reduce uncertainty

**A second, independent entrapment construction on the same spectra**, varying one factor
at a time:

1. a **fully-tryptic** search of the same combined database — isolates the `enzN`/`enzC` channel, which the feature analysis names as the dominant provenance signal;
2. a **non-homologous** entrapment proteome (e.g. shuffled-but-realistic sequences, or a prokaryotic proteome with conserved families removed) — separates the homology channel from the feature channel;
3. **peptide-shuffled instead of reversed** decoys — tests whether the terminus asymmetry is a Comet decoy-construction artefact.

If (1) or (3) removes most of the excess, the cause is decoy construction interacting with
the feature set, and a documented feature restriction is the right response. If the excess
survives all three, it is an intrinsic limitation of semi-supervised TDC and only
documentation should change.

## Direct answer to the question asked

Is the ~1.81% caused by percolator-rs being wrong, by the statistical assumptions not
holding, by the validation design, or by some combination?

**By the statistical assumptions not holding, amplified by the validation design.
percolator-rs is not wrong.**

- **percolator-rs being wrong: no.** The reference C++ implementation, run in a matched configuration on the identical PINs, gives 1.69% with R_ent/D_ent = 1.83 against percolator-rs's 1.79% and 1.94. No arithmetic or bookkeeping defect survived the checks.
- **Statistical assumptions not holding: yes — this is the origin.** Target-decoy exchangeability fails for the rescored score. The failure is absent from the raw search score (z = +1.25) and appears monotonically as semi-supervised training proceeds (z = +1.73 at 0 iterations → +6.2 at 10). The mechanism is that reversed decoys are distinguishable from real peptide sequences in the PIN feature space (AUC 0.521, driven by enzymatic-terminus features under a semi-tryptic search), and a discriminative rescorer exploits that difference in the tail where homolog matches sit.
- **Validation design: yes — it inflates the number by ~1.35×, and it overstates the precision.** The 1/f extrapolation invents native false discoveries that the homology channel does not produce; the unadjusted lower bound is 1.32%. The intended 1:1 opportunity ratio was not achieved (3.5:1). Reported uncertainty (±0.06%) is seed reproducibility, not sampling error; the honest interval is roughly [1.3%, 2.5%], and run-to-run values span 0.58%–2.09%.
- **A real effect underneath all of it: yes.** Stripped of every adjustment, more entrapment targets are accepted than entrapment decoys — 258 vs 133 at PSM level, 177 vs 123 distinct peptides. Both are certainly false. That is a genuine, if modest, anti-conservatism of the reported q-values, and it is shared with the reference implementation.
