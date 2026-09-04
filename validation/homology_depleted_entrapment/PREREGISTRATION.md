# Preregistration: pre-search homology-depleted entrapment

Frozen on 2026-08-31 after exact reproduction of the original seed-1 baseline
and before constructing or searching any intervention database.

## Question and causal contrast

Does removal, before searching, of entrapment proteins that can generate highly
native-homologous tryptic peptides reduce the excess of high-scoring false
entrapment targets over entrapment decoys after semi-supervised rescoring?

The causal comparison is the original database versus a homology-depleted
database versus three source- and protein-length-matched random-depletion
controls. The intervention is allowed to fail. No threshold or endpoint will be
changed after looking at intervention results.

## Frozen inputs and search model

- Six PXD032157 mzML runs from the canonical signal-present validation.
- Native target FASTA and plant entrapment records from the exact canonical
  `combined.fasta`.
- Comet 2019.01 rev. 5 embedded in Crux 4.0-b854e9b-2021-06-01.
- Deposited per-run parameters: trypsin, semi-tryptic search
  (`num_enzyme_termini=1`), two missed cleavages, MH+ 600--5000 Da,
  peptide length 1--63, 10 ppm precursor tolerance, static Cys +57.021464,
  I/L equivalence, concatenated automatic decoys (`decoy_search=1`), and five
  reported candidates per spectrum.
- Audited `percolator-rs` commit
  `e8d83d1c76e4cf651fdfcf22d98b0b499c35943a` and binary SHA-256
  `be9bf670bfd69df4dc3ba3b8be6c4c164acaf56a91f4a2819d115f49604b2c45`.
- Canonical rescoring command: `--canonical --no-select-c --seed 1
  --num-threads 1`; the canonical maximum iteration count is the program
  default (10). No production source will be changed or rebuilt during the
  experiment.

## Primary homology rule

The comparison universe is the canonical tryptic digest of the native and
entrapment proteins: cleavage after K/R except before P, zero to two missed
cleavages, unmodified MH+ in 600--5000 Da after the canonical static Cys
modification, peptide length 8--63, and no ambiguous residue. I and L are
canonicalized to L.

An entrapment peptide is a prohibited near-homolog when a native theoretical
peptide of the same length has:

- I/L-aware Hamming distance at most two; and
- I/L-aware sequence identity at least 85%.

Both conditions are required. Thus an 8--13 residue peptide can differ by at
most one residue, while a peptide of length 14 or greater can differ by at most
two. This is the single primary rule.

Rationale: the prior prospective lead was an exhaustive same-length,
I/L-collapsed distance-at-most-two enrichment, dominated by one- and
two-substitution actin/tubulin-like peptides. The 85% identity guard makes
"highly homologous" length-aware and prevents short, chance two-substitution
matches from defining the intervention. The theoretical peptide universe uses
canonical tryptic regions because these are coherent biological digestion
units, because most prior extreme-tail witnesses were fully tryptic, and
because the fully tryptic negative control showed that such regions remain
relevant even when semi-tryptic terminal asymmetry is removed. The actual
search remains canonical semi-tryptic.

No edit distance, mass-substitution list, alternate identity cutoff, alternate
minimum length, or post hoc homolog rule is a confirmatory endpoint. Other
similarity analyses may only be labeled exploratory.

## Unit of depletion

The unit is the whole entrapment protein. An entrapment protein is removed if it
generates at least one prohibited peptide under the primary rule. Native
proteins are never removed. Whole-protein deletion preserves a physically
coherent FASTA and avoids artificial protein termini, peptide records, or
post-search filtering. It is deliberately broader than deleting only the
matched peptide; the size-matched controls address the resulting opportunity
loss.

## Size-matched negative controls

Three random-depletion controls use frozen NumPy PCG64 seeds `130363`, `155921`,
and `196613`. Within each entrapment source proteome and fixed 50-residue
protein-length stratum (with a final >=2000 stratum), each control removes
exactly the same number of proteins as homology depletion. Sampling is from all
proteins in that stratum without using homology status. This matches source,
protein count, and protein-length opportunity while not targeting the most
homologous proteins. The exact theoretical/searchable peptide opportunity is
measured before search; any residual size mismatch is reported rather than
corrected after outcomes are seen.

The confirmatory size control is the mean of the three frozen controls. Each
seed is also reported separately. Percolator seed is held at 1 for the primary
comparison, so a database-control seed is never confused with a training seed.

## Databases and decoys

Five complete target FASTAs are built from the frozen inputs:

1. `original` (native plus all original entrapment proteins),
2. `homology_depleted`,
3. `size_control_130363`,
4. `size_control_155921`,
5. `size_control_196613`.

Comet regenerates its concatenated reversed decoys independently inside every
search (`decoy_search=1`). Old decoy matches and old PIN rows are never reused
as intervention results. Headers retain the unambiguous `ENT_` and `DECOY_`
conventions. Exact target and PIN hashes, target/entrapment collisions, and
I/L-collapsed peptide overlaps are reported.

## Outcomes

The primary endpoint is pooled `R_ent / D_ent` among seed-1 rescored PSMs with
reported `q < 0.01`, using strict inequality and the same six-run pooling as the
canonical 258/133 result. The primary causal contrast is:

`(homology_depleted - original) - (mean(size_controls) - original)`

on log ratio where all counts are nonzero; raw ratios and count differences are
also reported. A movement toward one that is materially larger than the
control mean supports the causal hypothesis. No pass/fail calibration cutoff
is imposed.

Predefined secondary q thresholds are `<0.001`, `<0.005`, `<0.01`, `<0.02`,
`<0.05`, and `<0.10`. At every threshold report accepted targets, `R_ent`,
`D_ent`, their ratio, the direct known-false fraction `R_ent/R`, and the
entrapment-adjusted FDP using the global pure-entrapment fraction among
non-mixed decoys. The direct fraction is the adjustment-free lower bound.

PEP bins are `[0,1e-4)`, `[1e-4,1e-3)`, `[1e-3,5e-3)`, `[5e-3,.01)`,
`[.01,.02)`, `[.02,.05)`, `[.05,.10)`, `[.10,.20)`, `[.20,.50)`, and
`[.50,1]`. For each report target PSM count, mean PEP, direct and adjusted
known-false fractions, entrapment target/decoy counts and ratio, residuals, and
Wilson uncertainty for the direct fraction. Pooled signed adjusted calibration
error uses all populated bins exactly as in the prior audit; direct `f=1`
calibration is reported beside it.

The prior bin-level relationship is retested over bins with positive predicted
PEP and nonzero entrapment counts using Pearson correlation between
`log(observed adjusted / predicted)` and `log(R_ent / D_ent)`, plus the median
absolute log-fold residual. Bins are not merged to improve the relationship.

## Raw-score, training, and feature diagnostics

Raw Comet XCorr is evaluated after the same spectrum-key competition used in
the prior audit. Frozen matched entrapment-decoy depths are 25, 50, 100, 133,
250, and 500; depth 133 is the headline comparison. This is computed before
interpreting rescored outcomes.

Training-depth diagnostics use `maxiter` 0, 1, 2, 3, and 10 with Percolator
seed 1 in every database condition. They report the q<0.01 internal-null ratio,
pooled PEP error, and accepted target count. The change from maxiter 0 to 10 is
the predefined amplification contrast.

The feature diagnostic compares all features with PIN copies lacking exactly
`enzN` and `enzC` in `original` and `homology_depleted`, at seed 1/maxiter 10.
It is diagnostic only and cannot select production features.

Training-seed reproducibility is checked with seeds 1, 2, and 3 for the
canonical maxiter-10 `original` and `homology_depleted` conditions. Seed spread
is not a biological confidence interval.

## Uncertainty and interpretation

- Exact seed-to-seed ranges describe algorithmic reproducibility only.
- Wilson intervals describe binomial PSM-level counting only and ignore
  clustering.
- A six-run cluster bootstrap (frozen seed `20260831`, 10,000 replicates)
  preserves all PSM multiplicity within LC-MS/MS runs and is the primary
  sampling-uncertainty sensitivity analysis.
- With only one dataset family, dataset-level generalization remains
  unquantified and is stated as a limitation.

The hypothesis is classified after analysis as STRONGLY SUPPORTED, SUPPORTED,
PARTIALLY SUPPORTED, NOT SUPPORTED, or INCONCLUSIVE. Strong support requires a
clear movement of the homology arm toward exchangeability, a larger change than
the control mean and each control's uncertainty permits, and concordant PEP/FDP
improvement. Similar changes in homology and random arms favor a search-space
size explanation. No change in the homology arm falsifies sufficiency of this
predefined channel. Partial changes are not upgraded by exploratory filtering.

Remaining extreme-tail false targets are characterized only after the frozen
analyses. No discovered pattern triggers a second confirmatory filter in this
experiment.

