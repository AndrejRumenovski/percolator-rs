# Protein-level calibration on the PrEST homology standard

> **These measurements predate the 2026-08-25 statistical repair** and describe commit `d83a7ba`,
> whose q-value and PEP estimators, cross-validation isolation and PIN feature selection were all
> subsequently found defective and replaced. They are kept as the record of what was measured then.
> For what the current implementation does and what it has been revalidated against, see
> [`../validation/REPAIR.md`](../validation/REPAIR.md).

## Why this benchmark

[PXD008425](https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD008425) was
constructed specifically to test protein inference under shared-peptide ambiguity. It contains two
pools of partially overlapping Protein Epitope Signature Tags (PrESTs): pool A has 192 sequences and
pool B has 191. Each pool contains one member of an overlapping pair, so a peptide can be shared by a
protein known to be present and its known-absent counterpart. The deposited database also contains
1,000 length-matched PrEST sequences absent from every sample.

The 12 Q Exactive runs comprise triplicates of pool A, pool B, A+B, and a blank containing only the
E. coli background. This supplies explicit protein-level present/absent labels, a strong null control,
and substantially more shared-peptide ambiguity than a conventional UPS spike-in. The standard and
its canonical group-level null hypothesis are described by
[The et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC6474350/): a reported group is considered
present if at least one member is present in that vial.

## Result

Calibration-only selection chose **α=0.1, β=0.0001, γ=0.001**, with peptide prior π=0.1. All
80 grid points converged within the 1,000-iteration selection cap. The selected model converged on
all held-out runs; the hardest selected run required 447 iterations.

At reported q ≤ 0.01 on the final test replicate:

| Method | A: accepted / absent / adjusted FDP (95% CI) | B | A+B | Blank: accepted / absent |
|---|---:|---:|---:|---:|
| Picked | 180 / 5 / 3.23% (1.38–7.36%) | 196 / 7 / 4.14% (2.02–8.34%) | 220 / 9 / 5.66% (3.00–10.50%) | 17 / 17 |
| Bayesian, fixed | 215 / 40 / 21.60% (16.22–28.26%) | 227 / 48 / 24.53% (18.95–31.23%) | 379 / 27 / 9.85% (6.83–14.06%) | 3 / 3 |
| Bayesian, selected | 176 / 3 / 1.98% (0.68–5.68%) | 185 / 5 / 3.14% (1.35–7.16%) | 321 / 3 / 1.29% (0.44–3.75%) | 3 / 3 |

Validation independently showed the same ordering. Selected Bayesian adjusted FDP was 2.61% for A,
5.04% for B, and 0.43% for A+B, versus 23.64%, 25.16%, and 9.51% with fixed defaults. On the final
test aggregate, selected Bayesian Brier score improved from 0.376 to 0.026 and 10-bin ECE from 0.523
to 0.023; ROC AUC was 0.983 selected, 0.980 fixed, and 0.989 picked. Group partitions differ, so
aggregate accepted-group counts are not sensitivity comparisons.

The conclusion is deliberately limited: **selection fixes most of the γ=0.5 failure, but does not
validate nominal 1% protein FDR**. The B test interval excludes 1%, and every method reports false
PrESTs in the held-out blank, where every PrEST is absent. Selected parameters are therefore recorded
for this search and retained as an explicit option, not installed as new global defaults. Fixed
defaults are strongly anti-conservative on this standard; picked inference also fails nominal 1%.

Across the four test runs, complete-process median wall times summed to 0.10 s picked, 0.20 s fixed
Bayesian, and 0.43 s selected Bayesian. Peak per-run RSS was 10.0, 12.7, and 12.4 MiB respectively.
The selected cost increase comes from the higher convergence cap and difficult B components, not
parameter-grid time. Full held-out rows are committed in
[`protein-calibration-results.tsv`](protein-calibration-results.tsv); generated reports retain every
threshold and composition row.

## Leakage-resistant design

Replicates have fixed roles before any inference is run:

| Replicate | Role |
|---:|---|
| 1 | calibration; select Bayesian α/β/γ |
| 2 | held-out validation |
| 3 | held-out final test |

The parameter grid is α ∈ {0.01, 0.05, 0.1, 0.3, 0.5}, β ∈ {0.0001, 0.001, 0.01, 0.05}, and
γ ∈ {0.001, 0.01, 0.1, 0.5}. Peptide prior π remains fixed at 0.1. Each combination is evaluated
only on the four replicate-1 runs. The deterministic selection objective is:

```text
mean Brier score + mean 10-bin ECE + 0.25 * (1 - mean ROC AUC)
```

Probability calibration is primary; the smaller ranking term prevents a flat prior from winning by
calibration alone. No replicate-2 or replicate-3 label or result participates in selection. Reports
keep the published fixed defaults (α=0.1, β=0.01, γ=0.5) alongside picked FDR and the selected model,
so parameter selection cannot erase the baseline.
Grid candidates get up to 1,000 belief-propagation iterations and are ineligible unless every
calibration run converges. The selected model retains that iteration cap on held-out runs; fixed
defaults retain the normal CLI default.

## Reproducible processing

`bench/protein_calibration/run.sh` performs the complete workflow:

1. Download all 12 Thermo RAW files and the three canonical FASTAs from PRIDE. Every deposited file
   is checked against the SHA-1 value returned by the PRIDE Archive API.
2. Convert RAW to centroided indexed mzML with the pinned self-contained Linux build of
   ThermoRawFileParser v2.0.0-dev, retaining MS2 spectra.
3. Concatenate the 1,383 target PrESTs and explicit paired reversed decoys. The script checks the
   expected A/B/random entry counts and writes a ground-truth table.
4. Search every run separately with pinned Sage 0.14.7: fully tryptic peptides, two missed cleavages,
   static carbamidomethyl cysteine, no variable modifications, 10 ppm precursor tolerance, and 20 ppm
   fragment tolerance. This follows the original benchmark's no-variable-modification and 10 ppm
   design while using Sage to produce a portable PIN directly.
5. Remove Sage's already-trained posterior-error feature before rescoring, canonically sort records,
   assign stable order-independent spectrum IDs, select Bayesian parameters on replicate 1, and run
   picked, fixed-Bayesian, and selected-Bayesian inference on all runs.
6. Write threshold calibration, ranking, probability calibration, composition, wall-time, and peak-RSS
   reports.

Conversion outputs are keyed by the converter checksum. Search and inference inputs are keyed by the
checksums of the database, Sage archive, search configuration, normalizer, and converter. A changed
pipeline therefore gets a new artifact directory instead of silently reusing stale PINs.

## Metrics and interpretation

For each run and method the report includes:

- known-present and known-absent groups at q ≤ 0.001, 0.005, 0.01, 0.02, 0.05, and 0.1;
- observed entrapment FDP (`known absent / accepted`), which is a direct lower bound because false
  matches landing on a sequence that happens to be present cannot be observed;
- search-space-adjusted FDP, dividing the observed false count by the vial-specific absent fraction of
  the 1,383 target sequences, with Wilson 95% intervals for the observed entrapment proportion scaled
  by the same fraction;
- tie-aware ROC AUC and normalized partial AUC through 5% FPR;
- Brier score and 10-bin expected calibration error for Bayesian protein posterior error
  probabilities (reported as NA for picked inference, whose shared-schema PEP field is the best
  peptide PEP rather than a protein-group posterior);
- group composition, including present/absent mixed groups and the 1,000-protein random entrapment;
- median complete-process wall time and peak RSS over three runs.

The count-based adjustment matches the benchmark's length-matched entrapment design but remains an
estimate: equal protein counts do not guarantee equal detectable-peptide search space. The raw
entrapment FDP and full composition table are retained so conclusions do not depend only on that
adjustment. Picked and Bayesian partitions can differ, so comparisons emphasize calibration and
ranking rather than assuming one-to-one group identity.
The intervals are descriptive: they condition on the count-based absent fraction and do not model
dependence among groups or uncertainty in effective peptide search space.

## Reproduce

```bash
bash bench/protein_calibration/run.sh
```

The workflow downloads roughly 12 GB of RAW data and stores all generated artifacts outside Git in
`$HOME/percolator_rs_out/protein-calibration`. Override this with `PROTEIN_CALIBRATION_OUT`; override
the Rust binary with `RS`; use `DOWNLOAD_JOBS` to control concurrent PRIDE downloads.
Set `REPEATS` to change the default three timing repetitions.

The final paths are printed by the script. Important files are `selected-params.json`,
`calibration-grid.tsv`, `report/thresholds.tsv`, `report/summary.tsv`, and
`report/composition.tsv`.
