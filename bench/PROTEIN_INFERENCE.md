# Picked-protein versus Bayesian inference

percolator-rs now provides two deliberately distinct protein-level models:

- `--protein-inference picked` (the default) retains the existing best-peptide ranking and
  picked target-decoy competition.
- `--protein-inference bayesian` implements the three-parameter noisy-OR model introduced by
  [Serang, MacCoss, and Noble](https://pmc.ncbi.nlm.nih.gov/articles/PMC2948606/). It reports
  marginal posterior probabilities for indistinguishable protein groups and derives q-values as
  the cumulative expected error from those posterior error probabilities.

These are not interchangeable q-value estimators. Picked-FDR is empirical target-decoy competition;
the Bayesian q-value depends on model fit and peptide PEP calibration. A larger Bayesian count is
therefore not evidence of higher accuracy without a known-mixture or entrapment evaluation.

## Model and implementation

For each protein, presence is Bernoulli with prior probability γ. Given `n` adjacent present
proteins, a peptide is emitted with probability:

```text
P(E=1 | n) = 1 - (1-beta) * (1-alpha)^n
```

The best PSM per peptide supplies a posterior error probability. Its likelihood contribution is
corrected against the peptide prior π following the original Fido formulation:

```text
P(D | n) ∝ ((1-PEP)/pi) * P(E=1 | n) + (PEP/(1-pi)) * P(E=0 | n)
```

Defaults are the published reasonable fixed values α=0.1, β=0.01, γ=0.5, and π=0.1. They are
overridable with `--protein-alpha`, `--protein-beta`, `--protein-gamma`, and
`--protein-peptide-prior`. `--protein-max-iter` controls inference iterations.

Proteins connected to exactly the same observed peptides are exchangeable. They are collapsed into
a count-valued variable with a binomial prior, preserving the probability that at least one group
member is present. Sum-product belief propagation is exact on tree-structured components. Cyclic
components use deterministic damped loopy belief propagation; the run log reports component counts,
iterations, and convergence. This is a Fido-style model, not the original Fido junction-tree or
graph-splitting implementation.

## Five-case benchmark

Same Ryzen 5 5600G host, local ext4 data, seed 1, canonical rescoring, and median of three complete
process runs. Timing includes the identical PSM rescoring stage, so it measures user-visible total
cost rather than isolated inference time. The PXD032157 case is the committed 12,000-PSM fixture,
not the full 8.6-million-PSM benchmark.

| Input | Picked groups / q<0.01 | Bayesian groups / q<0.01 | Wall, picked / Bayesian | Bayesian graph |
|---|---:|---:|---:|---:|
| PXD032157 fixture | 9,216 / 0 | 10,331 / 758 | 0.07 / 0.17 s | 9,216 components; 93 loopy |
| PXD007145 Tide | 17,164 / 4,037 | 17,164 / 3,687 | 0.45 / 0.68 s | 17,164; all trees |
| PXD020243 MSFragger | 5,849 / 190 | 5,849 / 163 | 0.04 / 0.08 s | 5,849; all trees |
| PXD060954 Sage | 6,061 / 1,598 | 6,149 / 1,542 | 0.23 / 0.33 s | 6,061; 23 loopy |
| Upstream yeast fixture | 15,105 / 460 | 15,814 / 215 | 0.10 / 0.22 s | 15,105; 42 loopy |

All five Bayesian runs converged in 37–46 iterations. Total-process cost was 1.4–2.4x wall time and
1.26–1.53x peak RSS relative to picked inference. On the four compact public extension cases,
Bayesian q<0.01 yield was 3.5–53% lower. It was dramatically higher on the redundant PXD032157
fixture, where many proteins share the same peptide evidence and the probability that at least one
member of a large group is present can be high under γ=0.5.

That PXD032157 reversal should be treated as a calibration warning, not a headline sensitivity gain.
The existing signal-present entrapment test already shows imperfect PSM-level calibration on this
search. A future protein-level known-mixture benchmark and dataset-specific α/β/γ selection are
needed before preferring one thresholded list as more accurate.

## Reproduction

First construct the compact public inputs, then run the comparison:

```bash
bash bench/multidataset/run.sh
bash bench/protein_inference.sh
```

Bulk outputs go to `$HOME/percolator_rs_out/protein-inference-benchmark` by default. Set
`PROTEIN_BENCH_OUT`, `MULTIDATASET_OUT`, `RS`, or `REPEATS` to override them. The recorded medians
are preserved in [`protein-inference-results.tsv`](protein-inference-results.tsv).
