# Frozen pre-search database characterization

This characterization was completed before any intervention Comet search. It
reads only target FASTAs and the frozen removal manifests. Complete length,
mass, enzymatic-terminus, and amino-acid distributions are in
`databases/presearch_summary.json` and
`databases/presearch_characterization.json`.

## Opportunity removed

| condition | entrapment proteins | residues | fully tryptic instances | unique fully tryptic (HLL) | searchable semi-tryptic instances | unique searchable (HLL) |
|---|---:|---:|---:|---:|---:|---:|
| original | 389,504 | 145,351,799 | 34,721,067 | 24,234,347 | 598,438,171 | 418,094,318 |
| homology-depleted | 189,363 | 49,874,647 | 10,736,954 | 8,403,939 | 201,087,250 | 156,459,902 |
| size control 130363 | 189,363 | 49,958,848 | 11,963,136 | 9,751,279 | 202,659,292 | 165,829,534 |
| size control 155921 | 189,363 | 49,954,424 | 11,961,153 | 9,728,091 | 202,581,927 | 165,724,726 |
| size control 196613 | 189,363 | 49,953,365 | 11,952,234 | 9,768,898 | 202,578,864 | 165,554,557 |

Homology depletion removes 200,141/389,504 proteins (51.38%), 65.69% of
entrapment residues, 66.40% of semi-tryptic peptide instances, and an estimated
62.58% of unique searchable sequences. Each random control removes exactly the
same source-by-50-aa-stratum protein counts. Relative to the control mean, the
homology arm has 0.75% fewer searchable peptide instances and 5.59% fewer
estimated unique searchable sequences. That residual mismatch was anticipated
in the preregistration and will not be corrected after observing outcomes; it
limits how sharply a homology effect can be separated from effective unique
search-space size.

The HLL estimates use `p=18` and have nominal relative standard error 0.203%.
Instance counts are exact for semi-tryptic peptides with one or two enzymatic
termini, zero to two missed cleavages, length 1--63, MH+ 600--5000 Da, static
Cys +57.021464, and no ambiguous residue.

## Similarity intervention check

| condition | retained proteins with a primary near-homolog | fraction of retained proteins | first-witness distance 0 / 1 / 2 | exact shared full-tryptic sequences |
|---|---:|---:|---:|---:|
| original | 200,141 | 51.38% | 6,698 / 189,792 / 3,651 | 7,857 |
| homology-depleted | 0 | 0% | 0 / 0 / 0 | 0 |
| size control 130363 | 81,391 | 42.98% | see JSON | 4,678 |
| size control 155921 | 81,187 | 42.87% | see JSON | 4,668 |
| size control 196613 | 81,183 | 42.87% | see JSON | 4,698 |

The intervention therefore eliminates the preregistered primary property,
whereas all three random controls retain about 81,000 qualifying proteins.
There are no duplicate accession headers and no I/L-canonical full-protein
native/entrapment collisions in any condition.

The rule is broad: 168,196/200,141 first qualifying witnesses are length 8 and
189,792 are one-substitution matches. This was observed only after the rule was
frozen. It is reported as a fixed design limitation and the threshold is not
changed. The controls are critical because short-peptide opportunity contributes
substantially to the removed set.

## Distribution matching

Compared with the mean of the three random controls, the homology-depleted
entrapment component has:

- 0.158% fewer residues;
- total-variation distance 0.0279 in searchable peptide length;
- total-variation distance 0.0277 in precursor-mass bins;
- total-variation distance 0.0158 in A--Z amino-acid composition; and
- 94.66% versus 94.10% one-enzymatic-terminus opportunity.

These differences are biological consequences of removing conserved proteins,
not post-result tuning. All condition-specific arrays and counts are preserved
machine-readably.

