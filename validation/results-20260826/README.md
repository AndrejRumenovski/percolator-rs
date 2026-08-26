# Machine-readable results, second repair (2026-08-26)

Raw outputs of the second repair's experiments. The prose that reads them is
[`../SECOND_REPAIR.md`](../SECOND_REPAIR.md); where the two differ, these files are authoritative.

Large per-run artifacts (result TSVs, stdout/stderr, full manifests) stay outside the repository
under `/run/media/andrej-rumenovski/New Volume/percolator_rs_repair2/`; `study-summaries.json`
carries their SHA-256 digests.

| file | produced by | what it holds |
|---|---|---|
| `competition-pre-repair.json`, `competition-repaired.json` | `validation/adversarial_competition.py` | nine permutations of a 200-spectrum exact-tie fixture, plus a four-way tie fixture and a second seed, per build. `permutation_invariant` and `distinct_winner_sets` are the verdict fields. |
| `cv-pre-repair.json`, `cv-repaired.json` | `validation/adversarial_cv.py` | fold-level and per-row held-out-label attacks on fixed-C, `--select-c` and `--ensemble`, per build. |
| `mutation-results.json` | `validation/mutation_test.py` | twelve reintroduced defects and the tests that caught each. |
| `entrapment-rootcause.json` | `validation/entrapment_rootcause.py` | raw-search-score and rescored entrapment accounting side by side, under three foreign-fraction assumptions. |
| `entrapment-rootcause-raw.json` | same | earlier run that also included `comet-out`, a duplicate of one of the six datasets; superseded, kept for the record. |
| `rootcause-maxiter{1,3,10}.json` | same | entrapment accounting at three semi-supervised training budgets. |
| `pep-entrapment-rust-seed{1..5}.json` | `validation/pep_entrapment.py` | PEP calibration bins over the six-dataset entrapment set, per seed. |
| `pep-by-dataset-*.json` | same | the same analysis per dataset at seed 1. |
| `study-summaries.json` | this directory's build step | the entrapment, complete-null and multi-seed manifests reduced to their calibration tables, counts, input digests and per-artifact SHA-256. |

Two arms are reported as **negative results** and should not be read as repairs: the entrapment
curve is statistically unchanged by this repair, and rebuilding the PEP estimator did not improve
its calibration.
