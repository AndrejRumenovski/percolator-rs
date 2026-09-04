# Scientific validation workspace

This directory contains adversarial, non-optimization validation of percolator-rs. Every study
preserves exact commands, software versions, seeds, inputs, output hashes, and machine-readable
results. Negative results are kept; nothing here is rewritten after a repair.

Read in this order:

1. [`IMPLEMENTATION_AUDIT.md`](IMPLEMENTATION_AUDIT.md) — the read-only source audit of commit
   `d83a7ba`, written before any experiment and before any change.
2. [`SCIENTIFIC_VALIDATION.md`](SCIENTIFIC_VALIDATION.md) — the adversarial validation that
   **rejected** that implementation, with its claim audit, limitations and verdicts.
3. [`COMPLETION_AUDIT.md`](COMPLETION_AUDIT.md) — requirement-by-requirement closeout of that study.
4. [`REPAIR.md`](REPAIR.md) — the first repair: root causes, corrections, tests, and the rerun of
   every predeclared experiment against a frozen build.
5. [`INDEPENDENT_AUDIT.md`](INDEPENDENT_AUDIT.md) + [`INDEPENDENT_AUDIT_RESULTS.json`](INDEPENDENT_AUDIT_RESULTS.json)
   — the independent adversarial audit that **rejected** the repaired build, with the attack scripts
   `adversarial_competition.py`, `adversarial_cv.py`, `adversarial_parser.py`,
   `adversarial_protein.rs` and `independent_stats_probe.rs`.
6. [`SECOND_REPAIR.md`](SECOND_REPAIR.md) — the response to that audit: competition tie
   handling, the PEP estimator, protein grouping/competition/posterior, the two leaking CV modes, a
   root-cause attribution of the entrapment failure that remains, and a twelve-mutation test of the
   suite itself.
7. [`POST_REPAIR_AUDIT.md`](POST_REPAIR_AUDIT.md) — the independent adversarial audit of
   the second repair. It preserves the successful repairs but rejects the build scientifically after
   new joined-input tie, PEP, and protein-construction failures and unchanged signal-present
   anti-conservatism.
8. [`ENTRAPMENT_ROOTCAUSE.md`](ENTRAPMENT_ROOTCAUSE.md) + [`ENTRAPMENT_ROOTCAUSE_RESULTS.json`](ENTRAPMENT_ROOTCAUSE_RESULTS.json)
   — why the signal-present entrapment experiment reports ~1.81% adjusted FDP at nominal q<0.01.
   Diagnosis only; no production code was changed.
9. [`PEP_ROOTCAUSE.md`](PEP_ROOTCAUSE.md) + [`PEP_ROOTCAUSE_RESULTS.json`](PEP_ROOTCAUSE_RESULTS.json)
   — why the reported PEPs remain optimistic after the exact-zero repair. Reconstructs
   the estimator from source, verifies every transformation against an independent
   greatest-convex-minorant oracle, measures calibration under a generative model where the
   assumptions hold, and compares against QVALITY 3.09 on byte-identical inputs. Diagnosis only;
   no production code was changed and nothing was fitted to the calibration outcome. Probes:
   `pep_rootcause_probe.rs`, `pep_rootcause_lib.py`, `pep_rootcause_synth.py`,
   `pep_rootcause_violation.py`, `pep_rootcause_mechanism.py`, `pep_rootcause_real.py`,
   `pep_rootcause_tables.py`, `pep_rootcause_controls.py`, `pep_rootcause_depth.py`,
   `pep_rootcause_qvality.py`, `pep_rootcause_uncertainty.py`.
10. [`homology_depleted_entrapment/`](homology_depleted_entrapment/) — **current.** A prospective,
    preregistered, pre-search intervention testing whether native-homologous entrapment proteins
    cause the high-confidence target/decoy exchangeability failure. It compares a freshly searched
    original database, a homology-depleted database, and three source/length-matched random controls;
    preserves all search, PIN, rescoring, analysis, uncertainty, and integrity artifacts; and finds
    the hypothesis supported for the q<0.01 failure in this dataset without recommending a production
    methodology change.

Documents 1–3 describe an implementation that no longer exists. They are the record of the failure
that document 4 responds to, and they are deliberately unchanged. Document 4's own claims are
corrected where document 5 falsified them and document 6 says so explicitly; it is not rewritten.

Attack scripts take a `--binary` argument so the same attack can be run against an old and a new
build. `mutation_test.py` runs in a throwaway git worktree and reverts every mutation;
`entrapment_rootcause.py` reimplements the entrapment accounting from the raw counts rather than
calling the production crate.

During the first repair, the canonical implementation *was* changed — that was the point — but no
experiment was altered after its results were seen. Two arms were added afterwards to test a
hypothesis the results raised (spectrum-level competition), and both are labelled as new experiments
in [`REPAIR.md`](REPAIR.md) §11 rather than folded into the predeclared ones. A third new experiment,
[`run_null_variants.py`](run_null_variants.py), exists purely to try to falsify the repaired
estimator on null constructions it was never checked against.

`psm_agreement.py` compares matching target and decoy PSM outputs at the row level. It accepts
either individual TSV files or run directories containing matching relative paths. Example:

```bash
python3 validation/psm_agreement.py \
  --rust /path/to/rust/results \
  --cpp /path/to/cpp/results \
  --output /path/to/agreement.json \
  --table /path/to/agreement-thresholds.tsv
```

The JSON is authoritative. Threshold membership is strict (`q < threshold`) and thresholds are
predeclared as 0.001, 0.005, 0.01, 0.02, 0.05, and 0.10.
