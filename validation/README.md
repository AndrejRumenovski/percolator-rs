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
4. [`REPAIR.md`](REPAIR.md) — **current.** The verified root causes, the corrections, the tests, the
   rerun of every predeclared experiment against a frozen build, and what still fails.

Documents 1–3 describe an implementation that no longer exists. They are the record of the failure
that document 4 responds to, and they are deliberately unchanged.

During the repair, the canonical implementation *was* changed — that was the point — but no
experiment was altered after its results were seen. Two arms were added afterwards to test a
hypothesis the results raised (spectrum-level competition), and both are labelled as new experiments
in [`REPAIR.md`](REPAIR.md) §11 rather than folded into the predeclared ones.

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
