# Scientific validation workspace

This directory contains adversarial, non-optimization validation of percolator-rs. The canonical
implementation is not changed to improve an experiment. Every study must preserve exact commands,
software versions, seeds, inputs, output hashes, and machine-readable results.

The implementation audit is in [`IMPLEMENTATION_AUDIT.md`](IMPLEMENTATION_AUDIT.md). It predates
methodological experiments and records failures rather than repairing them.

The consolidated results, claim audit, limitations, and final verdicts are in
[`SCIENTIFIC_VALIDATION.md`](SCIENTIFIC_VALIDATION.md).
The requirement-by-requirement closeout is in
[`COMPLETION_AUDIT.md`](COMPLETION_AUDIT.md).

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
