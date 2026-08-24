# Retention time, joint training, and intra-file threading

## Outcome

The advanced-feature measurements are now tied to deterministic PXD032157 input-selection rules,
input SHA-256 digests, a runnable driver, and three repeated runs. Earlier README figures did not
identify their input files and therefore could not be independently regenerated from the repository.
The current results replace those unverifiable figures.

At strict reported q<0.01, coarse ScanNr-based retention-time features change PSM yield by -4.81%,
+3.98%, and +2.00% on the three lexicographically first PIN files. This confirms the original
qualitative conclusion: the feature can help, but the ScanNr proxy is not consistently beneficial.

Joint training on the four smallest PIN files increases aggregate PSM yield from 1,524 to 1,606
(+82, +5.38%). Three files improve by 25, 30, and 43 PSMs; one loses 16. The shared model therefore
helps this prospectively defined small-file group, but not every constituent run.

On the largest PIN, three-fold threading reduces the fixed-weight median from 2.05 to 1.39 seconds.
Nine-thread legacy class-weight selection reduces its median from 3.93 to 1.85 seconds. Target PSM
and peptide files are byte-identical across the compared thread counts. These timings are medians
of three runs and include result-file writes.

Machine-readable results are in
[`advanced-feature-results.tsv`](advanced-feature-results.tsv); exact input sizes and hashes are in
[`advanced-feature-inputs.tsv`](advanced-feature-inputs.tsv). Complete repeated-run logs, timing
records, result files, and hashes from the 2026-08-24 run are under
`$HOME/percolator_rs_out/advanced-features`.

## Design

The driver uses rules fixed independently of treatment outcomes:

- retention time: the first three filenames under bytewise lexical ordering;
- joint training: the four smallest PIN files by byte size;
- threading: the single largest PIN file by byte size.

Every case uses `--canonical --seed 1`. Retention-time cases compare the default with
`--rt-features`. Joint cases compare four independent fits with one `--join` fit. Thread cases
compare `--num-threads 1` with 3 for fixed weights and 1 with 9 for `--select-c`. The driver requires
exactly 65 input PINs by default, records the selected inputs and their hashes, checks that all
repeats have identical counts and output hashes, and separately checks serial/parallel output
identity.

Run it with:

```bash
REPEATS=3 bash bench/advanced_features.sh
```

Override `ADVANCED_BENCH_INPUT`, `ADVANCED_BENCH_OUT`, or
`ADVANCED_BENCH_EXPECTED_FILES` for a different fixture. Timing is host-sensitive; counts, selected
input hashes, and result hashes are the reproducibility checks.

The experimental search-engine ensemble has no quantitative yield claim in this repository. It
requires multiple compatible searches of the same raw run, which the benchmark inputs do not
provide; combining unrelated runs would create false ScanNr agreement. Its namespacing,
cross-engine support counts, and fold grouping are instead covered by Rust unit tests. A biological
ensemble-yield claim should not be added until a matched multi-engine dataset and calibration study
are pinned. `tests/ensemble_regression.sh` additionally exercises the complete CLI, exact candidate
deduplication, protein-output rejection, and serial/parallel output identity using duplicate fixture
views; it is explicitly a structural test rather than biological evidence.
