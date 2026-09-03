# Behavior-preserving refactor record

This directory is the durable record for the architecture refactor.  The
scientific method is a constraint, not a refactor target: a change is accepted
only when its relevant outputs and adversarial observations match the frozen
baseline.

The files have distinct roles:

- `ARCHITECTURE.md` records the pre-refactor dependency map, risk classification,
  target boundaries, and intentionally excluded work.
- `RESULT.md` records the resulting boundaries, deliberate non-changes, and
  final acceptance and performance evidence.
- `freeze_baseline.py` builds and exercises one exact revision and writes a
  machine-readable artifact directory.
- `baseline/e8d83d1/` is the baseline captured before production source was
  changed.  It contains canonical output files, command records, adversarial
  evidence, and repeated benchmark measurements.

The baseline command is:

```bash
python3 refactor/freeze_baseline.py \
  --output refactor/baseline/e8d83d1 \
  --full-benchmarks --repeats 3
```

The harness refuses to overwrite an existing artifact directory.  It records
the dirty worktree because pre-existing validation research is deliberately
preserved; production files under `src/` and the build manifest must be clean
when the baseline starts.

For each coherent production slice:

1. state the boundary problem and its risk in the commit message;
2. make a mechanical or narrowly scoped extraction;
3. run the focused unit/integration/adversarial gate;
4. compare applicable TSV files byte-for-byte with `baseline/e8d83d1/outputs`;
5. run the full portable gate before committing;
6. benchmark only when a measured hot path or data layout was touched.

Known adverse scientific observations are baseline behavior too.  A refactor
must not silently “repair” one of them, because doing so would change the method
and is outside this task.

The repeatable acceptance command for the current checkout is:

```bash
python3 refactor/verify_baseline.py
```

It builds the release binary, runs the release test suite and six portable
shell gates, compares all fixed/selected/ensemble TSVs by size and SHA-256, and
reruns the frozen adversarial driver plus its standalone probes. Temporary
outputs are kept outside the worktree and removed after the check.
