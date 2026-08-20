# Benchmark dataset registry

`datasets.toml` is the single metadata registry for benchmark inputs. It is intentionally
separate from the existing benchmark scripts: adding a dataset only requires a new `[[datasets]]`
entry, and a future runner can consume the same registry without changing the rescoring code.

Validate it without downloading or running anything:

```bash
cargo run --quiet --bin validate-benchmark-manifests -- bench/datasets.toml
```

Each entry has these fields:

- Required: `id`, `source`, `organism`, `experiment_type`, `search_engine`, `pin_path`,
  `protein_level_evaluation`, and `notes`.
- Optional: `pride_accession`, `instrument`, `file_count`, `approximate_input_size`, and
  `preparation`.

`id` must be unique and use only letters, digits, `-`, and `_`. `pin_path` must end in `.pin` and
may contain `${UPPERCASE_ENV}` templates. For example, the committed PXD032157 entry uses
`${PERCOLATOR_BENCH_DATA}/PXD032157/**/*.pin`; set `PERCOLATOR_BENCH_DATA` to a local data root.
No benchmark data belongs in this repository.

The validator rejects unknown keys, missing required fields, unsupported schema versions, duplicate
IDs, invalid path templates, zero file counts, and empty text values with dataset-specific errors.
It intentionally does not expand paths or check whether files exist: that work belongs to the
future benchmark runner, so metadata validation remains portable and side-effect free.
