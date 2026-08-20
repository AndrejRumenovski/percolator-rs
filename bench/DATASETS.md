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
  `preparation`, and `reference_search_input`.

`id` must be unique and use only letters, digits, `-`, and `_`. `pin_path` must end in `.pin` and
may contain `${UPPERCASE_ENV}` templates. For example, the committed PXD032157 entry uses
`${PERCOLATOR_BENCH_DATA}/PXD032157/**/*.pin`; set `PERCOLATOR_BENCH_DATA` to a local data root.
No benchmark data belongs in this repository.

The validator rejects unknown keys, missing required fields, unsupported schema versions, duplicate
IDs, invalid path templates, zero file counts, and empty text values with dataset-specific errors.
It intentionally does not expand paths or check whether files exist: that work belongs to the
future benchmark runner, so metadata validation remains portable and side-effect free.

`reference_search_input` may be `concatenated` or `separate`. When set, the single-dataset runner
passes the corresponding `--search-input` option to C++ Percolator. Set it only when the input's
target/decoy interpretation is known; it is not a performance setting.

## Single-dataset runner

Build the Rust binary, then run both implementations over one manifest entry:

```bash
cargo build --release
cargo run --release --bin benchmark-dataset -- \
  --dataset PXD032157 \
  --output /local/benchmark-results \
  --percolator /path/to/percolator
```

The runner expands `pin_path`, verifies a known `file_count` before starting, processes every PIN
sequentially, and writes separate `rust/` and `cpp/` directories under
`OUTPUT/DATASET_ID/`. It writes `rust-summary.tsv`, `cpp-summary.tsv`, `per-file.tsv`, and
`failures.tsv`; every attempted file has a row, including non-zero exits. PSM, peptide, and protein
counts use the named `q-value` column and the strict threshold `q < 0.01`.
For safety, `OUTPUT/DATASET_ID` must not already exist; choose a fresh output directory per run so
stale result files cannot be reported as fresh output.

Use `--dry-run` with the same arguments to print the exact `/usr/bin/time` commands without running
or creating outputs. The runner uses seed 1 for both tools, Rust's `--canonical` profile, and
`--num-threads 1` for each implementation. It enables protein outputs only when
`protein_level_evaluation = true`. The algorithms' internal training implementations and defaults
are not forced to match beyond those public settings; the generated command preview is the
authoritative record of remaining CLI differences. Legacy PXD032157 scripts remain unchanged: use
`bench/regression.sh` for its existing parallel Rust-only performance gate.
