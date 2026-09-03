# Refactor result

The behavior-preserving architecture refactor is complete. The executable is
now a consumer of the reusable library, and process-level concerns no longer
own the scientific policies or output contracts.

## Resulting boundaries

```text
main.rs + cli.rs
  parse/validate argv, load inputs, manage profiling, print diagnostics
  -> pipeline.rs
       rescore -> select reported PSMs -> peptide statistics -> protein dispatch
       -> competition.rs   precursor competition and seeded exact ties
       -> peptide.rs       identity, representatives, q/PEP, mapping union
       -> percolator.rs    CV topology, selection, learner orchestration
            -> preprocessing.rs  fold-local normalization and RT matrices
            -> svm.rs / mlp.rs / simd.rs / stats.rs / rt.rs
       -> protein.rs / protein_bayes.rs
  -> output.rs
       byte-compatible PSM/peptide/protein/feature serialization
```

The public library boundary now exposes parsing, learning, statistics,
competition, peptide/protein inference, pipeline composition, and output
serialization. `preprocessing.rs` is deliberately crate-private: it exists to
make fold-local fitting cohesive, not to expose partially valid matrices as a
public API.

Fold construction, fold setup, legacy C selection, and nested selection remain
together in `percolator.rs`. They share partition-specific state and ordering;
splitting them further would make the held-out/training data flow less visible.
Likewise, compact `i8` labels were not replaced with wrapper types inside hot
loops, and the statistics, SVM/MLP mathematics, PIN layout/parser, and Bayesian
message schedule were not rewritten.

## Acceptance evidence

- `cargo test --release --all-targets --locked`: 138 tests pass (126 frozen
  tests plus 12 new CLI, peptide, and pipeline characterizations).
- `cargo fmt --all -- --check` passes.
- `cargo clippy --all-targets --all-features --locked -- -D warnings` passes.
- `RUSTDOCFLAGS="-D warnings" cargo doc --no-deps --all-features --locked`
  passes.
- All six portable regression scripts pass.
- Fixed serial, fixed parallel, selected-C, and ensemble TSV files match the
  frozen sizes and SHA-256 hashes exactly, including protein output where
  supported.
- The frozen adversarial summary and both independent Rust probes match. Known
  adverse observations remain adverse; no scientific repair is hidden here.

The final three-repeat performance comparison on PXD032157 also retained every
per-file output hash and discovery count:

| Benchmark | Frozen median | Refactored median | Ratio |
|---|---:|---:|---:|
| Largest PIN, 1 thread | 1.609 s | 1.589 s | 0.988 |
| Largest PIN, 3 threads | 0.882 s | 0.883 s | 1.001 |
| 65 files, sequential | 49.309 s | 49.202 s | 0.998 |
| 65 files, 4 processes | 15.241 s | 15.154 s | 0.994 |

Run the complete non-benchmark acceptance gate with:

```bash
python3 refactor/verify_baseline.py
```

The frozen baseline and its provenance remain under `baseline/e8d83d1/`.
