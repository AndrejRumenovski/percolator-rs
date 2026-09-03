# Frozen baseline

The production baseline is commit
`e8d83d1c76e4cf651fdfcf22d98b0b499c35943a`.  Before the baseline was
captured, the checkout had no changes under `src/`, `tests/`, `Cargo.toml`, or
`Cargo.lock`; the existing validation-report work is listed verbatim in
`baseline/e8d83d1/git-status.txt` and its tracked patch is preserved as
`pre-existing.patch`.

## Build and portable validation

- Command: `cargo build --release --locked`
- Binary: `target/release/percolator-rs`
- Size: 1,136,560 bytes
- SHA-256: `be9bf670bfd69df4dc3ba3b8be6c4c164acaf56a91f4a2819d115f49604b2c45`
- Toolchain: Rust/Cargo 1.97.0, LLVM 22.1.6
- `cargo test --release --all-targets --locked`: 126 tests passed across the
  binary and integration targets, with no failures
- All six portable shell gates passed: fixture regression, ensemble,
  selected-C, MLP, protein inference, and feature reporting

The committed `tests/fixtures/sample.pin` input is 3,101,218 bytes with SHA-256
`20a97c64394f93db705b0eb5dfd91566b76013e8f01d835153feccaac00f0737`.

## Canonical output freeze

`baseline/e8d83d1/outputs` contains the actual target and decoy PSM and peptide
TSVs for fixed-C, selected-C, and ensemble modes.  Fixed-C and selected-C also
contain picked target/decoy protein TSVs, including the required protein
`posterior_error_prob=NA` representation.

The fixed-C seed-1 outputs were byte-identical:

- across repeated serial executions;
- between `--num-threads 1` and `--num-threads 3`;
- for PSM, peptide, and protein files.

Every file's SHA-256 and byte count is recorded in `manifest.json`.  These files,
not newly generated expectations, are the byte-level refactor acceptance oracle.

## Adversarial observations

The fresh `final_repair_adversarial.py` rerun and both standalone final-repair
Rust probes are preserved under `baseline/e8d83d1`.  The following current
behavior is frozen:

- fixed-name joined file/row/label permutations are invariant;
- exact single-file tie winner identities, scores, q-values, and target PEPs are
  invariant, while decoy PEP presentation is not invariant;
- fixed-C, `--select-c`, and ensemble attacks detect no held-out-fold leakage;
- exact duplicate candidate rows create 101 false discoveries in the recorded
  adversarial construction;
- joined path aliases are not invariant;
- protein mapping in the competition construction changes with the seed;
- the independent q-value oracle passes 40,928 cases and 327,424 fast-path
  comparisons, while the accepted extreme `p` PEP-floor counterexample remains;
- protein grouping/tie probes pass their current graph and fairness checks,
  while the current picked pairing-key collision remains.

Negative observations are included because changing them belongs to a scientific
repair, not this architecture refactor.

## Preserved major scientific runs

The immediately preceding final scientific audit ran the expensive complete-null,
entrapment, multi-dataset, and protein-calibration workloads with the exact same
binary SHA-256 frozen above.  Their authoritative external manifests are recorded
in `manifest.json` with hashes.  Key observations are:

- complete null: 0/30 runs made any false discovery at each declared threshold;
- signal-present entrapment: mean adjusted FDP 0.018103981884256583 at nominal
  strict `q < 0.01` across five seeds;
- entrapment PEP weighted absolute and signed error:
  0.01868544626219112;
- all compact repeated outputs were byte-identical;
- picked-protein PEPs were all `NA`, Bayesian PEPs were numeric, and the then
  committed protein evaluator rejected the current `NA` schema.

This refactor reran the current end-to-end structural adversarial driver and
standalone arithmetic/graph probes.  The large empirical manifests are referenced
rather than duplicated because their recorded binary hash is identical.

## Fresh performance baseline

The 65-file PXD032157 workload comprises 2,295,401,156 input bytes.  Each number
below is the median of three fresh canonical fixed-C runs; all output hashes and
per-run yields are in `manifest.json`.

| Workload | Runs | Wall times (s) | Median (s) |
|---|---:|---|---:|
| Largest PIN, one thread | 3 | 1.765, 1.589, 1.609 | 1.609 |
| Largest PIN, three threads | 3 | 0.881, 0.883, 0.882 | 0.882 |
| Full workload, sequential | 3 | 53.803, 49.309, 49.301 | 49.309 |
| Full workload, file-level N=4 | 3 | 15.125, 15.256, 15.241 | 15.241 |

The full canonical workload produced 106,823 target PSMs and 35,886 target
peptides at strict q<0.01 on every repetition.  Performance comparisons should
use medians and investigate rather than compensate for a meaningful slowdown.

## Artifact identity

`baseline/e8d83d1/manifest.json` is the machine-readable index.  After the
clarifying annotation of pre-existing worktree state and external scientific
evidence, its SHA-256 is recorded when the baseline commit is made.
