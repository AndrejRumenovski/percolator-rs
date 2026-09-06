# percolator-rs

[![CI](https://github.com/AndrejRumenovski/percolator-rs/actions/workflows/ci.yml/badge.svg)](https://github.com/AndrejRumenovski/percolator-rs/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

A from-scratch Rust implementation of the
[Percolator](https://github.com/percolator/percolator) semi-supervised peptide-spectrum match (PSM)
rescoring workflow. The project focuses on deterministic execution, fold-isolated training,
adversarial validation, and high-throughput processing of Percolator input (`.pin`) files.

> [!WARNING]
> **Research status (2026-09-04):** this is experimental research software. The current build is
> not scientifically defensible as an unqualified source of calibrated PSM, PEP, or protein-level
> confidence. Its direct target-decoy q-value arithmetic is well tested, but calibration depends on
> assumptions that fail in available signal-present experiments, and several known implementation
> defects remain. Read [Scientific status](#scientific-status) before interpreting q-values or PEPs.

## Quick start

The project requires a stable Rust toolchain. Release builds target `x86-64-v3` by default (roughly
Haswell-era Intel or Zen-era AMD and newer):

```bash
cargo build --release --locked
```

Run the canonical linear-SVM workflow on one concatenated target-decoy PIN:

```bash
./target/release/percolator-rs \
  --seed 1 \
  --results-psms target.psms.tsv \
  --decoy-results-psms decoy.psms.tsv \
  --results-peptides target.peptides.tsv \
  --decoy-results-peptides decoy.peptides.tsv \
  input.pin
```

Output files are optional and are written only when their corresponding flags are supplied. Progress,
configuration, timing, and q<0.01 yield summaries are printed to standard error.

To tune the binary for the build host instead of the portable project baseline:

```bash
RUSTFLAGS="-C target-cpu=native" cargo build --release --locked
```

## Input and output contract

### PIN input

The supported statistical design is a **concatenated target-decoy search** with an equal-sized target
and decoy database. Separate target and decoy searches requiring mix-max post-processing are not
supported.

A PIN must be tab-delimited and contain:

- a row identifier in the first column (conventionally `SpecId`);
- `Label`, with target and decoy rows conventionally encoded as `1` and `-1`;
- `ScanNr` in the metadata block immediately after `Label`;
- at least one finite numeric feature;
- `Peptide`, followed optionally by one or more protein fields.

Metadata names are matched case-insensitively. The contiguous metadata prefix after `Label` may
contain `ScanNr`, `ExpMass`, `CalcMass`, `rt`, `retentiontime`, `FileName`, and `SpectraFile`.
The first unrecognized column starts the feature block; every column from there to `Peptide` must be
finite numeric data. Malformed labels, scans, masses, or features stop the run with a file, line,
column, and offending-value diagnostic.

The default PSM competition keeps one winner per `(source, ScanNr, ExpMass)` after rescoring. This
supports PINs that report multiple candidates for a precursor, subject to the duplicate-candidate
limitation described below. When `ExpMass` is absent, it is treated as zero.

### TSV output

Target and decoy PSM and peptide files share this schema:

```text
PSMId  score  q-value  posterior_error_prob  peptide  proteinIds
```

Protein files use:

```text
ProteinGroupId  q-value  posterior_error_prob  score  numPeptides  proteinIds
```

The actual files are tab-delimited. PSM and peptide files are sorted by descending score; protein
rows follow their inference method's deterministic ranking. Picked-protein inference does not
estimate a protein posterior, so its `posterior_error_prob` field is `NA`; Bayesian protein output
contains a numeric value.

Threshold counts reported by the program and validation tools use strict comparison (`q < 0.01`),
not `q <= 0.01`.

## How the canonical workflow works

1. **Parse the PIN.** Input is memory-mapped, metadata is separated from numeric features, and
   malformed or non-finite required values are rejected.
2. **Build three outer folds.** Candidates from the same `(source, ScanNr)` group remain together.
3. **Fit preprocessing inside each training partition.** Feature normalization and the best initial
   feature direction are learned without using held-out rows.
4. **Train semi-supervised models.** Each iteration scores the training partition, takes target rows
   below the training FDR threshold as positives, uses all decoys as negatives, and refits the model.
   The default learner is an L2-regularized squared-hinge linear SVM with fixed `Cpos=1` and `Cneg=4`.
5. **Score held-out rows.** Each PSM is scored by a model that did not train on it. Fold scores are
   standardized against their training-decoy distribution before being pooled.
6. **Compete candidates.** By default, only the highest rescored candidate for each precursor is
   reported. Exact score ties use a deterministic seed-dependent draw over the emitted tied rows.
7. **Recompute reported-list statistics.** The TDC+ estimate uses
   `(D + 1) * p / (1 - p) / max(T, 1)`, evaluates exact-score tie groups together, and applies a
   reverse cumulative minimum. `p` is `--null-target-win-prob`, defaulting to `0.5`.
8. **Estimate PEPs.** PEPs are derived from increments of the cumulative false-discovery estimate and
   monotonized with PAVA. They are deterministic target-list scores, not validated posterior
   probabilities.
9. **Collapse higher levels.** Peptide reporting keeps the best PSM per `(label, modified core
   peptide)`. Optional protein inference then groups proteins by identical observed peptide evidence
   and applies either picked target-decoy competition or a Fido-style Bayesian model.

The default run is deterministic for a fixed input identity, path/name, seed, options, thread mode,
and build. Serial and fold-parallel output is byte-identical in the regression suite.

## Command-line reference

The executable currently has a small hand-written parser rather than generated `--help`; invoking it
without an input prints a short usage and input-contract summary. Unknown option names are currently
ignored, so check the startup configuration line and requested output files carefully.

### Output options

| Option | Meaning |
|---|---|
| `--results-psms PATH`, `-m PATH` | Target PSM TSV |
| `--decoy-results-psms PATH`, `-M PATH` | Decoy PSM TSV |
| `--results-peptides PATH`, `-r PATH` | Target peptide TSV |
| `--decoy-results-peptides PATH`, `-B PATH` | Decoy peptide TSV |
| `--results-proteins PATH`, `-l PATH` | Target protein-group TSV; enables protein inference |
| `--decoy-results-proteins PATH`, `-L PATH` | Decoy protein-group TSV; enables protein inference |
| `--feature-report PATH` | Linear-SVM feature report |

### Core rescoring options

| Option | Default | Meaning |
|---|---:|---|
| `--seed N` | `1` | Fold and deterministic tie seed |
| `--rescore-model svm\|mlp`, `--model ...` | `svm` | Fold-local learner; `linear` and `neural` are aliases |
| `--maxiter N` | profile value | Semi-supervised iterations |
| `--subset-max-train N`, `-N N` | `0` | Maximum training rows per fold; `0` uses all rows |
| `--num-threads N` | `1` | `1` runs folds serially; values above `1` enable fold/grid parallelism |
| `--null-target-win-prob P` | `0.5` | Null probability that an incorrect target beats its decoy; must be in `(0,1)` |
| `--psm-competition` | on | Keep one rescored precursor winner |
| `--no-psm-competition` | — | Report all candidates; resulting q-values are not claimed as FDR estimates |

For `k` equivalent decoys per target, the intended null setting is `1 / (1 + k)`. This does not
repair candidate duplication or other violations of the competition model.

### Execution profiles

Explicit `--maxiter`, `--subset-max-train`, `--cpos`, `--cneg`, and C-selection flags override the
profile regardless of argument order.

| Profile | Training subset | Iterations | C selection | Intended use |
|---|---:|---:|---|---|
| `--fast` | 20,000 | 5 | off | Quick QA and development |
| `--balanced` | 40,000 | 10 | off | Reduced-cost exploratory runs |
| `--canonical` | all | 10 | off | Default full-sensitivity workflow |

The equivalent long form is `--profile fast|balanced|canonical`.

### Linear-SVM options

| Option | Default | Meaning |
|---|---:|---|
| `--cpos F` | `1` | Absolute positive-class weight |
| `--cneg F` | `4` | Absolute negative-class weight |
| `--select-c` | off | Nested per-outer-fold class-weight grid search |
| `--no-select-c` | on | Use fixed class weights |
| `--svm-tolerance F` | `1e-5` | Positive finite solver tolerance |
| `--auto-model`, `--nested-select` | off | Nested selection of SVM scale, class weights, feature count, and tolerance |
| `--no-auto-model` | on | Disable automatic model selection |

`--auto-model` supports only SVM and cannot be combined with `--select-c` or explicit
`--cpos`/`--cneg`. The feature report is also SVM-only. It records mean out-of-fold raw
coefficients, standardized effects, fold variability, label correlation, selection frequency, and
held-out permutation importance with fitted models held fixed. See
[`bench/AUTOMATIC_SELECTION.md`](bench/AUTOMATIC_SELECTION.md) for the selection study.

### Experimental MLP options

| Option | Default | Constraint |
|---|---:|---|
| `--mlp-hidden N` | `8` | `1..256` hidden units |
| `--mlp-epochs N` | `10` | `1..1000` epochs per semi-supervised iteration |
| `--mlp-learning-rate F` | `0.02` | Finite and positive |
| `--mlp-l2 F` | `0` | Finite and non-negative |

The MLP is a deterministic one-hidden-layer experimental learner using the same outer folds and
reported-list statistics as the SVM. It did not improve aggregate yield in the recorded evaluation;
the SVM remains the default. See [`bench/DEEP_LEARNING.md`](bench/DEEP_LEARNING.md).

### Multiple inputs and biological features

`--join file1.pin file2.pin ...` pools multiple compatible runs into one training problem. Inputs
must have the same feature layout. File and row argument order are canonicalized, but joined source
identity remains sensitive to lexical filenames; renaming otherwise identical inputs can change
folds and exact-tie draws.

`--ensemble ENGINE1=file1.pin ENGINE2=file2.pin ...` combines at least two search-engine views of the
same run. Engine names must be non-empty and unique. Feature spaces remain separate and two
label-free cross-engine agreement features are added. `--ensemble` and `--join` are mutually
exclusive, and protein inference is unavailable in ensemble mode.

`--rt-features` adds two experimental retention-time residual features. Sequence-based retention is
aligned to `ScanNr` as a within-run elution proxy, and the label-dependent alignment is refitted in
each outer training partition. A PIN `retentiontime` metadata column is not used as the observed RT
for this model.

### Protein inference

Protein inference runs only when a protein output path is requested.

| Option | Default | Meaning |
|---|---:|---|
| `--protein-inference picked\|bayesian` | `picked` | Protein method; `fido` aliases `bayesian` |
| `--protein-alpha F` | `0.1` | Bayesian peptide emission probability |
| `--protein-beta F` | `0.01` | Bayesian noise probability |
| `--protein-gamma F` | `0.5` | Bayesian protein-presence prior |
| `--protein-peptide-prior F` | `0.1` | Prior used to convert peptide PEPs to likelihood ratios |
| `--protein-max-iter N` | `200` | Bayesian message-passing iteration limit |

Bayesian inference is exact for tree-structured connected components and uses deterministic damped
loopy belief propagation for cyclic components. Both protein methods are experimental and currently
fail available protein-truth calibration requirements.

### Profiling build

Build with the optional instrumentation:

```bash
cargo build --release --locked --features profiling
```

This enables `--profile-json PATH`, `--profile-cpu PATH`, and `--profile-allocations`. These flags
are rejected by a normal build. Reproduction tooling and interpretation guidance are in
[`bench/RUNTIME_PROFILE.md`](bench/RUNTIME_PROFILE.md).

## Scientific status

The full validation record is intentionally cumulative: failed implementations and negative results
remain in the repository rather than being rewritten after repairs. Start with
[`validation/README.md`](validation/README.md). The most recent general audit is
[`validation/FINAL_REPAIR_SCIENTIFIC_AUDIT.md`](validation/FINAL_REPAIR_SCIENTIFIC_AUDIT.md), and the
latest causal experiment is the
[`homology-depleted entrapment report`](validation/homology_depleted_entrapment/FINAL_REPORT.md).

### What has strong evidence

- Direct reported-list TDC+ q-value arithmetic matched an independent oracle in 40,928 exhaustive
  cases, plus 327,424 optimized count/mask comparisons. Exact ties, strict thresholds, the `+1`
  safeguard, null probabilities `0.2`, `1/3`, `0.5`, and `0.8`, and reverse cumulative minima were
  covered.
- Fixed-C, nested C-selection, and ensemble scoring passed held-out-label, held-out-outlier, row-order,
  and fold-isolation attacks. Normalization, initial direction, RT alignment, and model selection are
  fitted inside the relevant training partition.
- Protein grouping by identical, class-separated peptide evidence passed adversarial graph and
  insertion-order tests in isolation.
- The current behavior-preserving architecture refactor passes 138 release tests, six portable
  regression scripts, frozen adversarial probes, and byte-for-byte output comparison against its
  recorded baseline.

### Known implementation defects

These are present in the current behavior baseline and are not repaired by the architectural
refactor:

1. Exact PSM ties are sampled over emitted rows. Duplicating one scientific candidate therefore
   changes its label's win probability; a minimized null fixture changed from 0 to 101 false q<0.01
   discoveries after target-row duplication.
2. Joined source numbering depends on lexical filenames. Accessing the same bytes through renamed
   symlinks can change folds, tie decisions, and discoveries.
3. Protein mappings are unioned only after PSM competition. A losing duplicate occurrence can carry
   a complementary protein mapping that is discarded before the union.
4. Picked-protein target/decoy group pairing serializes members with an unescaped `|`, so distinct
   member sets can collide.
5. At the accepted edge value `--null-target-win-prob 1e-15`, the fixed `1e-12` per-target PEP floor
   breaks the claimed relationship between total PEP mass and estimated false count.
6. Decoy PEP display values can change under row permutations even when winners, scores, q-values,
   and target PEPs do not. This is presentation-only because no decoy posterior claim is made.

The compact CLI also ignores unknown options, as noted in the command reference.

### Calibration evidence

The predefined complete-null experiment observed no rejections in 30 runs at thresholds from 0.001
through 0.10. That is a conservative result, but 30 dependent runs cannot establish calibration at
small nominal FDRs and do not exercise the duplicate-candidate counterexample.

In the original signal-present entrapment study, the current method's mean adjusted FDP at reported
`q < 0.01` was **1.8104%**, with above-nominal FDP at all six tested thresholds. Pooled PSM PEPs were
optimistic in every populated bin, with weighted signed and absolute calibration error of
`+0.018685`; 217 known-false PSMs had PEP below 0.001.

A prospective homology-depletion experiment subsequently supplied causal evidence for one important
failure mode. At `q < 0.01`, the entrapment-target/decoy ratio changed from **1.940** in a freshly
searched original database to **1.000** after removing proteins capable of highly native-homologous
peptides. Three size-matched random controls remained at 1.543–2.095 (mean 1.828), and adjusted FDP
fell from **1.687%** to **0.880%**. The near-homology hypothesis is classified as **supported**, not
strongly supported: global PEP error remained positive (`+0.01464`), some uncertainty intervals
included no improvement, and only one dataset family was tested. The experiment therefore does not
justify a production filter or statistical correction.

Protein confidence is weaker. On held-out PrEST A and B truth sets, picked-protein `q <= 0.01` had
raw known-absent FDP of **45.92%** and **48.08%**; predefined count-adjusted FDP was **53.32%** and
**55.78%**. Default Bayesian probabilities were also severely miscalibrated. Do not use either
protein mode as calibrated evidence.

### Historical C++ compatibility evidence

No C++ output is treated as a correctness oracle. In a historical comparison against C++ Percolator
3.09, mean PSM-count differences at `q < 0.01` were small on single-candidate Tide and Sage PINs and
less concordant on multi-candidate MSFragger and yeast inputs:

| Dataset | Rust − C++ PSMs | Discovery Jaccard | Score Spearman |
|---|---:|---:|---:|
| PXD007145, Tide | +0.8 | 0.9930 | 0.9988 |
| PXD060954, Sage | +15.2 | 0.9962 | 0.9955 |
| PXD020243, MSFragger | −2.8 | 0.9235 | 0.9627 |
| Upstream yeast fixture | +15.4 | 0.9201 | 0.9498 |

The C++ `--post-processing-tdc` setting used in that study does not make concatenated multi-candidate
reporting identical to this project's default precursor competition. Treat the table as compatibility
evidence, not calibration or superiority evidence. Commands, seeds, and caveats are in
[`bench/MULTI_DATASET.md`](bench/MULTI_DATASET.md) and
[`validation/SECOND_REPAIR.md`](validation/SECOND_REPAIR.md).

## Performance

The current authoritative profile used 65 Comet PINs from PXD032157: 8,639,746 PSMs in 2.295 GB.
Measurements were made on an AMD Ryzen 5 5600G (6 cores / 12 threads), Ubuntu 26.04, Rust 1.97.0,
and the project's `x86-64-v3` release target.

| Workload | Runs | Median wall time |
|---|---:|---:|
| Largest PIN, `--num-threads 1` | 5 | 1.616816 s |
| Largest PIN, fold-parallel mode | 5 | 0.890897 s |
| All 65 files, sequential | 3 | 49.619487 s |
| All 65 files, four concurrent processes | 3 | 15.482359 s |

Every full-corpus configuration produced 106,823 target PSMs and 35,886 target peptides at strict
`q < 0.01`. These are reproducibility baselines for a development dataset, not sensitivity or
accuracy estimates. The dataset was used during model development, and file-level yields are highly
skewed.

`--num-threads` uses a private Rayon pool for nested/selected-C modes, but the fixed-C path only
switches between serial and parallel execution of the three folds. Values above one therefore do not
provide more than three-way fold concurrency for the canonical model. Parallel folds also retain
three design matrices at once; prefer the default one-thread mode when processing many files with
external process-level concurrency.

Fresh profiling attributes 40.51% of sequential process time to q-value/count/mask work and 29.24%
inclusively to initial-direction selection. The next justified optimization target is exact-order
reuse and q-value sorting/scanning; solver bookkeeping and output formatting are no longer material
hotspots. See [`bench/RUNTIME_PROFILE.md`](bench/RUNTIME_PROFILE.md) for hashes, acquisition tools,
stage tables, overhead measurements, and the optimization decision.

Historical Rust-versus-C++ throughput and memory measurements are preserved in
[`bench/RESULTS.md`](bench/RESULTS.md) and [`bench/REPRODUCTION.md`](bench/REPRODUCTION.md). They use
specific hosts, commands, and post-processing modes and should not be generalized as current
cross-platform performance claims.

## Architecture

The executable is a thin process boundary over reusable library modules:

```text
main.rs + cli.rs
  -> pipeline.rs
       -> competition.rs
       -> peptide.rs
       -> percolator.rs
            -> preprocessing.rs
            -> svm.rs / mlp.rs / simd.rs / stats.rs / rt.rs
       -> protein.rs / protein_bayes.rs
  -> output.rs
```

- [`src/cli.rs`](src/cli.rs) parses and validates process options.
- [`src/pipeline.rs`](src/pipeline.rs) composes rescoring, reported-list selection, peptide scoring,
  and protein dispatch.
- [`src/percolator.rs`](src/percolator.rs) owns folds, semi-supervised learning, model selection, and
  score merging.
- [`src/competition.rs`](src/competition.rs), [`src/peptide.rs`](src/peptide.rs), and
  [`src/protein.rs`](src/protein.rs) isolate higher-level inference policies.
- [`src/stats.rs`](src/stats.rs) implements TDC q-values and PEP construction.
- [`src/output.rs`](src/output.rs) owns stable TSV serialization.

The completed behavior-preserving refactor, frozen outputs, risk map, and acceptance evidence are in
[`refactor/README.md`](refactor/README.md), [`refactor/ARCHITECTURE.md`](refactor/ARCHITECTURE.md), and
[`refactor/RESULT.md`](refactor/RESULT.md).

## Testing and development

Run the standard local checks with:

```bash
cargo fmt --all -- --check
cargo clippy --all-targets --all-features --locked -- -D warnings
cargo test --release --all-targets --locked
RUSTDOCFLAGS="-D warnings" cargo doc --no-deps --all-features --locked
```

Six portable shell gates cover canonical regression, the MLP, nested selection, feature reports,
ensembles, and protein inference:

```bash
bash tests/regression.sh
bash tests/model_regression.sh
bash tests/selection_regression.sh
bash tests/feature_report.sh
bash tests/ensemble_regression.sh
bash tests/protein_regression.sh
```

The comprehensive non-benchmark acceptance command builds the release binary, runs all Rust and
portable regression tests, compares fixed-C, selected-C, ensemble, and protein TSVs byte-for-byte
with the frozen baseline, and reruns the recorded adversarial probes:

```bash
python3 refactor/verify_baseline.py
```

GitHub Actions runs the release build/tests, all six shell gates, and benchmark-tool unit tests.
The full 2.3 GB performance gate is manual because it requires a self-hosted runner with PXD032157
and the C++ reference binary.

## Documentation map

- [`validation/README.md`](validation/README.md) — ordered scientific audit and repair history.
- [`validation/FINAL_REPAIR_SCIENTIFIC_AUDIT.md`](validation/FINAL_REPAIR_SCIENTIFIC_AUDIT.md) —
  current general adversarial verdict and minimized failures.
- [`validation/homology_depleted_entrapment/FINAL_REPORT.md`](validation/homology_depleted_entrapment/FINAL_REPORT.md)
  — latest preregistered causal validation.
- [`bench/REPRODUCTION.md`](bench/REPRODUCTION.md) — benchmark commands and result provenance.
- [`bench/RUNTIME_PROFILE.md`](bench/RUNTIME_PROFILE.md) — latest authoritative runtime profile.
- [`bench/ADVANCED_FEATURES.md`](bench/ADVANCED_FEATURES.md) — join, RT, threading, and protein feature
  evaluations.
- [`refactor/README.md`](refactor/README.md) — behavior-preserving architecture record and verifier.

## License

Licensed under the [MIT License](LICENSE).

## PRIDE Archive working cache

`percolator-rs pride` discovers public PRIDE projects, inspects storage costs, downloads and verifies selected files, and runs the existing analysis on validated PINs. The default large-data cache ceiling is 50 GB; ephemeral processing and `pride cache prune --all-evictable` reclaim recoverable data while preserving manifests, provenance and results.

See [PRIDE usage and storage guarantees](docs/PRIDE.md) and the [real-project demonstration](docs/PRIDE-demonstration.md). Start with `percolator-rs pride --help`.
