# percolator-rs

[![CI](https://github.com/AndrejRumenovski/percolator-rs/actions/workflows/ci.yml/badge.svg)](https://github.com/AndrejRumenovski/percolator-rs/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

**Semi-supervised peptide identification in Rust, with reproducible validation and performance analysis.**

`percolator-rs` is an independent Rust implementation of the
[Percolator](https://github.com/percolator/percolator) peptide-spectrum match (PSM) rescoring workflow.
It reads candidate peptide identifications and numeric search-engine features from Percolator input
(`.pin`) files, learns a linear support vector machine (SVM), and produces ranked PSM and peptide
tables with target-decoy q-values and posterior error probability (PEP) estimates. Optional modules
support protein inference, multiple search engines, and data acquisition from the PRIDE Archive.

The project combines algorithm implementation, systems engineering, and empirical evaluation.
Its central questions are whether the workflow can be implemented efficiently and reproducibly,
whether training remains isolated from held-out data, and under which experimental conditions its
reported confidence estimates agree with observed errors. The methodological foundation is
[Käll et al. (2007)](https://www.nature.com/articles/nmeth1113).

## Project objectives

- **Implement the rescoring workflow:** semi-supervised linear-SVM training, three-fold
  cross-validation, precursor competition, and peptide-level reporting.
- **Make computations reproducible:** seeded fold assignment and tie handling, explicit input and
  output contracts, and regression comparisons against recorded outputs.
- **Evaluate statistical behavior:** independent q-value oracles, tests for information leakage,
  entrapment experiments, and protein-truth benchmarks.
- **Measure computational cost:** memory-mapped input, vectorized numerical kernels, and recorded
  runtime and memory measurements on public proteomics data.
- **Preserve experimental provenance:** dataset manifests, checksums, parameter records, and
  reproducible analysis tools.

The repository includes the implementation, evaluation scripts, benchmark reports, and documented
limitations. [Scientific status](#scientific-status) distinguishes computational correctness from
empirical calibration; [Performance](#performance) describes a recorded benchmark and its conditions.

**Contents:** [Quick start](#quick-start) · [Input and output](#input-and-output-contract) ·
[Method](#how-the-canonical-workflow-works) · [Command-line reference](#command-line-reference) ·
[Scientific status](#scientific-status) · [Performance](#performance) ·
[Architecture](#architecture) · [Testing](#testing-and-development) ·
[PRIDE Archive](#pride-archive-integration) · [Documentation](#documentation-map) ·
[References](#references)

## Quick start

Run the following commands from the repository root with a stable Rust toolchain installed.
The default build configuration targets `x86-64-v3` (AVX2/FMA-capable x86-64 processors, including
Intel Haswell and AMD Zen or newer):

```bash
cargo build --release --locked
./target/release/percolator-rs --help
```

Run the canonical linear-SVM workflow on the included fixture:

```bash
mkdir -p target/readme-example
./target/release/percolator-rs \
  --seed 1 \
  --results-psms target/readme-example/target.psms.tsv \
  --decoy-results-psms target/readme-example/decoy.psms.tsv \
  --results-peptides target/readme-example/target.peptides.tsv \
  --decoy-results-peptides target/readme-example/decoy.peptides.tsv \
  tests/fixtures/sample.pin
```

This example writes four tab-delimited result files under `target/readme-example/`. The fixture
provides a small installation and workflow check. To analyze another dataset, replace
`tests/fixtures/sample.pin` with a compatible concatenated target-decoy PIN and choose output paths.

Output files are optional and are written only when their corresponding flags are supplied.
Configuration, progress, timing, and identification counts at `q < 0.01` are printed to standard error.

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

The default PSM competition keeps one winner per `(source, ScanNr, ExpMass)` after rescoring.
Repeated occurrences of a tied `(label, modified core peptide)` receive one candidate draw.
Training still uses the supplied rows, so duplicated observations can affect the fitted model.
When `ExpMass` is absent, it is treated as zero.

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

Q-values and PEPs are written with enough precision to recover the computed floating-point values,
preserving small probabilities and threshold membership. Scores use six decimal places.
Counts reported by the executable use strict comparison (`q < 0.01`). The historical PrEST
protein-calibration study uses `q <= 0.01`, as stated in its reports.

## How the canonical workflow works

Targets represent candidate matches to the biological sequence database; decoys provide negative
examples for training and target-decoy competition (TDC). The canonical profile uses all available
training rows and ten semi-supervised iterations.

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
   reported. Exact score ties use a deterministic seed-dependent draw over distinct tied candidates.
7. **Recompute reported-list statistics.** The TDC+ false discovery rate (FDR) estimate uses
   `(D + 1) * p / (1 - p) / max(T, 1)`, evaluates exact-score tie groups together, and applies a
   reverse cumulative minimum to obtain q-values. `T` and `D` are cumulative target and decoy
   counts; `p` is `--null-target-win-prob`, defaulting to `0.5`.
8. **Estimate PEPs.** PEPs are derived from increments of the cumulative false-discovery estimate and
   made monotonic with the pool-adjacent-violators algorithm (PAVA). Their empirical calibration is
   evaluated separately in [Scientific status](#scientific-status).
9. **Collapse higher levels.** Peptide reporting keeps the best PSM per `(label, modified core
   peptide)`. Optional protein inference then groups proteins by identical observed peptide evidence
   and applies either picked target-decoy competition or a Fido-style Bayesian model.

The default run is deterministic for a fixed input identity, path/name, seed, options, thread mode,
and build. Serial and fold-parallel output is byte-identical in the regression suite.

## Command-line reference

Use `percolator-rs --help` for the complete option list and `percolator-rs --version` for the version.
Unknown options, missing values, and invalid numeric settings produce explicit errors. Output paths
must be distinct from input files and other requested outputs. The separate `pride` subcommand
provides its own `--help`.

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
| `--maxiter N` | profile value | Semi-supervised iterations |
| `--subset-max-train N`, `-N N` | `0` | Maximum training rows per fold; `0` uses all rows |
| `--num-threads N` | `1` | `1` runs folds serially; values above `1` enable fold/grid parallelism |
| `--null-target-win-prob P` | `0.5` | Null probability that an incorrect target beats its decoy; must be in `(0,1)` |
| `--psm-competition` | on | Keep one rescored precursor winner |
| `--no-psm-competition` | — | Report all candidates; resulting q-values are not claimed as FDR estimates |

For `k` equivalent decoys per target, the intended null setting is `1 / (1 + k)`. This parameter
assumes the candidate-generation process satisfies the corresponding competition model.

### Execution profiles

Explicit `--maxiter`, `--subset-max-train`, `--cpos`, `--cneg`, and C-selection flags override the
profile regardless of argument order.

| Profile | Training subset | Iterations | C selection | Intended use |
|---|---:|---:|---|---|
| `--fast` | 20,000 | 5 | off | Quick QA and development |
| `--balanced` | 40,000 | 10 | off | Reduced-cost exploratory runs |
| `--canonical` | all | 10 | off | Default workflow using all training rows |

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

`--auto-model` cannot be combined with `--select-c` or explicit `--cpos`/`--cneg`. The feature
report records mean out-of-fold raw coefficients, standardized effects, fold variability,
label correlation, selection frequency, and held-out permutation importance with fitted models
held fixed. See
[`bench/AUTOMATIC_SELECTION.md`](bench/AUTOMATIC_SELECTION.md) for the selection study.

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
loopy belief propagation for cyclic components. Both protein methods remain experimental;
their protein-truth calibration results are summarized in [Scientific status](#scientific-status).

### Profiling build

Build with the optional instrumentation:

```bash
cargo build --release --locked --features profiling
```

This enables `--profile-json PATH`, `--profile-cpu PATH`, and `--profile-allocations`. These flags
are rejected by a normal build. Reproduction tooling and interpretation guidance are in
[`bench/RUNTIME_PROFILE.md`](bench/RUNTIME_PROFILE.md).

## Scientific status

Evaluation covers numerical correctness, training isolation, reproducibility, and empirical
calibration. These are distinct claims: agreement with a q-value formula establishes correct
arithmetic, while calibrated confidence also depends on the candidate-generation process and
target-decoy assumptions.

The cumulative [validation record](validation/README.md) retains both successful tests and negative
results. The [general scientific audit](validation/FINAL_REPAIR_SCIENTIFIC_AUDIT.md) documents
implementation-level findings, and the
[homology-depleted entrapment report](validation/homology_depleted_entrapment/FINAL_REPORT.md)
examines a biological source of calibration error. Results below refer to the revisions and
experimental designs recorded in those reports.

### Verified computational properties

- Direct reported-list TDC+ q-value arithmetic matched an independent oracle in 40,928 exhaustive
  cases, plus 327,424 optimized count/mask comparisons. Exact ties, strict thresholds, the `+1`
  safeguard, null probabilities `0.2`, `1/3`, `0.5`, and `0.8`, and reverse cumulative minima were
  covered.
- Fixed-C, nested C-selection, and ensemble scoring passed held-out-label, held-out-outlier, row-order,
  and fold-isolation attacks. Normalization, initial direction, RT alignment, and model selection are
  fitted inside the relevant training partition.
- Protein grouping by identical, class-separated peptide evidence passed adversarial graph and
  insertion-order tests in isolation.
- The recorded architecture refactor passed the release test suite, five portable regression
  scripts, frozen adversarial probes, and byte-for-byte output comparison against its recorded
  baseline.

### Current implementation checks

The [readiness review](validation/READINESS_REVIEW.md) records repairs and verification against
the audit's minimized counterexamples. Regression tests now cover duplicate-invariant candidate
tie draws, preservation of protein mappings from losing occurrences of a reported peptide,
collision-free protein-group pairing, PEP mass at very small null probabilities, and decoy PEP
display under row permutations. An independent q-value oracle is part of the regular test suite.

CLI tests cover invalid arguments, incompatible joined feature layouts, conflicting output paths,
and write failures. Joined per-file counts use the same competed PSM statistics as the result tables.

### Remaining limitations

- Joined source numbering depends on lexical filenames. Accessing identical data through renamed
  files or symlinks can change folds and tie decisions; preserve source names when reproducing a run.
- Candidate deduplication applies to exact reporting ties. Repeated input observations can still
  change training weights and learned scores.
- Numerical correctness and reproducibility do not establish empirical calibration. PSM PEPs and
  protein-level confidence remain subject to the experimental limitations below; the available
  studies do not establish generalization to an untouched biological dataset.

### Empirical calibration

False discovery proportion (FDP) is the observed fraction of accepted identifications that are false
under an experiment's truth definition. The studies below compare observed or adjusted FDP with
the reported q-value threshold and assess PEPs against known-false identifications.

The predefined complete-null experiment observed no rejections in 30 runs at thresholds from 0.001
through 0.10. That is a conservative result, but 30 dependent runs cannot establish calibration at
small nominal FDRs and do not exercise the duplicate-candidate counterexample.

In the original signal-present entrapment study, the evaluated method's mean adjusted FDP at reported
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

On held-out PrEST A and B truth sets, picked-protein `q <= 0.01` had raw known-absent FDP of
**45.92%** and **48.08%**; predefined count-adjusted FDP was **53.32%** and
**55.78%**. Default Bayesian probabilities also showed substantial calibration error. These results
leave calibrated PSM PEPs and protein-level confidence as open research problems for this
implementation.

### Historical C++ compatibility evidence

Comparison with C++ Percolator assesses compatibility independently of the truth-based calibration
studies. In a recorded comparison against version 3.09, mean PSM-count differences at `q < 0.01`
were small on single-candidate Tide and Sage PINs, with lower discovery agreement on multi-candidate
MSFragger and yeast inputs:

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

The recorded **2026-09-03 baseline** used 65 Comet PINs from PXD032157: 8,639,746 PSMs in 2.295 GB.
Measurements were made on an AMD Ryzen 5 5600G (6 cores / 12 threads), Ubuntu 26.04, Rust 1.97.0,
and the project's `x86-64-v3` release target. The linked
[runtime report](bench/RUNTIME_PROFILE.md) records the measured revision and binary hashes;
these timings describe that baseline rather than every subsequent build.

| Workload | Runs | Median wall time |
|---|---:|---:|
| Largest PIN, `--num-threads 1` | 5 | 1.616816 s |
| Largest PIN, fold-parallel mode | 5 | 0.890897 s |
| All 65 files, sequential | 3 | 49.619487 s |
| All 65 files, four concurrent processes | 3 | 15.482359 s |

In that baseline, every full-corpus configuration produced 106,823 target PSMs and 35,886 target peptides at strict
`q < 0.01`. These are reproducibility baselines for a development dataset, not sensitivity or
accuracy estimates. The dataset was used during model development, and file-level yields are highly
skewed.

The [2026-09-21 readiness check](validation/READINESS_REVIEW.md) completed the same 65-file workload
with four concurrent processes in **15.69 s**, reporting **106,823 target PSMs** and **35,885 target
peptides** at `q < 0.01`. This is a single local gate run. The one-peptide change was traced to
the revised ordering of distinct exact-tie candidates; it is documented in the review.

`--num-threads` uses a private Rayon pool for nested/selected-C modes, but the fixed-C path only
switches between serial and parallel execution of the three folds. Values above one therefore do not
provide more than three-way fold concurrency for the canonical model. Parallel folds also retain
three design matrices at once; prefer the default one-thread mode when processing many files with
external process-level concurrency.

That profile attributed 40.51% of sequential process time to q-value/count/mask work and 29.24%
inclusively to initial-direction selection, motivating investigation of score-order reuse and
q-value sorting/scanning. See [`bench/RUNTIME_PROFILE.md`](bench/RUNTIME_PROFILE.md) for acquisition
tools, stage tables, and instrumentation overhead measurements.

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
            -> svm.rs / simd.rs / stats.rs / rt.rs
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

Five portable shell gates cover canonical regression, nested selection, feature reports, ensembles,
and protein inference:

```bash
bash tests/regression.sh
bash tests/selection_regression.sh
bash tests/feature_report.sh
bash tests/ensemble_regression.sh
bash tests/protein_regression.sh
```

Run the complete current check suite, including both ordinary and profiling-feature builds,
documentation checks, the five shell gates, and Python evaluation-tool tests:

```bash
bash scripts/check.sh
```

GitHub Actions runs the same script. With the 65 PXD032157 PINs available under `data/PXD032157/`,
`bash scripts/check.sh --full-benchmark` also runs the Rust performance gate. The separate manual
CI benchmark job additionally compares with the C++ reference binary.

`refactor/verify_baseline.py` preserves the historical refactor comparison. Its frozen outputs and
known-failure assertions predate the correctness repairs and probability-output precision changes;
it is an archival comparison, not the current acceptance gate.

## PRIDE Archive integration

`percolator-rs pride` discovers public PRIDE projects, inspects storage requirements, downloads and
verifies selected files, and runs the rescoring pipeline on validated PINs. The default large-data
cache ceiling is 50 GB. Ephemeral processing and `pride cache prune --all-evictable` reclaim
downloaded data while retaining manifests, provenance, and results.

```bash
./target/release/percolator-rs pride --help
```

See [PRIDE usage and storage behavior](docs/PRIDE.md) and the
[real-project demonstration](docs/PRIDE-demonstration.md) for a complete acquisition, analysis,
and cleanup example.

## Documentation map

- [`validation/README.md`](validation/README.md) — ordered scientific audit and repair history.
- [`validation/READINESS_REVIEW.md`](validation/READINESS_REVIEW.md) — current repairs, executed checks,
  and remaining calibration limitations.
- [`validation/FINAL_REPAIR_SCIENTIFIC_AUDIT.md`](validation/FINAL_REPAIR_SCIENTIFIC_AUDIT.md) —
  historical general adversarial verdict and minimized failures.
- [`validation/homology_depleted_entrapment/FINAL_REPORT.md`](validation/homology_depleted_entrapment/FINAL_REPORT.md)
  — latest preregistered causal validation.
- [`bench/REPRODUCTION.md`](bench/REPRODUCTION.md) — benchmark commands and result provenance.
- [`bench/RUNTIME_PROFILE.md`](bench/RUNTIME_PROFILE.md) — 2026-09-03 runtime baseline and profiling method.
- [`bench/ADVANCED_FEATURES.md`](bench/ADVANCED_FEATURES.md) — join, RT, threading, and protein feature
  evaluations.
- [`refactor/README.md`](refactor/README.md) — behavior-preserving architecture record and verifier.
- [`docs/PRIDE.md`](docs/PRIDE.md) — public dataset acquisition, provenance, and cache management.

## References

Käll, L., Canterbury, J. D., Weston, J., Noble, W. S., and MacCoss, M. J. (2007).
Semi-supervised learning for peptide identification from shotgun proteomics datasets.
*Nature Methods*, **4**, 923–925. [doi:10.1038/nmeth1113](https://www.nature.com/articles/nmeth1113).

The [upstream Percolator project](https://github.com/percolator/percolator) provides the reference
implementation of the original method. `percolator-rs` is an independent implementation with its
own statistical reporting choices and evaluation record. When describing experiments using this
repository, record the Git revision, input identities, random seed, build configuration, and
complete command line alongside the results.

## License

Licensed under the [MIT License](LICENSE).
