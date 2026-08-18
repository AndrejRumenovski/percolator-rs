# Small-MLP rescoring benchmark

## Outcome

The optional neural scorer does **not** improve identification yield on the current benchmarks. It
is retained as an experimental opt-in model; the linear SVM remains the default. Peptide-sequence
embeddings are deferred because the tabular MLP has not established that extra model complexity is
useful.

At reported q<0.01 on the 65-file PXD032157 benchmark, the MLP reports 105,521 PSMs and 36,447
peptides versus 107,046 and 37,469 for the SVM: -1.42% PSMs and -2.73% peptides. The MLP takes
120.8 seconds at four-file concurrency versus 27.1 seconds for the SVM (4.46x slower). It wins on
21 files, loses on 42, and ties on two at PSM level; peptide counts improve on 22, decline on 40,
and tie on three.

The four independent extension cases also do not show an aggregate gain:

| case | SVM PSM / peptide | MLP PSM / peptide | MLP delta |
|---|---:|---:|---:|
| PXD020243, MSFragger | 1,554 / 1,177 | 1,525 / 1,179 | -29 / +2 |
| PXD060954, Sage | 26,614 / 11,433 | 26,642 / 11,408 | +28 / -25 |
| Hogrebe, Tide | 29,264 / 20,614 | 29,264 / 20,562 | 0 / -52 |
| Percolator yeast | 1,126 / 903 | 1,079 / 871 | -47 / -32 |
| **aggregate** | **58,558 / 34,127** | **58,510 / 34,020** | **-48 / -107** |

The six-run signal-present entrapment check reaches the same practical conclusion. At reported
q<=0.01 the MLP accepts 19,149 PSMs at an adjusted entrapment FDP of 2.45% (95% CI 2.21-2.71%),
while the SVM accepts 19,666 at 2.78% (2.53-3.06%). The MLP is slightly less anti-conservative but
finds 517 fewer PSMs, and neither model validates nominal 1% FDR.

Machine-readable results are in [`model-comparison-results.tsv`](model-comparison-results.tsv).

## Model and fair-comparison design

`--rescore-model mlp` replaces only the fold-local learner. Both models share:

- feature parsing and global z-score normalization;
- the seeded three-fold assignment;
- the best-single-feature initial direction;
- iterative selection of targets below the training-fold q-value cutoff and all decoys as
  negatives;
- the same optional training subset, class-weight grid, and iteration count;
- out-of-fold scoring, target-decoy q-values, isotonic PEPs, and peptide rollup.

The network is deliberately small: one eight-unit tanh hidden layer plus a trainable linear skip
connection. The skip starts at the same best-feature direction as the SVM, while hidden output
weights start at zero, so the initial neural score exactly equals the linear initialization. Each
semi-supervised iteration performs ten deterministic mini-batch Adam epochs at learning rate 0.02.
The fixed architecture contains only a few hundred parameters for normal PIN feature counts and has
no external runtime dependency.

Defaults were frozen using only the committed 12,000-PSM development fixture. That fixture is not
included in the evaluation totals above. Seeds 1-5 all improved its PSM count relative to the SVM,
which made it a useful optimizer-stability check but did not predict generalization.

## Usage

```bash
cargo build --release
target/release/percolator-rs --canonical --seed 1 \
  --rescore-model mlp --results-psms target.psms.tsv input.pin
```

The SVM remains the default. Neural parameters can be overridden with `--mlp-hidden`,
`--mlp-epochs`, `--mlp-learning-rate`, and `--mlp-l2`. `--select-c`, `--cpos`, and `--cneg` apply
to either learner.

Reproduce the 65-file or alternate-directory yield comparison with:

```bash
bash bench/model_comparison.sh
MODEL_BENCH_INPUT=/path/to/pins MODEL_BENCH_OUT=/path/to/output \
  bash bench/model_comparison.sh
```

After generating the entrapment searches with `bench/entrapment/run.sh`, reproduce the calibration
comparison with:

```bash
bash bench/model_entrapment.sh
```

## Decision on sequence embeddings

Sequence embeddings would add a second information source, so this negative result does not prove
they cannot help. They also introduce substantially more parameters, preprocessing, and validation
risk. The next neural experiment should therefore start only with a predeclared evaluation plan and
an independent development split, and should compare a sequence-aware model with both this MLP and
the SVM under the same out-of-fold/FDR path. The tabular MLP result alone does not justify making
that investment or changing the default scorer.
