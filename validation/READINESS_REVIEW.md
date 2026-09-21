# Implementation and scientific readiness review

Date: 2026-09-21. Base revision: `a2e6061` plus the working-tree changes identified in
[the machine-readable results](readiness-review-results.json).

The reviewed implementation passes its software checks and the local public-data performance gate.
Several reproducible correctness and usability defects were repaired. Empirical confidence
calibration remains an open research problem: the real-data reruns continue to show excess false
discoveries at a nominal 1% threshold. This review supports demonstration and further evaluation
of the implementation; it does not establish calibrated confidence across biological datasets.

## Repairs

| Area | Previous behavior | Reviewed behavior |
|---|---|---|
| CLI | Unknown flags and several malformed numeric values silently selected defaults; `--help` failed | Help and version commands succeed; unknown flags, missing values, non-finite class weights, and invalid thread counts fail clearly |
| Output handling | Buffered write failures could be ignored in ordinary builds; other output errors aborted the process | Every writer explicitly flushes and propagates errors; the CLI reports the affected path and exits unsuccessfully |
| Output identities | Inputs or another result table could be overwritten through a shared output path | Conflicting paths, including symlink aliases and Unix hard links, are rejected before analysis |
| Joined features | Equal-width feature blocks with different names or order could be pooled | Joined inputs require the same feature names in the same order |
| Joined counts | Per-file summaries used pre-competition statistics | Counts use the competed, reported PSM list and its recomputed q-values |
| Exact candidate ties | Repeating one candidate gave it additional chances to win | The draw uses distinct `(label, modified core peptide)` candidates within a precursor |
| Protein evidence | A losing PSM occurrence could discard a complementary protein mapping | A reported peptide retains the union of its mappings across all input occurrences of the same label and core sequence |
| Protein pairing | Unescaped `|` separators could conflate distinct protein-member sets | Structured member lists identify groups; tie seeds encode each member separately |
| Tiny PEPs | A fixed `1e-12` floor could add artificial false-discovery mass | The lower bound is the smallest positive representable `f64`; ordinary floating-point precision limits still apply |
| Decoy PEP display | Mixed-label exact ties could display different values after row permutation | Complete score-tie groups are processed together |
| Probability output | Six decimal places could round a q-value across a cutoff or turn a small PEP into zero | Q-values and PEPs round-trip through their text representation; scores retain six decimal places |
| Protein evaluation | The calibration reader rejected the picked method's `NA` posterior field | Missing posteriors remain missing; ranking and threshold metrics work without inventing a protein posterior |

These changes intentionally alter some historical outputs. In particular, candidate ordering is
defined by modified core sequence, protein-pair tie seeds use an unambiguous representation, and
probabilities have greater output precision. The frozen refactor artifacts retain the old behavior
and known counterexamples; they are not overwritten.

## Executed software checks

The common local and CI entry point is:

```bash
bash scripts/check.sh
```

Verification covered:

- 212 Rust tests in each of the ordinary and `profiling` configurations, including the independent
  oracle promoted into `tests/qvalue_oracle.rs`; the optional live-network test is excluded from
  these offline totals.
- 40,928 independent q-value oracle cases, 327,424 count/mask comparisons, and six exact threshold
  boundary checks in each configuration.
- Five portable regression scripts covering the canonical workflow, nested selection, feature
  reports, ensembles, and protein inference.
- Ten Python tests covering Sage normalization, protein truth metrics and schema validation,
  and PSM agreement calculations.
- Formatting, all-target/all-feature Clippy with warnings denied, and Rust documentation with
  warnings denied.
- The separately invoked live PRIDE metadata smoke test passed against the public archive.
  The ordinary CI suite remains independent of network availability.

The existing adversarial driver was also rerun against a frozen normal release binary. It found:

- One identical complete result across single-file row permutations, including decoy PEP fields.
- Invariance to joined file/row permutations when filenames are preserved.
- Zero additional false discoveries from the minimized 94-fold target-row duplication attack;
  the historical implementation created 101.
- No held-out-label or held-out-outlier score leakage in fixed-C, selected-C, or ensemble tests.
- Both complementary protein mappings retained at each tested seed, with insertion-order invariance.

The historical driver's `seed_changes_protein_mapping` field compares complete protein results,
including q-values and seed-dependent picked competition. It remains true even though the two
specific missing mappings are now both retained. The focused mapping checks and regression test
assess that repair directly.

## Public-data performance check

The local gate processed all 65 PXD032157 Comet PINs, containing 8,639,746 PSMs, with four concurrent
processes. One warm run completed in **15.69 seconds**, with sampled aggregate peak RSS of
**606,732 kB**, below the existing 45-second and 1.5-GiB limits. This is a local gate measurement,
not a speed comparison between versions.

At strict `q < 0.01`, the run reported **106,823 target PSMs** and **35,885 target peptides**.
The PSM total matches the prior baseline; the peptide total is one lower. A per-file comparison
against the hash-verified archived September 11 binary localized that difference to
`28May2015-QE-HF-Anopheles-38-MAGs-P-3rd-02-comet.pin`. Two exact target-candidate ties select different
peptides under the core-sequence ordering. The high-confidence alternative was already represented
elsewhere, reducing the distinct accepted-peptide count by one. This is an explained reporting
change, not evidence of improved or reduced biological accuracy.

Reproduce the full software and local performance checks with:

```bash
bash scripts/check.sh --full-benchmark
```

The 65 input files must already be available under `data/PXD032157/`.

## Empirical accuracy reruns

The existing experimental design, thresholds, five seeds, and selected protein parameters were
retained. The reviewed normal release binary was frozen and identified by SHA-256 in the results
JSON. The rerun covered six entrapment inputs, four compact search-engine datasets with repeated
runs, and twelve PrEST samples evaluated with three protein methods. All repeated compact-dataset
outputs were byte-identical, and the repaired protein evaluator accepted the current `NA` schema.

| Evaluation | Threshold | Observed result |
|---|---|---:|
| Signal-present entrapment, five-seed mean adjusted FDP | `q < 0.01` | **1.8104%** |
| Pooled PSM PEP calibration, observed minus predicted error | All populated bins | **+0.018678** |
| Known-false PSMs assigned PEP below 0.001 | `PEP < 0.001` | **217** |
| Held-out PrEST A, picked proteins: raw / adjusted FDP | `q <= 0.01` | **45.92% / 53.32%** |
| Held-out PrEST B, picked proteins: raw / adjusted FDP | `q <= 0.01` | **48.08% / 55.78%** |

Default Bayesian protein probabilities also remain poorly calibrated: the held-out split has
Brier score **0.22444** and ten-bin expected calibration error **0.39445**. Previously selected
parameters transfer better on this benchmark, but this does not establish broader generalization.
Complete threshold summaries and per-sample protein counts are in
[the results JSON](readiness-review-results.json).

The implementation repairs therefore do not resolve the biological exchangeability and confidence
calibration problems documented by the earlier studies. No correction was fitted to these rerun
outcomes, and no production homology filter was introduced.

## Remaining scope

- Joined source identifiers still depend on lexical filenames. The alias experiment changed 332
  PSM records when the same bytes were accessed under different names. Preserve source names and
  paths in experiment records. Replacing this identity model requires a separate design that also
  preserves fold isolation and distinguishes genuinely different runs.
- Tie deduplication does not make model training invariant to repeated observations. Training uses
  the rows supplied by the input producer.
- The available biological datasets have already been used during development. An untouched
  evaluation set is still needed for claims of biological generalization.
- These checks exercise the documented workflows and targeted failure cases; they do not prove
  correctness for every possible input or establish empirical FDR control from arithmetic alone.

## Evidence and reproduction

[readiness-review-results.json](readiness-review-results.json) records the base revision, source
hashes, normal-release binary hash, statistical results, and hashes of the detailed local evidence.
The reviewed working tree includes numerical optimizations that predated this review; their
presence is recorded rather than attributed to the repairs above.

Detailed command logs and output tables are retained locally under `target/readiness-review/`.
The existing adversarial and empirical drivers can be rerun without editing their historical
results:

```bash
python3 validation/final_repair_adversarial.py \
  --binary target/release/percolator-rs \
  --work-dir "$(pwd)/target/new-adversarial-review" \
  --json "$(pwd)/target/new-adversarial-review.json"
python3 validation/final_repair_empirical.py \
  --binary target/release/percolator-rs \
  --output target/new-empirical-review
```

Use a new output directory for each run. The empirical driver requires the locally prepared data
at the locations declared in that script; it does not download or regenerate search results.
