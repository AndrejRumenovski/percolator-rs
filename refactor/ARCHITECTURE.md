# Architecture and risk audit

## Frozen scope

This audit describes commit `e8d83d1c76e4cf651fdfcf22d98b0b499c35943a`
before production refactoring.  The checkout also contains pre-existing,
uncommitted validation reports and probes; `git status` in the baseline artifact
is authoritative.  There were no pre-existing changes under `src/`, `tests/`,
`Cargo.toml`, or `Cargo.lock`.

This is a behavior-preserving architecture change.  It intentionally does not
address the scientific and validation limitations in
`validation/FINAL_REPAIR_SCIENTIFIC_AUDIT.md`, including candidate-multiplicity
bias, path-alias sensitivity, the picked-protein key collision, and PEP
calibration.  Those observations are frozen so an architectural extraction
cannot accidentally disguise a methodology change as cleanup.

## Dependency map before refactoring

```text
CLI parsing and profiles (main.rs)
  -> PIN parse / join / ensemble merge (pin.rs)
  -> optional RT feature reservation (rt.rs)
  -> semi-supervised orchestration (percolator.rs)
       -> fold-local RT residuals (rt.rs)
       -> normalization and design matrices (percolator.rs)
       -> spectrum-grouped folds and nested splits (percolator.rs)
       -> initial direction and positive selection (percolator.rs -> stats.rs)
       -> linear SVM or MLP (svm.rs / mlp.rs -> simd.rs)
       -> fold score calibration and held-out score merge (percolator.rs)
       -> reported q-values and PEPs (stats.rs)
  -> PSM competition / ensemble deduplication (main.rs -> tiebreak.rs)
  -> reported-list q-values and PEPs (main.rs -> stats.rs)
  -> peptide best-PSM selection (main.rs -> stats.rs)
  -> complete peptide-to-protein mapping union (main.rs)
  -> picked or Bayesian protein inference (protein.rs / protein_bayes.rs)
  -> PSM, peptide, protein, and feature-report serialization (main.rs)
  -> benchmark and validation hooks (profile.rs plus CLI subprocess drivers)
```

`src/lib.rs` exposes only benchmark manifest/result/comparison support.  The
scientific graph above is compiled as private modules of `src/main.rs`.
Consequently, integration tests cannot call the statistical pipeline directly,
and standalone audit probes use `#[path = "../src/..."]` source inclusion.

## Responsibility and coupling inventory

| Area | Current responsibility, inputs, outputs, and mutable state | Risk and coupling |
|---|---|---|
| `main.rs::parse_args` and `Args` | Reads process argv, applies profiles and overrides, validates combinations, mutates `percolator::Params`, and exits on errors. | Medium. Policy is mixed with process control; malformed option values can become defaults/NaN before validation. No statistical arithmetic, but flags select scientific paths. |
| `pin::Dataset` | Structure-of-arrays PSM representation: row-major features, raw `i8` labels, strings, source ids/names, and ensemble marker. | Scientifically and performance critical. Row positions index every parallel vector. `source` controls folds and competition identity; strings control peptide/protein identity. Layout is tuned for hot loops. |
| `pin::parse` | Memory-mapped TSV parsing, header/feature discovery, numeric validation, allocations, and source metadata. Produces `Dataset`; local buffers and mapped input are mutable/read-only respectively. | High and performance critical. Parsing defines accepted scientific input and silently changing feature boundaries changes models. Measured at about 12% of the historical full workload. |
| `pin::canonicalize_rows`, `compare_parts`, `merge` | Canonicalizes each joined source, orders named parts, assigns numeric sources, then concatenates all columns. | Scientifically critical. Controls accumulation order, folds, tie keys, and joined permutation invariance. File names remain part of identity by validated current behavior. |
| `pin::merge_ensemble` | Builds namespaced engine feature blocks, exact agreement features, and combined identity columns. | Scientifically critical/high. Label-free agreement prevents leakage; engine/source semantics differ from joined input. Must not be generalized with `merge`. |
| `rt::{augment, Alignment::residuals}` | Reserves two feature columns, derives sequence/scan series, and fits source-specific target-label-dependent alignments inside supplied rows. | Scientifically critical. Placement inside each fit partition is a fold-isolation invariant. `augment` also changes the feature matrix layout. |
| `percolator` parameter/config types | Own model, CV, training, selection, null-probability, RT, and thread settings. | High. CLI policy and learning configuration are tightly coupled in one broad mutable struct. |
| `percolator::{fit_normalization, transform_matrix, build_matrix_fit}` | Fits mean/variance on explicit rows, materializes all transformed rows plus bias, optionally substitutes fold RT values. | Scientifically critical and performance sensitive. Floating-point traversal order and fit-row isolation are observable. |
| `percolator::{assign_dataset_folds, outer_fold_assignments, inner_splits}` | Groups spectrum identities, deterministically balances them into folds, and constructs nested partitions. | Scientifically critical. Row/source identity and seed determine leakage isolation. |
| `percolator::initial_direction` | Tests both orientations of every feature using the training TDC heuristic and chooses the first best yield. | Scientifically and performance critical. It performs many q-value scans; tie/grid order and floating-point values affect the complete training trajectory. |
| `percolator::{FoldModel, train_fold}` | Iteratively scores training rows, selects confident target positives, subsamples, packs rows, and updates SVM/MLP state. Mutable model and reusable workspaces live here. | Scientifically/performance critical. SVM objective evaluation is the largest measured hot path. Do not rewrite loops or allocation order for style. |
| `percolator::{training_null_calibration, standardized_heldout_scores}` | Fits a training-decoy location/scale and maps held-out scores to a fold-comparable scale. | Scientifically critical. Using held-out rows or altering accumulation changes cross-fold ranking. |
| `percolator::{cv_scores, select_c_for_fold, cv_scores_with_selected_c, nested_cv_scores, run}` | Builds outer folds, performs fixed-C or nested selection, runs fold models serially or with Rayon, merges held-out scores, and computes full-list statistics. | Scientifically critical/high. Several modes share orchestration but have intentionally distinct selection semantics. `run` is the main scientific entry point but is private to the binary. |
| `svm::{Problem, Workspace, train}` and `simd` | Squared-hinge objective, active set, gradient/Hessian, Cholesky solve, line search, and specialized dot/AXPY kernels. | Performance critical and scientifically sensitive. Accumulation order affects margins and convergence. `f_and_active` was about 25% of the historical workload. Intentionally left structurally intact. |
| `mlp::Network` | Alternative seeded nonlinear fold learner with its own parameters and mutable optimizer state. | High but not part of the canonical validated profile. It shares CV/statistics contracts with SVM. |
| `stats::Tdc` | Represents reported versus training estimator configuration using public scalar fields and a boolean. | Scientifically critical. The reported/training distinction is explicit via constructors but invalid scalar combinations remain representable inside the crate. |
| `stats::{qvalues*, target_*}` | Deterministic descending score ordering, exact numeric tie groups, TDC+ scan, reverse cumulative minimum, and reusable fast count/mask paths. | Scientifically and performance critical. Sorting was about 15% and all q-value work about 22% historically. Fast paths must remain equivalent to materialization. |
| `stats::{peps, peps_from_competition_into, pava_non_decreasing}` | Derives target PEP increments from competition counts and applies isotonic pooling. | Scientifically critical. Tie grouping, mass transfer, bounds/floor, and target/decoy presentation are observable current behavior. |
| `main.rs::competition_winners` | Groups rows by `(source, scan, mass bits)`, finds exact best-score ties, canonically orders tied rows, and performs seeded draws. | Scientifically critical. It is hidden in CLI code even though it defines which hypotheses reach reported statistics. Candidate multiplicity and path behavior are frozen known limitations. |
| ensemble output dedup block in `main` | Keeps one best row for exact `(scan, label, peptide)` candidates when competition is disabled. | High. Similar shape to PSM competition but different scientific semantics; must remain separate. |
| peptide block in `main` | Defines core peptide identity, keeps the first best-scoring representative by class, preserves input index order, and recalculates peptide q/PEP. | Scientifically critical/performance sensitive. Representative identity, exact `>=` tie behavior, class separation, and order affect output and protein inference. |
| `main.rs::protein_entries` | Unions protein associations over all reported occurrences of `(label, core peptide)` while retaining representative peptide scores/PEPs. | Scientifically critical. This is the repaired `M(p)` invariant, but it sits between CLI row materialization and inference rather than in an inference boundary. |
| `protein::{infer, picked_fdr}` | Constructs protein evidence sets, class-aware equivalence groups, best-peptide group scores, target/decoy picking, and q-values. Mutates `picked/qval` fields. | Scientifically critical. Hash iteration is explicitly neutralized before tie-sensitive operations. Pairing serialization has a frozen known collision; do not repair here. |
| `protein_bayes::infer` | Constructs a peptide/protein factor graph, exact or loopy inference, posterior/q-value assignment, and diagnostics. | Scientifically critical/high but not the canonical default. Its entry type is an unlabelled tuple shared with picked inference. |
| `main.rs::{Row, write_results, write_proteins, write_feature_report}` | Borrows/materializes output rows, preserves a 96-byte sort-layout compatibility constraint, sorts, formats fixed precision, buffers writes, and profiles serialization. | Medium to high and performance critical. Formatting is about 11% historically; output bytes and equal-score ordering are acceptance criteria. |
| `profile.rs` | Feature-gated allocator/event/CPU instrumentation and JSON output. | Low scientific risk, medium architectural coupling because hot modules contain many conditional hooks. |
| benchmark modules/binaries | Dataset registry, process orchestration, result schema, and comparisons. | Low scientific risk. Already correctly exposed from the library and largely separate from rescoring. |
| `tests/` and `validation/` | Unit tests nested in binary modules, CLI integration attacks, shell gates, standalone source-inclusion probes, and historical research artifacts. | Medium. Provenance is strong but organization exposes the missing library boundary; older probes can intentionally fail after later repairs. |

## Duplicated or mixed logic

- The PSM, peptide, and protein paths each build score/label vectors and call the
  same statistics API, but their identity and competition rules differ.  Only
  the mechanical statistic invocation is common; the inference rules must not
  be collapsed into a generic “level” abstraction.
- Core peptide parsing is duplicated conceptually in `main.rs` and RT sequence
  parsing.  Their purposes differ: inference identity retains modifications,
  while RT prediction strips non-amino-acid syntax.  They should be separately
  named rather than mechanically unified.
- Target/decoy is repeatedly expressed as `i8 > 0` / `< 0`, while protein class
  is independently inferred from accession prefixes.  A boundary enum can make
  decisions explicit without replacing compact label arrays in hot loops.
- Fold setup, selected-C setup, and nested model selection repeat partition,
  seed, train, score, and merge mechanics.  Consolidation is safe only after
  explicit training/held-out partition types exist; the selection topology must
  remain visible.
- Output creation and error exits are interleaved with science in `main`, making
  it difficult to test materialization without filesystem/process behavior.

## Implicit invariants to make explicit

1. `Dataset` parallel vectors have exactly `n_psm` rows and `features` has
   `n_psm * n_feat` values.
2. PIN labels are target or decoy only; learning hot paths currently infer class
   from sign.
3. Source id plus scan plus normalized mass bits defines a precursor, except an
   ensemble deliberately erases source from that identity.
4. Fold fitting may read only training rows; outer held-out rows cannot affect
   normalization, RT alignment, initial direction, model selection, training,
   or fold score calibration.
5. Reported q-values/PEPs are recalculated after any competition or
   deduplication and use the reported TDC constructor.
6. Exact numeric equality defines statistical score ties; total ordering exists
   only to make sorting deterministic.
7. Peptide identity is `(target/decoy class, modified core peptide)` and best-PSM
   ties retain the earlier canonical row under the current `>=` rule.
8. For peptide `p`, protein mapping is the union over every reported occurrence;
   no representative row owns the mapping.
9. Protein groups are equal only when both their evidence sets and their
   target/decoy classes are equal.  Strict subsets remain distinct.
10. Picked-protein inference has no protein posterior; its PEP is unavailable
    and serializes as `NA`.
11. Seeded tie behavior is content/identity driven within the currently
    validated identity contract, and serial/parallel execution is byte-identical.
12. Output precision, headers, row ordering, and the row-layout-dependent
    unstable-sort behavior are frozen compatibility constraints.

## Risk classification

| Classification | Code |
|---|---|
| Low risk | benchmark schemas/comparison, help text, module declarations, validation harness documentation |
| Medium risk | CLI parser extraction, error plumbing, feature-report serialization, profiling ownership |
| High risk | PIN parsing/merge, output materialization/order, peptide identity, Bayesian inference, shared configuration changes |
| Scientifically critical | competition/ties, TDC/q-values/PEPs, folds/CV, normalization, RT fit, initial direction, model selection, peptide mapping, protein grouping and picking |
| Performance critical | PIN row parser; normalization matrix construction; SVM active-set/objective/Hessian loops; score sorting and q scans; initial-direction scans; result formatting; peptide dedup; protein inference on enabled runs |

## Target architecture

The target uses cohesive flat modules rather than a wholesale directory rewrite:

```text
lib.rs
  domain.rs             small boundary types and checked identities only
  pin.rs                PIN parsing plus joined/ensemble input construction
  preprocessing.rs      fold-local normalization and RT feature materialization
  learning/
    mod.rs              public run/config/result boundary
    folds.rs            deterministic partitions and isolated fold setup
    selection.rs        legacy selected-C and nested selection topology
    train.rs            existing semi-supervised learner orchestration
  svm.rs, mlp.rs, simd.rs
  statistics/
    mod.rs              TDC configuration, q-values, PEPs
    competition.rs      PSM competition only
    tiebreak.rs         seeded identity draw primitive
  peptide.rs            peptide identity, representative selection, statistics,
                        and complete protein-association union
  protein.rs            staged picked inference
  protein_bayes.rs      separate probabilistic inference
  output.rs             output rows and byte-compatible serialization
  pipeline.rs           science-level PSM -> peptide -> protein composition
  benchmark_*           existing benchmark support

main.rs
  cli.rs                argv/profile parsing and validation
  load -> pipeline -> write -> diagnostics only
```

This is a direction, not permission for a big-bang rewrite.  Existing files may
remain flat when moving them adds no clarity.  In particular, the optimized SVM,
SIMD, parser, q-value/PEP loops, and Bayesian message passing stay intact unless
a later boundary extraction can be proven byte-identical and benchmark-neutral.

Useful explicit types are deliberately limited:

- `TargetDecoy` at parsing, competition, peptide, protein, and output boundaries;
  keep `Vec<i8>` internally for the validated hot layout.
- `FoldId`/partition views around fold construction, if they prevent mixing
  training and held-out rows without changing stored `u8` assignments.
- a peptide identity key containing class plus core sequence, owned only where
  lifetime or representative-row coupling otherwise leaks through.
- `Option<Pep>` or an equivalent output-facing availability type for protein
  posterior values; do not wrap every `f64` score/q-value in hot code.

## Incremental slices

1. Expose existing scientific modules through the library and make the binary a
   consumer.  No functions move and no loops change.  This removes source-text
   inclusion as a prerequisite for independent tests.
2. Extract CLI configuration from `main.rs`.  Preserve argv behavior and process
   exit codes with characterization tests.
3. Extract byte-compatible output serialization.  Retain the 96-byte row layout,
   unstable sort, formatting implementation, and profiling counters exactly.
4. Extract PSM competition into a statistics-adjacent module.  Move existing
   tests with it and add public pure-input tests; do not change tie identity.
5. Extract peptide inference/materialization and the protein-entry union into
   `peptide.rs`.  Freeze representative and mapping behavior with direct tests.
6. Reduce `main` to orchestration and move reusable science composition into a
   pipeline boundary.
7. Only then evaluate a mechanical split of preprocessing/folds from
   `percolator.rs`.  This slice touches critical floating-point and CV code, so
   it requires all CV attacks, canonical byte comparisons, and current
   benchmarks.  Stop if the split makes partition flow less obvious.
8. Organize validation entry points and fixture helpers without rewriting or
   deleting historical reports.

## Areas intentionally not refactored

- SVM/MLP mathematics, iteration order, line search, data packing, and SIMD.
- q-value/PEP mathematics, score sorting, tie grouping, PAVA, and numeric floors.
- PIN storage layout and parser hot loop.
- joined-input filename identity, candidate multiplicity behavior, and current
  competition identity semantics.
- protein pairing-key semantics or inference methodology.
- Bayesian factor graph/message schedule.
- research claims, historical expected results, and historical probes.

These exclusions keep architecture work from becoming an unrequested scientific
repair or performance project.
