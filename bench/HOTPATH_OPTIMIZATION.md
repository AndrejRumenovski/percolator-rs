# Exact q-value and initial-direction optimization

Starting revision: `03c754f41911dfb717eb6e00c7505f6a01aa7af5`.
Historical reference: 49.619 s sequential, 15.482 s N=4. Acceptance uses
alternating fresh before/after runs on this host, not historical timings alone.
Artifacts: `/tmp/percolator-opt-FoGSSG`.

## Result

Four exact optimizations are retained; three unsuccessful experiments are
reverted. The final N=4 median is **13.399 s**. The **<12 s N=4 stretch target
was not reached**.

| Workload | Fresh before (s) | Final (s) | Improvement |
|---|---:|---:|---:|
| Largest PIN, one thread (5 pairs) | 1.617 | 1.354 | 16.28% |
| Full corpus, sequential (3 pairs) | 50.149 | 41.184 | 17.88% |
| Full corpus, N=2 (3 pairs) | 27.480 | 22.877 | 16.75% |
| Full corpus, N=4 (3 pairs) | 15.723 | 13.399 | 14.78% |
| Full corpus, N=6 (3 pairs) | 12.160 | 10.630 | 12.58% |

The supplied historical reference was 49.619 s sequential and 15.482 s N=4.
Relative to it, the final times improve by 17.00% and 13.45%, respectively.
Fresh paired comparisons above establish that the improvement persists under
the current host load. Every final pair improves, at every tested concurrency.

Sequential pairs: 50.149 → 41.325, 50.150 → 41.184, 50.070 → 41.162 s.
N=4 pairs: 15.781 → 13.457, 15.686 → 13.376, 15.723 → 13.399 s.
All 260 corpus TSV size/hash records agree on every pass, including the rejected
experiments; direct frozen-fixture byte comparisons also agree. Yield remains
**106,823 target PSMs and 35,886 target peptides at strict q < 0.01**.

Final median peak RSS is essentially unchanged: 213.877 → 213.848 MB sequential
and 820.584 → 821.010 MB at N=4 (decimal MB). Final runs have zero major faults;
swap remains unused. N=4 scaling is 3.07× after optimization versus 3.19× before,
consistent with the remaining scoring/training/normalization costs limiting
scaling more than the reduced q-value work. N=2/N=6 are diagnostics; the primary
acceptance metric remains N=4.

The compact evidence is in [`hotpath-optimization-results.json`](hotpath-optimization-results.json).
The proofs, individual acceptance decisions, and reproduction command follow.

### Final profile

Percent of summed in-process elapsed time in the separate stage-timed build;
these inclusive rows overlap and must not be added together.

| Operation | Fresh reference sequential | Final sequential | Final N=4 |
|---|---:|---:|---:|
| All q-value/count/mask work | 40.51% | 27.15% | 24.17% |
| Initial direction | 29.24% | 18.82% | 18.98% |
| Q-value score sorting | 24.37% | 18.04% | 15.96% |
| Reverse-rank sorting | 5.14% | 0% | 0% |
| Tie/count scans | 10.14% | 8.14% | 7.06% |

The final sequential profile records 10.491 s of q-value work, including
6.971 s of score sorting and 3.148 s of tie/count scans. Initial direction takes
7.271 s. There are 6,085 q-value/count/mask passes, down from 10,335; 4,095 reverse
passes have been fused and 155 unchanged-model passes reused. All 1,950 training
iterations remain. The final N=4/serial aggregate stage-time ratios are 1.10 for
sorting, 1.08 for scans, 1.24 for normalization, 1.46 for model scoring, and 1.33
for SVM training. These are timing ratios, not measured cache-miss ratios.

The stage-timed build's full wall times were 40.254 s sequential and 13.134 s
N=4, about 2% below the normal-build medians. Profiles therefore describe that
build's hotspot distribution; the acceptance gains above come exclusively from
the repeated normal-build comparisons.

Across the campaign, 88 complete corpus runs and 70 single-file runs passed
per-file output identity checks. The aggregate corpus output-manifest SHA-256 is
`830e208e80ef71b85cc78a790560704c9ba646264ff10f7aacbcf5ac29d9e7be`;
its canonical serialization and all binary hashes are recorded in the JSON.

## 1. Reverse-rank elimination — equivalence argument

The first orientation rejects non-finite scores before sorting. For finite
scores, negation reverses the order of numeric tie groups. Reversing the existing
descending permutation therefore visits the same groups as sorting negated
scores or descending dense ranks. Signed zeros remain one numeric tie group.
Within-group permutation cannot affect cumulative counts: only additions of one
are performed, and FDP is evaluated after the whole group. Both implementations
use the same floating-point counts, formula, and strict comparison. Thus the
last accepted group and target count are identical. The feature/orientation
first-maximum selection rule is unchanged.

Adversarial coverage enumerates labels across empty through nine-row inputs,
finite extrema, subnormals, signed zeros, reported/training estimators, unequal
opportunity ratios, NaN/infinite thresholds, exact q thresholds and their next
representable values. The old public dense-rank API remains available, but the
production initial-direction path no longer calls it.

Paired full-corpus times (seconds), before → after:

| Repetition | Sequential | N=4 |
|---|---|---|
| 1 | 55.018 → 47.041 | 15.690 → 15.990 |
| 2 | 51.966 → 47.348 | 16.318 → 15.319 |
| 3 | 50.210 → 46.393 | 15.547 → 14.778 |
| Median | 51.966 → 47.041 (9.48%) | 15.690 → 15.319 (2.37%) |

The first before sequential run incurred 65 major faults; subsequent runs had
none. Both warm sequential pairs also improve. N=4 improves in two of three
pairs, with substantial host variation. Every corpus output size/hash agrees;
all frozen fixture files additionally pass direct byte comparison in fixed
serial, fixed parallel, selected-C, and ensemble modes, including proteins.

## 2. Proposed feature-order reuse — equivalence argument

Raw feature sorting requires no labels or fitted statistics. Removing held-out
rows from that permutation leaves the training feature values sorted. Each
fold still fits its own mean and positive standard deviation on training rows.
The existing rounded subtraction followed by positive division is numerically
monotone on finite results, although it can merge originally distinct values.
Consequently, counts must use equality of the actual normalized scores, never
cached raw-score ranks. Signed zeros remain a single numeric group. Changing
held-out values can only change where removed rows occur in the global order;
the remaining groups and their counts are unaffected. Label changes cannot
affect the cache at all. The selected feature/orientation uses the existing
first-maximum rule. Label-dependent RT residuals must use fold-local ordering.

The intended cache holds only row indices, with no labels, counts, q-values, or
fold-fitted values. A bounded index representation needs an explicit size guard
and an uncached fallback. These conditions are prerequisites, not permission
to assume raw and normalized tie groups are identical.

After step 1, the full sequential profile measured q-value/count/mask work at
36.63%, score sorting at 24.72%, tie/count scans at 11.06%, and initial direction
at 23.44%. Reverse-rank sorts disappeared. This supports trying step 2 next.

The step-2 experiment caches `u32` row permutations for the fixed-C outer CV
pass, falling back to the existing path with RT or more than `u32::MAX` rows.
Both cached and uncached paths are compared under held-out feature/label attacks
and a constructed normalization-rounding collision, including threshold edges.

**Rejected and reverted.** Largest-file five-repeat medians: 1.496 → 1.516 s.
N=4 pairs: 14.607 → 15.095, 14.741 → 15.090, 14.710 → 15.202 s;
median regression 2.62%, with all three pairs worse. Peak N=4 RSS rises from
820–827 MB to 907–923 MB. All outputs agree. The failed source snapshots and
profiles are retained in the artifact directory; only the normalization-order
adversarial property remains in the test suite.

## 3. Unchanged-model reuse — equivalence argument

Inside `train_fold`, the design matrix, row sequence, labels, TDC parameters,
and threshold remain fixed. If every SVM weight has identical `to_bits()` before
and after training, the next deterministic scoring call must return the same
score bits. Therefore its score order, q-value decisions, and positive mask are
unchanged. Retaining scores and mask is exact. This does not use tolerance-based
equality. Training, positive selection, subsampling, RNG draws, and iteration
count must still run exactly as before. A changed weight invalidates reuse;
MLP training is conservatively always invalidated. If no training occurs, the
model and retained results remain valid without a weight comparison.

**Accepted.** N=4 pairs: 15.033 → 14.963, 15.068 → 14.927,
15.096 → 14.948 s (median improvement 0.79%, all three pairs improve).
Largest-file five-repeat medians are flat at 1.517 s. Frozen bytes and all
corpus hashes match. Tests compare forced recomputation against reuse, including
final score bits and RNG state, with subsampling, ordinary/zero updates,
immediate convergence, empty positives, signed-zero initialization, and MLP.

## 4. Finite numeric comparator — equivalence argument

Every production caller of `sort_score_order` first asserts that all scores
are finite. On that domain the old comparator always selects its `(false,
false)` arm, which is exactly `b.total_cmp(&a)`. Removing the per-comparison NaN
branches therefore preserves every comparator result, including the distinct
ordering of signed zeros. Index type, input permutation, sorting algorithm,
and tie grouping remain unchanged. Non-finite inputs still fail at the existing
entry-point assertions. Exact index-permutation tests against the old comparator
will cover finite bit patterns, duplicates, signed zeros, and sort-size edges.

**Accepted.** Largest-file five-repeat median: 1.516 → 1.415 s (6.67%).
N=4 pairs: 14.838 → 14.175, 14.918 → 14.073, 14.961 → 14.119 s;
median improvement 5.36%, all three pairs improve. Frozen bytes and all corpus
hashes match; exact legacy permutation and non-finite rejection tests pass.

Step 3's sequential profile recorded 155 reused training passes out of 1,950
iterations. Q-value/count/mask work was 35.66%, sorting 23.85%, and tie/count
scans 10.91%, which supported the finite comparator experiment.

## 5. Fused orientation scans — equivalence argument

At a descending tie group, the reverse-orientation cumulative counts are the
whole-list counts minus the counts strictly before that group. For lists no
longer than 2^53−1, every count, safeguard addition, and subtraction is an exactly
representable integer in f64, so these counts have exactly the bits produced by
the old reverse walk. Larger lists retain the two-walk fallback. FDP evaluation
and its strict comparison are unchanged. Forward orientation still takes the
last qualifying group. Reverse cumulative target counts decrease along this
walk, so the first qualifying reverse group gives the largest accepted reverse
target set; subsequent groups cannot increase its count. Both orientations use
the actual numeric score ties and the original feature/orientation choice order.

**Accepted.** Largest-file five-repeat median: 1.415 → 1.355 s (4.27%).
N=4 pairs: 14.054 → 13.666, 14.174 → 13.736, 14.128 → 13.658 s;
median improvement 3.27%, all three pairs improve. All frozen bytes and corpus
hashes agree. Sequential tie/count time falls from 4.795 to 3.143 s. The new
profile is q-values 26.96%, sorting 17.91%, tie/count scans 8.08%, initial
direction 18.74% (previously 29.66%, 16.95%, 11.63%, 21.93%).

## 6. Training-mask specialization — equivalence argument

For the training estimator with `pi0 == 1`, adding only decoys cannot decrease
raw FDP: the nonnegative opportunity ratio and the target denominator are
unchanged, and rounding/clamping preserve monotonicity. A qualifying decoy-only
group therefore cannot admit a target that was not already admitted by the
preceding target-containing group. Skipping its FDP evaluation leaves every
target mask byte unchanged. Mixed tie groups still include every target and
decoy before evaluation. Other TDC configurations keep the original scan.
No division is replaced by a threshold-times-count comparison.

**Rejected and reverted.** Largest-file medians are flat at 1.354–1.355 s.
N=4 pairs: 13.644 → 13.733, 13.686 → 13.710, 13.717 → 13.673 s;
median regression 0.18%, two of three pairs worse. Output checks pass. The
exhaustive mask/threshold/decoy-tail test remains; the specialization does not.

## 7. N=4 memory/cache experiment — equivalence argument

Sort `(score, index)` entries contiguously for count/mask-only operations, then
recover the index permutation. The comparator is the same finite `total_cmp`.
Element layout can change within-tie permutations, but numeric groups and their
cumulative target/decoy counts remain identical, so count and mask results do
too. Reported q-value/PEP paths must retain the existing index sort because
decoy PEP presentation can depend on within-tie permutation. Reusing a mask
workspace in a reported calculation still resets and sorts its index vector.

This experiment trades an additional reusable 16-byte-per-row scratch buffer
for fewer dependent score loads during sorting. It will be accepted only if
N=4 measurements justify the added memory. Hardware cache counters are denied
(`perf_event_paranoid=4`); cache-miss rates cannot be claimed. After step 5,
N=4/serial aggregate stage-time ratios were 1.11 for sorting, 1.10 for tie scans,
1.28 for normalization, 1.50 for model scoring, and 1.37 for training. These
show that the largest remaining scaling penalties are outside q-value scans.

**Rejected and reverted.** Largest-file medians: 1.354 → 1.355 s. N=4 pairs:
13.542 → 14.346, 13.648 → 14.270, 13.688 → 14.265 s; median regression 4.56%,
all three pairs worse. Peak RSS increases from 820–822 MB to 831–832 MB;
minor faults rise from about 2.927 million to 2.983 million per corpus.
All outputs agree, including the explicit reported-PEP workspace/tie test.
The extra scratch storage and changed sort layout are not retained.

## Final validation and reproduction

`python3 refactor/verify_baseline.py` passes the complete release all-targets
suite, all six portable shell gates, frozen-output comparisons, the structural
adversarial driver, and the independent arithmetic and protein probes.
`cargo fmt --all -- --check`, all-target/all-feature Clippy with warnings denied,
and all-feature rustdoc with warnings denied also pass. The only lint repair was
a profiling-only equivalent byte-count expression.

The retained implementation changes only `src/stats.rs` and
`src/percolator.rs`. Raw cross-fold feature caches, mask specialization, and
contiguous score/index scratch storage have been removed. Public statistics APIs
remain available. PEP formulas, protein inference, CV topology, all statistical
thresholds/formulas, and identification yield are unchanged.

Use the paired runner with separately built binaries; it refuses to overwrite
artifacts, checks the frozen fixture by direct byte comparison, and compares
every per-file TSV size/hash after each untimed hashing interval. For example:

```bash
python3 bench/optimize_hotpaths.py \
  --before /path/to/baseline --after /path/to/candidate \
  --profile /path/to/profiling-candidate \
  --artifacts /tmp/new-hotpath-experiment \
  --concurrency 1 2 4 6 --profile-concurrency 1 4 --repeats 3
```

Single-file screens use `--single --concurrency 1 --repeats 5`. Normal binaries
determine acceptance. Stage profiles are separate full-corpus runs after the
timing repetitions. The input order and 20 ms process sampling interval match
the fresh-runtime collector. Compilation and test runs do not overlap timings.
Environment: Ryzen 5 5600G (6 cores/12 threads, 16 MiB L3), 30.69 GiB RAM,
unused 8 GiB swap, Linux 7.0.0-30-generic, Rust 1.97.0/LLVM 22.1.6, thin LTO,
one codegen unit, `x86-64-v3`; corpus inputs on ntfs3 and output on `/tmp` tmpfs.
