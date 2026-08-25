# Validation requirement completion audit

This checklist maps the requested adversarial validation phases to preserved evidence. “Complete”
means the question was answered; it does not imply a positive scientific result. Where legitimate
ground truth was unavailable, completion means that the absence was established and no substitute
was mislabeled as truth.

| Requested phase | Status | Evidence |
|---|---|---|
| 1. Read-only implementation audit | Complete | [`IMPLEMENTATION_AUDIT.md`](IMPLEMENTATION_AUDIT.md), including every named component and one-PSM trace |
| Leakage influence inventory | Complete | audit data-availability and leakage tables; default, `--select-c`, `--auto-model`, and RT-feature paths separated |
| 2. C++ PSM agreement | Complete | [`psm_agreement.py`](psm_agreement.py); five-seed PSM correlations, overlaps, rank/q distributions, exclusives, and examples |
| 3. Multi-seed reproducibility | Complete | seeds 1--5 for Rust and C++ on four compact datasets; all runs retained in v3 manifest |
| 4. Pure-null validation | Complete; failed scientifically | three source PINs x ten exact-balance relabelings, six thresholds, source-stratified Wilson intervals |
| 5. Negative controls | Complete | exact-label null plus all-features-tied fixture; interpretation and dependence limitations documented |
| 6. Ground truth | Complete with explicit availability limit | no complete PSM truth found; entrapment treated as partial truth; PrEST used only for protein present/absent evidence |
| 7. Multiple datasets | Complete for compatibility/agreement; calibration unavailable on compact cases | Comet, Tide, MSFragger, Sage, and legacy fixture across diverse organisms/instruments; PrEST protein standard |
| 8. Threshold calibration | Complete where labels support it | 0.001, 0.005, 0.01, 0.02, 0.05, 0.10 on null and five-seed entrapment studies |
| 9. Ablation | Complete to the evidence available | one-factor fixed-score q-estimator variants plus recorded weighting, selector, nested, and learner comparisons; precision caveat explicit |
| 10. Statistical analysis | Complete | effect sizes, distributions, paired overlap, variability, descriptive intervals, and reasons for avoiding significance tests |
| 11. Edge cases | Complete | small, imbalance, ties, duplicates, malformed, missing, NaN, extreme, constant, and protein-mapping fixtures |
| 12. Claim audit | Complete | every major README scientific claim classified; README and stale benchmark wording corrected |
| Experiment provenance | Complete for authoritative runs | dataset/input hashes, commit, C++ version, argv, seeds, parameters, platform, outputs, and evaluator hashes in manifests |
| Failed experiments | Complete | missing C++ library, fail-closed null run, smoke test, and extra entrapment PIN retained and disclosed |
| Final 16-section report | Complete | [`SCIENTIFIC_VALIDATION.md`](SCIENTIFIC_VALIDATION.md), separate verdicts and final claim boundary |

## Non-negotiable-rule audit

- No Rust production or statistical methodology was changed.
- No seed, source PIN, dataset, or failed run was removed.
- Thresholds were fixed at the requested six values.
- C++ was forced to the same concatenated search-input mode in new matched comparisons; the old
  auto/mix-max headline mismatch remains disclosed.
- No identification count is called accuracy or sensitivity.
- No target-decoy estimate is called ground truth.
- Runtime was recorded only as provenance/secondary information and was not optimized.
- Negative results are the headline conclusion, not hidden in limitations.

## Verification

The final workspace passed `cargo test --release`, all six portable shell regression gates, the
Sage-normalizer test, four PrEST report tests, both PSM-agreement tests, Python byte-compilation,
Markdown relative-link checks, and `git diff --check`.
