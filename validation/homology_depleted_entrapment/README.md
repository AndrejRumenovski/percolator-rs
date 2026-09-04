# Homology-depleted entrapment experiment

This directory contains a prospective, pre-search causal intervention on the
foreign-proteome entrapment database used by the current Percolator validation.
It does not modify production `percolator-rs` code or statistical methodology.

The experiment is ordered deliberately:

1. `baseline_reproduction/` freezes the audited inputs and reproduces the
   canonical seed-1 result before intervention construction.
2. `PREREGISTRATION.md` freezes the homology rule, depletion unit, controls,
   endpoints, seeds, and interpretation before intervention searches.
3. `construct_databases.py` constructs all target FASTAs from the frozen source
   database and writes pre-search characterizations.
4. `run_pipeline.sh` reruns Comet from the six raw mzML files, allowing Comet to
   regenerate reversed decoys for every condition, then runs the unmodified
   audited `percolator-rs` binary.
5. `analyze.py` computes raw-search and rescored endpoints without filtering
   PSMs from one condition into another.

Large generated FASTAs, PINs, pepXML files, and PSM tables are kept below this
directory rather than overwriting any earlier validation evidence.

## Result and artifact map

The confirmatory conclusion is **SUPPORTED**: pre-search homology depletion
restored the primary q<0.01 internal-null ratio from 258/133 = 1.940 to 102/102
= 1.000, while three size-matched controls remained at 1.543--2.095. The
complete causal interpretation and limitations are in `FINAL_REPORT.md`.

- `METHODOLOGY.md`: design, definitions, exact commands, and uncertainty model.
- `PREREGISTRATION.md`: immutable prospective rule and endpoints.
- `PRESEARCH_CHARACTERIZATION.md`: result-blind database comparison.
- `analysis/rescored_results.json`: machine-readable confirmatory analysis.
- `FINAL_SUMMARY.json`: compact machine-readable conclusions and headline
  endpoints.
- `analysis/raw_xcorr.json`: raw search-engine exchangeability.
- `analysis/exploratory_remaining.json`: explicitly exploratory residual tail.
- `tables/`: machine-readable endpoint, calibration, dose, seed, and ablation
  tables.
- `figures/`: PNG and SVG diagnostic figures.
- `search_manifest.json`: hashes and integrity checks for all 30 fresh searches.
- `ARTIFACT_MANIFEST.json` and `SHA256SUMS.txt`: complete integrity inventory.

No production `percolator-rs` source or statistical methodology was changed.
