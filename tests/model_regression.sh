#!/usr/bin/env bash
# Portable regression gate for the optional MLP scorer. It verifies that the
# model uses out-of-fold inference deterministically in serial and parallel.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

BIN=target/release/percolator-rs
FIX=tests/fixtures/sample.pin
[ -x "$BIN" ] || cargo build --release
TMP_MODEL_REGRESSION="$(mktemp -d)"
trap 'rm -rf "$TMP_MODEL_REGRESSION"' EXIT

run_model() {
  local name=$1 threads=$2
  "$BIN" --canonical --seed 1 --rescore-model mlp --num-threads "$threads" \
    --results-psms "$TMP_MODEL_REGRESSION/$name.target.psms.tsv" \
    --decoy-results-psms "$TMP_MODEL_REGRESSION/$name.decoy.psms.tsv" \
    --results-peptides "$TMP_MODEL_REGRESSION/$name.target.peptides.tsv" \
    --decoy-results-peptides "$TMP_MODEL_REGRESSION/$name.decoy.peptides.tsv" \
    "$FIX" 2>"$TMP_MODEL_REGRESSION/$name.log"
}

run_model serial 1
run_model parallel 3

for kind in target.psms decoy.psms target.peptides decoy.peptides; do
  cmp "$TMP_MODEL_REGRESSION/serial.$kind.tsv" "$TMP_MODEL_REGRESSION/parallel.$kind.tsv"
done

psm=$(sed -n 's/.*target PSMs q<0.01: \([0-9]*\).*/\1/p' "$TMP_MODEL_REGRESSION/serial.log")
# q<0.05 at peptide level: on a fixture this small the corrected estimator's best
# statement is 1/38, so a q<0.01 peptide gate would assert 0 for any build.
peptide=$(awk -F'\t' 'NR>1 && $3<0.05' "$TMP_MODEL_REGRESSION/serial.target.peptides.tsv" | wc -l)

echo "== percolator-rs MLP regression gate =="
[ "$psm" -eq 130 ] || { echo "  FAIL  PSM q<0.01: $psm (expected 130)"; exit 1; }
[ "$peptide" -eq 38 ] || { echo "  FAIL  peptide q<0.01: $peptide (expected 38)"; exit 1; }
echo "  PASS  PSM q<0.01           $psm"
echo "  PASS  peptide q<0.05       $peptide"
echo "  PASS  serial/parallel      byte-identical outputs"
echo "ALL CHECKS PASSED"
