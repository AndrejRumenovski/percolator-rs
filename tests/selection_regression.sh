#!/usr/bin/env bash
# Regression gate for leakage-free nested SVM model selection.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

BIN=target/release/percolator-rs
FIX=tests/fixtures/sample.pin
[ -x "$BIN" ] || cargo build --release
TMP_SELECTION_REGRESSION="$(mktemp -d)"
trap 'rm -rf "$TMP_SELECTION_REGRESSION"' EXIT

run_selection() {
  local name=$1 threads=$2
  "$BIN" --canonical --seed 1 --auto-model --num-threads "$threads" \
    --results-psms "$TMP_SELECTION_REGRESSION/$name.target.psms.tsv" \
    --decoy-results-psms "$TMP_SELECTION_REGRESSION/$name.decoy.psms.tsv" \
    --results-peptides "$TMP_SELECTION_REGRESSION/$name.target.peptides.tsv" \
    --decoy-results-peptides "$TMP_SELECTION_REGRESSION/$name.decoy.peptides.tsv" \
    "$FIX" 2>"$TMP_SELECTION_REGRESSION/$name.log"
}

run_selection serial 1
run_selection parallel 3

for kind in target.psms decoy.psms target.peptides decoy.peptides; do
  cmp "$TMP_SELECTION_REGRESSION/serial.$kind.tsv" \
    "$TMP_SELECTION_REGRESSION/parallel.$kind.tsv"
done
cmp <(sed -n '/^  fold /p' "$TMP_SELECTION_REGRESSION/serial.log") \
  <(sed -n '/^  fold /p' "$TMP_SELECTION_REGRESSION/parallel.log")

if "$BIN" --auto-model --select-c "$FIX" >/dev/null 2>&1; then
  echo "FAIL: --auto-model accepted conflicting --select-c"
  exit 1
fi
if "$BIN" --auto-model --rescore-model mlp "$FIX" >/dev/null 2>&1; then
  echo "FAIL: --auto-model accepted unsupported MLP learner"
  exit 1
fi

folds=$(sed -n '/^  fold /p' "$TMP_SELECTION_REGRESSION/serial.log" | wc -l)
psm=$(sed -n 's/.*target PSMs q<0.01: \([0-9]*\).*/\1/p' "$TMP_SELECTION_REGRESSION/serial.log")
# q<0.05 at peptide level: see tests/expected.env.
peptide=$(awk -F'\t' 'NR>1 && $3<0.05' "$TMP_SELECTION_REGRESSION/serial.target.peptides.tsv" | wc -l)

echo "== percolator-rs nested-selection regression gate =="
[ "$folds" -eq 3 ] || { echo "  FAIL  selected outer folds: $folds"; exit 1; }
[ "$psm" -eq 117 ] || { echo "  FAIL  PSM q<0.01: $psm (expected 117)"; exit 1; }
[ "$peptide" -eq 43 ] || { echo "  FAIL  peptide q<0.01: $peptide (expected 43)"; exit 1; }
echo "  PASS  selected outer folds  $folds"
echo "  PASS  PSM q<0.01           $psm"
echo "  PASS  peptide q<0.05       $peptide"
echo "  PASS  serial/parallel      byte-identical choices and outputs"
echo "ALL CHECKS PASSED"
