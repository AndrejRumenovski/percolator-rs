#!/usr/bin/env bash
# Compare SVM and MLP reported q-values on the existing signal-present
# entrapment searches. Run bench/entrapment/run.sh first to create the PINs.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
BIN="$ROOT/target/release/percolator-rs"
SEARCH_ROOT="${ENTRAPMENT_WORK:-$HOME/percolator_rs_out/entrapment}"
OUT="${MODEL_ENTRAPMENT_OUT:-$HOME/percolator_rs_out/model-comparison-entrapment}"
JOBS="${MODEL_BENCH_JOBS:-4}"

cd "$ROOT"
case "$OUT" in
  ""|/|"$HOME") echo "FAIL: unsafe MODEL_ENTRAPMENT_OUT: $OUT"; exit 2 ;;
esac
[ -f "$SEARCH_ROOT/database-stats.txt" ] || {
  echo "FAIL: missing $SEARCH_ROOT/database-stats.txt; run bench/entrapment/run.sh first"
  exit 2
}
mapfile -t pins < <(find "$SEARCH_ROOT" -mindepth 2 -maxdepth 2 -type f -name comet.pin | sort)
[ "${#pins[@]}" -eq 6 ] || { echo "FAIL: expected 6 entrapment PINs, found ${#pins[@]}"; exit 2; }
cargo build --release
rm -rf "$OUT"
mkdir -p "$OUT"

export BIN OUT
run_entrapment_model() {
  local task=$1 model pin sample
  model=${task%%::*}
  pin=${task#*::}
  sample=$(basename "$(dirname "$pin")")
  "$BIN" --canonical --seed 1 --rescore-model "$model" \
    --results-psms "$OUT/$model.$sample.target.psms.tsv" \
    --decoy-results-psms "$OUT/$model.$sample.decoy.psms.tsv" \
    "$pin" 2>"$OUT/$model.$sample.log"
}
export -f run_entrapment_model

tasks=()
for model in svm mlp; do
  for pin in "${pins[@]}"; do
    tasks+=("$model::$pin")
  done
done
printf '%s\n' "${tasks[@]}" | xargs -P "$JOBS" -I{} bash -c 'run_entrapment_model "$1"' _ {}

fraction=$(sed -n 's/^entrapment_fraction=//p' "$SEARCH_ROOT/database-stats.txt")
mapfile -t targets < <(find "$OUT" -maxdepth 1 -type f -name '*.target.psms.tsv' | sort)
python3 "$ROOT/bench/entrapment/report.py" --entrapment-fraction "$fraction" \
  --output "$OUT/report.tsv" "${targets[@]}"
echo "Report: $OUT/report.tsv"
