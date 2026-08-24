#!/usr/bin/env bash
# Compare fixed and leakage-free nested SVM selection on existing entrapment PINs.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
BIN="$ROOT/target/release/percolator-rs"
SEARCH_ROOT="${ENTRAPMENT_WORK:-$HOME/percolator_rs_out/entrapment}"
OUT="${SELECTION_ENTRAPMENT_OUT:-$HOME/percolator_rs_out/selection-comparison-entrapment}"
JOBS="${SELECTION_BENCH_JOBS:-4}"

cd "$ROOT"
case "$OUT" in
  ""|/|"$HOME") echo "FAIL: unsafe SELECTION_ENTRAPMENT_OUT: $OUT"; exit 2 ;;
esac
[ -f "$SEARCH_ROOT/database-stats.txt" ] || {
  echo "FAIL: missing $SEARCH_ROOT/database-stats.txt; run bench/entrapment/run.sh first"
  exit 2
}
# Crux can leave its generic default output directory beside the six explicitly
# named sample directories.  It is not a seventh sample.
mapfile -t pins < <(find "$SEARCH_ROOT" -mindepth 2 -maxdepth 2 -type f -name comet.pin \
  ! -path "$SEARCH_ROOT/comet-out/*" | sort)
[ "${#pins[@]}" -eq 6 ] || { echo "FAIL: expected 6 entrapment PINs, found ${#pins[@]}"; exit 2; }
cargo build --release
rm -rf "$OUT"
mkdir -p "$OUT"

export BIN OUT
run_entrapment_selection() {
  local task=$1 method pin sample selection_flag=()
  method=${task%%::*}
  pin=${task#*::}
  sample=$(basename "$(dirname "$pin")")
  [ "$method" = nested ] && selection_flag=(--auto-model)
  "$BIN" --canonical --seed 1 --rescore-model svm "${selection_flag[@]}" \
    --results-psms "$OUT/$method.$sample.target.psms.tsv" \
    --decoy-results-psms "$OUT/$method.$sample.decoy.psms.tsv" \
    "$pin" 2>"$OUT/$method.$sample.log"
}
export -f run_entrapment_selection

tasks=()
for method in fixed nested; do
  for pin in "${pins[@]}"; do
    tasks+=("$method::$pin")
  done
done
printf '%s\n' "${tasks[@]}" | xargs -P "$JOBS" -I{} bash -c 'run_entrapment_selection "$1"' _ {}

fraction=$(sed -n 's/^entrapment_fraction=//p' "$SEARCH_ROOT/database-stats.txt")
mapfile -t targets < <(find "$OUT" -maxdepth 1 -type f -name '*.target.psms.tsv' | sort)
python3 "$ROOT/bench/entrapment/report.py" --entrapment-fraction "$fraction" \
  --output "$OUT/report.tsv" "${targets[@]}"
echo "Report: $OUT/report.tsv"
