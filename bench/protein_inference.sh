#!/usr/bin/env bash
# Compare picked-protein FDR with Fido-style Bayesian inference on five PIN schemas.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
OUT="${PROTEIN_BENCH_OUT:-$HOME/percolator_rs_out/protein-inference-benchmark}"
INPUT="${MULTIDATASET_OUT:-$HOME/percolator_rs_out/multidataset}/inputs"
BIN="${RS:-$ROOT/target/release/percolator-rs}"
REPEATS="${REPEATS:-3}"

[[ -x "$BIN" ]] || cargo build --release --manifest-path "$ROOT/Cargo.toml"
mkdir -p "$OUT"

datasets=(PXD032157_fixture hogrebe_tide PXD020243_msfragger PXD060954_sage percolator_yeast)
pins=(
  "$ROOT/tests/fixtures/sample.pin"
  "$INPUT/hogrebe_tide.pin"
  "$INPUT/PXD020243_msfragger.pin"
  "$INPUT/PXD060954_sage.pin"
  "$INPUT/percolator_yeast.pin"
)

for pin in "${pins[@]}"; do
  if [[ ! -f "$pin" ]]; then
    echo "missing $pin; run bash bench/multidataset/run.sh first" >&2
    exit 2
  fi
done

median_column() {
  local file="$1" column="$2"
  cut -f"$column" "$file" | sort -n |
    awk '{value[NR]=$1} END { print value[int((NR+1)/2)] }'
}

count_q01() {
  awk -F'\t' 'NR > 1 && $2+0 < 0.01 { count++ } END { print count+0 }' "$1"
}

printf 'dataset\tmethod\twall_seconds\tpeak_rss_kb\tprotein_groups\ttarget_q_lt_0.01\tcomponents\tloopy_components\tbp_iterations\tconverged\n' \
  > "$OUT/results.tsv"

for index in "${!datasets[@]}"; do
  dataset="${datasets[$index]}"
  pin="${pins[$index]}"
  for method in picked bayesian; do
    run="$OUT/$dataset/$method"
    mkdir -p "$run"
    : > "$run/times.tsv"
    for repeat in $(seq 1 "$REPEATS"); do
      /usr/bin/time -f '%e\t%M' -o "$run/time.$repeat.tsv" \
        "$BIN" --canonical --seed 1 --protein-inference "$method" \
          --results-proteins "$run/target.tsv" \
          --decoy-results-proteins "$run/decoy.tsv" \
          "$pin" >"$run/stdout.log" 2>"$run/stderr.log"
      cat "$run/time.$repeat.tsv" >> "$run/times.tsv"
    done

    wall="$(median_column "$run/times.tsv" 1)"
    rss="$(median_column "$run/times.tsv" 2)"
    groups="$(sed -n 's/^protein groups: \([0-9]*\).*/\1/p' "$run/stderr.log")"
    q01="$(count_q01 "$run/target.tsv")"
    if [[ "$method" == bayesian ]]; then
      diagnostics="$(sed -n \
        's/.*components: \([0-9]*\) (.* \([0-9]*\) loopy); BP iterations: \([0-9]*\), converged: \(.*\)$/\1\t\2\t\3\t\4/p' \
        "$run/stderr.log")"
    else
      diagnostics=$'-\t-\t-\t-'
    fi
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
      "$dataset" "$method" "$wall" "$rss" "$groups" "$q01" "$diagnostics" \
      >> "$OUT/results.tsv"
  done
done

column -t -s $'\t' "$OUT/results.tsv" 2>/dev/null || cat "$OUT/results.tsv"
echo "outputs: $OUT"
