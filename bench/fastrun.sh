#!/usr/bin/env bash
# Trimmed C++ Percolator run over PXD032157. This is a speed/yield tradeoff,
# not the canonical reference configuration.
# Usage: bash bench/fastrun.sh [file concurrency]
set -euo pipefail

N="${1:-5}"
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
PERC="${PERC:-$HOME/opt/percolator-root/usr/bin/percolator}"
LIBS="${PERC_LIBS:-$HOME/opt/perc-libs}"
IN="${CPP_FAST_INPUT:-$ROOT/data/PXD032157}"
OUT="${CPP_FAST_OUT:-$HOME/percolator_fast_out/PXD032157}"
EXPECTED_FILES="${CPP_FAST_EXPECTED_FILES:-65}"
export LC_ALL=C

[[ "$N" =~ ^[1-9][0-9]*$ ]] || { echo "FAIL: concurrency must be positive"; exit 2; }
[[ "$EXPECTED_FILES" =~ ^[1-9][0-9]*$ ]] || { echo "FAIL: expected file count must be positive"; exit 2; }
[ -d "$IN" ] || { echo "FAIL: input directory not found: $IN"; exit 2; }
[ -x "$PERC" ] || { echo "FAIL: C++ Percolator not executable: $PERC"; exit 2; }
case "$OUT" in ""|/|"$HOME") echo "FAIL: unsafe CPP_FAST_OUT: $OUT"; exit 2 ;; esac
export LD_LIBRARY_PATH="$LIBS:${LD_LIBRARY_PATH:-}"

mapfile -t size_order < <(find "$IN" -maxdepth 1 -type f -name '*.pin' -printf '%s\t%p\n' | sort -rn | cut -f2-)
[ "${#size_order[@]}" -eq "$EXPECTED_FILES" ] || {
  echo "FAIL: expected $EXPECTED_FILES PINs, found ${#size_order[@]}"; exit 2;
}
order=()
left=0
right=$((${#size_order[@]}-1))
while [ "$left" -le "$right" ]; do
  order+=("${size_order[$left]}")
  [ "$left" -eq "$right" ] || order+=("${size_order[$right]}")
  left=$((left+1))
  right=$((right-1))
done

rm -rf "$OUT"
mkdir -p "$OUT"
peak_file="$OUT/peak-rss-kb.txt"
stop_file="$OUT/.monitor-stop"
echo 0 >"$peak_file"
(
  peak=0
  while [ ! -f "$stop_file" ]; do
    current=$(ps --no-headers -o rss -C percolator 2>/dev/null |
      awk '{sum+=$1} END{print sum+0}' || true)
    if [ "${current:-0}" -gt "$peak" ]; then
      peak=$current
      echo "$peak" >"$peak_file"
    fi
    sleep 0.1
  done
) &
monitor=$!
stop_monitor() {
  touch "$stop_file"
  if [ -n "${monitor:-}" ]; then wait "$monitor" 2>/dev/null || true; monitor=""; fi
}
trap stop_monitor EXIT

export PERC OUT
run_fast_file() {
  local pin=$1 stem destination
  stem=$(basename "$pin" .pin)
  destination="$OUT/$stem"
  mkdir -p "$destination"
  "$PERC" --seed 1 --num-threads 1 --subset-max-train 20000 --maxiter 5 \
    --results-psms "$destination/target.psms.tsv" \
    --decoy-results-psms "$destination/decoy.psms.tsv" \
    --results-peptides "$destination/target.peptides.tsv" \
    --decoy-results-peptides "$destination/decoy.peptides.tsv" \
    "$pin" >"$destination/stdout.log" 2>"$destination/stderr.log"
}
export -f run_fast_file

start=$(date +%s.%N)
printf '%s\n' "${order[@]}" | xargs -P "$N" -I{} bash -c 'run_fast_file "$1"' _ {}
end=$(date +%s.%N)
stop_monitor
trap - EXIT

count_q_values() {
  awk -F'\t' '
    FNR == 1 {q=0; for (i=1; i<=NF; i++) if ($i=="q-value") q=i; next}
    q && $q+0 < 0.01 {count++}
    END {print count+0}
  ' "$@"
}
valid=$(find "$OUT" -mindepth 2 -maxdepth 2 -type f -name target.psms.tsv -size +0c | wc -l)
[ "$valid" -eq "$EXPECTED_FILES" ] || { echo "FAIL: only $valid/$EXPECTED_FILES result files are valid"; exit 1; }
psm=$(count_q_values "$OUT"/*/target.psms.tsv)
peptide=$(count_q_values "$OUT"/*/target.peptides.tsv)
peak=$(cat "$peak_file")
awk -v n="$N" -v start="$start" -v end="$end" -v peak="$peak" -v valid="$valid" \
    -v psm="$psm" -v peptide="$peptide" '
  BEGIN {
    print "implementation\tconcurrency\twall_seconds\tpeak_rss_kb\tvalid_files\tpsm_q_lt_0.01\tpeptide_q_lt_0.01"
    printf "cpp-fast\t%d\t%.3f\t%d\t%d\t%d\t%d\n",n,end-start,peak,valid,psm,peptide
  }
' | tee "$OUT/summary.tsv"
