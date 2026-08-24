#!/usr/bin/env bash
# Batch-run canonical C++ Percolator over the PXD032157 benchmark.
# Usage: bash bench/run_cpp.sh [file-level concurrency]
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
IN="${CPP_BENCH_INPUT:-$ROOT/data/PXD032157}"
OUT="${CPP_BENCH_OUT:-$HOME/percolator_rs_out/cpp-canonical}"
PERC="${PERC:-$HOME/opt/percolator-root/usr/bin/percolator}"
LIBS="${PERC_LIBS:-$HOME/opt/perc-libs}"
N="${1:-4}"
EXPECTED_FILES="${CPP_BENCH_EXPECTED_FILES:-65}"

[[ -d "$IN" ]] || { echo "FAIL: input directory not found: $IN"; exit 2; }
[[ -x "$PERC" ]] || { echo "FAIL: C++ Percolator not executable: $PERC"; exit 2; }
case "$OUT" in
  ""|/|"$HOME") echo "FAIL: unsafe CPP_BENCH_OUT: $OUT"; exit 2 ;;
esac
[[ "$N" =~ ^[1-9][0-9]*$ ]] || { echo "FAIL: concurrency must be a positive integer"; exit 2; }
[[ "$EXPECTED_FILES" =~ ^[1-9][0-9]*$ ]] || {
  echo "FAIL: CPP_BENCH_EXPECTED_FILES must be a positive integer"
  exit 2
}

export LD_LIBRARY_PATH="$LIBS:${LD_LIBRARY_PATH:-}"
rm -rf "$OUT"
mkdir -p "$OUT"

mapfile -t inputs < <(find "$IN" -type f -name '*.pin' -printf '%s\t%p\n' | sort -rn | cut -f2-)
[[ "${#inputs[@]}" -eq "$EXPECTED_FILES" ]] || {
  echo "FAIL: expected $EXPECTED_FILES PINs, found ${#inputs[@]}"
  exit 2
}

peak_file="$OUT/peak-rss-kb.txt"
stop_file="$OUT/.monitor-stop"
echo 0 > "$peak_file"
(
  peak=0
  while [[ ! -f "$stop_file" ]]; do
    current=$(
      ps --no-headers -o rss -C percolator 2>/dev/null |
        awk '{sum += $1} END {print sum + 0}' || true
    )
    if (( current > peak )); then
      peak=$current
      echo "$peak" > "$peak_file"
    fi
    sleep 0.1
  done
) &
monitor=$!
stop_monitor() {
  touch "$stop_file"
  if [[ -n "${monitor:-}" ]]; then
    wait "$monitor" 2>/dev/null || true
    monitor=""
  fi
}
trap stop_monitor EXIT

export PERC OUT
run_cpp_file() {
  local pin=$1 stem destination
  stem=$(basename "$pin" .pin)
  destination="$OUT/$stem"
  mkdir -p "$destination"
  "$PERC" --seed 1 --num-threads 1 \
    --results-psms "$destination/target.psms.tsv" \
    --decoy-results-psms "$destination/decoy.psms.tsv" \
    --results-peptides "$destination/target.peptides.tsv" \
    --decoy-results-peptides "$destination/decoy.peptides.tsv" \
    "$pin" >"$destination/stdout.log" 2>"$destination/stderr.log"
}
export -f run_cpp_file

start=$(date +%s.%N)
printf '%s\n' "${inputs[@]}" | xargs -P "$N" -I{} bash -c 'run_cpp_file "$1"' _ {}
end=$(date +%s.%N)
stop_monitor
trap - EXIT

count_q_values() {
  awk -F '\t' '
    FNR == 1 {
      q = 0
      for (i = 1; i <= NF; i++) if ($i == "q-value") q = i
      next
    }
    q && ($q + 0) < 0.01 { count++ }
    END { print count + 0 }
  ' "$@"
}

psm=$(count_q_values "$OUT"/*/target.psms.tsv)
peptide=$(count_q_values "$OUT"/*/target.peptides.tsv)
valid=$(find "$OUT" -mindepth 2 -maxdepth 2 -type f -name target.psms.tsv -size +0c | wc -l)
peak=$(cat "$peak_file")
awk -v start="$start" -v end="$end" -v peak="$peak" -v valid="$valid" \
  -v psm="$psm" -v peptide="$peptide" -v n="$N" '
  BEGIN {
    printf "implementation\tconcurrency\twall_seconds\tpeak_rss_kb\tvalid_files\tpsm_q_lt_0.01\tpeptide_q_lt_0.01\n"
    printf "cpp\t%d\t%.3f\t%d\t%d\t%d\t%d\n", n, end-start, peak, valid, psm, peptide
  }
' | tee "$OUT/summary.tsv"
