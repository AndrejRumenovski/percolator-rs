#!/usr/bin/env bash
# Compare fixed SVM defaults with leakage-free nested automatic selection.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
BIN="$ROOT/target/release/percolator-rs"
INPUT_ROOT="${SELECTION_BENCH_INPUT:-$ROOT/data/PXD032157}"
SELECTION_BENCH_OUT="${SELECTION_BENCH_OUT:-$HOME/percolator_rs_out/selection-comparison}"
JOBS="${SELECTION_BENCH_JOBS:-4}"

cd "$ROOT"
[ -d "$INPUT_ROOT" ] || { echo "FAIL: input directory not found: $INPUT_ROOT"; exit 2; }
case "$SELECTION_BENCH_OUT" in
  ""|/|"$HOME") echo "FAIL: unsafe SELECTION_BENCH_OUT: $SELECTION_BENCH_OUT"; exit 2 ;;
esac
cargo build --release
mkdir -p "$SELECTION_BENCH_OUT"
mapfile -t inputs < <(find "$INPUT_ROOT" -type f -name '*.pin' -printf '%s\t%p\n' | sort -rn | cut -f2-)
[ "${#inputs[@]}" -gt 0 ] || { echo "FAIL: no PIN files under $INPUT_ROOT"; exit 2; }

run_method() {
  local method=$1 out="$SELECTION_BENCH_OUT/$1" start end valid psm peptide
  rm -rf "$out"
  mkdir -p "$out"
  export BIN out method
  run_selection_file() {
    local pin=$1 stem destination selection_flag=()
    stem=$(basename "$pin" .pin)
    destination="$out/$stem"
    mkdir -p "$destination"
    [ "$method" = nested ] && selection_flag=(--auto-model)
    "$BIN" --canonical --seed 1 --rescore-model svm "${selection_flag[@]}" \
      --results-psms "$destination/target.psms.tsv" \
      --decoy-results-psms "$destination/decoy.psms.tsv" \
      --results-peptides "$destination/target.peptides.tsv" \
      --decoy-results-peptides "$destination/decoy.peptides.tsv" \
      "$pin" 2>"$destination/log"
  }
  export -f run_selection_file

  start=$(date +%s.%N)
  printf '%s\n' "${inputs[@]}" | xargs -P "$JOBS" -I{} bash -c 'run_selection_file "$1"' _ {}
  end=$(date +%s.%N)
  valid=$(rg -l 'target PSMs q<0.01:' "$out"/*/log | wc -l)
  psm=$(rg -o 'target PSMs q<0.01: [0-9]+' "$out"/*/log | awk '{sum+=$NF} END{print sum+0}')
  peptide=$(rg -o 'target peptides q<0.01: [0-9]+' "$out"/*/log | awk '{sum+=$NF} END{print sum+0}')
  awk -v method="$method" -v files="$valid" -v psm="$psm" -v peptide="$peptide" \
      -v start="$start" -v end="$end" \
      'BEGIN{printf "%s\t%d\t%d\t%d\t%.3f\n",method,files,psm,peptide,end-start}' \
      >>"$SELECTION_BENCH_OUT/summary.tsv"
}

printf 'method\tfiles\tpsm_q_lt_0.01\tpeptide_q_lt_0.01\twall_seconds\n' \
  >"$SELECTION_BENCH_OUT/summary.tsv"
run_method fixed
run_method nested

printf 'sample\tfixed_psm\tnested_psm\tdelta_psm\tfixed_peptide\tnested_peptide\tdelta_peptide\n' \
  >"$SELECTION_BENCH_OUT/per-file.tsv"
for fixed_log in "$SELECTION_BENCH_OUT"/fixed/*/log; do
  sample=$(basename "$(dirname "$fixed_log")")
  nested_log="$SELECTION_BENCH_OUT/nested/$sample/log"
  fixed_psm=$(sed -n 's/.*target PSMs q<0.01: \([0-9]*\).*/\1/p' "$fixed_log")
  nested_psm=$(sed -n 's/.*target PSMs q<0.01: \([0-9]*\).*/\1/p' "$nested_log")
  fixed_peptide=$(sed -n 's/.*target peptides q<0.01: \([0-9]*\).*/\1/p' "$fixed_log")
  nested_peptide=$(sed -n 's/.*target peptides q<0.01: \([0-9]*\).*/\1/p' "$nested_log")
  printf '%s\t%d\t%d\t%d\t%d\t%d\t%d\n' "$sample" \
    "$fixed_psm" "$nested_psm" "$((nested_psm-fixed_psm))" \
    "$fixed_peptide" "$nested_peptide" "$((nested_peptide-fixed_peptide))" \
    >>"$SELECTION_BENCH_OUT/per-file.tsv"
done

cat "$SELECTION_BENCH_OUT/summary.tsv"
echo "Per-file results: $SELECTION_BENCH_OUT/per-file.tsv"
