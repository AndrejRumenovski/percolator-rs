#!/usr/bin/env bash
# Fair SVM-versus-MLP comparison on the same PIN files. Both models use the
# canonical profile, seed, three folds, semi-supervised labels, and FDR code.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
BIN="$ROOT/target/release/percolator-rs"
INPUT_ROOT="${MODEL_BENCH_INPUT:-$ROOT/data/PXD032157}"
MODEL_BENCH_OUT="${MODEL_BENCH_OUT:-$HOME/percolator_rs_out/model-comparison}"
N="${MODEL_BENCH_JOBS:-4}"

cd "$ROOT"
[ -d "$INPUT_ROOT" ] || { echo "FAIL: input directory not found: $INPUT_ROOT"; exit 2; }
case "$MODEL_BENCH_OUT" in
  ""|/|"$HOME") echo "FAIL: unsafe MODEL_BENCH_OUT: $MODEL_BENCH_OUT"; exit 2 ;;
esac
cargo build --release
mkdir -p "$MODEL_BENCH_OUT"

mapfile -t inputs < <(find "$INPUT_ROOT" -type f -name '*.pin' -printf '%s\t%p\n' | sort -rn | cut -f2-)
[ "${#inputs[@]}" -gt 0 ] || { echo "FAIL: no PIN files under $INPUT_ROOT"; exit 2; }

run_model() {
  local model=$1 out="$MODEL_BENCH_OUT/$1" peak_file="$MODEL_BENCH_OUT/$1.peak"
  local stop_file="$MODEL_BENCH_OUT/$1.stop" start end monitor peak valid psm peptide
  rm -rf "$out"
  mkdir -p "$out"
  rm -f "$stop_file"
  echo 0 >"$peak_file"

  (
    local peak_rss=0 current
    while [ ! -f "$stop_file" ]; do
      current=$(ps --no-headers -o rss -C percolator-rs 2>/dev/null |
        awk '{sum+=$1} END{print sum+0}' || true)
      if [ "${current:-0}" -gt "$peak_rss" ]; then
        peak_rss=$current
        echo "$peak_rss" >"$peak_file"
      fi
      # Small SVM inputs can finish between 100 ms samples.
      sleep 0.02
    done
  ) &
  monitor=$!

  export BIN out model
  run_model_file() {
    local pin=$1 stem destination
    stem=$(basename "$pin" .pin)
    destination="$out/$stem"
    mkdir -p "$destination"
    "$BIN" --canonical --seed 1 --rescore-model "$model" \
      --results-psms "$destination/target.psms.tsv" \
      --decoy-results-psms "$destination/decoy.psms.tsv" \
      --results-peptides "$destination/target.peptides.tsv" \
      --decoy-results-peptides "$destination/decoy.peptides.tsv" \
      "$pin" 2>"$destination/log"
  }
  export -f run_model_file

  start=$(date +%s.%N)
  printf '%s\n' "${inputs[@]}" | xargs -P "$N" -I{} bash -c 'run_model_file "$1"' _ {}
  end=$(date +%s.%N)
  touch "$stop_file"
  wait "$monitor" 2>/dev/null || true

  peak=$(cat "$peak_file")
  valid=$(rg -l 'target PSMs q<0.01:' "$out"/*/log | wc -l)
  psm=$(rg -o 'target PSMs q<0.01: [0-9]+' "$out"/*/log | awk '{sum+=$NF} END{print sum+0}')
  peptide=$(rg -o 'target peptides q<0.01: [0-9]+' "$out"/*/log | awk '{sum+=$NF} END{print sum+0}')
  awk -v model="$model" -v files="$valid" -v psm="$psm" -v peptide="$peptide" \
      -v start="$start" -v end="$end" -v peak="$peak" \
      'BEGIN{printf "%s\t%d\t%d\t%d\t%.3f\t%d\n",model,files,psm,peptide,end-start,peak}' \
      >>"$MODEL_BENCH_OUT/summary.tsv"
}

printf 'model\tfiles\tpsm_q_lt_0.01\tpeptide_q_lt_0.01\twall_seconds\tpeak_rss_kb\n' \
  >"$MODEL_BENCH_OUT/summary.tsv"
run_model svm
run_model mlp

printf 'sample\tsvm_psm\tmlp_psm\tdelta_psm\tsvm_peptide\tmlp_peptide\tdelta_peptide\n' \
  >"$MODEL_BENCH_OUT/per-file.tsv"
for svm_log in "$MODEL_BENCH_OUT"/svm/*/log; do
  sample=$(basename "$(dirname "$svm_log")")
  mlp_log="$MODEL_BENCH_OUT/mlp/$sample/log"
  svm_psm=$(sed -n 's/.*target PSMs q<0.01: \([0-9]*\).*/\1/p' "$svm_log")
  mlp_psm=$(sed -n 's/.*target PSMs q<0.01: \([0-9]*\).*/\1/p' "$mlp_log")
  svm_peptide=$(sed -n 's/.*target peptides q<0.01: \([0-9]*\).*/\1/p' "$svm_log")
  mlp_peptide=$(sed -n 's/.*target peptides q<0.01: \([0-9]*\).*/\1/p' "$mlp_log")
  printf '%s\t%d\t%d\t%d\t%d\t%d\t%d\n' "$sample" \
    "$svm_psm" "$mlp_psm" "$((mlp_psm-svm_psm))" \
    "$svm_peptide" "$mlp_peptide" "$((mlp_peptide-svm_peptide))" \
    >>"$MODEL_BENCH_OUT/per-file.tsv"
done

cat "$MODEL_BENCH_OUT/summary.tsv"
echo "Per-file results: $MODEL_BENCH_OUT/per-file.tsv"
