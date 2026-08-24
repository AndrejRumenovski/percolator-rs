#!/usr/bin/env bash
# Deterministic benchmarks for the README's retention-time, pooled-training,
# and intra-file threading claims on PXD032157.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
BIN="${ADVANCED_BENCH_BIN:-$ROOT/target/release/percolator-rs}"
INPUT_ROOT="${ADVANCED_BENCH_INPUT:-$ROOT/data/PXD032157}"
OUT="${ADVANCED_BENCH_OUT:-$HOME/percolator_rs_out/advanced-features}"
REPEATS="${REPEATS:-3}"
EXPECTED_FILES="${ADVANCED_BENCH_EXPECTED_FILES:-65}"
export LC_ALL=C

cd "$ROOT"
[ -d "$INPUT_ROOT" ] || { echo "FAIL: input directory not found: $INPUT_ROOT"; exit 2; }
case "$OUT" in
  ""|/|"$HOME") echo "FAIL: unsafe ADVANCED_BENCH_OUT: $OUT"; exit 2 ;;
esac
[[ "$REPEATS" =~ ^[1-9][0-9]*$ ]] || { echo "FAIL: REPEATS must be positive"; exit 2; }
[[ "$EXPECTED_FILES" =~ ^[1-9][0-9]*$ ]] || { echo "FAIL: ADVANCED_BENCH_EXPECTED_FILES must be positive"; exit 2; }
[ -x "$BIN" ] || cargo build --release

mapfile -d '' -t lexical_inputs < <(find "$INPUT_ROOT" -maxdepth 1 -type f -name '*.pin' -print0 | sort -z)
mapfile -t size_inputs < <(find "$INPUT_ROOT" -maxdepth 1 -type f -name '*.pin' -printf '%s\t%p\n' | sort -n | cut -f2-)
[ "${#lexical_inputs[@]}" -eq "$EXPECTED_FILES" ] || {
  echo "FAIL: expected $EXPECTED_FILES PIN files, found ${#lexical_inputs[@]}"; exit 2;
}

rm -rf "$OUT/runs"
mkdir -p "$OUT/runs"
printf 'case\trepeat\twall_seconds\tpeak_rss_kb\tpsm_q_lt_0.01\tpeptide_q_lt_0.01\tpsm_sha256\tpeptide_sha256\n' >"$OUT/raw.tsv"

median_column() {
  local case_name=$1 column=$2
  awk -F'\t' -v c="$case_name" -v column="$column" '$1 == c {print $column}' "$OUT/raw.tsv" |
    sort -n | awk '{value[NR]=$1} END {if (NR) print value[int((NR+1)/2)]}'
}

case_value() {
  local case_name=$1 column=$2
  awk -F'\t' -v c="$case_name" -v column="$column" '$1 == c {print $column; exit}' "$OUT/raw.tsv"
}

run_case() {
  local case_name=$1
  shift
  local destination="$OUT/runs/$case_name" repeat wall rss psm peptide psm_sha peptide_sha
  mkdir -p "$destination"
  for repeat in $(seq 1 "$REPEATS"); do
    /usr/bin/time -f '%e\t%M' -o "$destination/time.$repeat.tsv" \
      "$BIN" --canonical --seed 1 "$@" \
        --results-psms "$destination/target.psms.$repeat.tsv" \
        --decoy-results-psms "$destination/decoy.psms.$repeat.tsv" \
        --results-peptides "$destination/target.peptides.$repeat.tsv" \
        --decoy-results-peptides "$destination/decoy.peptides.$repeat.tsv" \
        >"$destination/stdout.$repeat.log" 2>"$destination/stderr.$repeat.log"
    IFS=$'\t' read -r wall rss <"$destination/time.$repeat.tsv"
    psm=$(sed -n 's/.*target PSMs q<0.01: \([0-9]*\).*/\1/p' "$destination/stderr.$repeat.log")
    peptide=$(sed -n 's/.*target peptides q<0.01: \([0-9]*\).*/\1/p' "$destination/stderr.$repeat.log")
    [ -n "$psm" ] && [ -n "$peptide" ] || { echo "FAIL: missing yield for $case_name repeat $repeat"; exit 1; }
    psm_sha=$(sha256sum "$destination/target.psms.$repeat.tsv" | cut -d' ' -f1)
    peptide_sha=$(sha256sum "$destination/target.peptides.$repeat.tsv" | cut -d' ' -f1)
    printf '%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\n' \
      "$case_name" "$repeat" "$wall" "$rss" "$psm" "$peptide" "$psm_sha" "$peptide_sha" >>"$OUT/raw.tsv"
  done
  [ "$(awk -F'\t' -v c="$case_name" '$1 == c {print $5 FS $6 FS $7 FS $8}' "$OUT/raw.tsv" | sort -u | wc -l)" -eq 1 ] || {
    echo "FAIL: nondeterministic result for $case_name"; exit 1;
  }
}

printf 'selection_rule\tordinal\tbytes\tsha256\tfile\n' >"$OUT/selected-inputs.tsv"
for ordinal in 0 1 2; do
  pin=${lexical_inputs[$ordinal]}
  printf 'lexicographically_first_three\t%d\t%s\t%s\t%s\n' "$((ordinal+1))" \
    "$(stat -c %s "$pin")" "$(sha256sum "$pin" | cut -d' ' -f1)" "$(basename "$pin")" >>"$OUT/selected-inputs.tsv"
done
for ordinal in 0 1 2 3; do
  pin=${size_inputs[$ordinal]}
  printf 'four_smallest_by_bytes\t%d\t%s\t%s\t%s\n' "$((ordinal+1))" \
    "$(stat -c %s "$pin")" "$(sha256sum "$pin" | cut -d' ' -f1)" "$(basename "$pin")" >>"$OUT/selected-inputs.tsv"
done
largest=${size_inputs[$((${#size_inputs[@]}-1))]}
printf 'largest_by_bytes\t1\t%s\t%s\t%s\n' \
  "$(stat -c %s "$largest")" "$(sha256sum "$largest" | cut -d' ' -f1)" "$(basename "$largest")" >>"$OUT/selected-inputs.tsv"

# Retention time: a prospective, filename-only sample rather than cherry-picked outcomes.
for ordinal in 0 1 2; do
  pin=${lexical_inputs[$ordinal]}
  run_case "rt$((ordinal+1))-baseline" "$pin"
  run_case "rt$((ordinal+1))-features" --rt-features "$pin"
done
printf 'ordinal\tfile\tbaseline_psm\trt_psm\tdelta_psm\tdelta_percent\tbaseline_peptide\trt_peptide\tbaseline_wall_seconds\trt_wall_seconds\n' >"$OUT/retention-time.tsv"
for ordinal in 0 1 2; do
  pin=${lexical_inputs[$ordinal]}
  baseline="rt$((ordinal+1))-baseline"
  features="rt$((ordinal+1))-features"
  bpsm=$(case_value "$baseline" 5); rpsm=$(case_value "$features" 5)
  bpep=$(case_value "$baseline" 6); rpep=$(case_value "$features" 6)
  awk -v ordinal="$((ordinal+1))" -v file="$(basename "$pin")" \
      -v bpsm="$bpsm" -v rpsm="$rpsm" -v bpep="$bpep" -v rpep="$rpep" \
      -v bwall="$(median_column "$baseline" 3)" -v rwall="$(median_column "$features" 3)" \
      'BEGIN {pct=bpsm ? 100*(rpsm-bpsm)/bpsm : 0; printf "%d\t%s\t%d\t%d\t%+d\t%+.3f\t%d\t%d\t%.3f\t%.3f\n", ordinal,file,bpsm,rpsm,rpsm-bpsm,pct,bpep,rpep,bwall,rwall}' \
      >>"$OUT/retention-time.tsv"
done

# Joint training: "small" is defined mechanically as the four smallest PINs.
join_inputs=("${size_inputs[@]:0:4}")
for ordinal in 0 1 2 3; do
  run_case "join$((ordinal+1))-standalone" "${join_inputs[$ordinal]}"
done
run_case join-pooled --join "${join_inputs[@]}"
printf 'ordinal\tfile\tbytes\tstandalone_psm\tpooled_psm\tdelta_psm\n' >"$OUT/joint-training.tsv"
join_log="$OUT/runs/join-pooled/stderr.1.log"
standalone_total=0
pooled_source_total=0
for ordinal in 0 1 2 3; do
  pin=${join_inputs[$ordinal]}
  file=$(basename "$pin")
  standalone=$(case_value "join$((ordinal+1))-standalone" 5)
  pooled=$(awk -v prefix="  [$file] " 'index($0, prefix) == 1 {print $NF}' "$join_log")
  [ -n "$pooled" ] || { echo "FAIL: no pooled per-file yield for $file"; exit 1; }
  printf '%d\t%s\t%s\t%d\t%d\t%+d\n' "$((ordinal+1))" "$file" "$(stat -c %s "$pin")" \
    "$standalone" "$pooled" "$((pooled-standalone))" >>"$OUT/joint-training.tsv"
  standalone_total=$((standalone_total+standalone))
  pooled_source_total=$((pooled_source_total+pooled))
done
pooled_total=$(case_value join-pooled 5)
[ "$pooled_source_total" -eq "$pooled_total" ] || {
  echo "FAIL: pooled per-source total $pooled_source_total differs from reported total $pooled_total"; exit 1;
}
printf 'aggregate\tall_four\tNA\t%d\t%d\t%+d\n' "$standalone_total" "$pooled_total" \
  "$((pooled_total-standalone_total))" >>"$OUT/joint-training.tsv"

# Intra-file scaling: largest input, fixed outputs, and identical result hashes.
run_case threads-fixed-1 --num-threads 1 "$largest"
run_case threads-fixed-3 --num-threads 3 "$largest"
run_case threads-select-1 --select-c --num-threads 1 "$largest"
run_case threads-select-9 --select-c --num-threads 9 "$largest"
printf 'mode\tthreads\tfile\tmedian_wall_seconds\tmedian_peak_rss_kb\tpsm_q_lt_0.01\tpeptide_q_lt_0.01\tpsm_sha256\tpeptide_sha256\n' >"$OUT/threading.tsv"
for spec in 'fixed 1 threads-fixed-1' 'fixed 3 threads-fixed-3' 'select-c 1 threads-select-1' 'select-c 9 threads-select-9'; do
  read -r mode threads case_name <<<"$spec"
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$mode" "$threads" "$(basename "$largest")" \
    "$(median_column "$case_name" 3)" "$(median_column "$case_name" 4)" \
    "$(case_value "$case_name" 5)" "$(case_value "$case_name" 6)" \
    "$(case_value "$case_name" 7)" "$(case_value "$case_name" 8)" >>"$OUT/threading.tsv"
done
[ "$(awk -F'\t' 'NR>1 && $1=="fixed" {print $8 FS $9}' "$OUT/threading.tsv" | sort -u | wc -l)" -eq 1 ] || {
  echo 'FAIL: fixed-weight output differs by thread count'; exit 1;
}
[ "$(awk -F'\t' 'NR>1 && $1=="select-c" {print $8 FS $9}' "$OUT/threading.tsv" | sort -u | wc -l)" -eq 1 ] || {
  echo 'FAIL: select-C output differs by thread count'; exit 1;
}

{
  uname -a
  lscpu | sed -n 's/^Model name:[[:space:]]*/CPU: /p'
  printf 'percolator-rs commit: %s\n' "$(git rev-parse HEAD)"
  rustc --version
} >"$OUT/environment.txt"

cat "$OUT/retention-time.tsv"
cat "$OUT/joint-training.tsv"
cat "$OUT/threading.tsv"
echo "outputs: $OUT"
