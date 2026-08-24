#!/usr/bin/env bash
# Reproduce protein-level calibration on the PrEST homology standard PXD008425.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
HERE="$ROOT/bench/protein_calibration"
WORK="${PROTEIN_CALIBRATION_OUT:-$HOME/percolator_rs_out/protein-calibration}"
BASE="https://ftp.pride.ebi.ac.uk/pride/data/archive/2017/12/PXD008425"
RAW="$WORK/raw"
MZML_ROOT="$WORK/mzml"
TOOLS="$WORK/tools"
INPUT="$WORK/input"
RESULTS="$WORK/results"
RS="${RS:-$ROOT/target/release/percolator-rs}"
DOWNLOAD_JOBS="${DOWNLOAD_JOBS:-4}"
REPEATS="${REPEATS:-3}"

TRFP_ARCHIVE="$TOOLS/ThermoRawFileParser-v.2.0.0-dev-linux.zip"
TRFP_DIR="$TOOLS/ThermoRawFileParser-v.2.0.0-dev-linux"
TRFP="$TRFP_DIR/ThermoRawFileParser"
SAGE_ARCHIVE="$TOOLS/sage-v0.14.7-x86_64-unknown-linux-gnu.tar.gz"
SAGE_DIR="$TOOLS/sage-0.14.7"
SAGE="$SAGE_DIR/sage-v0.14.7-x86_64-unknown-linux-gnu/sage"

mkdir -p "$RAW" "$MZML_ROOT" "$TOOLS" "$INPUT" "$RESULTS"
export LC_ALL=C

expected_hash() {
  local manifest="$1" name="$2"
  awk -v name="$name" '$2 == name { print $1 }' "$manifest"
}

fetch_sha1() {
  local name="$1" destination="$2" expected
  expected="$(expected_hash "$HERE/sources.sha1" "$name")"
  if [[ -s "$destination" ]] && [[ "$(sha1sum "$destination" | cut -d' ' -f1)" == "$expected" ]]; then
    return
  fi
  echo "download $name"
  curl --fail --location --continue-at - --silent --show-error \
    --retry 5 --retry-all-errors -o "$destination.part" "$BASE/$name"
  printf '%s  %s\n' "$expected" "$destination.part" | sha1sum --check --status
  mv -f "$destination.part" "$destination"
}

fetch_sha256() {
  local name="$1" url="$2" destination="$3" expected
  expected="$(expected_hash "$HERE/tools.sha256" "$name")"
  if [[ -s "$destination" ]] && [[ "$(sha256sum "$destination" | cut -d' ' -f1)" == "$expected" ]]; then
    return
  fi
  echo "download $name"
  curl --fail --location --continue-at - --silent --show-error \
    --retry 5 --retry-all-errors -o "$destination.part" "$url"
  printf '%s  %s\n' "$expected" "$destination.part" | sha256sum --check --status
  mv -f "$destination.part" "$destination"
}

raw_files=(
  mixtureArep1.raw mixtureArep2.raw mixtureArep3.raw
  mixtureBrep1.raw mixtureBrep2.raw mixtureBrep3.raw
  mixtureABrep1.raw mixtureABrep2.raw mixtureABrep3.raw
  blankrep1.raw blankrep2.raw blankrep3.raw
)
for name in prest_pool_a.fasta prest_pool_b.fasta prest_1000_random.fasta; do
  fetch_sha1 "$name" "$INPUT/$name"
done
export BASE HERE RAW
export -f expected_hash fetch_sha1
printf '%s\0' "${raw_files[@]}" | xargs -0 -P "$DOWNLOAD_JOBS" -n 1 \
  bash -euo pipefail -c 'fetch_sha1 "$1" "$RAW/$1"' _

fetch_sha256 \
  ThermoRawFileParser-v.2.0.0-dev-linux.zip \
  https://github.com/CompOmics/ThermoRawFileParser/releases/download/v.2.0.0-dev/ThermoRawFileParser-v.2.0.0-dev-linux.zip \
  "$TRFP_ARCHIVE"
fetch_sha256 \
  sage-v0.14.7-x86_64-unknown-linux-gnu.tar.gz \
  https://github.com/lazear/sage/releases/download/v0.14.7/sage-v0.14.7-x86_64-unknown-linux-gnu.tar.gz \
  "$SAGE_ARCHIVE"

if [[ ! -x "$TRFP" ]]; then
  mkdir -p "$TRFP_DIR"
  unzip -q -o "$TRFP_ARCHIVE" -d "$TRFP_DIR"
  chmod +x "$TRFP"
fi
if [[ ! -x "$SAGE" ]]; then
  mkdir -p "$SAGE_DIR"
  tar -xzf "$SAGE_ARCHIVE" -C "$SAGE_DIR"
fi

python3 "$HERE/build_database.py" \
  --pool-a "$INPUT/prest_pool_a.fasta" \
  --pool-b "$INPUT/prest_pool_b.fasta" \
  --random "$INPUT/prest_1000_random.fasta" \
  --database "$INPUT/prest-target-decoy.fasta" \
  --truth "$INPUT/ground-truth.tsv" \
  > "$INPUT/database-stats.txt"

converter_key="$(sha256sum "$TRFP_ARCHIVE" | cut -c1-12)"
MZML="$MZML_ROOT/$converter_key"
mkdir -p "$MZML"
for raw in "${raw_files[@]}"; do
  stem="${raw%.raw}"
  if [[ ! -s "$MZML/$stem.mzML" ]]; then
    echo "convert $raw"
    "$TRFP" -i="$RAW/$raw" -o="$MZML" -f=2 -L=2 \
      >"$MZML/$stem.convert.stdout" 2>"$MZML/$stem.convert.stderr"
  fi
done

pipeline_key="$({
  sha256sum "$HERE/sage-prest.json" "$HERE/normalize_sage_pin.py" \
    "$INPUT/prest-target-decoy.fasta" "$SAGE_ARCHIVE"
  printf 'converter=%s\n' "$converter_key"
} | sha256sum | cut -c1-16)"
SEARCH="$WORK/search/$pipeline_key"
RUNS="$WORK/runs/$pipeline_key"
mkdir -p "$SEARCH" "$RUNS"

declare -A vial_for=(
  [mixtureA]=A [mixtureB]=B [mixtureAB]=AB [blank]=BLANK
)
printf 'sample\tvial\treplicate\tsplit\tpin\n' > "$WORK/manifest.tsv"
for prefix in mixtureA mixtureB mixtureAB blank; do
  for replicate in 1 2 3; do
    sample="${prefix}rep${replicate}"
    search="$SEARCH/$sample"
    pin="$INPUT/$pipeline_key/$sample.pin"
    mkdir -p "$search" "$(dirname "$pin")"
    if [[ ! -s "$pin" ]]; then
      echo "search $sample"
      "$SAGE" "$HERE/sage-prest.json" \
        --fasta "$INPUT/prest-target-decoy.fasta" \
        --output_directory "$search" --write-pin \
        --disable-telemetry-i-dont-want-to-improve-sage \
        "$MZML/$sample.mzML" \
        >"$search/stdout.log" 2>"$search/stderr.log"
      python3 "$HERE/normalize_sage_pin.py" "$search/results.sage.pin" "$pin"
    fi
    case "$replicate" in
      1) split=calibration ;;
      2) split=validation ;;
      3) split=test ;;
    esac
    printf '%s\t%s\t%s\t%s\t%s\n' \
      "$sample" "${vial_for[$prefix]}" "$replicate" "$split" "$pin" \
      >> "$WORK/manifest.tsv"
  done
done

[[ -x "$RS" ]] || cargo build --release --manifest-path "$ROOT/Cargo.toml"
SELECTION="$RESULTS/$pipeline_key/selection"
python3 "$HERE/select_params.py" \
  --binary "$RS" --truth "$INPUT/ground-truth.tsv" \
  --manifest "$WORK/manifest.tsv" --output-dir "$SELECTION"
# This file is generated by select_params.py and contains only three numeric assignments.
# shellcheck disable=SC1090
source "$SELECTION/selected-params.env"

run_method() {
  local sample="$1" pin="$2" method="$3" run="$RUNS/$sample/$method"
  local repeat wall rss
  local -a protein_args=()
  case "$method" in
    picked)
      protein_args=(--protein-inference picked)
      ;;
    bayes-fixed)
      protein_args=(--protein-inference bayesian)
      ;;
    bayes-selected)
      protein_args=(
        --protein-inference bayesian
        --protein-alpha "$PROTEIN_ALPHA"
        --protein-beta "$PROTEIN_BETA"
        --protein-gamma "$PROTEIN_GAMMA"
        --protein-max-iter "$PROTEIN_MAX_ITER"
      )
      ;;
  esac
  mkdir -p "$run"
  : > "$run/times.tsv"
  for repeat in $(seq 1 "$REPEATS"); do
    /usr/bin/time -f '%e\t%M' -o "$run/time.$repeat.tsv" \
      "$RS" --canonical --seed 1 "${protein_args[@]}" \
        --results-proteins "$run/target.tsv" \
        --decoy-results-proteins "$run/decoy.tsv" \
        "$pin" >"$run/stdout.log" 2>"$run/stderr.log"
    cat "$run/time.$repeat.tsv" >> "$run/times.tsv"
  done
  wall="$(cut -f1 "$run/times.tsv" | sort -n | awk '{a[NR]=$1} END{print a[int((NR+1)/2)]}')"
  rss="$(cut -f2 "$run/times.tsv" | sort -n | awk '{a[NR]=$1} END{print a[int((NR+1)/2)]}')"
  printf '%s\t%s\n' "$wall" "$rss" > "$run/time.tsv"
}

while IFS=$'\t' read -r sample _vial _replicate _split pin; do
  [[ "$sample" == sample ]] && continue
  for method in picked bayes-fixed bayes-selected; do
    echo "infer $sample $method"
    run_method "$sample" "$pin" "$method"
  done
done < "$WORK/manifest.tsv"

REPORT="$RESULTS/$pipeline_key/report"
python3 "$HERE/report.py" \
  --truth "$INPUT/ground-truth.tsv" --manifest "$WORK/manifest.tsv" \
  --runs "$RUNS" --output-dir "$REPORT"

(cd "$INPUT/$pipeline_key" && sha256sum ./*.pin) > "$RESULTS/$pipeline_key/generated-pins.sha256"
{
  uname -a
  lscpu | sed -n 's/^Model name:[[:space:]]*/CPU: /p'
  printf 'percolator-rs commit: %s\n' "$(git -C "$ROOT" rev-parse HEAD)"
  "$TRFP" --version 2>&1
  "$SAGE" --version 2>&1
  printf 'pipeline key: %s\n' "$pipeline_key"
  printf 'converter key: %s\n' "$converter_key"
} > "$RESULTS/$pipeline_key/environment.txt"

awk -F'\t' 'NR == 1 || ($1 == "ALL" && ($3 == "validation" || $3 == "test"))' \
  "$REPORT/summary.tsv" | column -t -s $'\t' 2>/dev/null || true
echo "selected parameters: $SELECTION/selected-params.json"
echo "report: $REPORT"
