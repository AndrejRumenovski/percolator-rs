#!/usr/bin/env bash
# Reproduce the compact four-input extension to the existing PXD032157 benchmark.
# Large/generated artifacts stay on the local ext4 home volume, not in the repo.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
HERE="$ROOT/bench/multidataset"
WORK="${MULTIDATASET_OUT:-$HOME/percolator_rs_out/multidataset}"
INPUT="$WORK/inputs"
RUNS="$WORK/runs"
REPEATS="${REPEATS:-3}"
RS="${RS:-$ROOT/target/release/percolator-rs}"
PERC="${PERC:-$HOME/opt/percolator-root/usr/bin/percolator}"
SAGE_DIR="$WORK/tools/sage-0.14.7"
SAGE="$SAGE_DIR/sage-v0.14.7-x86_64-unknown-linux-gnu/sage"

mkdir -p "$INPUT" "$RUNS" "$SAGE_DIR"
export LC_ALL=C
export LD_LIBRARY_PATH="$HOME/opt/perc-libs:${LD_LIBRARY_PATH:-}"

expected_hash() {
  awk -v name="$1" '$2 == name { print $1 }' "$HERE/sources.sha256"
}

fetch() {
  local name="$1" url="$2" path expected
  path="$INPUT/$name"
  expected="$(expected_hash "$name")"
  if [[ -f "$path" ]] && [[ "$(sha256sum "$path" | cut -d' ' -f1)" == "$expected" ]]; then
    return
  fi
  echo "download $name"
  curl -fL --retry 5 --retry-all-errors -o "$path.download" "$url"
  echo "$expected  $path.download" | sha256sum --check --status
  mv -f "$path.download" "$path"
}

fetch hogrebe_tide.pin \
  https://raw.githubusercontent.com/wfondrie/mokapot/5bad097eabe528e17a6e8fb11f4d20cb5376ebb5/data/phospho_rep1.pin
fetch percolator_yeast.pin \
  https://raw.githubusercontent.com/percolator/percolator/8c7412ea0e556dddbc2dfa26c0641f49c948378e/data/percolator/tab/percolatorTab
fetch PXD020243_msfragger.pepXML \
  https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/10/PXD020243/MSB32231WmutBand_01.pepXML
fetch PXD060954_A1.mgf \
  https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/09/PXD060954/A1.mgf
fetch PXD060954_A_Breve_cont_decoy.fasta \
  https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/09/PXD060954/Amuciniphila_BAA835_BreveJCM1192_cont_decoy.fasta

python3 "$HERE/pepxml_to_pin.py" \
  "$INPUT/PXD020243_msfragger.pepXML" "$INPUT/PXD020243_msfragger.pin"

if [[ ! -x "$SAGE" ]]; then
  archive="$WORK/tools/sage-v0.14.7.tar.gz"
  expected="$(expected_hash sage-v0.14.7.tar.gz)"
  curl -fL --retry 5 --retry-all-errors -o "$archive.download" \
    https://github.com/lazear/sage/releases/download/v0.14.7/sage-v0.14.7-x86_64-unknown-linux-gnu.tar.gz
  echo "$expected  $archive.download" | sha256sum --check --status
  mv -f "$archive.download" "$archive"
  tar -xzf "$archive" -C "$SAGE_DIR"
fi

mkdir -p "$WORK/sage-pxd060954"
"$SAGE" "$HERE/sage-pxd060954.json" \
  --fasta "$INPUT/PXD060954_A_Breve_cont_decoy.fasta" \
  --output_directory "$WORK/sage-pxd060954" --write-pin \
  --disable-telemetry-i-dont-want-to-improve-sage \
  "$INPUT/PXD060954_A1.mgf" \
  >"$WORK/sage-pxd060954/stdout.log" 2>"$WORK/sage-pxd060954/stderr.log"
python3 "$HERE/normalize_sage_pin.py" \
  "$WORK/sage-pxd060954/results.sage.pin" "$INPUT/PXD060954_sage.pin"

(cd "$INPUT" && sha256sum --check "$HERE/generated.sha256")

[[ -x "$RS" ]] || cargo build --release --manifest-path "$ROOT/Cargo.toml"
[[ -x "$PERC" ]] || { echo "C++ Percolator not found: $PERC" >&2; exit 2; }

count_q01() {
  awk -F'\t' 'NR == 1 { for (i=1; i<=NF; i++) if ($i == "q-value") q=i; next }
    q && $q+0 <= 0.01 { n++ } END { print n+0 }' "$1"
}

benchmark() {
  local dataset="$1" implementation="$2" binary="$3" pin="$4"
  local out="$RUNS/$dataset/$implementation" wall rss psm peptides repeat
  local -a extra=()
  # Use direct target-decoy post-processing to match percolator-rs. C++'s
  # count-based auto-detector misclassifies the target-heavy MSFragger and Sage
  # concatenated searches as separate, which selects mix-max and emits warnings;
  # the legacy fixture no longer retains enough metadata to resolve its mode.
  [[ "$implementation" == "cpp" ]] && extra=(--search-input concatenated)
  mkdir -p "$out"
  : > "$out/times.tsv"
  for repeat in $(seq 1 "$REPEATS"); do
    /usr/bin/time -f '%e\t%M' -o "$out/time.$repeat.tsv" \
      "$binary" --seed 1 "${extra[@]}" \
        --results-psms "$out/target.psms.tsv" \
        --decoy-results-psms "$out/decoy.psms.tsv" \
        --results-peptides "$out/target.peptides.tsv" \
        --decoy-results-peptides "$out/decoy.peptides.tsv" \
        "$pin" >"$out/stdout.log" 2>"$out/stderr.log"
    cat "$out/time.$repeat.tsv" >> "$out/times.tsv"
  done
  wall="$(cut -f1 "$out/times.tsv" | sort -n | awk '{a[NR]=$1} END{print a[int((NR+1)/2)]}')"
  rss="$(cut -f2 "$out/times.tsv" | sort -n | awk '{a[NR]=$1} END{print a[int((NR+1)/2)]}')"
  psm="$(count_q01 "$out/target.psms.tsv")"
  peptides="$(count_q01 "$out/target.peptides.tsv")"
  printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$dataset" "$implementation" "$wall" "$rss" "$psm" "$peptides" \
    >> "$WORK/results.tsv"
}

printf 'dataset\timplementation\twall_seconds\tpeak_rss_kb\tpsm_q_le_0.01\tpeptide_q_le_0.01\n' \
  > "$WORK/results.tsv"
datasets=(hogrebe_tide PXD020243_msfragger PXD060954_sage percolator_yeast)
for dataset in "${datasets[@]}"; do
  pin="$INPUT/$dataset.pin"
  benchmark "$dataset" rust "$RS" "$pin"
  benchmark "$dataset" cpp "$PERC" "$pin"
done

(cd "$INPUT" && sha256sum \
  PXD020243_msfragger.pin PXD060954_sage.pin) \
  > "$WORK/generated.sha256"
{
  uname -a
  lscpu | sed -n 's/^Model name:[[:space:]]*/CPU: /p'
  printf 'percolator-rs commit: %s\n' "$(git -C "$ROOT" rev-parse HEAD)"
  printf 'C++ %s\n' "$("$PERC" -h 2>&1 | sed -n '1p')"
  "$SAGE" --version 2>&1 || true
} > "$WORK/environment.txt"

column -t -s $'\t' "$WORK/results.tsv" 2>/dev/null || cat "$WORK/results.tsv"
echo "outputs: $WORK"
