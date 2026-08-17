#!/usr/bin/env bash
# End-to-end signal-present calibration using foreign-proteome entrapment.
# Downloads are ~6.3 GB and are kept off the repository's NTFS volume.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
WORK="${ENTRAPMENT_WORK:-$HOME/percolator_rs_out/entrapment}"
BASE="https://ftp.pride.ebi.ac.uk/pride/data/archive/2022/04/PXD032157"
CRUX="${CRUX:-$WORK/crux-4.0.Linux.x86_64/bin/crux}"
PERC="${PERC:-$HOME/opt/percolator-root/usr/bin/percolator}"
RS="$ROOT/target/release/percolator-rs"
mkdir -p "$WORK"

stems=(
  "28May2015-QE-HF-Anopheles-22-atrium-S-24H-2nd-01"
  "28May2015-QE-HF-Anopheles-23-atrium-P-24H-2nd-01"
  "28May2015-QE-HF-Anopheles-38-MAGs-P-3rd-02"
  "09Dec2015-QEHF1-Anopheles-5-atrium-P-12hpm-3rd-01"
  "22Oct2014-Anopheles-8-MAGs-S-01"
  "9March2015-29-MAGs-pellet-2ndRep-14N-male-02"
)
declare -A mzml_sha1=(
  [28May2015-QE-HF-Anopheles-22-atrium-S-24H-2nd-01]=c0a4bb0b8128fad1f605a0292c4a0ea8492d597a
  [28May2015-QE-HF-Anopheles-23-atrium-P-24H-2nd-01]=2d2c0fbdc4cbe40587192d0431d2a6f08d98eff6
  [28May2015-QE-HF-Anopheles-38-MAGs-P-3rd-02]=c896949013bb57ab9fd75924032a29fc25e1720a
  [09Dec2015-QEHF1-Anopheles-5-atrium-P-12hpm-3rd-01]=6858f12eb939e065872844647e4e5f8e9cadc7b1
  [22Oct2014-Anopheles-8-MAGs-S-01]=d0be7aa253d84067fa50af2d09872799b344c0f7
  [9March2015-29-MAGs-pellet-2ndRep-14N-male-02]=97c0bf16dde3c778f593b0987f0ed33a8c432150
)

fetch() {
  local url=$1 out=$2 part="$2.part"
  [ -s "$out" ] && return
  curl -fL -C - --retry 5 "$url" -o "$part"
  mv "$part" "$out"
}

fetch "$BASE/AnogambiaeVB54_AnocoluzziiMaliVB54_SacCereUP2021_RNAseq2016-3FT_contam.fasta" "$WORK/native.fasta"
fetch "https://noble.gs.washington.edu/crux-downloads/crux-4.0/crux-4.0.Linux.x86_64.zip" "$WORK/crux.zip"
[ -x "$CRUX" ] || unzip -q "$WORK/crux.zip" -d "$WORK"
printf '%s  %s\n' 'befa6c4bd5614a923e0cf59c010119036577a97c' "$WORK/native.fasta" | sha1sum -c -

for stem in "${stems[@]}"; do
  fetch "$BASE/$stem.mzML" "$WORK/$stem.mzML"
  fetch "$BASE/$stem-comet.params.txt" "$WORK/$stem-comet.params.txt"
  printf '%s  %s\n' "${mzml_sha1[$stem]}" "$WORK/$stem.mzML" | sha1sum -c -
done

proteomes=(
  "arabidopsis UP000006548/UP000006548_3702"
  "rice UP000059680/UP000059680_39947"
  "maize UP000007305/UP000007305_4577"
  "soybean UP000008827/UP000008827_3847"
  "wheat UP000019116/UP000019116_4565"
  "tomato UP000004994/UP000004994_4081"
  "potato UP000011115/UP000011115_4113"
  "grape UP000009183/UP000009183_29760"
  "apple UP000290289/UP000290289_3750"
)
foreign=()
for item in "${proteomes[@]}"; do
  read -r name proteome_path <<<"$item"
  path="$WORK/$name.fasta.gz"
  fetch "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/$proteome_path.fasta.gz" "$path"
  gzip -t "$path"
  foreign+=("$path")
done
(cd "$WORK" && sha256sum -c "$ROOT/bench/entrapment/proteomes.sha256")

python3 "$ROOT/bench/entrapment/build_database.py" \
  "$WORK/native.fasta" "$WORK/combined.fasta" "${foreign[@]}" | tee "$WORK/database-stats.txt"
fraction=$(sed -n 's/^entrapment_fraction=//p' "$WORK/database-stats.txt")

cargo build --release --manifest-path "$ROOT/Cargo.toml"
export LD_LIBRARY_PATH="$HOME/opt/perc-libs:${LD_LIBRARY_PATH:-}"
report_inputs=()
for stem in "${stems[@]}"; do
  out="$WORK/comet-$stem"
  mkdir -p "$out"
  "$CRUX" comet --parameter-file "$WORK/$stem-comet.params.txt" \
    "$WORK/$stem.mzML" "$WORK/combined.fasta" --output-dir "$out"
  pin="$out/comet.pin"
  [ -s "$pin" ] || { echo "FAIL: Comet did not produce $pin"; exit 1; }

  "$RS" --canonical --seed 1 --results-psms "$WORK/rs.$stem.target.psms.tsv" \
    --decoy-results-psms "$WORK/rs.$stem.decoy.psms.tsv" "$pin" 2>"$WORK/rs.$stem.log"
  "$PERC" --seed 1 --num-threads 1 --results-psms "$WORK/cpp.$stem.target.psms.tsv" \
    --decoy-results-psms "$WORK/cpp.$stem.decoy.psms.tsv" "$pin" \
    >"$WORK/cpp.$stem.stdout" 2>"$WORK/cpp.$stem.log"
  report_inputs+=("$WORK/rs.$stem.target.psms.tsv" "$WORK/cpp.$stem.target.psms.tsv")
done

python3 "$ROOT/bench/entrapment/report.py" --entrapment-fraction "$fraction" \
  --output "$WORK/report.tsv" "${report_inputs[@]}"
echo "Artifacts: $WORK"
