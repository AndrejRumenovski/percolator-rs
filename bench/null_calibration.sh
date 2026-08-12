#!/usr/bin/env bash
# Null-calibration check: are the reported q-values honest, or does the
# class-weight grid search manufacture signal out of noise?
#
# Construction: keep ONLY the decoy rows of a real .pin and randomly relabel half
# of them as targets. Both classes are then drawn from the same null distribution,
# so every "target" identification is false by construction. A calibrated method
# must report ~0 targets at q<0.01; anything substantially above 0 is the
# anti-conservative bias of the pipeline, measured directly.
#
# Rows are split by a seeded RNG rather than by row order: the input holds ~5 Comet
# ranks per scan and rank correlates strongly with score, so an alternating split
# would stack the better-scoring ranks into one class and fake a real signal.
#
# Usage: bash bench/null_calibration.sh [n_files] [seed]
set -u
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT" || exit 2
BIN="$ROOT/target/release/percolator-rs"
IN="$ROOT/data/PXD032157"
N="${1:-3}"
SEED="${2:-42}"
# Keep scratch off the NTFS mount and off the quota-limited scratchpad.
WORK="${NULLCAL_WORK:-$HOME/percolator_rs_out/nullcal}"

[ -x "$BIN" ] || cargo build --release || { echo "FAIL: build"; exit 2; }
[ -d "$IN" ] || { echo "SKIP: dataset $IN not present"; exit 0; }
mkdir -p "$WORK" || exit 2

printf '%-52s %10s %10s %10s\n' "file" "null-PSMs" "selectC" "fixedC"
printf '%-52s %10s %10s %10s\n' "----" "---------" "-------" "------"

> "$WORK/.files"
find "$IN" -name '*.pin' -print0 | sort -z | head -z -n "$N" > "$WORK/.files"
while IFS= read -r -d '' f; do
  b=$(basename "$f" .pin)
  null="$WORK/$b.null.pin"

  # Header, then decoy rows only, with Label reassigned +1/-1 at random.
  awk -F'\t' -v OFS='\t' -v seed="$SEED" '
    NR==1 { print; next }
    NR==2 && $1 ~ /^DefaultDirection/ { next }
    {
      if (!init) { srand(seed); init = 1 }
      if ($2 + 0 < 0) { $2 = (rand() < 0.5) ? 1 : -1; print }
    }' "$f" > "$null"

  tot=$(awk -F'\t' 'NR>1{n++} END{print n+0}' "$null")
  pos=$(awk -F'\t' 'NR>1 && $2+0>0{n++} END{print n+0}' "$null")

  sel=$("$BIN" --canonical --seed 1 "$null" 2>&1 | sed -n 's/.*target PSMs q<0.01: \([0-9]*\).*/\1/p')
  fix=$("$BIN" --canonical --seed 1 --no-select-c "$null" 2>&1 | sed -n 's/.*target PSMs q<0.01: \([0-9]*\).*/\1/p')

  printf '%-52s %10s %10s %10s\n' "${b:0:52}" "$tot ($pos+)" "${sel:-ERR}" "${fix:-ERR}"
  rm -f "$null"
done < "$WORK/.files"
rm -f "$WORK/.files"

echo
echo "Every reported identification above is FALSE by construction."
echo "Calibrated behaviour is ~0; the counts are the realized false-discovery load at q<0.01."
