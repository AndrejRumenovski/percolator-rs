#!/usr/bin/env bash
# Real-data picked-protein regression (requires the local single-organism F_3.pin).
# This complements tests/protein_regression.sh: the committed synthetic fixture
# checks exact invariants in hosted CI, while this checks that picked FDR improves
# sensitivity on a conventional protein-rich biological search.
set -u
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT" || exit 2

BIN="$ROOT/target/release/percolator-rs"
PIN="${PROTEIN_PIN:-$ROOT/data/F_3.pin}"
[ -x "$BIN" ] || cargo build --release || { echo "FAIL: build"; exit 2; }
[ -f "$PIN" ] || { echo "SKIP: single-organism PIN $PIN not present"; exit 0; }

WORK="$(mktemp -d "${TMPDIR:-/tmp}/percolator-protein-real.XXXXXX")" || exit 2
trap 'rm -rf "$WORK"' EXIT
LOG="$WORK/run.log"

"$BIN" --canonical --seed 1 --num-threads 3 \
  --results-proteins "$WORK/target.tsv" \
  --decoy-results-proteins "$WORK/decoy.tsv" \
  "$PIN" >/dev/null 2>"$LOG" || { cat "$LOG"; exit 1; }

line=$(grep 'target proteins q<0.01:' "$LOG" | tail -1)
groups=$(printf '%s\n' "$line" | sed -n 's/.*protein groups: \([0-9]*\).*/\1/p')
picked_entries=$(printf '%s\n' "$line" | sed -n 's/.*picked entries: \([0-9]*\).*/\1/p')
picked=$(printf '%s\n' "$line" | sed -n 's/.*q<0.01: \([0-9]*\) (picked-FDR).*/\1/p')
classic=$(printf '%s\n' "$line" | sed -n 's/.*vs \([0-9]*\) (classic).*/\1/p')

echo "== real single-organism picked-protein check ($(basename "$PIN")) =="
printf '  protein groups: %s; picked entries: %s\n' "${groups:-ERR}" "${picked_entries:-ERR}"
printf '  q<0.01: %s picked vs %s classic\n' "${picked:-ERR}" "${classic:-ERR}"

if [ -z "$picked" ] || [ -z "$classic" ] || [ "$picked" -le "$classic" ]; then
  echo "FAIL: picked-FDR did not strictly improve on classic TDA"
  cat "$LOG"
  exit 1
fi

echo "PASS: picked-FDR adds $((picked - classic)) protein groups at q<0.01"
