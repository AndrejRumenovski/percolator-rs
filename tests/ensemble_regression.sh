#!/usr/bin/env bash
# End-to-end ensemble smoke test using two identical engine views. This is a
# structural/determinism test, not a biological ensemble-yield benchmark.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
BIN="$ROOT/target/release/percolator-rs"
FIX="$ROOT/tests/fixtures/sample.pin"
WORK=$(mktemp -d "${TMPDIR:-/tmp}/percolator-rs-ensemble.XXXXXX")
trap 'rm -rf "$WORK"' EXIT

[ -x "$BIN" ] || cargo build --release --manifest-path "$ROOT/Cargo.toml"
[ -f "$FIX" ] || { echo "FAIL: fixture missing: $FIX"; exit 2; }

run_ensemble() {
  local name=$1 threads=$2 destination="$WORK/$1"
  mkdir -p "$destination"
  "$BIN" --canonical --seed 1 --num-threads "$threads" --ensemble \
    "comet=$FIX" "tide=$FIX" \
    --results-psms "$destination/target.psms.tsv" \
    --decoy-results-psms "$destination/decoy.psms.tsv" \
    --results-peptides "$destination/target.peptides.tsv" \
    --decoy-results-peptides "$destination/decoy.peptides.tsv" \
    >"$destination/stdout.log" 2>"$destination/stderr.log"
}

run_ensemble serial 1
run_ensemble parallel 3
for file in target.psms.tsv decoy.psms.tsv target.peptides.tsv decoy.peptides.tsv; do
  cmp "$WORK/serial/$file" "$WORK/parallel/$file"
done

# The two inputs are identical, so exact candidate deduplication must reduce the
# emitted PSM set to the number of distinct (ScanNr, Label, Peptide) candidates.
expected=$(awk -F'\t' '
  NR==1 {for(i=1;i<=NF;i++){if($i=="ScanNr") scan=i; if($i=="Label") label=i; if($i=="Peptide") peptide=i} next}
  $1 !~ /^DefaultDirection/ {seen[$scan SUBSEP $label SUBSEP $peptide]=1}
  END {for(key in seen) n++; print n+0}
' "$FIX")
observed=$(( $(wc -l <"$WORK/serial/target.psms.tsv") + $(wc -l <"$WORK/serial/decoy.psms.tsv") - 2 ))
[ "$observed" -eq "$expected" ] || {
  echo "FAIL: ensemble emitted $observed candidates, expected $expected after deduplication"; exit 1;
}
awk -F'\t' 'FNR>1 && $1 !~ /^(comet|tide):/ {bad=1} END{exit bad}' \
  "$WORK/serial/target.psms.tsv" "$WORK/serial/decoy.psms.tsv"
rg -q 'ensemble from 2 engines' "$WORK/serial/stderr.log"

if "$BIN" --ensemble "comet=$FIX" "tide=$FIX" \
    --results-proteins "$WORK/proteins.tsv" >/dev/null 2>&1; then
  echo 'FAIL: ensemble mode accepted protein inference output'
  exit 1
fi

echo "PASS: ensemble candidate deduplication and serial/parallel outputs are exact ($observed candidates)"
