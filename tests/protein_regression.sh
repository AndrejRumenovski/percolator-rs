#!/usr/bin/env bash
# Synthetic protein-inference regression gate: generates a PIN with a large block
# of unpaired target-only seed proteins plus matched target/decoy protein pairs.
# The seed block makes the semi-supervised training converge reproducibly; the paired
# block then yields a deterministic picked-FDR sensitivity gain over classic TDA.
set -u
if [ ! -f Cargo.toml ] || [ ! -d tests ]; then
  cd "$(dirname "$0")/.." || exit 2
fi

BIN=target/release/percolator-rs
if [ ! -x "$BIN" ]; then
  echo "[build] $BIN missing, building release..."
  cargo build --release || { echo "FAIL: build"; exit 2; }
fi

tmpdir=$(mktemp -d)
trap 'rm -rf "$tmpdir"' EXIT
pin="$tmpdir/synthetic_proteins.pin"
tprot="$tmpdir/targets.proteins.tsv"
dprot="$tmpdir/decoys.proteins.tsv"
err="$tmpdir/run.err"

{
  printf 'SpecId\tLabel\tScanNr\tXcorr\tPeptide\tProteins\n'
  scan=1
  for i in $(seq 1 260); do
    score=$((5000 - i))
    printf 'seed_t%03d\t1\t%d\t%d.0\tK.SEEDT%03d.R\tSEEDPROT_%03d\n' "$i" "$scan" "$score" "$i" "$i"
    scan=$((scan + 1))
  done
  for i in $(seq 1 110); do
    target=$((2000 - i * 2))
    decoy=$((target - 2))
    printf 't%03d_a\t1\t%d\t%d.0\tK.TGT%03dA.R\tPROT_%03d\n' "$i" "$scan" "$target" "$i" "$i"
    scan=$((scan + 1))
    printf 't%03d_b\t1\t%d\t%d.8\tK.TGT%03dB.R\tPROT_%03d\n' "$i" "$scan" "$((target - 1))" "$i" "$i"
    scan=$((scan + 1))
    printf 'd%03d_a\t-1\t%d\t%d.5\tK.DEC%03dA.R\tDECOY_PROT_%03d\n' "$i" "$scan" "$decoy" "$i" "$i"
    scan=$((scan + 1))
    printf 'd%03d_b\t-1\t%d\t%d.3\tK.DEC%03dB.R\tDECOY_PROT_%03d\n' "$i" "$scan" "$((decoy - 1))" "$i" "$i"
    scan=$((scan + 1))
  done
} >"$pin"

"$BIN" --canonical --seed 1 --results-proteins "$tprot" --decoy-results-proteins "$dprot" "$pin" \
  >/dev/null 2>"$err" || { echo "FAIL: run errored"; cat "$err"; exit 1; }

groups=$(grep -oP 'protein groups: \K[0-9]+' "$err")
picked_entries=$(grep -oP 'picked entries: \K[0-9]+' "$err")
picked_q01=$(grep -oP 'target proteins q<0.01: \K[0-9]+' "$err")
classic_q01=$(grep -oP 'target proteins q<0.01: [0-9]+ \(picked-FDR\) vs \K[0-9]+' "$err")

fail=0
assert_eq() { # name value expected
  if [ "${2:-}" = "$3" ]; then
    printf '  PASS  %-22s %s\n' "$1" "$2"
  else
    printf '  FAIL  %-22s %s (expected %s)\n' "$1" "${2:-missing}" "$3"
    fail=1
  fi
}

echo "== percolator-rs protein regression gate (synthetic picked-FDR fixture) =="
assert_eq "protein groups" "${groups:-}" 480
assert_eq "picked entries" "${picked_entries:-}" 370
assert_eq "picked q<0.01" "${picked_q01:-}" 276
assert_eq "classic q<0.01" "${classic_q01:-}" 263
assert_eq "target rows" "$(wc -l < "$tprot")" 339
assert_eq "decoy rows" "$(wc -l < "$dprot")" 33

np1=$(awk 'NR > 1 && $5 == 1 { c++ } END { print c + 0 }' "$tprot")
np2=$(awk 'NR > 1 && $5 == 2 { c++ } END { print c + 0 }' "$tprot")
assert_eq "numPeptides=1 rows" "$np1" 260
assert_eq "numPeptides=2 rows" "$np2" 78

if [ "${picked_q01:-0}" -le "${classic_q01:-0}" ]; then
  echo "  FAIL  picked>classic          picked-FDR should strictly beat classic on this fixture"
  fail=1
else
  echo "  PASS  picked>classic          picked-FDR strictly beats classic on this fixture"
fi

if [ "$fail" -ne 0 ]; then
  echo "REGRESSION FAILED"
  exit 1
fi
echo "ALL CHECKS PASSED"
