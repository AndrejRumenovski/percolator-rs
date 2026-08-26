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
bprot="$tmpdir/bayesian.targets.tsv"
bdprot="$tmpdir/bayesian.decoys.tsv"
err="$tmpdir/run.err"
berr="$tmpdir/bayesian.err"

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
"$BIN" --canonical --seed 1 --protein-inference bayesian \
  --results-proteins "$bprot" --decoy-results-proteins "$bdprot" "$pin" \
  >/dev/null 2>"$berr" || { echo "FAIL: Bayesian run errored"; cat "$berr"; exit 1; }

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
assert_eq "picked q<0.01" "${picked_q01:-}" 267
assert_eq "classic q<0.01" "${classic_q01:-}" 262
assert_eq "target rows" "$(wc -l < "$tprot")" 336
assert_eq "decoy rows" "$(wc -l < "$dprot")" 36

np1=$(awk 'NR > 1 && $5 == 1 { c++ } END { print c + 0 }' "$tprot")
np2=$(awk 'NR > 1 && $5 == 2 { c++ } END { print c + 0 }' "$tprot")
assert_eq "numPeptides=1 rows" "$np1" 260
assert_eq "numPeptides=2 rows" "$np2" 75

assert_eq "Bayesian converged" "$(grep -oP 'converged: \K(true|false)' "$berr")" true
assert_eq "Bayesian header" "$(head -1 "$bprot")" \
  $'ProteinGroupId\tq-value\tposterior_error_prob\tscore\tnumPeptides\tproteinIds'
if ! awk -F'\t' 'NR > 1 {
  if ($2 < 0 || $2 > 1 || $3 < 0 || $3 > 1 || $4 < 0 || $4 > 1) {
    print "out-of-range probability at row " NR > "/dev/stderr"; exit 1
  }
  if (($3 + $4 - 1)^2 > 1e-8) {
    print "PEP + score != 1 at row " NR > "/dev/stderr"; exit 1
  }
  if ($2 + 1e-9 < previous) {
    print "q-value decreased at row " NR ": " previous " -> " $2 > "/dev/stderr"; exit 1
  }
  previous=$2
}' "$bprot"; then
  echo "  FAIL  Bayesian probabilities invalid or q-values non-monotone"
  fail=1
else
  echo "  PASS  Bayesian probabilities bounded and q-values monotone"
fi

# Picked-protein FDR estimates a cumulative error rate over protein groups and
# no protein-level posterior, so the column must read NA rather than carry a
# peptide-level PEP under a protein-level name. Requiring picked >= classic here
# would encode a sensitivity claim this fixture was built to produce; the counts
# above are recorded instead, in either direction.
if awk -F'\t' 'NR > 1 && $3 != "NA" { exit 1 }' "$tprot" &&
   awk -F'\t' 'NR > 1 && $3 != "NA" { exit 1 }' "$dprot"; then
  echo "  PASS  picked protein PEP      NA (no protein-level posterior is estimated)"
else
  echo "  FAIL  picked protein PEP      a value was reported where none is estimated"
  fail=1
fi

if [ "$fail" -ne 0 ]; then
  echo "REGRESSION FAILED"
  exit 1
fi
echo "ALL CHECKS PASSED"
