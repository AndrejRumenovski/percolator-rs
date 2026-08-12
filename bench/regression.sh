#!/usr/bin/env bash
# Full-dataset performance + accuracy gate (requires the PXD032157 data locally;
# intended for a self-hosted / nightly runner, not hosted CI).
# Runs percolator-rs --canonical at N=4 over all 65 files and asserts wall/RSS
# budgets and aggregate q<0.01 yield within +/-1% of recorded reference.
set -u
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT" || exit 2
source tests/expected.env

BIN="$ROOT/target/release/percolator-rs"
IN="$ROOT/data/PXD032157"
[ -x "$BIN" ] || cargo build --release || { echo "FAIL: build"; exit 2; }
[ -d "$IN" ] || { echo "SKIP: dataset $IN not present (run on a machine with the data)"; exit 0; }

N=4
OUT="${TMPDIR:-/tmp}/perc_rs_reg"; rm -rf "$OUT"; mkdir -p "$OUT"
PEAKF="${TMPDIR:-/tmp}/reg_peak"; echo 0 >"$PEAKF"; STOP="${TMPDIR:-/tmp}/reg_stop"; rm -f "$STOP"
( peak=0; while [ ! -f "$STOP" ]; do
    s=$(ps --no-headers -o rss -C percolator-rs 2>/dev/null | awk '{t+=$1} END{print t+0}')
    [ "${s:-0}" -gt "$peak" ] && { peak=$s; echo "$peak" >"$PEAKF"; }; sleep 0.1
  done ) & SAMP=$!
export BIN OUT
run_one(){ local f="$1" b d; b=$(basename "$f" .pin); d="$OUT/$b"; mkdir -p "$d"
  "$BIN" --canonical --seed 1 --results-psms "$d/t" --results-peptides "$d/p" "$f" 2>"$d/log"; }
export -f run_one
ORDER=$(find "$IN" -name '*.pin' -printf '%s\t%p\n' | sort -rn | cut -f2- | awk '
  {a[NR]=$0} END{ i=1;j=NR; while(i<=j){ if(i==j){print a[i]} else {print a[i]; print a[j]} i++; j-- } }')
t0=$(date +%s.%N)
printf '%s\n' "$ORDER" | xargs -P "$N" -I{} bash -c 'run_one "$1"' _ {}
t1=$(date +%s.%N)
touch "$STOP"; wait "$SAMP" 2>/dev/null

wall=$(echo "$t1 - $t0" | bc); peak=$(cat "$PEAKF")
psm=$(grep -hoP 'target PSMs q<0.01: \K[0-9]+' "$OUT"/*/log | awk '{s+=$1} END{print s+0}')
pep=$(grep -hoP 'target peptides q<0.01: \K[0-9]+' "$OUT"/*/log | awk '{s+=$1} END{print s+0}')
ok=$(grep -l "target PSMs" "$OUT"/*/log | wc -l)

fail=0
aw(){ awk -v v="$2" -v e="$3" -v t="$4" -v n="$1" 'BEGIN{lo=e*(1-t/100);hi=e*(1+t/100);
  if(v+0>=lo&&v+0<=hi){printf"  PASS  %-18s %s (ref %s +/-%s%%)\n",n,v,e,t}
  else{printf"  FAIL  %-18s %s (ref %s allowed %.1f..%.1f)\n",n,v,e,lo,hi;exit 1}}'||return 1; }
le(){ awk -v v="$2" -v b="$3" -v n="$1" -v u="$4" 'BEGIN{
  if(v+0<=b+0){printf"  PASS  %-18s %s%s (budget %s%s)\n",n,v,u,b,u}
  else{printf"  FAIL  %-18s %s%s (budget %s%s)\n",n,v,u,b,u;exit 1}}'||return 1; }

echo "== percolator-rs full-dataset gate (canonical, N=$N, 65 files) =="
le "valid runs"     "$ok"  "65" ""  && [ "$ok" -eq 65 ] || { echo "  FAIL  only $ok/65 valid"; fail=1; }
aw "PSM q<0.01"     "$psm" "$FULL_PSM_Q01" "$YIELD_TOLERANCE_PCT" || fail=1
aw "peptide q<0.01" "$pep" "$FULL_PEP_Q01" "$YIELD_TOLERANCE_PCT" || fail=1
le "wall time"      "$wall" "$FULL_TIME_BUDGET_S" "s"  || fail=1
le "peak RSS"       "$peak" "$FULL_MEM_BUDGET_KB" "kB" || fail=1
rm -rf "$OUT"
[ "$fail" -eq 0 ] && { echo "ALL CHECKS PASSED"; exit 0; } || { echo "REGRESSION FAILED"; exit 1; }
