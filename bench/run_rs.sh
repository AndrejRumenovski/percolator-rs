#!/usr/bin/env bash
# Batch-run percolator-rs over all 65 PXD032157 files with a chosen execution profile.
#   usage: bash run_rs.sh [profile] [N]
#   profile = canonical (default) | balanced | fast ;  N = file-level concurrency (default 4)
set -u
PROFILE="${1:-canonical}"; N="${2:-4}"
ROOT="/run/media/andrej-rumenovski/New Volume/Code/percolator-rs"
BIN="$ROOT/target/release/percolator-rs"
IN="$ROOT/data/PXD032157"
export OUT="$HOME/percolator_rs_out/$PROFILE"; rm -rf "$OUT"; mkdir -p "$OUT"
PEAKF="/tmp/rs_run_peak"; echo 0 >"$PEAKF"; STOP="/tmp/rs_run_stop"; rm -f "$STOP"

( peak=0
  while [ ! -f "$STOP" ]; do
    s=$(ps --no-headers -o rss -C percolator-rs 2>/dev/null | awk '{t+=$1} END{print t+0}')
    [ "${s:-0}" -gt "$peak" ] && { peak=$s; echo "$peak" >"$PEAKF"; }
    sleep 0.1
  done ) & SAMP=$!

export BIN OUT PROFILE
run_one(){ local f="$1" b d; b=$(basename "$f" .pin); d="$OUT/$b"; mkdir -p "$d"
  "$BIN" "--$PROFILE" --seed 1 \
    --results-psms "$d/target.psms.tsv" --decoy-results-psms "$d/decoy.psms.tsv" \
    --results-peptides "$d/target.peptides.tsv" --decoy-results-peptides "$d/decoy.peptides.tsv" \
    "$f" 2>"$d/log"; }
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
printf 'profile=%-9s N=%s  WALL=%.1fs  PEAK_RSS=%.2fGiB  valid=%s/65  PSMq01=%s  pepq01=%s\n  out=%s\n' \
  "$PROFILE" "$N" "$wall" "$(echo "$peak/1048576"|bc -l)" "$ok" "$psm" "$pep" "$OUT"
