#!/usr/bin/env bash
# Fast 65-file percolator run: all of PXD032157 in <60s at minimal RAM.
# Default N=5 (~49s, ~1.2GiB). Use N=4 for min RAM (~59s, ~0.9GiB).
#   usage: bash fastrun.sh [N]
set -u
N="${1:-5}"
export LD_LIBRARY_PATH="$HOME/opt/perc-libs:${LD_LIBRARY_PATH:-}"
export PERC="$HOME/opt/percolator-root/usr/bin/percolator"
IN="/run/media/andrej-rumenovski/New Volume/Code/percolator-rs/data/PXD032157"
export OUT="$HOME/percolator_fast_out/PXD032157"   # local ext4 — NOT the slow ntfs-3g drive
rm -rf "$OUT"; mkdir -p "$OUT"
PEAKF="/tmp/perc_peak"; echo 0 >"$PEAKF"; STOP="/tmp/perc_stop"; rm -f "$STOP"

# size-interleaved order keeps concurrent total file-size (=> peak RAM) steady
ORDER=$(find "$IN" -name '*.pin' -printf '%s\t%p\n' | sort -rn | cut -f2- | awk '
  {a[NR]=$0} END{ i=1; j=NR; while(i<=j){ if(i==j){print a[i]} else {print a[i]; print a[j]} i++; j-- } }')

# peak-RSS sampler (throttled — do NOT remove the sleep)
( peak=0
  while [ ! -f "$STOP" ]; do
    s=$(ps --no-headers -o rss -C percolator 2>/dev/null | awk '{t+=$1} END{print t+0}')
    [ "${s:-0}" -gt "$peak" ] && { peak=$s; echo "$peak" >"$PEAKF"; }
    sleep 0.2
  done ) & SAMP=$!

run_one(){ local f="$1" b d; b=$(basename "$f" .pin); d="$OUT/$b"; mkdir -p "$d"
  "$PERC" --seed 1 --num-threads 1 --subset-max-train 20000 --maxiter 5 \
    --results-psms "$d/target.psms.tsv" --decoy-results-psms "$d/decoy.psms.tsv" \
    --results-peptides "$d/target.peptides.tsv" --decoy-results-peptides "$d/decoy.peptides.tsv" \
    "$f" >/dev/null 2>"$d/log"; }
export -f run_one

t0=$(date +%s.%N)
printf '%s\n' "$ORDER" | xargs -P "$N" -I{} bash -c 'run_one "$1"' _ {}
t1=$(date +%s.%N)
touch "$STOP"; wait "$SAMP" 2>/dev/null

wall=$(echo "$t1 - $t0" | bc); peak=$(cat "$PEAKF")
ok=$(grep -l "Processing took" "$OUT"/*/log 2>/dev/null | wc -l)
printf 'N=%s  WALL=%.1fs  PEAK_RSS=%.2fGiB  valid=%s/65  out=%s\n' \
  "$N" "$wall" "$(echo "$peak/1048576"|bc -l)" "$ok" "$OUT"
