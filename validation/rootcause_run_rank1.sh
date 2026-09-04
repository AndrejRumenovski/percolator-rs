#!/usr/bin/env bash
set -uo pipefail
SP="/tmp/claude-1000/-run-media-andrej-rumenovski-New-Volume-Code-percolator-rs/9d5d4d63-9a67-4e88-ab50-b68d02980cac/scratchpad"
BIN="/run/media/andrej-rumenovski/New Volume/Code/percolator-rs/target/release/percolator-rs"
mkdir -p "$SP/rank1_out"
for seed in 1 2 3; do
for p in "$SP"/rank1/*.pin; do
  n=$(basename "$p" .pin)
  o="$SP/rank1_out/seed-$seed/$n"; mkdir -p "$o"
  [ -s "$o/target.tsv" ] && continue
  "$BIN" --canonical --no-select-c --seed $seed --num-threads 1 \
    --results-psms "$o/target.tsv" --decoy-results-psms "$o/decoy.tsv" "$p" \
    > "$o/stdout.log" 2> "$o/stderr.log"
  echo "rank1 seed=$seed done $n rc=$? $(date -Is)" >> "$SP/rank1_out/progress.log"
done
done
echo "RANK1 ALL DONE $(date -Is)" >> "$SP/rank1_out/progress.log"
