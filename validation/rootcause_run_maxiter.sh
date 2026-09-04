#!/usr/bin/env bash
set -uo pipefail
SP="/tmp/claude-1000/-run-media-andrej-rumenovski-New-Volume-Code-percolator-rs/9d5d4d63-9a67-4e88-ab50-b68d02980cac/scratchpad"
BIN="/run/media/andrej-rumenovski/New Volume/Code/percolator-rs/target/release/percolator-rs"
ROOT="$HOME/percolator_rs_out/entrapment"
for mi in 0 1 2 3 5 10; do
for seed in 1 2 3; do
for d in "$ROOT"/comet-2* "$ROOT"/comet-0* "$ROOT"/comet-9*; do
  n=$(basename "$d"); [ "$n" = "comet-out" ] && continue
  o="$SP/maxiter/mi-$mi/seed-$seed/$n"; mkdir -p "$o"
  [ -s "$o/target.tsv" ] && continue
  "$BIN" --canonical --no-select-c --maxiter $mi --seed $seed --num-threads 1 \
    --results-psms "$o/target.tsv" --decoy-results-psms "$o/decoy.tsv" \
    "$d/comet.pin" > "$o/stdout.log" 2> "$o/stderr.log"
  echo "mi=$mi seed=$seed $n rc=$? $(date -Is)" >> "$SP/maxiter/progress.log"
done; done; done
echo "MAXITER ALL DONE $(date -Is)" >> "$SP/maxiter/progress.log"
