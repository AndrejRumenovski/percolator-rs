#!/usr/bin/env bash
set -uo pipefail
SP="/tmp/claude-1000/-run-media-andrej-rumenovski-New-Volume-Code-percolator-rs/9d5d4d63-9a67-4e88-ab50-b68d02980cac/scratchpad"
export LD_LIBRARY_PATH="$HOME/opt/perc-libs:${LD_LIBRARY_PATH:-}"
PERC="$HOME/opt/percolator-root/usr/bin/percolator"
ROOT="$HOME/percolator_rs_out/entrapment"
for d in "$ROOT"/comet-*/; do
  n=$(basename "$d")
  [ -s "$SP/cpp_tdc/$n.target.tsv" ] && continue
  /usr/bin/time -v "$PERC" --post-processing-tdc --seed 1 --num-threads 1 \
    --results-psms "$SP/cpp_tdc/$n.target.tsv" \
    --decoy-results-psms "$SP/cpp_tdc/$n.decoy.tsv" \
    "$d/comet.pin" > "$SP/cpp_tdc/$n.stdout" 2> "$SP/cpp_tdc/$n.stderr"
  echo "done $n rc=$? $(date -Is)" >> "$SP/cpp_tdc/progress.log"
done
echo "ALL DONE $(date -Is)" >> "$SP/cpp_tdc/progress.log"
