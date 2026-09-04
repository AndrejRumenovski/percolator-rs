#!/usr/bin/env bash
# Exact, resumable execution of the preregistered search and rescoring pipeline.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
EXP="$ROOT/validation/homology_depleted_entrapment"
SOURCE="/home/andrej-rumenovski/percolator_rs_out/entrapment"
CRUX="$SOURCE/crux-4.0.Linux.x86_64/bin/crux"
RS="$ROOT/target/release/percolator-rs"

conditions=(original homology_depleted size_control_130363 size_control_155921 size_control_196613)
stems=(
  28May2015-QE-HF-Anopheles-22-atrium-S-24H-2nd-01
  28May2015-QE-HF-Anopheles-23-atrium-P-24H-2nd-01
  28May2015-QE-HF-Anopheles-38-MAGs-P-3rd-02
  09Dec2015-QEHF1-Anopheles-5-atrium-P-12hpm-3rd-01
  22Oct2014-Anopheles-8-MAGs-S-01
  9March2015-29-MAGs-pellet-2ndRep-14N-male-02
)

search_all() {
  mkdir -p "$EXP/searches"
  for condition in "${conditions[@]}"; do
    fasta="$EXP/databases/$condition.fasta"
    for stem in "${stems[@]}"; do
      out="$EXP/searches/$condition/comet-$stem"
      mkdir -p "$out"
      if [ -s "$out/comet.pin" ] && [ -s "$out/comet.log.txt" ] && \
         grep -q 'Return Code:0' "$out/comet.log.txt"; then
        continue
      fi
      /usr/bin/time -v -o "$out/time.txt" \
        "$CRUX" comet --parameter-file "$EXP/baseline_reproduction/search_parameters/$stem-comet.params.txt" \
        "$SOURCE/$stem.mzML" "$fasta" --output-dir "$out" \
        >"$out/command.stdout.log" 2>"$out/command.stderr.log"
      test -s "$out/comet.pin"
    done
  done
}

run_rs() {
  local pin=$1 destination=$2 seed=$3 maxiter=$4
  mkdir -p "$destination"
  if [ -s "$destination/target.tsv" ] && [ -s "$destination/decoy.tsv" ]; then
    return
  fi
  args=(--canonical --no-select-c --seed "$seed" --num-threads 1)
  if [ "$maxiter" != canonical ]; then
    args+=(--maxiter "$maxiter")
  fi
  "$RS" "${args[@]}" --results-psms "$destination/target.tsv" \
    --decoy-results-psms "$destination/decoy.tsv" "$pin" \
    >"$destination/stdout.log" 2>"$destination/stderr.log"
}

rescore_primary() {
  for condition in "${conditions[@]}"; do
    for stem in "${stems[@]}"; do
      pin="$EXP/searches/$condition/comet-$stem/comet.pin"
      test -s "$pin"
      run_rs "$pin" "$EXP/percolator/primary/$condition/seed-1/comet-$stem" 1 canonical
    done
  done
}

rescore_dose() {
  for maxiter in 0 1 2 3; do
    for condition in "${conditions[@]}"; do
      for stem in "${stems[@]}"; do
        pin="$EXP/searches/$condition/comet-$stem/comet.pin"
        run_rs "$pin" "$EXP/percolator/dose/mi-$maxiter/$condition/seed-1/comet-$stem" 1 "$maxiter"
      done
    done
  done
}

rescore_seeds() {
  for condition in original homology_depleted; do
    for seed in 2 3; do
      for stem in "${stems[@]}"; do
        pin="$EXP/searches/$condition/comet-$stem/comet.pin"
        run_rs "$pin" "$EXP/percolator/seeds/$condition/seed-$seed/comet-$stem" "$seed" canonical
      done
    done
  done
}

prepare_enz_ablation() {
  python3 "$ROOT/validation/pep_rootcause_experiments.py" prepare-ablation \
    --input-root "$EXP/searches/original" --output-root "$EXP/ablation_pins/original" \
    --drop enzN enzC --manifest "$EXP/ablation_pins/original.manifest.json"
  python3 "$ROOT/validation/pep_rootcause_experiments.py" prepare-ablation \
    --input-root "$EXP/searches/homology_depleted" --output-root "$EXP/ablation_pins/homology_depleted" \
    --drop enzN enzC --manifest "$EXP/ablation_pins/homology_depleted.manifest.json"
}

rescore_enz_ablation() {
  prepare_enz_ablation
  for condition in original homology_depleted; do
    for stem in "${stems[@]}"; do
      pin="$EXP/ablation_pins/$condition/comet-$stem/comet.pin"
      run_rs "$pin" "$EXP/percolator/enz_ablation/$condition/seed-1/comet-$stem" 1 canonical
    done
  done
}

case "${1:-}" in
  search) search_all ;;
  primary) rescore_primary ;;
  dose) rescore_dose ;;
  seeds) rescore_seeds ;;
  enz-ablation) rescore_enz_ablation ;;
  *) echo "usage: $0 search|primary|dose|seeds|enz-ablation" >&2; exit 2 ;;
esac
