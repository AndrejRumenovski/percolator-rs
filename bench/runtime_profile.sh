#!/usr/bin/env bash
# Reproducible, measurement-only runtime profile for the PXD032157 workload.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
INPUT="${RUNTIME_PROFILE_INPUT:-$ROOT/data/PXD032157}"
STAMP="$(date +%Y%m%d-%H%M%S)"
ARTIFACTS="${RUNTIME_PROFILE_OUT:-$HOME/percolator_rs_out/runtime-profile-$STAMP}"
SINGLE_REPEATS="${RUNTIME_PROFILE_SINGLE_REPEATS:-3}"
FULL_REPEATS="${RUNTIME_PROFILE_FULL_REPEATS:-1}"
EXPECTED_FILES="${RUNTIME_PROFILE_EXPECTED_FILES:-65}"
export LC_ALL=C

[[ "$SINGLE_REPEATS" =~ ^[1-9][0-9]*$ ]] || { echo "invalid single repeat count" >&2; exit 2; }
[[ "$FULL_REPEATS" =~ ^[1-9][0-9]*$ ]] || { echo "invalid full repeat count" >&2; exit 2; }
[ -d "$INPUT" ] || { echo "missing input directory: $INPUT" >&2; exit 2; }
if [ -e "$ARTIFACTS" ]; then
  echo "refusing to overwrite existing profiling directory: $ARTIFACTS" >&2
  exit 2
fi
mkdir -p "$ARTIFACTS/bin" "$ARTIFACTS/build" "$ARTIFACTS/profiles" "$ARTIFACTS/outputs"

mapfile -t SIZE_ORDER < <(
  find "$INPUT" -maxdepth 1 -type f -name '*.pin' -printf '%s\t%p\n' |
    sort -rn | cut -f2-
)
[ "${#SIZE_ORDER[@]}" -eq "$EXPECTED_FILES" ] || {
  echo "expected $EXPECTED_FILES PIN files, found ${#SIZE_ORDER[@]}" >&2
  exit 2
}
LARGE_PIN="${SIZE_ORDER[0]}"

# Alternate large/small files to balance batches in the same manner as bench/run_rs.sh.
ORDER=()
left=0
right=$((${#SIZE_ORDER[@]} - 1))
while [ "$left" -le "$right" ]; do
  ORDER+=("${SIZE_ORDER[$left]}")
  [ "$left" -eq "$right" ] || ORDER+=("${SIZE_ORDER[$right]}")
  left=$((left + 1))
  right=$((right - 1))
done
printf '%s\n' "${ORDER[@]}" >"$ARTIFACTS/input-order.txt"

rustc -Vv >"$ARTIFACTS/build/rustc.txt"
cargo -V >"$ARTIFACTS/build/cargo.txt"
uname -a >"$ARTIFACTS/build/uname.txt"
lscpu >"$ARTIFACTS/build/lscpu.txt"
git -C "$ROOT" rev-parse HEAD >"$ARTIFACTS/build/git-head.txt"
git -C "$ROOT" diff >"$ARTIFACTS/build/instrumentation.patch"

# Record why perf was not selected if this host disallows performance counters.
set +e
perf stat -e task-clock -- true >"$ARTIFACTS/build/perf-probe.txt" 2>&1
perf_status=$?
set -e
printf '%s\n' "$perf_status" >"$ARTIFACTS/build/perf-probe-exit-code.txt"

NORMAL_TARGET="$ARTIFACTS/build/target-normal"
PROFILE_TARGET="$ARTIFACTS/build/target-profile"
CPU_TARGET="$ARTIFACTS/build/target-cpu"
cargo build --release --manifest-path "$ROOT/Cargo.toml" --target-dir "$NORMAL_TARGET"
cp "$NORMAL_TARGET/release/percolator-rs" "$ARTIFACTS/bin/percolator-rs-normal"
cargo build --release --features profiling --manifest-path "$ROOT/Cargo.toml" \
  --target-dir "$PROFILE_TARGET"
cp "$PROFILE_TARGET/release/percolator-rs" "$ARTIFACTS/bin/percolator-rs-profile"
RUSTFLAGS='-C force-frame-pointers=yes -C debuginfo=1' \
  cargo build --release --features profiling --manifest-path "$ROOT/Cargo.toml" \
  --target-dir "$CPU_TARGET"
cp "$CPU_TARGET/release/percolator-rs" "$ARTIFACTS/bin/percolator-rs-cpu"

NORMAL_BIN="$ARTIFACTS/bin/percolator-rs-normal"
PROFILE_BIN="$ARTIFACTS/bin/percolator-rs-profile"
CPU_BIN="$ARTIFACTS/bin/percolator-rs-cpu"
TIMINGS="$ARTIFACTS/timings.tsv"
printf 'configuration\tbuild\trepetition\twall_ns\tprocesses\tintra_file_threads\tcpu_sampling\n' >"$TIMINGS"

now_ns() { date +%s%N; }

run_file() {
  local bin=$1 pin=$2 threads=$3 destination=$4 json=${5:-} cpu=${6:-} proteins=${7:-0} allocations=${8:-0}
  mkdir -p "$destination"
  local command=(
    "$bin" --canonical --seed 1 --num-threads "$threads"
    --results-psms "$destination/target.psms.tsv"
    --decoy-results-psms "$destination/decoy.psms.tsv"
    --results-peptides "$destination/target.peptides.tsv"
    --decoy-results-peptides "$destination/decoy.peptides.tsv"
  )
  if [ "$proteins" -eq 1 ]; then
    command+=(
      --results-proteins "$destination/target.proteins.tsv"
      --decoy-results-proteins "$destination/decoy.proteins.tsv"
    )
  fi
  [ -z "$json" ] || command+=(--profile-json "$json")
  [ -z "$cpu" ] || command+=(--profile-cpu "$cpu")
  [ "$allocations" -eq 0 ] || command+=(--profile-allocations)
  command+=("$pin")
  "${command[@]}" >"$destination/stdout.log" 2>"$destination/stderr.log"
}

timed_single() {
  local config=$1 build=$2 repetition=$3 bin=$4 threads=$5 profile=$6 cpu=$7 proteins=${8:-0} allocations=${9:-0}
  local tag="${config}_${build}_r${repetition}"
  local destination="$ARTIFACTS/outputs/$tag"
  local profile_dir="$ARTIFACTS/profiles/$tag"
  local json='' cpu_prefix=''
  mkdir -p "$profile_dir"
  if [ "$profile" -eq 1 ]; then json="$profile_dir/profile.json"; fi
  if [ "$cpu" -eq 1 ]; then cpu_prefix="$profile_dir/cpu"; fi
  local start end
  start=$(now_ns)
  run_file "$bin" "$LARGE_PIN" "$threads" "$destination" "$json" "$cpu_prefix" "$proteins" "$allocations"
  end=$(now_ns)
  printf '%s\t%s\t%s\t%s\t1\t%s\t%s\n' \
    "$config" "$build" "$repetition" "$((end - start))" "$threads" "$cpu" >>"$TIMINGS"
}

run_batch_file() {
  local pin=$1 stem destination json='' cpu_prefix=''
  stem=$(basename "$pin" .pin)
  destination="$BATCH_OUTPUT/$stem"
  if [ "$BATCH_PROFILE" -eq 1 ]; then json="$BATCH_PROFILES/$stem.json"; fi
  if [ "$BATCH_CPU" -eq 1 ]; then cpu_prefix="$BATCH_PROFILES/$stem.cpu"; fi
  run_file "$BATCH_BIN" "$pin" 1 "$destination" "$json" "$cpu_prefix" 0
}
export -f run_file run_batch_file

timed_batch() {
  local config=$1 build=$2 repetition=$3 bin=$4 concurrency=$5 profile=$6 cpu=$7
  local tag="${config}_${build}_r${repetition}"
  export BATCH_BIN="$bin"
  export BATCH_OUTPUT="$ARTIFACTS/outputs/$tag"
  export BATCH_PROFILES="$ARTIFACTS/profiles/$tag"
  export BATCH_PROFILE="$profile"
  export BATCH_CPU="$cpu"
  mkdir -p "$BATCH_OUTPUT" "$BATCH_PROFILES"
  local start end
  start=$(now_ns)
  if [ "$concurrency" -eq 1 ]; then
    local pin
    for pin in "${ORDER[@]}"; do run_batch_file "$pin"; done
  else
    printf '%s\0' "${ORDER[@]}" |
      xargs -0 -P "$concurrency" -I{} bash -c 'run_batch_file "$1"' _ '{}'
  fi
  end=$(now_ns)
  printf '%s\t%s\t%s\t%s\t%s\t1\t%s\n' \
    "$config" "$build" "$repetition" "$((end - start))" "$EXPECTED_FILES" "$cpu" >>"$TIMINGS"
}

for repetition in $(seq 1 "$SINGLE_REPEATS"); do
  if [ $((repetition % 2)) -eq 1 ]; then
    timed_single single_file_t1 normal "$repetition" "$NORMAL_BIN" 1 0 0
    timed_single single_file_t1 instrumented "$repetition" "$PROFILE_BIN" 1 1 0
    timed_single single_file_t3 normal "$repetition" "$NORMAL_BIN" 3 0 0
    timed_single single_file_t3 instrumented "$repetition" "$PROFILE_BIN" 3 1 0
  else
    timed_single single_file_t1 instrumented "$repetition" "$PROFILE_BIN" 1 1 0
    timed_single single_file_t1 normal "$repetition" "$NORMAL_BIN" 1 0 0
    timed_single single_file_t3 instrumented "$repetition" "$PROFILE_BIN" 3 1 0
    timed_single single_file_t3 normal "$repetition" "$NORMAL_BIN" 3 0 0
  fi
done

# The extra single-file run covers the conditional protein-inference stage.
timed_single single_file_t1_protein instrumented 1 "$PROFILE_BIN" 1 1 0 1
timed_single single_file_t1_allocations instrumented 1 "$PROFILE_BIN" 1 1 0 0 1
timed_single single_file_t1_cpu cpu 1 "$CPU_BIN" 1 1 1

for repetition in $(seq 1 "$FULL_REPEATS"); do
  timed_batch full_sequential normal "$repetition" "$NORMAL_BIN" 1 0 0
  timed_batch full_sequential instrumented "$repetition" "$PROFILE_BIN" 1 1 0
  timed_batch full_n4 normal "$repetition" "$NORMAL_BIN" 4 0 0
  timed_batch full_n4 instrumented "$repetition" "$PROFILE_BIN" 4 1 0
done
timed_batch full_n4_cpu cpu 1 "$CPU_BIN" 4 1 1

python3 "$ROOT/bench/runtime_profile_report.py" \
  --artifacts "$ARTIFACTS" \
  --json "$ARTIFACTS/runtime-profile.json" \
  --markdown "$ARTIFACTS/RUNTIME_PROFILE.md"

echo "runtime profile complete: $ARTIFACTS"
