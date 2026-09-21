#!/usr/bin/env bash
# Repeatable local/CI checks. Public-data performance runs are opt-in.
set -euo pipefail
cd "$(dirname "$0")/.."

case "${1:-}" in
  ""|--full-benchmark) ;;
  *) echo "usage: bash scripts/check.sh [--full-benchmark]" >&2; exit 2 ;;
esac
if [ "$#" -gt 1 ]; then
  echo "usage: bash scripts/check.sh [--full-benchmark]" >&2
  exit 2
fi

cargo fmt --all -- --check
cargo clippy --all-targets --all-features --locked -- -D warnings
cargo test --release --all-targets --locked
cargo test --release --all-targets --all-features --locked
RUSTDOCFLAGS="-D warnings" cargo doc --no-deps --all-features --locked

# Restore the ordinary executable after the profiling-feature test build.
cargo build --release --locked
for gate in regression selection_regression feature_report ensemble_regression protein_regression; do
  bash "tests/${gate}.sh"
done
python3 -m unittest -v bench/multidataset/test_normalize_sage_pin.py
(cd bench/protein_calibration && python3 -m unittest -v test_report.py)
python3 -m unittest -v validation.test_psm_agreement

if [ "${1:-}" = --full-benchmark ]; then
  if [ ! -d data/PXD032157 ]; then
    echo "Full benchmark requires local PIN files under data/PXD032157." >&2
    exit 2
  fi
  bash bench/regression.sh
fi
