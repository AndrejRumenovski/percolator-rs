#!/usr/bin/env bash
# Checks the linear-SVM explanation report is complete, deterministic, and
# rejects a nonlinear model whose feature weights would be misleading.
set -u
if [ ! -f Cargo.toml ] || [ ! -d tests ]; then
  cd "$(dirname "$0")/.." || exit 2
fi

BIN=target/release/percolator-rs
[ -x "$BIN" ] || { cargo build --release || exit 2; }
FIX=tests/fixtures/sample.pin
FIRST=/tmp/percolator_feature_report_first.tsv
SECOND=/tmp/percolator_feature_report_second.tsv
ERR=/tmp/percolator_feature_report_mlp.err

"$BIN" --canonical --seed 1 --feature-report "$FIRST" "$FIX" >/dev/null
"$BIN" --canonical --seed 1 --feature-report "$SECOND" "$FIX" >/dev/null

if ! cmp -s "$FIRST" "$SECOND"; then
  echo "FAIL: feature reports differ between identical seeded runs"
  exit 1
fi
count=$(awk 'BEGIN { count=0 } !/^#/ && $1 != "feature_index" { count++ } END { print count }' "$FIRST")
baseline=$(awk -F= '/^# baseline_target_psms_q<0.01=/{ print $2 }' "$FIRST")
if [ "$count" != 21 ] || [ "$baseline" != 132 ]; then
  echo "FAIL: expected 21 feature rows and q<0.01 baseline 132; got $count and ${baseline:-missing}"
  exit 1
fi
if "$BIN" --rescore-model mlp --feature-report "$FIRST" "$FIX" >/dev/null 2>"$ERR"; then
  echo "FAIL: nonlinear MLP feature report unexpectedly succeeded"
  exit 1
fi
if ! grep -q 'feature-report currently supports only' "$ERR"; then
  echo "FAIL: missing clear MLP feature-report error"
  exit 1
fi

echo "== percolator-rs feature-report regression gate =="
echo "  PASS  feature rows             $count"
echo "  PASS  baseline PSM q<0.01      $baseline"
echo "  PASS  deterministic report     byte-identical"
echo "  PASS  MLP rejection            clear unsupported-model error"
echo "ALL CHECKS PASSED"
