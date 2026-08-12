#!/usr/bin/env bash
# CI regression gate (no external data needed): runs percolator-rs on the committed
# fixture and asserts (a) q<0.01 yield within +/-TOL% of recorded reference, and
# (b) wall time + peak RSS within budget. Exits non-zero on any failure.
set -u
cd "$(dirname "$0")/.." || exit 2
source tests/expected.env

BIN=target/release/percolator-rs
if [ ! -x "$BIN" ]; then
  echo "[build] $BIN missing, building release..."
  cargo build --release || { echo "FAIL: build"; exit 2; }
fi

FIX=tests/fixtures/sample.pin
[ -f "$FIX" ] || { echo "FAIL: missing fixture $FIX"; exit 2; }

# time + RSS via GNU time if available, else fall back to the tool's own timer.
TIMEV=""
if command -v /usr/bin/time >/dev/null 2>&1; then TIMEV="/usr/bin/time -v -o /tmp/reg_time.txt"; fi

err=/tmp/reg_run.err
$TIMEV "$BIN" --canonical --seed 1 "$FIX" >/dev/null 2>"$err" || { echo "FAIL: run errored"; cat "$err"; exit 1; }

psm=$(grep -oP 'target PSMs q<0.01: \K[0-9]+' "$err")
pep=$(grep -oP 'target peptides q<0.01: \K[0-9]+' "$err")
wall=$(grep -oP 'or \K[0-9.]+(?= seconds)' "$err" | tail -1)   # from tool's own line
rss=""
if [ -f /tmp/reg_time.txt ]; then
  rss=$(awk -F': ' '/Maximum resident/{print $2}' /tmp/reg_time.txt)
  gnu_wall=$(awk -F': ' '/Elapsed \(wall/{split($2,a,":"); print (length(a)==3? a[1]*3600+a[2]*60+a[3] : a[1]*60+a[2])}' /tmp/reg_time.txt)
  [ -n "$gnu_wall" ] && wall="$gnu_wall"
fi

# ---- assertions ----
fail=0
assert_within() { # name value expected tol_pct
  awk -v v="$2" -v e="$3" -v t="$4" -v n="$1" 'BEGIN{
    lo=e*(1-t/100); hi=e*(1+t/100);
    if (v+0>=lo && v+0<=hi) { printf "  PASS  %-22s %s (ref %s, +/-%s%%)\n", n, v, e, t }
    else { printf "  FAIL  %-22s %s (ref %s, allowed %.2f..%.2f)\n", n, v, e, lo, hi; exit 1 }
  }' || return 1
}
assert_leq() { # name value budget unit
  awk -v v="$2" -v b="$3" -v n="$1" -v u="$4" 'BEGIN{
    if (v+0<=b+0) { printf "  PASS  %-22s %s%s (budget %s%s)\n", n, v, u, b, u }
    else { printf "  FAIL  %-22s %s%s (budget %s%s)\n", n, v, u, b, u; exit 1 }
  }' || return 1
}

echo "== percolator-rs regression gate (fixture: $(basename "$FIX")) =="
assert_within "PSM q<0.01"     "${psm:-0}" "$FIXTURE_PSM_Q01" "$YIELD_TOLERANCE_PCT" || fail=1
assert_within "peptide q<0.01" "${pep:-0}" "$FIXTURE_PEP_Q01" "$YIELD_TOLERANCE_PCT" || fail=1
assert_leq    "wall time"      "${wall:-999}" "$FIXTURE_TIME_BUDGET_S" "s" || fail=1
[ -n "$rss" ] && { assert_leq "peak RSS"  "$rss" "$FIXTURE_MEM_BUDGET_KB" "kB" || fail=1; }

if [ "$fail" -ne 0 ]; then
  echo "REGRESSION FAILED"
  exit 1
fi
echo "ALL CHECKS PASSED"
