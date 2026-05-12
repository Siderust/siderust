#!/usr/bin/env bash
# check_dynamics_coverage.sh — CI gate for src/astro/dynamics/ line coverage.
#
# Runs cargo llvm-cov for the dynamics unit-test suite, then fails if any
# .rs file under src/astro/dynamics/ has line coverage below the threshold.
#
# Usage (from the siderust crate directory):
#   bash scripts/check_dynamics_coverage.sh
#
# Exit codes:
#   0  all dynamics files meet the coverage threshold
#   1  one or more files are below the threshold (table printed to stderr)
#   2  prerequisite missing (cargo-llvm-cov not installed)

set -euo pipefail

THRESHOLD=90.0
COV_JSON="target/dyn-cov.json"
FILTER="src/astro/dynamics/"

# ---------------------------------------------------------------------------
# Prerequisites
# ---------------------------------------------------------------------------
if ! cargo llvm-cov --version >/dev/null 2>&1; then
    echo "error: cargo-llvm-cov is not installed." >&2
    echo "       Install with: cargo install cargo-llvm-cov" >&2
    exit 2
fi

# ---------------------------------------------------------------------------
# Generate coverage report
# ---------------------------------------------------------------------------
echo "Running cargo llvm-cov (dynamics unit tests)…"
cargo llvm-cov --lib --json --output-path "${COV_JSON}" -- astro::dynamics

# ---------------------------------------------------------------------------
# Parse and gate
# ---------------------------------------------------------------------------
python3 - "${COV_JSON}" "${FILTER}" "${THRESHOLD}" <<'PYEOF'
import json
import sys

cov_json   = sys.argv[1]
path_token = sys.argv[2]
threshold  = float(sys.argv[3])

with open(cov_json) as fh:
    data = json.load(fh)

all_files = data.get("data", [{}])[0].get("files", [])
dyn_files = [f for f in all_files if path_token in f.get("filename", "")]

if not dyn_files:
    print(f"warning: no files matched '{path_token}' in {cov_json}", file=sys.stderr)
    sys.exit(0)

failing = []
for f in sorted(dyn_files, key=lambda x: x["filename"]):
    short = f["filename"].split(path_token)[-1]
    lm    = f.get("summary", {}).get("lines", {})
    pct   = lm.get("percent", 0.0)
    count = lm.get("count", 0)
    cov   = lm.get("covered", 0)
    if pct < threshold:
        failing.append((short, cov, count, pct))

if failing:
    print(f"\nFAIL: {len(failing)} file(s) below {threshold}% line coverage:\n",
          file=sys.stderr)
    print(f"  {'File':<60} {'Covered':>9}  {'%':>7}", file=sys.stderr)
    print(f"  {'-'*60} {'-'*9}  {'-'*7}", file=sys.stderr)
    for short, cov, cnt, pct in failing:
        print(f"  {short:<60} {cov:>4}/{cnt:<4}   {pct:>6.1f}%", file=sys.stderr)
    print(file=sys.stderr)
    sys.exit(1)

print(f"OK: all {len(dyn_files)} dynamics files >= {threshold}% line coverage.")
PYEOF
