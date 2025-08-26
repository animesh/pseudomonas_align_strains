#!/usr/bin/env bash
# Wrapper to run scripts/make_rects.py safely.
# By default it writes to pa_comparison.rects.auto.* to avoid overwriting a hand-edited pa_comparison.rects.gp

set -euo pipefail
SCRIPT_DIR=$(dirname "$(readlink -f "$0")")
REPO_ROOT=$(dirname "$SCRIPT_DIR")

python3 "$SCRIPT_DIR/make_rects.py" --out-gp "$REPO_ROOT/pa_comparison.rects.auto.gp" \
  --out-png "$REPO_ROOT/pa_comparison.rects.auto.png" --export-svg --export-pdf

echo "Wrote pa_comparison.rects.auto.* (gp/png/svg/pdf). To overwrite canonical files run: make rects-force or run with --out-gp pa_comparison.rects.gp" 
