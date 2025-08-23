#!/usr/bin/env bash
set -euo pipefail

# Usage: scripts/make_plot_png.sh [prefix]
# Default prefix is pa_comparison

PREFIX=${1:-pa_comparison}
WD=$(pwd)

# Ensure required commands exist
command -v mummerplot >/dev/null 2>&1 || { echo "mummerplot not found in PATH" >&2; exit 2; }
command -v gnuplot >/dev/null 2>&1 || { echo "gnuplot not found in PATH" >&2; exit 2; }

FILTERED="${PREFIX}.filtered.delta"
DELTA="${PREFIX}.delta"
GPFILE="${PREFIX}.gp"
OUTPNG="${PREFIX}.png"

# Prefer filtered delta
if [ -s "$FILTERED" ]; then
  echo "Using filtered delta: $FILTERED"
  mummerplot --fat --layout --filter -p "$PREFIX" "$FILTERED"
elif [ -s "$DELTA" ]; then
  echo "Using delta: $DELTA"
  mummerplot --fat --layout -p "$PREFIX" "$DELTA"
else
  echo "No delta file found: $FILTERED or $DELTA" >&2
  exit 3
fi

if [ ! -f "$GPFILE" ]; then
  echo "$GPFILE not created by mummerplot" >&2
  exit 4
fi

# Create a non-interactive version of the .gp file suitable for PNG output
TMPGP=$(mktemp --suffix=.gp)
# Replace the first 'set terminal' with a pngcairo terminal; remove interactive lines and unsupported mouse options
awk '
  BEGIN{repl=0}
  /^set terminal/ && repl==0 { print "set terminal pngcairo size 1600,1200 enhanced font \"Arial,10\""; print "set output \"'"$OUTPNG"'\""; repl=1; next }
  /^pause -1/ { next }
  /^print / { next }
  /^set mouse/ { next }
  /^set mouse / { next }
  /^set mouse mouseformat/ { next }
  /^set mouse clipboardformat/ { next }
  { print }
' "$GPFILE" > "$TMPGP"

# Run gnuplot non-interactively to write PNG
echo "Running gnuplot to generate $OUTPNG (this may overwrite existing file)"
gnuplot "$TMPGP"

if [ -f "$OUTPNG" ]; then
  echo "Created $OUTPNG"
  ls -lh "$OUTPNG"
else
  echo "Failed to produce $OUTPNG" >&2
  rm -f "$TMPGP"
  exit 5
fi

# Cleanup
rm -f "$TMPGP"

echo "Done"
