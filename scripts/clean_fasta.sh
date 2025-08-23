#!/usr/bin/env bash
set -euo pipefail

if [ "$#" -lt 1 ]; then
  echo "Usage: $0 <fasta1> [fasta2 ...]" >&2
  exit 2
fi

for f in "$@"; do
  if [ ! -f "$f" ]; then
    echo "File not found: $f" >&2
    continue
  fi
  out="${f%.fasta}.fixed.fasta"
  # Remove everything after the first space in header lines, and strip all whitespace from sequence lines
  awk '/^>/{sub(/ .*/,"",
$0)} {gsub(/[ \t\r]/,""); print}' "$f" > "$out"
  echo "Wrote $out"
done
