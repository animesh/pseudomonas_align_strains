#!/usr/bin/env python3
"""Count nucleotides in FASTA files.

Usage:
  python3 scripts/count_nucleotides.py [file1.fasta file2.fasta ...]

If no filenames are provided the script will default to: ATCC.fasta ST.fasta
It prints per-file total bases and a breakdown for A C G T N and other characters.
"""
import sys
from collections import Counter
from pathlib import Path

def count_nucleotides(path: Path):
    counts = Counter()
    total = 0
    seq_count = 0
    try:
        with path.open('r') as fh:
            for line in fh:
                if line.startswith('>'):
                    seq_count += 1
                    continue
                seq = line.strip().upper()
                if not seq:
                    continue
                for ch in seq:
                    counts[ch] += 1
                    total += 1
    except FileNotFoundError:
        return None, None, None
    return total, counts, seq_count

def pretty_print(path: Path, total, counts: Counter, seq_count=None):
    print(f"File: {path}")
    if total is None:
        print("  NOT FOUND")
        return
    print(f"  Total bases: {total}")
    if seq_count is not None:
        print(f"  Sequences: {seq_count}")
    bases = ['A','C','G','T','N']
    for b in bases:
        print(f"    {b}: {counts.get(b,0)}")
    # show any other characters (e.g., ambiguous IUPAC codes) if present
    others = {k:v for k,v in counts.items() if k not in bases}
    if others:
        print("    Other characters:")
        for k,v in sorted(others.items()):
            print(f"      {k}: {v}")

def main(args):
    if len(args) >= 1:
        files = args
    else:
        files = ['ATCC.fasta', 'ST.fasta']

    grand_total = 0
    for fname in files:
        path = Path(fname)
        total, counts, seq_count = count_nucleotides(path)
        if total is None:
            pretty_print(path, None, None)
            continue
        pretty_print(path, total, counts, seq_count=seq_count)
        grand_total += total

    print('\nGrand total bases (all files processed):', grand_total)

if __name__ == '__main__':
    main(sys.argv[1:])
