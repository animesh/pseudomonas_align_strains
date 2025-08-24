#!/usr/bin/env python3
"""Extract unique (non-aligned) segments from two FASTA files using a show-coords file.

Usage:
  python3 scripts/extract_unique_segments.py --ref ATCC.fasta --qry ST.fasta \ 
        --coords pa_comparison.coords --min-len 200

The script parses the coords file produced by `show-coords -rcl` and removes aligned
intervals from each genome. Remaining segments (>= min length) are written to
`<prefix>.unique.fasta` (defaults to ATCC.unique.fasta and ST.unique.fasta).

If `--coords` is not present or cannot be parsed, the script will exit with an error.
"""
import argparse
from pathlib import Path
from collections import defaultdict
import sys


def parse_coords(coords_path: Path):
    """Parse a show-coords file and return aligned intervals per file.

    Returns: dict with keys 'ref' and 'qry', each a dict mapping sequence name -> list of (start0, end0)
    """
    ref_intervals = defaultdict(list)
    qry_intervals = defaultdict(list)

    with coords_path.open('r') as fh:
        for line in fh:
            line = line.rstrip('\n')
            if not line.strip():
                continue
            # Skip header lines that don't start with a digit or '>'
            parts = line.split()
            # find first four integer columns in the line
            ints = []
            for p in parts[:12]:
                try:
                    ints.append(int(p))
                except ValueError:
                    continue
            if len(ints) < 4:
                continue
            # show-coords arrangement: ref_start ref_end qry_start qry_end ... then ref and qry names at end
            ref_s, ref_e, qry_s, qry_e = ints[0:4]
            # names are usually the last two columns
            if len(parts) >= 2:
                ref_name = parts[-1]
                qry_name = parts[-2]
            else:
                # fallback generic names
                ref_name = 'ref'
                qry_name = 'qry'

            # convert to 0-based half-open
            ref_intervals[ref_name].append((ref_s - 1, ref_e))
            qry_intervals[qry_name].append((qry_s - 1, qry_e))

    return {'ref': ref_intervals, 'qry': qry_intervals}


def merge_intervals(intervals):
    if not intervals:
        return []
    intervals = sorted(intervals, key=lambda x: x[0])
    merged = [intervals[0]]
    for s, e in intervals[1:]:
        last_s, last_e = merged[-1]
        if s <= last_e:
            merged[-1] = (last_s, max(last_e, e))
        else:
            merged.append((s, e))
    return merged


def fasta_iter(path: Path):
    """Yield (header, sequence) from a FASTA file (simple, memory-friendly)."""
    header = None
    seq_lines = []
    with path.open('r') as fh:
        for line in fh:
            if line.startswith('>'):
                if header is not None:
                    yield header, ''.join(seq_lines)
                header = line[1:].strip().split()[0]
                seq_lines = []
            else:
                seq_lines.append(line.strip())
        if header is not None:
            yield header, ''.join(seq_lines)


def subtract_intervals(seq_len, intervals):
    """Given seq length and a list of merged intervals (0-based half-open),
    return list of (start, end) intervals that are NOT covered by the intervals.
    """
    if not intervals:
        return [(0, seq_len)]
    res = []
    cur = 0
    for s, e in intervals:
        if cur < s:
            res.append((cur, s))
        cur = max(cur, e)
    if cur < seq_len:
        res.append((cur, seq_len))
    return res


def write_fasta_segments(out_path: Path, seq_id, seq, segments):
    with out_path.open('a') as out:
        for i, (s, e) in enumerate(segments, start=1):
            seg = seq[s:e]
            hdr = f">{seq_id}_seg{i} {s+1}-{e}\n"
            out.write(hdr)
            # wrap at 80
            for i in range(0, len(seg), 80):
                out.write(seg[i:i+80] + '\n')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--ref', required=True, help='Reference FASTA (e.g., ATCC.fasta)')
    p.add_argument('--qry', required=True, help='Query FASTA (e.g., ST.fasta)')
    p.add_argument('--coords', required=True, help='show-coords file (output of show-coords -rcl)')
    p.add_argument('--min-len', type=int, default=200, help='Minimum unique segment length to keep')
    p.add_argument('--out-prefix', default=None, help='Output prefix (defaults to REF/ QRY base names)')
    args = p.parse_args()

    ref_path = Path(args.ref)
    qry_path = Path(args.qry)
    coords_path = Path(args.coords)

    if not coords_path.exists():
        print(f"Coords file not found: {coords_path}", file=sys.stderr)
        sys.exit(2)

    parsed = parse_coords(coords_path)
    ref_intervals = parsed['ref']
    qry_intervals = parsed['qry']

    # prepare output files
    out_ref = Path(args.out_prefix + '.ATCC.unique.fasta') if args.out_prefix else ref_path.with_suffix('.unique.fasta')
    out_qry = Path(args.out_prefix + '.ST.unique.fasta') if args.out_prefix else qry_path.with_suffix('.unique.fasta')

    # clear outputs
    for pth in (out_ref, out_qry):
        if pth.exists():
            pth.unlink()

    # process reference
    ref_written = 0
    for seq_id, seq in fasta_iter(ref_path):
        name = seq_id
        intervals = ref_intervals.get(name, [])
        merged = merge_intervals(intervals)
        unique = subtract_intervals(len(seq), merged)
        # filter by min length
        unique = [(s, e) for s, e in unique if (e - s) >= args.min_len]
        if unique:
            write_fasta_segments(out_ref, seq_id, seq, unique)
            ref_written += sum((e - s) for s, e in unique)

    # process query
    qry_written = 0
    for seq_id, seq in fasta_iter(qry_path):
        name = seq_id
        intervals = qry_intervals.get(name, [])
        merged = merge_intervals(intervals)
        unique = subtract_intervals(len(seq), merged)
        unique = [(s, e) for s, e in unique if (e - s) >= args.min_len]
        if unique:
            write_fasta_segments(out_qry, seq_id, seq, unique)
            qry_written += sum((e - s) for s, e in unique)

    print(f"Wrote unique segments: {out_ref} ({ref_written} bp total) and {out_qry} ({qry_written} bp total)")


if __name__ == '__main__':
    main()
