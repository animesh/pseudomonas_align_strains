#!/usr/bin/env python3
"""Produce matched and unmatched segments for reference and query FASTA files
based on a `show-coords -rcl` file.

Outputs (by default, using pa_comparison.coords):
  ATCC.matched.tsv, ATCC.matched.fasta
  ATCC.unmatched.tsv, ATCC.unmatched.fasta
  ST.matched.tsv, ST.matched.fasta
  ST.unmatched.tsv, ST.unmatched.fasta

Each TSV has columns: seq_id, start(1-based), end(inclusive), length
"""
import sys
from pathlib import Path
from collections import defaultdict
import argparse


def parse_coords(coords_path: Path):
    ref_intervals = defaultdict(list)
    qry_intervals = defaultdict(list)
    with coords_path.open() as fh:
        for line in fh:
            line = line.rstrip('\n')
            if not line.strip():
                continue
            parts = line.split()
            # try to find four integers near the start (S1 E1 S2 E2)
            ints = []
            for p in parts[:12]:
                try:
                    ints.append(int(p))
                except ValueError:
                    continue
            if len(ints) < 4:
                continue
            s1, e1, s2, e2 = ints[0:4]
            # assume last two tokens are ref_name and qry_name (show-coords typical)
            if len(parts) >= 2:
                ref_name = parts[-2]
                qry_name = parts[-1]
            else:
                continue
            # convert to 0-based half-open
            ref_intervals[ref_name].append((s1 - 1, e1))
            qry_intervals[qry_name].append((s2 - 1, e2))
    return ref_intervals, qry_intervals


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
    header = None
    seq_lines = []
    with path.open() as fh:
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


def write_tsv_and_fasta(out_prefix: Path, seq_id: str, seq: str, segments, mode):
    # mode: 'matched' or 'unmatched'
    tsv_path = out_prefix.with_suffix('.' + mode + '.tsv')
    fasta_path = out_prefix.with_suffix('.' + mode + '.fasta')
    with tsv_path.open('a') as tsv, fasta_path.open('a') as f:
        for i, (s, e) in enumerate(segments, start=1):
            length = e - s
            # positions 1-based inclusive
            start1 = s + 1
            end1 = e
            tsv.write(f"{seq_id}\t{start1}\t{end1}\t{length}\n")
            header = f">{seq_id}.{mode}.seg{i} {start1}-{end1}\n"
            f.write(header)
            # wrap sequence lines
            seg = seq[s:e]
            for j in range(0, len(seg), 80):
                f.write(seg[j:j+80] + '\n')


def process(fasta_path: Path, intervals_map, out_prefix: Path, min_len: int = 1, genome_label: str = 'REF'):
    # clear existing output files for this prefix
    for suffix in ('.matched.tsv', '.matched.fasta', '.unmatched.tsv', '.unmatched.fasta'):
        p = out_prefix.with_suffix(suffix)
        if p.exists():
            p.unlink()

    total_matched = 0
    total_unmatched = 0
    seq_count = 0
    for seq_id, seq in fasta_iter(fasta_path):
        seq_count += 1
        seq_len = len(seq)
        intervals = intervals_map.get(seq_id, [])
        merged = merge_intervals(intervals)
        matched = [(s, e) for s, e in merged if (e - s) >= min_len]
        unmatched_raw = subtract_intervals(seq_len, merged)
        unmatched = [(s, e) for s, e in unmatched_raw if (e - s) >= min_len]
        if matched:
            write_tsv_and_fasta(out_prefix, seq_id, seq, matched, 'matched')
            total_matched += sum(e - s for s, e in matched)
        if unmatched:
            write_tsv_and_fasta(out_prefix, seq_id, seq, unmatched, 'unmatched')
            total_unmatched += sum(e - s for s, e in unmatched)
    return seq_count, total_matched, total_unmatched


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--ref', default='ATCC.fasta', help='Reference FASTA')
    p.add_argument('--qry', default='ST.fasta', help='Query FASTA')
    p.add_argument('--coords', default='pa_comparison.coords', help='show-coords file')
    p.add_argument('--min-len', type=int, default=1, help='minimum segment length to output')
    args = p.parse_args()

    coords = Path(args.coords)
    if not coords.exists():
        print(f"Coords file not found: {coords}", file=sys.stderr)
        sys.exit(2)

    ref_intervals, qry_intervals = parse_coords(coords)

    ref_fa = Path(args.ref)
    qry_fa = Path(args.qry)

    ref_prefix = Path('ATCC')
    qry_prefix = Path('ST')

    r_seq_count, r_mat, r_unm = process(ref_fa, ref_intervals, ref_prefix, min_len=args.min_len, genome_label='REF')
    q_seq_count, q_mat, q_unm = process(qry_fa, qry_intervals, qry_prefix, min_len=args.min_len, genome_label='QRY')

    print(f"Reference: {r_seq_count} sequences; matched bp: {r_mat}; unmatched bp: {r_unm}")
    print(f"Query: {q_seq_count} sequences; matched bp: {q_mat}; unmatched bp: {q_unm}")


if __name__ == '__main__':
    main()
