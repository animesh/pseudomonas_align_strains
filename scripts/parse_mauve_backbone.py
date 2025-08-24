#!/usr/bin/env python3
"""Parse Mauve backbone and write per-genome matched/unmatched TSV + FASTA.

Usage:
  python3 scripts/parse_mauve_backbone.py --backbone alignment.mauve.backbone \
      --ref ATCC.fasta --qry ST.fasta --out-prefix mauve

Outputs:
  <prefix>.ATCC.matched.fasta / .tsv
  <prefix>.ATCC.unmatched.fasta / .tsv
  <prefix>.ST.matched.fasta / .tsv
  <prefix>.ST.unmatched.fasta / .tsv

This script assumes backbone coordinates are 1-based positions into a
concatenated representation of each FASTA file (the order of contigs is the
same as in the FASTA). Negative coordinates indicate reverse orientation; the
absolute value is used to map to positions.
"""
import argparse
from pathlib import Path
from collections import defaultdict
import sys


def read_backbone(path: Path):
    rows = []
    with path.open() as fh:
        header = fh.readline()
        for line in fh:
            if not line.strip():
                continue
            parts = line.split()
            if len(parts) < 4:
                continue
            a, b, c, d = parts[:4]
            try:
                a = int(a); b = int(b); c = int(c); d = int(d)
            except ValueError:
                continue
            rows.append((a, b, c, d))
    return rows


def build_contig_index(fasta: Path):
    contigs = []
    cur = 1
    header = None
    seq_len = 0
    with fasta.open() as fh:
        for line in fh:
            if line.startswith('>'):
                if header is not None:
                    contigs.append((header, seq_len, cur))
                    cur += seq_len
                header = line[1:].strip().split()[0]
                seq_len = 0
            else:
                seq_len += len(line.strip())
        if header is not None:
            contigs.append((header, seq_len, cur))
    # contigs: list of (name, length, cum_start) where cum_start is 1-based
    return contigs


def map_pos_to_contig(contigs, pos):
    # pos is 1-based absolute position (use abs for negatives)
    if pos <= 0:
        return None
    for name, length, start in contigs:
        end = start + length - 1
        if start <= pos <= end:
            local_start = pos - start + 1
            return name, local_start
    return None


def add_interval(intervals_map, contig_name, s_abs, e_abs, contigs):
    # s_abs and e_abs are absolute positions (1-based). Map to contig and add local interval
    if s_abs <= 0 or e_abs <= 0:
        return
    # ensure s <= e
    s, e = min(s_abs, e_abs), max(s_abs, e_abs)
    # find contig for s
    s_map = map_pos_to_contig(contigs, s)
    e_map = map_pos_to_contig(contigs, e)
    if s_map is None or e_map is None:
        return
    s_contig, s_loc = s_map
    e_contig, e_loc = e_map
    if s_contig != e_contig:
        # split across contigs: add to start contig until its end, then to intermediate, then to end contig
        # naive split: iterate contigs in order and assign overlapping portions
        started = False
        for name, length, start in contigs:
            end = start + length - 1
            if end < s:
                continue
            if start > e:
                break
            seg_s = max(s, start)
            seg_e = min(e, end)
            local_s = seg_s - start + 1
            local_e = seg_e - start + 1
            intervals_map[name].append((local_s, local_e))
    else:
        intervals_map[s_contig].append((s_loc, e_loc))


def merge_intervals(intervals):
    if not intervals:
        return []
    intervals = sorted(intervals, key=lambda x: x[0])
    merged = [intervals[0]]
    for s, e in intervals[1:]:
        last_s, last_e = merged[-1]
        if s <= last_e + 1:
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


def write_segments_for_genome(contigs, intervals_map, fasta_path: Path, out_prefix: Path, genome_label):
    # Build a map from contig -> sequence
    seq_map = {name: seq for name, seq in fasta_iter(fasta_path)}
    matched_tsv = out_prefix.with_suffix('.' + genome_label + '.mauve.matched.tsv')
    matched_fa = out_prefix.with_suffix('.' + genome_label + '.mauve.matched.fasta')
    unmatched_tsv = out_prefix.with_suffix('.' + genome_label + '.mauve.unmatched.tsv')
    unmatched_fa = out_prefix.with_suffix('.' + genome_label + '.mauve.unmatched.fasta')
    # clear files
    for p in (matched_tsv, matched_fa, unmatched_tsv, unmatched_fa):
        if p.exists():
            p.unlink()

    total_matched = 0
    total_unmatched = 0
    for name, length, start in contigs:
        seq = seq_map.get(name, '')
        ints = intervals_map.get(name, [])
        merged = merge_intervals(ints)
        # write matched
        if merged:
            with matched_tsv.open('a') as mt, matched_fa.open('a') as mf:
                for i, (s, e) in enumerate(merged, start=1):
                    mt.write(f"{name}\t{s}\t{e}\t{e-s+1}\n")
                    mf.write(f">{name}.mauve.match.{i} {s}-{e}\n")
                    mf.write(seq[s-1:e] + '\n')
                    total_matched += (e - s + 1)
        # compute unmatched by subtracting merged from full
        unmatched = []
        cur = 1
        for s, e in merged:
            if cur < s:
                unmatched.append((cur, s-1))
            cur = e + 1
        if cur <= length:
            unmatched.append((cur, length))
        if unmatched:
            with unmatched_tsv.open('a') as ut, unmatched_fa.open('a') as uf:
                for i, (s, e) in enumerate(unmatched, start=1):
                    ut.write(f"{name}\t{s}\t{e}\t{e-s+1}\n")
                    uf.write(f">{name}.mauve.unmatch.{i} {s}-{e}\n")
                    uf.write(seq[s-1:e] + '\n')
                    total_unmatched += (e - s + 1)

    return total_matched, total_unmatched


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--backbone', default='alignment.mauve.backbone')
    p.add_argument('--ref', default='ATCC.fasta')
    p.add_argument('--qry', default='ST.fasta')
    p.add_argument('--out-prefix', default='mauve')
    args = p.parse_args()

    backbone = Path(args.backbone)
    if not backbone.exists():
        print('Backbone file not found', file=sys.stderr)
        sys.exit(2)

    rows = read_backbone(backbone)
    ref_contigs = build_contig_index(Path(args.ref))
    qry_contigs = build_contig_index(Path(args.qry))

    # intervals per contig
    ref_intervals = defaultdict(list)
    qry_intervals = defaultdict(list)

    for a, b, c, d in rows:
        # for ref (seq0)
        if abs(a) > 0 and abs(b) > 0:
            add_interval(ref_intervals, None, abs(a), abs(b), ref_contigs)
        # for qry (seq1)
        if abs(c) > 0 and abs(d) > 0:
            add_interval(qry_intervals, None, abs(c), abs(d), qry_contigs)

    out_prefix = Path(args.out_prefix)
    r_mat, r_unm = write_segments_for_genome(ref_contigs, ref_intervals, Path(args.ref), out_prefix, 'ATCC')
    q_mat, q_unm = write_segments_for_genome(qry_contigs, qry_intervals, Path(args.qry), out_prefix, 'ST')

    print(f"Wrote mauve segments. ATCC matched: {r_mat} bp, unmatched: {r_unm} bp")
    print(f"ST matched: {q_mat} bp, unmatched: {q_unm} bp")


if __name__ == '__main__':
    main()
