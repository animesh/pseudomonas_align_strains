#!/usr/bin/env python3
"""Compute ATCC x ST intersection rectangles for blaVIM2 matches.

Reads:
 - pa_comparison.coords (ATCC vs ST)
 - bla_ST.coords (ST vs blaVIM2)

Writes:
 - pa_comparison.rects.gp (a gnuplot script that draws precise rects)

Run from repository root: python3 scripts/make_rects.py
"""
import math
from pathlib import Path
import argparse
import logging
import subprocess

# Default repo root (two levels up from this script: repo_root/scripts/..)
ROOT = Path(__file__).resolve().parents[1]



def parse_coords(fn, reverse_tags=False):
    """Parse a nuccmer-style .coords file into a list of segments.

    Returns list of dicts: {at_s, at_e, st_s, st_e, at_tag, st_tag}
    For files where the reference/query order differs, set reverse_tags accordingly.
    """
    segs = []
    with open(fn) as f:
        started = False
        for line in f:
            if not started:
                if line.strip().startswith('==='):
                    started = True
                continue
            line = line.rstrip('\n')
            if not line.strip():
                continue
            # skip header separators
            if line.lstrip().startswith('='):
                continue
            if '|' not in line:
                continue
            parts = [p.strip() for p in line.split('|')]
            if len(parts) < 4:
                continue
            # parts[0]: ATCC S1 E1
            # parts[1]: ST   S2 E2
            # parts[-1]: TAGS (two tokens)
            try:
                at_s, at_e = [int(x) for x in parts[0].split()[:2]]
                st_s, st_e = [int(x) for x in parts[1].split()[:2]]
            except Exception:
                continue
            tags = parts[-1].split()
            if len(tags) >= 2:
                if not reverse_tags:
                    at_tag = tags[0]
                    st_tag = tags[1]
                else:
                    # some .coords have the query/ref reversed in tags
                    at_tag = tags[1]
                    st_tag = tags[0]
            else:
                at_tag = st_tag = ''
            segs.append({
                'at_s': at_s,
                'at_e': at_e,
                'st_s': st_s,
                'st_e': st_e,
                'at_tag': at_tag,
                'st_tag': st_tag,
            })
    return segs


def parse_bla_st(fn):
    """Parse bla_ST.coords where the coordinates are ST | BLA.

    Returns list of dicts: {st_start, st_end, st_tag}
    """
    matches = []
    with open(fn) as f:
        started = False
        for line in f:
            if not started:
                if line.strip().startswith('==='):
                    started = True
                continue
            if not line.strip():
                continue
            if '|' not in line:
                continue
            parts = [p.strip() for p in line.split('|')]
            if len(parts) < 3:
                continue
            try:
                st_s, st_e = [int(x) for x in parts[0].split()[:2]]
            except Exception:
                continue
            tags = parts[-1].split()
            st_tag = tags[0] if tags else ''
            matches.append({'st_start': st_s, 'st_end': st_e, 'st_tag': st_tag})
    return matches


def parse_unqry(fn):
    """Parse pa_comparison.unqry into ordered contig offsets.

    Returns dict contig -> (start_global, length) where start_global is 1-based global offset
    assuming concatenation order in the file.
    """
    offsets = {}
    cur = 1
    with open(fn) as f:
        for line in f:
            parts = line.strip().split()
            if not parts:
                continue
            name = parts[0]
            # last numeric field is length
            try:
                length = int(parts[-1])
            except Exception:
                continue
            offsets[name] = (cur, length)
            cur += length
    return offsets


def map_overlap(pa_seg, bla_st_range):
    """Given one pa_comparison segment and a bla ST-range, compute the AT overlap interval
    corresponding to the overlapping portion (if any) and return (at_lo, at_hi) as floats.
    Returns None if no overlap.
    Handles reversed orientations.
    """
    st_a = pa_seg['st_s']
    st_b = pa_seg['st_e']
    st_min = min(st_a, st_b)
    st_max = max(st_a, st_b)
    bla_min = min(bla_st_range['st_start'], bla_st_range['st_end'])
    bla_max = max(bla_st_range['st_start'], bla_st_range['st_end'])
    ov_s = max(st_min, bla_min)
    ov_e = min(st_max, bla_max)
    if ov_e < ov_s:
        return None

    st_len = st_max - st_min
    if st_len == 0:
        return None

    # fraction along the pa segment for overlap endpoints
    if st_a <= st_b:
        frac0 = (ov_s - st_min) / st_len
        frac1 = (ov_e - st_min) / st_len
    else:
        # reversed
        frac0 = (st_max - ov_s) / st_len
        frac1 = (st_max - ov_e) / st_len

    at_a = pa_seg['at_s']
    at_b = pa_seg['at_e']
    at_min = min(at_a, at_b)
    at_max = max(at_a, at_b)
    at_len = at_max - at_min
    if at_len == 0:
        # degenerate -> single coordinate
        at0 = at_min
        at1 = at_max
    else:
        if at_a <= at_b:
            at0 = at_min + frac0 * at_len
            at1 = at_min + frac1 * at_len
        else:
            at0 = at_max - frac0 * at_len
            at1 = at_max - frac1 * at_len

    at_lo = min(at0, at1)
    at_hi = max(at0, at1)
    return (at_lo, at_hi, ov_s, ov_e)


def map_overlap_global(pa_seg, bla_match, unqry):
    """Try to map a bla ST match (which is on a specific ST contig) to the given pa_seg
    by converting both to global ST coordinates using `unqry`. If the pa_seg's ST contig
    lacks unqry info, return None. If there is a global overlap, compute the corresponding
    local ST coordinates within the pa_seg and call map_overlap() to get AT interval.
    Returns (at_lo, at_hi, global_ov_s, global_ov_e) or None.
    """
    st_tag_match = bla_match['st_tag']
    st_tag_seg = pa_seg['st_tag']
    if st_tag_seg not in unqry or st_tag_match not in unqry:
        return None

    seg_offset, seg_len = unqry[st_tag_seg]
    match_offset, match_len = unqry[st_tag_match]

    seg_local_min = min(pa_seg['st_s'], pa_seg['st_e'])
    seg_local_max = max(pa_seg['st_s'], pa_seg['st_e'])
    seg_global_s = seg_offset + seg_local_min - 1
    seg_global_e = seg_offset + seg_local_max - 1

    match_local_min = min(bla_match['st_start'], bla_match['st_end'])
    match_local_max = max(bla_match['st_start'], bla_match['st_end'])
    match_global_s = match_offset + match_local_min - 1
    match_global_e = match_offset + match_local_max - 1

    ov_s = max(seg_global_s, match_global_s)
    ov_e = min(seg_global_e, match_global_e)
    if ov_e < ov_s:
        return None

    # convert overlap back to local coordinates on the pa_seg's contig
    ov_local_s = ov_s - seg_offset + 1
    ov_local_e = ov_e - seg_offset + 1
    local_bla_range = {'st_start': ov_local_s, 'st_end': ov_local_e}
    mapped = map_overlap(pa_seg, local_bla_range)
    if not mapped:
        return None
    at_lo, at_hi, ov_local_s2, ov_local_e2 = mapped
    # return AT interval and global overlap
    return (at_lo, at_hi, ov_s, ov_e)


def main(pa_coords_path=None):
    parser = argparse.ArgumentParser(description='Compute rects for bla matches and emit a gnuplot script')
    parser.add_argument('--pa-coords', dest='pa_coords', default=str(ROOT / 'pa_comparison.coords'), help='path to pa_comparison.coords')
    parser.add_argument('--bla-st', dest='bla_st', default=str(ROOT / 'bla_ST.coords'), help='path to bla_ST.coords')
    parser.add_argument('--unqry', dest='unqry', default=str(ROOT / 'pa_comparison.unqry'), help='path to pa_comparison.unqry')
    parser.add_argument('--out-gp', dest='out_gp', default=str(ROOT / 'pa_comparison.rects.gp'), help='output gnuplot script path')
    parser.add_argument('--out-png', dest='out_png', default=str(ROOT / 'pa_comparison.rects.png'), help='output PNG path referenced in gnuplot script')
    parser.add_argument('--reverse-tags', dest='reverse_tags', action='store_true', help='treat tags as reversed in pa coords')
    parser.add_argument('--export-svg', dest='export_svg', action='store_true', help='also emit SVG via gnuplot')
    parser.add_argument('--export-pdf', dest='export_pdf', action='store_true', help='also emit PDF via gnuplot')
    args = parser.parse_args(pa_coords_path.split() if isinstance(pa_coords_path, str) and pa_coords_path else None)

    pa_path = Path(args.pa_coords)
    bla_path = Path(args.bla_st)
    gp_out = Path(args.out_gp)
    png_out = Path(args.out_png)
    unqry_path = Path(args.unqry)

    logging.basicConfig(level=logging.INFO)
    log = logging.getLogger('make_rects')

    pa_segs = parse_coords(pa_path, reverse_tags=args.reverse_tags)
    bla_matches = parse_bla_st(bla_path)
    unqry = parse_unqry(unqry_path) if unqry_path.exists() else {}

    rects = []
    for midx, match in enumerate(bla_matches, start=1):
        st_tag = match['st_tag']
        # find pa segments for the same ST contig
        overlaps = []
        for seg in pa_segs:
            if seg['st_tag'] != st_tag:
                continue
            mapped = map_overlap(seg, match)
            if mapped:
                at_lo, at_hi, ov_s, ov_e = mapped
                overlaps.append((at_lo, at_hi, ov_s, ov_e, seg))

        if not overlaps:
            # fallback: if the ST contig exists in unqry, compute global ST positions and
            # create a rectangle spanning the full AT range (1..max AT)
            if st_tag in unqry:
                st_offset, st_len = unqry[st_tag]
                st_global_lo = st_offset + min(match['st_start'], match['st_end']) - 1
                st_global_hi = st_offset + max(match['st_start'], match['st_end']) - 1
                # full AT range
                at_min = 1
                at_max = max( (s['at_e'] for s in pa_segs), default=1 )
                rects.append({'idx': midx, 'st_tag': st_tag, 'at_lo': at_min, 'at_hi': at_max, 'st_lo': st_global_lo, 'st_hi': st_global_hi, 'overlaps': []})
                print(f"Fallback: creating full-AT rectangle for bla match {midx} on {st_tag} {match['st_start']}-{match['st_end']} mapped to ST global {st_global_lo}-{st_global_hi}")
                continue
            print(f"No pa_comparison overlaps found for bla match {midx} on {st_tag} {match['st_start']}-{match['st_end']}")
            continue

        # union AT ranges across overlaps
        at_min = min(x[0] for x in overlaps)
        at_max = max(x[1] for x in overlaps)
        # ST rectangle should be the original bla match range
        st_lo = min(match['st_start'], match['st_end'])
        st_hi = max(match['st_start'], match['st_end'])
        rects.append({'idx': midx, 'st_tag': st_tag, 'at_lo': at_min, 'at_hi': at_max, 'st_lo': st_lo, 'st_hi': st_hi, 'overlaps': overlaps})

    # build gnuplot script (reuse plotting settings from annotated3.gp but change output)
    header = []
    annotated = ROOT / 'pa_comparison.annotated3.gp'
    if annotated.exists():
        with open(annotated) as f:
            for line in f:
                # replace output and png filename
                if line.strip().startswith('set output'):
                    header.append(f"set output '{png_out}'\n")
                elif line.strip().startswith('set terminal'):
                    header.append(line)
                else:
                    # copy until objects/rects area (we'll append our rects)
                    header.append(line)
    else:
        # minimal header
        header = [
            'set terminal pngcairo size 1600,1200 enhanced font "Courier,10"\n',
            f"set output '{png_out}'\n",
            'set size 1,1\nset grid\nunset key\nset border 10\nset tics scale 0\nset xlabel "CP011857.1"\nset ylabel "QRY (ST coordinates)"\nset format "%.0f"\nset xrange [1:6833187]\nset yrange [1:7550891]\n',
        ]

    # start writing gp
    gp_out.parent.mkdir(parents=True, exist_ok=True)
    with open(gp_out, 'w') as out:
        out.writelines(header)
        out.write('\n# Precise ATCC x ST rectangles computed from pa_comparison.coords and bla_ST.coords\n')
        for i, r in enumerate(rects, start=1):
            at_lo_i = int(math.floor(r['at_lo']))
            at_hi_i = int(math.ceil(r['at_hi']))
            st_lo_i = int(r['st_lo'])
            st_hi_i = int(r['st_hi'])
            out.write(f"# rect {i}: ST {r['st_tag']} {st_lo_i}-{st_hi_i} -> ATCC {at_lo_i}-{at_hi_i}\n")
            out.write(f"set object {100+i} rect from {at_lo_i},{st_lo_i} to {at_hi_i},{st_hi_i} fc rgb \"green\" fillstyle solid 0.25 noborder\n")
        out.write('\n# Plot the existing forward/reverse mummerplot data\n')
        out.write(f"plot \\\n+ '{ROOT / 'pa_comparison.fplot'}' title 'FWD' with lp ls 1, \\\n+ '{ROOT / 'pa_comparison.rplot'}' title 'REV' with lp ls 2\n")
    print(f"Wrote {gp_out} with {len(rects)} rectangle(s). Output PNG will be {png_out}")

    # optional vector exports
    def write_vector_gp(src_gp: Path, dest_gp: Path, terminal_line: str, out_path: Path):
        """Create a copy of src_gp replacing terminal and output lines to produce vector output.

        If the original gp lacks set terminal/set output lines, the function will prepend them.
        """
        with open(src_gp, 'r') as f:
            lines = f.readlines()

        had_terminal = False
        had_output = False
        new_lines = []
        for line in lines:
            if line.strip().startswith('set terminal') and not had_terminal:
                new_lines.append(terminal_line + '\n')
                had_terminal = True
                continue
            if line.strip().startswith('set output') and not had_output:
                new_lines.append(f"set output '{out_path}'\n")
                had_output = True
                continue
            new_lines.append(line)

        if not had_terminal:
            new_lines.insert(0, terminal_line + '\n')
        if not had_output:
            new_lines.insert(1 if had_terminal else 1, f"set output '{out_path}'\n")

        dest_gp.parent.mkdir(parents=True, exist_ok=True)
        with open(dest_gp, 'w') as f:
            f.writelines(new_lines)

    if args.export_svg:
        svg_out = gp_out.with_suffix('.svg')
        svg_gp = gp_out.with_suffix('.svg.gp')
        svg_term = 'set terminal svg size 1600,1200 enhanced font "Courier,10"'
        write_vector_gp(gp_out, svg_gp, svg_term, svg_out)
        print(f"Wrote SVG gnuplot script {svg_gp} -> will render to {svg_out}")
        rc = subprocess.call(['gnuplot', str(svg_gp)])
        print('gnuplot svg rc=', rc)

    if args.export_pdf:
        pdf_out = gp_out.with_suffix('.pdf')
        pdf_gp = gp_out.with_suffix('.pdf.gp')
        pdf_term = 'set terminal pdfcairo enhanced font "Courier,10"'
        write_vector_gp(gp_out, pdf_gp, pdf_term, pdf_out)
        print(f"Wrote PDF gnuplot script {pdf_gp} -> will render to {pdf_out}")
        rc = subprocess.call(['gnuplot', str(pdf_gp)])
        print('gnuplot pdf rc=', rc)


if __name__ == '__main__':
    import sys
    # allow passing the raw argument string for compatibility with previous usage
    arg = sys.argv[1] if len(sys.argv) > 1 else None
    main(arg)
