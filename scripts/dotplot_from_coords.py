#!/usr/bin/env python3
"""
Generate a dotplot from a tab-delimited coords file produced by the aligner.
Usage:
  python3 scripts/dotplot_from_coords.py <coords.tab> [--out out.png] [--min-len N]

Output: PNG file with reference positions on X and query positions on Y. Each alignment
is drawn as a line from (S1,S2) to (E1,E2). Color indicates percent identity.
"""
import sys
from pathlib import Path
import argparse
import re

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


def parse_coords_tab(path):
    text = Path(path).read_text(encoding='utf-8', errors='surrogateescape')
    lines = text.splitlines()
    # find header
    hidx = None
    for i, l in enumerate(lines[:40]):
        if '% IDY' in l or '% ID' in l or '%ID' in l or '[% IDY]' in l:
            hidx = i
            break
    if hidx is None:
        for i, l in enumerate(lines[:60]):
            if 'LEN' in l and 'ID' in l:
                hidx = i
                break
    if hidx is None:
        raise SystemExit('Could not find header in coords.tab')
    header = [h.strip() for h in lines[hidx].split('\t') if h.strip()]
    expected = len(header)

    rows = []
    i = hidx + 1
    while i < len(lines):
        cur = lines[i]
        parts = cur.split('\t')
        j = i + 1
        # if wrapped TAGS or long fields, merge next lines until we have expected columns
        while len(parts) < expected and j < len(lines):
            cur = cur + ' ' + lines[j]
            parts = cur.split('\t')
            j += 1
        if len(parts) >= expected:
            rows.append(parts[:expected])
            i = j
        else:
            # skip malformed
            i += 1

    # find columns
    cols = {h.lower(): idx for idx, h in enumerate(header)}
    # possible names
    def find(names):
        for n in names:
            for h, idx in cols.items():
                if n in h:
                    return idx
        return None

    s1_idx = find(['s1', '[s1]'])
    e1_idx = find(['e1', '[e1]'])
    s2_idx = find(['s2', '[s2]'])
    e2_idx = find(['e2', '[e2]'])
    pid_idx = find(['% idy', '% id', '%idy', '%id'])

    # fallback: typical order seen in header
    if None in (s1_idx, e1_idx, s2_idx, e2_idx, pid_idx):
        # try by positions: expected header in example is [S1][E1][S2][E2][LEN 1][LEN 2][% IDY]...
        # so S1=0,E1=1,S2=2,E2=3, PID=6
        s1_idx = s1_idx if s1_idx is not None else 0
        e1_idx = e1_idx if e1_idx is not None else 1
        s2_idx = s2_idx if s2_idx is not None else 2
        e2_idx = e2_idx if e2_idx is not None else 3
        pid_idx = pid_idx if pid_idx is not None else 6

    out = []
    for parts in rows:
        try:
            s1 = int(re.sub('[^0-9\-]', '', parts[s1_idx]))
            e1 = int(re.sub('[^0-9\-]', '', parts[e1_idx]))
            s2 = int(re.sub('[^0-9\-]', '', parts[s2_idx]))
            e2 = int(re.sub('[^0-9\-]', '', parts[e2_idx]))
            pid = float(re.sub('[^0-9\.\-]', '', parts[pid_idx]))
        except Exception:
            continue
        out.append((s1, e1, s2, e2, pid))
    return out


def plot_dotplot(alignments, outpath, min_len=0, figsize=(10,10), mode='lines'):
    # alignments: list of (s1,e1,s2,e2,pid)
    lines = []
    pids = []
    lens = []
    for s1,e1,s2,e2,pid in alignments:
        L = abs(e1 - s1) + 1
        if L < min_len:
            continue
        lines.append(((s1, s2), (e1, e2)))
        pids.append(pid)
        lens.append(L)
    if not lines:
        raise SystemExit('No alignments pass min_len filter')
    pids = np.array(pids)
    lens = np.array(lens)

    # set up plot
    fig, ax = plt.subplots(figsize=figsize)

    cmap = plt.get_cmap('viridis')
    norm = matplotlib.colors.Normalize(vmin=max(0, pids.min()), vmax=min(100, pids.max()))
    xs = []
    ys = []
    ws = []
    for (x0,y0),(x1,y1),pid,L in zip([l[0] for l in lines],[l[1] for l in lines],pids,lens):
        # center point for point/hexbin modes
        cx = (x0 + x1) / 2.0
        cy = (y0 + y1) / 2.0
        xs.append(cx)
        ys.append(cy)
        ws.append((pid, L))

    if mode == 'lines':
        for (x0,y0),(x1,y1),pid in zip([l[0] for l in lines],[l[1] for l in lines],pids):
            ax.plot([x0,x1],[y0,y1], color=cmap(norm(pid)), linewidth=1, alpha=0.8)
    elif mode == 'points':
        # scatter by center position, size by length, color by pid
        colors = [cmap(norm(pid)) for pid, L in ws]
        sizes = [max(1, min(200, L/50)) for pid, L in ws]
        ax.scatter(xs, ys, c=colors, s=sizes, alpha=0.8, edgecolors='none')
    elif mode == 'hexbin':
        # hexbin density colored by counts; use counts weighted by length
        hb = ax.hexbin(xs, ys, C=[L for pid, L in ws], reduce_C_function=np.sum, gridsize=200, cmap='inferno')
        cbar_h = fig.colorbar(hb, ax=ax)
        cbar_h.set_label('sum alignment length (bp)')
    else:
        raise SystemExit(f'Unknown mode: {mode}')

    ax.set_xlabel('Reference position (ATCC)')
    ax.set_ylabel('Query position (ST)')
    ax.set_title('Dotplot: ATCC vs ST (lines colored by % identity)')
    # colorbar
    if mode in ('lines','points'):
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        cbar = fig.colorbar(sm, ax=ax)
        cbar.set_label('% identity')

    ax.set_aspect('auto')
    plt.tight_layout()
    plt.savefig(outpath, dpi=300)
    print('Wrote', outpath)


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('coords', help='coords.tab file')
    p.add_argument('--out', default='ATCC_vs_ST.dotplot.png')
    p.add_argument('--min-len', type=int, default=0, help='minimum alignment length to plot')
    p.add_argument('--figsize', type=float, nargs=2, default=(10,10))
    p.add_argument('--mode', choices=['lines', 'points', 'hexbin'], default='lines', help='plot mode')
    args = p.parse_args()
    aln = parse_coords_tab(args.coords)
    plot_dotplot(aln, args.out, min_len=args.min_len, figsize=tuple(args.figsize), mode=args.mode)
