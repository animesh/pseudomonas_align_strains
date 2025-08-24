#!/usr/bin/env python3
"""Plot Mauve backbone as a simple dotplot (headless PNG).

Reads `alignment.mauve.backbone` with four columns per line:
 ref_left ref_right qry_left qry_right
Negative coordinates indicate reverse orientation.

Produces `alignment.mauve.backbone.png` in the current directory.
"""
import sys
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def read_backbone(path):
    pts = []
    with open(path, 'rt') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) < 4:
                continue
            try:
                a1 = int(parts[0])
                a2 = int(parts[1])
                b1 = int(parts[2])
                b2 = int(parts[3])
            except ValueError:
                continue
            # compute midpoints and orientation
            arev = (a1 < 0) or (a2 < 0)
            brev = (b1 < 0) or (b2 < 0)
            amid = (abs(a1) + abs(a2)) / 2.0
            bmid = (abs(b1) + abs(b2)) / 2.0
            pts.append((amid, bmid, arev or brev))
    return pts


def plot(pts, out):
    if not pts:
        print('No backbone points found in file', file=sys.stderr)
        return 2
    x = [p[0] for p in pts]
    y = [p[1] for p in pts]
    color = ['tab:orange' if p[2] else 'tab:blue' for p in pts]

    plt.figure(figsize=(10, 8))
    plt.scatter(x, y, c=color, s=8, edgecolors='none')
    plt.xlabel('Reference (A) coordinate')
    plt.ylabel('Query (B) coordinate')
    plt.title('Mauve backbone dotplot (midpoints)')
    # legend proxy
    import matplotlib.patches as mpatches
    blue = mpatches.Patch(color='tab:blue', label='same orientation')
    orange = mpatches.Patch(color='tab:orange', label='reverse orientation')
    plt.legend(handles=[blue, orange])
    plt.grid(alpha=0.2)
    plt.tight_layout()
    plt.savefig(out, dpi=150)
    print('Wrote', out)
    return 0


def main():
    bk = Path('alignment.mauve.backbone')
    if not bk.exists():
        print('alignment.mauve.backbone not found in current directory', file=sys.stderr)
        return 1
    pts = read_backbone(bk)
    return plot(pts, 'alignment.mauve.backbone.png')


if __name__ == '__main__':
    raise SystemExit(main())
