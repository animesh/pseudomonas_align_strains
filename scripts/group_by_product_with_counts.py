#!/usr/bin/env python3
from pathlib import Path
from collections import defaultdict

def group_by_product(inpath, outpath):
    p = Path(inpath)
    text = p.read_text(encoding='utf-8', errors='surrogateescape').splitlines()
    header_line = None
    for i, l in enumerate(text[:10]):
        if l.startswith('#') and 'Sequence' in l:
            header_line = l.lstrip('#')
            start = i
            break
    if header_line is None:
        header_line = text[0]
        start = 0
    fieldnames = [f.strip() for f in header_line.split('\t')]
    rows = []
    for l in text[start+1:]:
        if not l.strip():
            continue
        parts = l.split('\t')
        if len(parts) < len(fieldnames):
            parts += [''] * (len(fieldnames) - len(parts))
        row = dict(zip(fieldnames, parts))
        rows.append(row)
    groups = defaultdict(list)
    for r in rows:
        prod = r.get('Product','')
        groups[prod].append(r)
    out_header = ['Product','Matches'] + [c for c in fieldnames if c!='Product']
    out_lines = []
    for prod, items in sorted(groups.items(), key=lambda x: (-len(x[1]), x[0])):
        matches = len(items)
        agg = {}
        for col in fieldnames:
            vals = []
            for it in items:
                v = it.get(col,'')
                if v and v not in vals:
                    vals.append(v)
            agg[col] = ';'.join(vals)
        rowout = [prod, str(matches)] + [agg[c] for c in fieldnames if c!='Product']
        out_lines.append('\t'.join(rowout))
    outp = Path(outpath)
    outp.parent.mkdir(parents=True, exist_ok=True)
    with outp.open('w', encoding='utf-8') as fh:
        fh.write('\t'.join(out_header) + '\n')
        for l in out_lines:
            fh.write(l + '\n')

if __name__ == '__main__':
    group_by_product('ATCC.tsv','ATCC_grouped_by_product_with_counts.tsv')
    group_by_product('ST.tsv','ST_grouped_by_product_with_counts.tsv')
