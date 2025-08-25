#!/usr/bin/env python3
from pathlib import Path

lines = Path('combined_grouped_with_counts.tsv').read_text(encoding='utf-8', errors='surrogateescape').splitlines()
header = lines[0].split('\t')
with open('combined_grouped_with_counts_marked.tsv','w',encoding='utf-8') as fh:
    fh.write('\t'.join(header + ['Source']) + '\n')
    for ln in lines[1:]:
        if not ln.strip():
            continue
        parts = ln.split('\t')
        row = dict(zip(header, parts))
        a = int(row.get('ATCC_Matches','0') or 0)
        s = int(row.get('ST_Matches','0') or 0)
        if a>0 and s>0:
            src='both'
        elif a>0:
            src='ATCC_only'
        elif s>0:
            src='ST_only'
        else:
            src='none'
        fh.write('\t'.join(parts + [src]) + '\n')
