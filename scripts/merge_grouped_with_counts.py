#!/usr/bin/env python3
from pathlib import Path

def read_grouped(path):
    text = Path(path).read_text(encoding='utf-8', errors='surrogateescape').splitlines()
    header = text[0].split('\t')
    rows = {}
    for ln in text[1:]:
        if not ln.strip():
            continue
        parts = ln.split('\t')
        if len(parts) < len(header):
            parts += [''] * (len(header)-len(parts))
        d = dict(zip(header, parts))
        prod = d.get('Product','')
        rows[prod] = {'Matches': int(d.get('Matches','0') or 0), 'data': d}
    return header, rows

ha, ra = read_grouped('ATCC_grouped_by_product_with_counts.tsv')
hs, rs = read_grouped('ST_grouped_by_product_with_counts.tsv')
cols_a = [c for c in ha if c not in ('Product','Matches')]
cols_s = [c for c in hs if c not in ('Product','Matches')]
all_products = sorted(set(list(ra.keys()) + list(rs.keys())))
outcols = ['Product','ATCC_Matches','ST_Matches','Total_Matches'] + [f'ATCC_{c}' for c in cols_a] + [f'ST_{c}' for c in cols_s]
with open('combined_grouped_with_counts.tsv','w',encoding='utf-8') as fh:
    fh.write('\t'.join(outcols) + '\n')
    for prod in all_products:
        a = ra.get(prod)
        s = rs.get(prod)
        a_matches = a['Matches'] if a else 0
        s_matches = s['Matches'] if s else 0
        total = a_matches + s_matches
        row = [prod, str(a_matches), str(s_matches), str(total)]
        for c in cols_a:
            row.append(a['data'].get(c,'') if a else '')
        for c in cols_s:
            row.append(s['data'].get(c,'') if s else '')
        fh.write('\t'.join(row) + '\n')
