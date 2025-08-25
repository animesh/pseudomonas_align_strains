#!/usr/bin/env python3
from pathlib import Path
from matplotlib_venn import venn2
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

lines = Path('combined_grouped_with_counts_marked.tsv').read_text(encoding='utf-8', errors='surrogateescape').splitlines()
header = lines[0].split('\t')
src_idx = header.index('Source')
counts = {'ATCC_only':0,'ST_only':0,'both':0}
for ln in lines[1:]:
    if not ln.strip():
        continue
    parts = ln.split('\t')
    src = parts[src_idx]
    if src in counts:
        counts[src]+=1
set1 = counts['ATCC_only'] + counts['both']
set2 = counts['ST_only'] + counts['both']
inter = counts['both']
plt.figure(figsize=(6,6))
venn = venn2(subsets=(set1-inter, set2-inter, inter), set_labels=('ATCC','ST'))
if venn.get_label_by_id('10'):
    venn.get_label_by_id('10').set_text(str(counts['ATCC_only']))
if venn.get_label_by_id('01'):
    venn.get_label_by_id('01').set_text(str(counts['ST_only']))
if venn.get_label_by_id('11'):
    venn.get_label_by_id('11').set_text(str(counts['both']))
plt.title('Product overlap: ATCC vs ST')
plt.savefig('venn_proportional_ATCC_ST.png', dpi=200, bbox_inches='tight')
