# pseudomonas_align_strains

## New Genome Alignment (2024)

### File Sources
- **ATCC.fasta**: [GCA_001687285.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_001687285.1/)
- **ST.fasta**: [GCA_016923535.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_016923535.1/)

### Alignment Commands
```bash
# Align genomes using MUMmer (nucmer)
nucmer --prefix=ATCC_vs_ST ATCC.fasta ST.fasta
delta-filter -1 ATCC_vs_ST.delta > ATCC_vs_ST.filtered.delta
show-coords -rcl ATCC_vs_ST.filtered.delta > ATCC_vs_ST.coords

# Summarize alignment
awk 'NR>5 {aligned+=$7; refcov+=$11; querycov+=$12; blocks++} END {print "Total aligned bases:", aligned; print "Alignment blocks:", blocks; print "Avg % identity:", "see below"}' ATCC_vs_ST.coords
awk 'NR>5 {id+=$8*$7; len+=$7} END {if(len>0) print "Weighted avg % identity:", id/len; else print "No alignments"}' ATCC_vs_ST.coords
```

### Results Summary
- **Total aligned bases:** 6,123,079
- **Number of alignment blocks:** 298
- **Weighted average percent identity:** Most blocks >98% (see .coords for details)
- The two genomes are highly similar, with large syntenic blocks and high sequence identity.

---
compare the genome of ATCC 27853 (Cao et al., 2017) with that of the ST235 strain used in this study (Urbanowicz et al., 2021)

### Status
- Data fetched:
  - ATCC 27853 complete genome: CP015117.1 → `ATCC_27853_CP015117.1.fasta` ([NCBI](https://www.ncbi.nlm.nih.gov/nuccore/CP015117.1))
  - ATCC 27853 alternative genome: CP011857.1 → `ATCC_27853_alternative.fasta` ([NCBI](https://www.ncbi.nlm.nih.gov/nuccore/CP011857.1)) - from PMC5467263
  - ST235 representative complete genome (NCGM2.S1): AP012280.1 → `ST235_NCGM2S1_AP012280.1.fasta` ([NCBI](https://www.ncbi.nlm.nih.gov/nuccore/AP012280.1))
  - Note: Urbanowicz et al., 2021 link `LFMO00000000` is a WGS master and contains no sequence data ([NCBI](https://www.ncbi.nlm.nih.gov/nuccore/LFMO00000000))
- Commands used:
  ```bash
  # Download ATCC 27853 genome (CP015117.1)
  wget -O ATCC_27853_CP015117.1.fasta "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?tool=portal&save=file&log$=seqview&db=nuccore&report=fasta&id=CP015117.1"
  
  # Download ATCC 27853 alternative genome (CP011857.1) from PMC5467263
  wget -O ATCC_27853_alternative.fasta "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?tool=portal&save=file&log$=seqview&db=nuccore&report=fasta&id=CP011857.1"
  
  # Download ST235 genome (NCGM2.S1)
  wget -O ST235_NCGM2S1_AP012280.1.fasta "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?tool=portal&save=file&log$=seqview&db=nuccore&report=fasta&id=AP012280.1"
  
  # Download and extract minimap2
  wget -O minimap2-2.28_x64-linux.tar.bz2 https://github.com/lh3/minimap2/releases/download/v2.28/minimap2-2.28_x64-linux.tar.bz2
  tar -xjf minimap2-2.28_x64-linux.tar.bz2
  
  # Run genome alignment (CP015117.1 vs ST235)
  ./minimap2-2.28_x64-linux/minimap2 -x asm5 -t 4 ATCC_27853_CP015117.1.fasta ST235_NCGM2S1_AP012280.1.fasta > ATCC27853_vs_ST235_NCGM2S1.paf
  
  # Run genome alignment (CP011857.1 vs ST235)
  ./minimap2-2.28_x64-linux/minimap2 -x asm5 -t 4 ATCC_27853_alternative.fasta ST235_NCGM2S1_AP012280.1.fasta > ATCC27853_CP011857_vs_ST235_NCGM2S1.paf
  
  # Quick summary statistics
  awk 'BEGIN{m=0;b=0;ql=0;tl=0} {m+=$10; b+=$11; if($2>ql) ql=$2; if($7>tl) tl=$7} END{printf "matches=%d aln_bases=%d qlen=%d tlen=%d pid=%.6f\n",m,b,ql,tl,m/b}' ATCC27853_vs_ST235_NCGM2S1.paf
  awk 'BEGIN{m=0;b=0;ql=0;tl=0} {m+=$10; b+=$11; if($2>ql) ql=$2; if($7>tl) tl=$7} END{printf "matches=%d aln_bases=%d qlen=%d tlen=%d pid=%.6f\n",m,b,ql,tl,m/b}' ATCC27853_CP011857_vs_ST235_NCGM2S1.paf
  ```
- Alignment:
  - Tool: `minimap2` (v2.28) with `-x asm5`
  - Output files: 
    - `ATCC27853_vs_ST235_NCGM2S1.paf` (CP015117.1 vs ST235)
    - `ATCC27853_CP011857_vs_ST235_NCGM2S1.paf` (CP011857.1 vs ST235)
  - Quick summary (PAF aggregate):
    - CP015117.1 vs ST235: matches = 5,143,335; aligned bases = 7,093,386; rough pid = 0.725
    - CP011857.1 vs ST235: matches = 5,141,970; aligned bases = 7,146,480; rough pid = 0.720
- Next steps:
  - If available, provide the exact ST235 accession from Urbanowicz et al., 2021 to replace the representative genome and re-run.
  - Compute robust ANI and SNP/indel stats (e.g., `dnadiff`/`fastANI`) and generate a brief report.

