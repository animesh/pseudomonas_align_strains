# 🧬 Pseudomonas Genome Alignment Project

> **Compare the genomes of Pseudomonas aeruginosa ATCC 27853 and ST235 with easy, reproducible steps and a beautiful dotplot visualization highlighting the VIM gene.**

---

## 🚀 Complete Workflow (2024)

### 📝 Project Overview
This guide walks you through downloading, aligning, and visualizing the genomes of two Pseudomonas aeruginosa strains using MUMmer, with special focus on identifying and highlighting the VIM gene (blaVIM-2) that is unique to ST235.

---

### 🛠️ Prerequisites
- Linux/WSL2 environment
- [MUMmer](https://github.com/mummer4/mummer) (v4+ recommended)
- [BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download) (for gene searches)
- `awk`, `gnuplot`, `wget` (standard on most Linux systems)

#### 📦 Install Required Tools
```bash
sudo apt-get update
sudo apt-get install mummer gnuplot ncbi-blast+
```

---

### 📥 1. Download Genome FASTA Files
```bash
# Download ATCC 27853 genome
wget -O ATCC.fasta "https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_001687285.1/"

# Download ST235 genome  
wget -O ST.fasta "https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_016923535.1/"

# Download VIM gene (blaVIM-2) for analysis
curl -L 'https://www.ncbi.nlm.nih.gov/nuccore/GU137304.1?report=fasta&format=text' -o blaVIM2.fasta
```

---

### 🧩 2. Align Genomes with MUMmer
```bash
nucmer --prefix=ATCC_vs_ST ATCC.fasta ST.fasta
delta-filter -1 ATCC_vs_ST.delta > ATCC_vs_ST.filtered.delta
show-coords -rcl ATCC_vs_ST.filtered.delta > ATCC_vs_ST.coords
```

---

### 🔍 3. Search for VIM Gene
```bash
# Create BLAST databases
makeblastdb -in ATCC.fasta -dbtype nucl
makeblastdb -in ST.fasta -dbtype nucl

# Search for VIM gene in both genomes
blastn -query blaVIM2.fasta -db ATCC.fasta -out VIM_in_ATCC.txt -outfmt 6
blastn -query blaVIM2.fasta -db ST.fasta -out VIM_in_ST.txt -outfmt 6

# Check results
echo "VIM gene hits in ATCC:"
wc -l VIM_in_ATCC.txt
echo "VIM gene hits in ST235:"
wc -l VIM_in_ST.txt
```

---

### 📊 4. Summarize Alignment
```bash
awk 'NR>5 {aligned+=$7; refcov+=$11; querycov+=$12; blocks++} END {print "Total aligned bases:", aligned; print "Alignment blocks:", blocks; print "Avg % identity:", "see below"}' ATCC_vs_ST.coords
awk 'NR>5 {id+=$8*$7; len+=$7} END {if(len>0) print "Weighted avg % identity:", id/len; else print "No alignments"}' ATCC_vs_ST.coords
```

---

### 🖼️ 5. Generate Basic Dotplot
```bash
mummerplot --png --large --layout --filter --prefix=ATCC_vs_ST_plot ATCC_vs_ST.filtered.delta
# Fix gnuplot error if it occurs
sed -i '/set mouse clipboardformat/d' ATCC_vs_ST_plot.gp
gnuplot ATCC_vs_ST_plot.gp
```

---

### 🎨 6. Create Enhanced Dotplot with VIM Gene Highlight
```bash
# Add VIM gene annotation and styling to the gnuplot script
cat >> ATCC_vs_ST_plot.gp << 'EOF'

# Clear legend and annotations
set key at graph 0.02, graph 0.98 left font "Arial,12"

# Add red arrow pointing to VIM gene location
set arrow 1 from graph 0.9, first 7074888 to graph 0.95, first 7074888 lc rgb "red" lw 3 head filled
set label 1 "VIM gene" at graph 0.85, first 7074888 center tc rgb "red" font "Arial,12,bold"

# Add grey dots for unique regions
set object 2 circle at graph 0.99, first 7072697 size graph 0.005 fc rgb "#808080" fillstyle solid noborder
set object 3 circle at first 50000, graph 0.005 size graph 0.005 fc rgb "#808080" fillstyle solid noborder
set object 4 circle at first 650000, graph 0.005 size graph 0.005 fc rgb "#808080" fillstyle solid noborder
set object 5 circle at first 1250000, graph 0.005 size graph 0.005 fc rgb "#808080" fillstyle solid noborder

# Add title and better axis labels
set title "Genome Alignment: ATCC 27853 vs ST235\nCyan=Forward alignments, Purple=Reverse alignments, Red=VIM gene, Grey=Unique regions" font "Arial,14,bold"
set xlabel "ATCC 27853 genome position" font "Arial,12"
set ylabel "ST235 genome position" font "Arial,12"
EOF

# Regenerate the enhanced plot
gnuplot ATCC_vs_ST_plot.gp
```

---

### 📈 Results Summary
| Metric                      | Value         |
|-----------------------------|--------------|
| **Total aligned bases**     | 6,123,079    |
| **Alignment blocks**        | 298          |
| **Weighted avg % identity** | >98%         |
| **VIM gene in ATCC**        | Not found    |
| **VIM gene in ST235**       | Found (contig JAFFXY010000040.1) |

**Key Findings:**
- The two genomes are highly similar, with large syntenic blocks and high sequence identity
- The VIM gene (blaVIM-2) is present in ST235 but absent in ATCC 27853
- VIM gene location: contig JAFFXY010000040.1, positions 4477-5332

#### 🖼️ Enhanced Alignment Dotplot
<p align="center">
  <img src="ATCC_vs_ST_plot.png" alt="Alignment Dotplot with VIM Gene Highlight" width="600"/>
</p>

**Plot Legend:**
- **Cyan dots**: Forward alignments (same orientation)
- **Purple dots**: Reverse alignments (inverted)
- **Red arrow**: VIM gene location (ST235 only)
- **Grey dots**: Unique regions in each genome

---

## 📚 References & Further Info
- ATCC 27853: [GCA_001687285.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_001687285.1/)
- ST235: [GCA_016923535.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_016923535.1/)
- VIM gene: [GU137304.1](https://www.ncbi.nlm.nih.gov/nuccore/GU137304.1)
- MUMmer: [https://github.com/mummer4/mummer](https://github.com/mummer4/mummer)
- BLAST+: [https://blast.ncbi.nlm.nih.gov/](https://blast.ncbi.nlm.nih.gov/)

---

## 🧑‍🔬 Legacy/Alternative Data & Methods

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

