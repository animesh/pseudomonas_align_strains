
Step 1: Download Genome Sequences
bash
# Create working directory
mkdir pa_comparison && cd pa_comparison

# Download ATCC 27853 genome
wget -O atcc27853.fasta "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?tool=portal&save=file&log$=seqview&db=nuccore&report=fasta&id=CP009001.1"

# For ST235 strain, you'll need to choose one. Popular options:
# ST235 strain VRFPA04 (complete genome)
wget -O st235.fasta "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?tool=portal&save=file&log$=seqview&db=nuccore&report=fasta&id=CP003149.1"
Step 2: Install Required Tools
bash
# Using conda (recommended)
conda create -n genome_comparison
conda activate genome_comparison

# Install tools
conda install -c bioconda mummer samtools prokka roary blast
conda install -c conda-forge matplotlib seaborn

# Or using apt/yum for Ubuntu/CentOS
sudo apt install mummer samtools
Step 3: Genome Alignment with NUCmer
bash
# Create NUCmer alignment
nucmer --prefix=pa_comparison atcc27853.fasta st235.fasta

# Generate alignment coordinates
show-coords -rcl pa_comparison.delta > pa_comparison.coords

# Generate SNP/indel summary
show-snps -Clr pa_comparison.delta > pa_comparison.snps

# Generate alignment statistics
dnadiff -p pa_comparison atcc27853.fasta st235.fasta
Step 4: Visualization
bash
# Generate dot plot
mummerplot --fat --layout --filter -p pa_comparison pa_comparison.delta

# Convert to PNG (requires gnuplot)
sudo apt install gnuplot
gnuplot pa_comparison.gp
Step 5: Detailed Analysis
bash
# Extract unaligned regions (potential insertions/deletions)
show-diff pa_comparison.delta > pa_comparison.diff

# Get alignment coverage statistics
show-tiling pa_comparison.delta > pa_comparison.tiling
Step 6: Gene-Level Comparison (Optional)
bash
# Annotate both genomes
sudo apt install prokka
prokka --outdir atcc27853_annotation --prefix atcc27853 atcc27853.fasta
prokka --outdir st235_annotation --prefix st235 st235.fasta

# Compare gene content using BLAST
makeblastdb -in atcc27853_annotation/atcc27853.ffn -dbtype nucl -out atcc27853_genes
blastn -query st235_annotation/st235.ffn -db atcc27853_genes -outfmt 6 -out gene_comparison.blast
Interpreting Key Output Files
pa_comparison.report - Overall alignment statistics
pa_comparison.snps - SNPs and small indels
pa_comparison.coords - Alignment coordinates
pa_comparison.[png/ps] - Dot plot visualization
Quick Analysis Script
Here's a Python script to summarize key differences:

python
#!/usr/bin/env python3
import pandas as pd

# Read SNPs
snps = pd.read_csv('pa_comparison.snps', sep='\t', skiprows=5, 
                   names=['P1', 'SUB', 'P2', 'BUFF', 'DIST', 'LEN_R', 'LEN_Q', 'FRM', 'TAGS'])

print(f"Total SNPs/indels: {len(snps)}")
print(f"SNPs: {len(snps[snps['SUB'] != '.'])}")
print(f"Indels: {len(snps[snps['SUB'] == '.'])}")

# Read alignment stats
with open('pa_comparison.report', 'r') as f:
    content = f.read()
    print("\nAlignment Summary:")
    print(content)
Would you like me to help you with:

Selecting a specific ST235 strain genome?
Setting up any of these tools?
Interpreting the results once you have them?
Adding resistance gene or virulence factor analysis?






