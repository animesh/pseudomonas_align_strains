#!/usr/bin/env bash
set -euo pipefail
# Helper to BLAST unique segment FASTA files.
# Usage:
#   scripts/blast_unique.sh ATCC.unique.fasta ST.unique.fasta
#
# By default this will run remote BLAST (requires internet) using blastn -task blastn.
# For large jobs or privacy, create a local BLAST DB and replace the remote call with -db <localdb>.

FA_FILES=("${@:-ATCC.unique.fasta ST.unique.fasta}")

for f in "${FA_FILES[@]}"; do
    if [ ! -s "$f" ]; then
        echo "File not found or empty: $f" >&2
        continue
    fi
    out="${f%.*}.blast.tsv"
    echo "Running remote BLAST for $f -> $out"
    # Remote blast; adjust -evalue and -max_target_seqs as needed
    blastn -query "$f" -db nt -remote -outfmt '6 qseqid sseqid pident length evalue bitscore stitle' -max_target_seqs 5 -evalue 1e-5 -out "$out"
    echo "Wrote $out"
done

echo "Done. If you have a local nt/nt-like DB, replace '-db nt -remote' with '-db /path/to/db' and remove -remote." 
