#!/usr/bin/env bash
set -euo pipefail
# Reproducible script to install and run Bakta annotation used in this project.
# Adjust /path/to/db-light.tar.xz to your DB archive location if different.

# Optional: create & activate conda env (recommended)
conda create -n bakta python=3.10 -y || true
# shellcheck disable=SC1091
source $(conda info --base)/etc/profile.d/conda.sh
conda activate bakta

# install bakta
python3 -m pip install --upgrade pip
python3 -m pip install bakta

# extract the provided Bakta DB bundle (change path if needed)
mkdir -p bakta_db
tar -xJf /mnt/z/Download/db-light.tar.xz -C bakta_db

# Fix C++ runtime if needed
conda install -n bakta libgcc-ng libstdcxx-ng -y || true

# Run Bakta annotation
bakta annotate --db-dir bakta_db --output bakta_ATCC --prefix ATCC ATCC.fasta
bakta annotate --db-dir bakta_db --output bakta_ST --prefix ST ST.fasta

# Copy TSV outputs to repo root for downstream scripts
cp bakta_ATCC/ATCC.tsv ./ATCC.tsv
cp bakta_ST/ST.tsv ./ST.tsv

# Done
echo "Bakta runs completed; TSVs copied to ATCC.tsv and ST.tsv"
