#!/usr/bin/env bash
set -e

# Install dependencies
pip install -r ../requirements.txt
pip install primer3-py

# Paths
FASTA_DIR="../outputs/fastas"
PARAMS="../src/params.json"
OUT_DIR="../outputs/primers"

mkdir -p "$OUT_DIR"

# Run primer design for each CDS FASTA
for fasta in "$FASTA_DIR"/*_cds.fasta; do
  python ../src/design_primers.py "$fasta" "$PARAMS" "$OUT_DIR"
done
