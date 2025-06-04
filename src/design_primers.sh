#!/usr/bin/env bash
set -euo pipefail

# Install dependencies
pip install -r ../requirements.txt
pip install primer3-py

# Paths
FASTA_DIR="../outputs/fastas"
PARAMS="../src/params.json"
OUT_DIR="../outputs/primers"

if [[ ! -d "$FASTA_DIR" ]]; then
  echo "Error: FASTA directory '$FASTA_DIR' does not exist" >&2
  exit 1
fi

shopt -s nullglob
fastas=("$FASTA_DIR"/*_cds.fasta)
if (( ${#fastas[@]} == 0 )); then
  echo "Error: no FASTA files found in $FASTA_DIR" >&2
  exit 1
fi

mkdir -p "$OUT_DIR"

# Run primer design for each CDS FASTA
for fasta in "${fastas[@]}"; do
  python ../src/design_primers.py "$fasta" "$PARAMS" "$OUT_DIR"
done
