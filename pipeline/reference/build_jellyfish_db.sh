#!/bin/bash
# Build a genome-wide Jellyfish k-mer count database from the hg38 reference.
# Output is reused by ffind.py for all genes — only needs to be run once.
#
# Usage: sbatch scripts/build_jellyfish_db.sh
#
#SBATCH --job-name=jf_build
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --output=logs/jf_build_%j.out

set -euo pipefail

BASE="${COASM_BASE:-/nfs/roberts/project/pi_hc878/jst72}"
REF="$BASE/reference/hg38_orig_reference.fasta"
OUT="$BASE/reference/hg38_orig_reference.k15.jf"

export PATH="/home/jst72/.conda/envs/biotools/bin:$PATH"

if [[ ! -f "$REF" ]]; then
    echo "ERROR: reference not found: $REF"
    exit 1
fi

if [[ -f "$OUT" ]]; then
    echo "Database already exists: $OUT — skipping"
    exit 0
fi

echo "Building k=15 Jellyfish database from $(basename "$REF")..."
echo "Threads: $SLURM_CPUS_PER_TASK  Hash size: 3G"

jellyfish count \
    -m 15 \
    -s 3G \
    -t "$SLURM_CPUS_PER_TASK" \
    -o "$OUT" \
    "$REF"

echo "Done: $OUT ($(du -sh "$OUT" | cut -f1))"
