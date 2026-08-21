#!/bin/bash
# build_k21_CYP2D6.sh
# Build a k=21 Jellyfish DB for hg38, then extract genome-unique 21-mers
# for the CYP2D6 locus using ffind.py.
#
# Outputs:
#   reference/hg38_orig_reference.k21.jf       (k=21 Jellyfish DB, ~8-10G)
#   reference/hg38CYP2D6.unique_k21.txt        (one 21-mer per line)
#
# Skip-if-present on both outputs.
#
#SBATCH --job-name=build_k21_CYP
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=logs/build_k21_CYP2D6_%j.out

set -euo pipefail
export PATH="$HOME/.conda/envs/biotools/bin:$PATH"

BASE="${COASM_BASE:-/nfs/roberts/project/pi_hc878/jst72}"
REF="$BASE/reference/hg38_orig_reference.fasta"
JF_DB="$BASE/reference/hg38_orig_reference.k21.jf"
KMERS_OUT="$BASE/reference/hg38CYP2D6.unique_k21.txt"
FFIND="$BASE/scripts/ffind.py"

# CYP2D6 gene body from genes.txt (col 3-4; ffind.py finds the flank)
CHROM=chr22
GENE_START=42126499
GENE_END=42130810

echo "=== Build k=21 Jellyfish DB + CYP2D6 unique 21-mers ==="
echo "Started: $(date)"
echo ""

# --- Step 1: build k=21 Jellyfish DB ---
if [[ -f "$JF_DB" && -s "$JF_DB" ]]; then
    echo "SKIP DB build — $JF_DB already exists ($(du -sh "$JF_DB" | cut -f1))"
else
    echo "--- Building k=21 Jellyfish DB ---"
    echo "Reference: $REF"
    jellyfish count -m 21 -s 4G -t 8 -o "$JF_DB" "$REF"
    echo "  Done: $(du -sh "$JF_DB" | cut -f1)"
fi
echo ""

# --- Step 2: run ffind.py for CYP2D6 at k=21 ---
if [[ -f "$KMERS_OUT" && -s "$KMERS_OUT" ]]; then
    echo "SKIP ffind — $KMERS_OUT already exists ($(wc -l < "$KMERS_OUT") k-mers)"
else
    echo "--- Running ffind.py for CYP2D6 (k=21, t=1, m=10000) ---"
    echo "Region: $CHROM:$GENE_START-$GENE_END"
    python3 "$FFIND" "$REF" "$CHROM" "$GENE_START" "$GENE_END" \
        -k 21 -t 1 -m 10000 \
        --db "$JF_DB" \
        --output-kmers "$KMERS_OUT"
    echo "  => $(wc -l < "$KMERS_OUT") unique 21-mers written to $KMERS_OUT"
fi

echo ""
echo "Done: $(date)"
