#!/bin/bash
# Build genome-unique k=21 k-mer files for the 7 remaining genes
# (CYP2D6 already done — reference/hg38CYP2D6.unique_k21.txt).
#
# Requires: reference/hg38_orig_reference.k21.jf (22 GB, already built)
# Outputs:  reference/hg38<GENE>.unique_k21.txt  for each gene
#
# Submit:
#   sbatch --array=0-6 build_k21_all_genes.sh
#   (or omit --array to let SBATCH directive below control it)
#
#SBATCH --job-name=build_k21_genes
#SBATCH --array=0-6
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --output=logs/build_k21_genes_%A_%a.out

set -euo pipefail
export PATH="$HOME/.conda/envs/biotools/bin:$PATH"

BASE="${COASM_BASE:-/nfs/roberts/project/pi_hc878/jst72}"
REF="$BASE/reference/hg38_orig_reference.fasta"
JF_DB="$BASE/reference/hg38_orig_reference.k21.jf"
FFIND="$BASE/scripts/ffind.py"

# Gene body coordinates from genes.txt (col 3-4); CYP2D6 excluded (already done)
genes=(  TERT      TCF3      SLC6A3    LMF1    IGH          KIR          HLA        )
chroms=( chr5      chr19     chr5      chr16   chr14        chr19        chr6       )
starts=( 1253147   1609290   1392794   853634  105586437    54816468     25726063   )
ends=(   1295068   1652615   1445440   981318  106879844    54830778     33400644   )

GENE="${genes[$SLURM_ARRAY_TASK_ID]}"
CHROM="${chroms[$SLURM_ARRAY_TASK_ID]}"
START="${starts[$SLURM_ARRAY_TASK_ID]}"
END="${ends[$SLURM_ARRAY_TASK_ID]}"
KMERS_OUT="$BASE/reference/hg38${GENE}.unique_k21.txt"

echo "=== Build k=21 unique k-mers: $GENE ==="
echo "Region: $CHROM:$START-$END"
echo "Output: $KMERS_OUT"
echo "Started: $(date)"
echo ""

if [[ ! -f "$JF_DB" || ! -s "$JF_DB" ]]; then
    echo "ERROR: k=21 Jellyfish DB not found: $JF_DB"
    echo "Run build_k21_CYP2D6.sh first to build the DB."
    exit 1
fi
echo "DB: $JF_DB ($(du -sh "$JF_DB" | cut -f1))"
echo ""

if [[ -f "$KMERS_OUT" && -s "$KMERS_OUT" ]]; then
    echo "SKIP — $KMERS_OUT already exists ($(wc -l < "$KMERS_OUT") k-mers)"
    exit 0
fi

echo "--- Running ffind.py (k=21, t=1, m=10000) ---"
python3 "$FFIND" "$REF" "$CHROM" "$START" "$END" \
    -k 21 -t 1 -m 10000 \
    --db "$JF_DB" \
    --output-kmers "$KMERS_OUT"

echo ""
echo "=> $(wc -l < "$KMERS_OUT") unique 21-mers written to $KMERS_OUT"
echo "Done: $(date)"
