#!/bin/bash
# Submit pedigree alignment jobs for all (or one) gene(s) in genes.txt.
# Each gene becomes a SLURM job array with 11 tasks (one per sample).
# All genes run in parallel; within each gene all 11 samples run in parallel.
#
# Usage:
#   bash submit_align_pedigree.sh           # submit all genes in genes.txt
#   bash submit_align_pedigree.sh CYP2D6    # submit one specific gene

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
GENES_FILE="$SCRIPT_DIR/../genes.txt"
WORKER="$SCRIPT_DIR/align_pedigree_gene.sh"
BASE="${COASM_BASE:-/nfs/roberts/project/pi_hc878/jst72}"
N_SAMPLES=11
FILTER="${1:-}"

if [[ ! -f "$GENES_FILE" ]]; then
    echo "ERROR: genes.txt not found at $GENES_FILE"
    exit 1
fi

while read -r gene chrom gene_start gene_end flank_start flank_end; do
    # skip comments and blank lines
    [[ -z "$gene" || "$gene" =~ ^# ]] && continue
    # if a gene filter was given, skip non-matching genes
    [[ -n "$FILTER" && "$gene" != "$FILTER" ]] && continue

    ref="$BASE/reference/hg38${gene}.fa"
    if [[ ! -f "$ref" ]]; then
        echo "SKIP $gene: reference not found ($ref)"
        echo "  Run: bash scripts/extract_region.sh $gene ${chrom}:${flank_start}-${flank_end}"
        continue
    fi

    echo "Submitting pedigree alignment for $gene ($N_SAMPLES samples in parallel)"
    jobid=$(sbatch \
        --array=0-$((N_SAMPLES - 1)) \
        --job-name="ped_${gene}" \
        --output="$SCRIPT_DIR/logs/ped_${gene}_%A_%a.out" \
        --export=GENE="$gene" \
        "$WORKER" | awk '{print $NF}')
    echo "  Job array ID: $jobid  (logs: scripts/logs/ped_${gene}_${jobid}_*.out)"
done < "$GENES_FILE"
