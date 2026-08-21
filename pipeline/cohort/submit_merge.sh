#!/bin/bash
# Submit cohort merge jobs for all genes × datasets (16 total).
# Each job: MAPQ-filter all haplotype BAMs for one gene+dataset, merge into
# complete_subsample.bam in {pedigree,hprc}_results/hg38<GENE>/.
#
# Usage:
#   bash submit_merge.sh                   # all 8 genes × both datasets
#   bash submit_merge.sh CYP2D6            # one gene × both datasets
#   bash submit_merge.sh CYP2D6 pedigree   # one gene × one dataset

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
GENES_FILE="$SCRIPT_DIR/../genes.txt"
WORKER="$SCRIPT_DIR/merge_cohort.sh"
BASE="${COASM_BASE:-/nfs/roberts/project/pi_hc878/jst72}"
GENE_FILTER="${1:-}"
DATASET_FILTER="${2:-}"

if [[ ! -f "$GENES_FILE" ]]; then
    echo "ERROR: genes.txt not found at $GENES_FILE"
    exit 1
fi

while read -r gene chrom gene_start gene_end flank_start flank_end; do
    [[ -z "$gene" || "$gene" =~ ^# ]] && continue
    [[ -n "$GENE_FILTER" && "$gene" != "$GENE_FILTER" ]] && continue

    for dataset in pedigree hprc; do
        [[ -n "$DATASET_FILTER" && "$dataset" != "$DATASET_FILTER" ]] && continue

        results_dir="$BASE/${dataset}_results/hg38${gene}"
        final="$BASE/cohort_bams/${dataset}_${gene}_complete.bam"

        if [[ -f "$final" && -f "${final}.bai" ]]; then
            echo "SKIP $dataset $gene: complete_subsample.bam already exists"
            continue
        fi

        if [[ ! -d "$results_dir" ]]; then
            echo "SKIP $dataset $gene: results directory not found ($results_dir)"
            continue
        fi

        echo "Submitting merge: $dataset $gene"
        jobid=$(sbatch \
            --job-name="merge_${dataset}_${gene}" \
            --output="$SCRIPT_DIR/logs/merge_${dataset}_${gene}_%j.out" \
            --export=GENE="$gene",DATASET="$dataset" \
            "$WORKER" | awk '{print $NF}')
        echo "  Job ID: $jobid  (log: scripts/logs/merge_${dataset}_${gene}_${jobid}.out)"
    done
done < "$GENES_FILE"
