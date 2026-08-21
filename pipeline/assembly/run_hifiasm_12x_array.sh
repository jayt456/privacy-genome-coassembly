#!/bin/bash
# Run hifiasm on the remaining 10 pedigree 12x samples (NA12878 already done).
# One SLURM array task per sample.
#
#SBATCH --job-name=hifiasm_12x
#SBATCH --array=0-9
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --output=logs/hifiasm_12x_%a_%j.out

set -euo pipefail

module load hifiasm/0.25.0-GCCcore-13.3.0

samples=(NA12877 NA12879 NA12881 NA12882 NA12885 NA12886 NA12889 NA12890 NA12891 NA12892)
SAMPLE="${samples[$SLURM_ARRAY_TASK_ID]}"

BASE="${COASM_BASE:-/nfs/roberts/project/pi_hc878/jst72}"
READS="$BASE/pedigree_reads/wholegenome_12x/${SAMPLE}/${SAMPLE}_12x.fastq.gz"
OUTDIR="$BASE/pedigree_results/hifiasm_12x/${SAMPLE}"
PREFIX="${OUTDIR}/${SAMPLE}_12x"

echo "=== hifiasm: $SAMPLE (task $SLURM_ARRAY_TASK_ID) ==="
echo "Input:  $READS"
echo "Output: $PREFIX.*"

[[ ! -f "$READS" ]] && { echo "ERROR: reads not found: $READS"; exit 1; }

if [[ -f "${PREFIX}.bp.hap1.p_ctg.gfa" && -f "${PREFIX}.bp.hap2.p_ctg.gfa" ]]; then
    echo "Output GFAs already exist — skipping hifiasm"
else
    mkdir -p "$OUTDIR"
    hifiasm -o "$PREFIX" -t 32 "$READS"
fi

echo ""
echo "=== Converting GFAs to FASTA ==="
for gfa in "${PREFIX}".bp.hap1.p_ctg.gfa \
            "${PREFIX}".bp.hap2.p_ctg.gfa \
            "${PREFIX}".bp.p_ctg.gfa; do
    [[ ! -f "$gfa" ]] && { echo "WARNING: $gfa not found, skipping"; continue; }
    fa="${gfa%.gfa}.fa"
    if [[ -f "$fa" ]]; then
        echo "  $fa already exists, skipping"
    else
        awk '/^S/{print ">"$2; print $3}' "$gfa" > "$fa"
        echo "  $(basename $fa): $(grep -c '^>' "$fa") contigs"
    fi
done

echo ""
echo "=== Assembly stats ==="
for fa in "${PREFIX}".bp.hap1.p_ctg.fa \
           "${PREFIX}".bp.hap2.p_ctg.fa \
           "${PREFIX}".bp.p_ctg.fa; do
    [[ ! -f "$fa" ]] && continue
    ctgs=$(grep -c '^>' "$fa")
    total_bp=$(grep -v '^>' "$fa" | tr -d '\n' | wc -c)
    echo "  $(basename $fa): $ctgs contigs, $total_bp bp"
done

echo ""
echo "Done: $SAMPLE"
