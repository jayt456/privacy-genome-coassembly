#!/bin/bash
# Run hifiasm on NA12878 12x whole-genome HiFi reads.
# Output: pedigree_results/hifiasm_12x/NA12878/
# Produces hap1/hap2 phased contig GFAs + FASTAs.
#
#SBATCH --job-name=hifiasm_NA12878
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --output=logs/hifiasm_NA12878_%j.out

set -euo pipefail

module load hifiasm/0.25.0-GCCcore-13.3.0

SAMPLE=NA12878
BASE="${COASM_BASE:-/nfs/roberts/project/pi_hc878/jst72}"
READS="$BASE/pedigree_reads/wholegenome_12x/${SAMPLE}/${SAMPLE}_12x.fastq.gz"
OUTDIR="$BASE/pedigree_results/hifiasm_12x/${SAMPLE}"
PREFIX="${OUTDIR}/${SAMPLE}_12x"

if [[ -f "${PREFIX}.bp.hap1.p_ctg.gfa" && -f "${PREFIX}.bp.hap2.p_ctg.gfa" ]]; then
    echo "Output GFAs already exist — skipping hifiasm"
else
    mkdir -p "$OUTDIR"
    echo "=== hifiasm: $SAMPLE (12x) ==="
    echo "Input:  $READS"
    echo "Output: $PREFIX.*"
    echo "Threads: 32  Memory: 64G"
    echo ""

    hifiasm -o "$PREFIX" -t 32 "$READS"
fi

# Convert all primary contig GFAs to FASTA
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
        echo "  $gfa -> $fa"
        awk '/^S/{print ">"$2; print $3}' "$gfa" > "$fa"
        echo "  $(grep -c '^>' "$fa") contigs"
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
    echo "  $(basename $fa): $ctgs contigs, $total_bp bp total"
done

echo ""
echo "Done. Outputs in $OUTDIR"
