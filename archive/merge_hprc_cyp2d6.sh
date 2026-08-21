#!/bin/bash
# MAPQ-filter all 94 HPRC CYP2D6 BAMs (47 samples x paternal/maternal),
# add read groups, merge per sample, then merge all into complete_subsample.bam.
# Output: hprc_results/hg38CYP2D6/complete_subsample.bam
#
#SBATCH --job-name=hprc_merge
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=logs/hprc_merge_%j.out

set -euo pipefail

module load awscli
module load SAMtools/1.21-GCC-13.3.0

MAPQ=21
THREADS=8
BASE="${COASM_BASE:-/nfs/roberts/project/pi_hc878/jst72}"
RESULTS="$BASE/hprc_results/hg38CYP2D6"
FINAL="$RESULTS/complete_subsample.bam"

if [[ -f "$FINAL" && -f "${FINAL}.bai" ]]; then
    echo "Output already exists: $FINAL — skipping"
    exit 0
fi

samples=(
  HG002 HG005 HG00438 HG00621 HG00673 HG00733 HG00735 HG00741
  HG01071 HG01106 HG01109 HG01123 HG01175 HG01243 HG01258 HG01358
  HG01361 HG01891 HG01928 HG01952 HG01978 HG02055 HG02080 HG02109
  HG02145 HG02148 HG02257 HG02486 HG02559 HG02572 HG02622 HG02630
  HG02717 HG02723 HG02818 HG02886 HG03098 HG03453 HG03486 HG03492
  HG03516 HG03540 HG03579 NA18906 NA19240 NA20129 NA21309
)

combined_bams=()

for s in "${samples[@]}"; do
    echo "=== $s ==="
    hap_bams=()

    for hap in paternal maternal; do
        in_bam="$RESULTS/$s/${s}_${hap}_to_hg38CYP2D6.bam"
        if [[ ! -f "$in_bam" ]]; then
            echo "  WARNING: $in_bam not found, skipping"
            continue
        fi

        sub_bam="$RESULTS/$s/${s}_${hap}_subsample.bam"
        echo "  Filtering $hap (MAPQ >= $MAPQ) -> $sub_bam"

        samtools view -@ "$THREADS" -bhq "$MAPQ" "$in_bam" \
            | samtools addreplacerg -O bam \
                -r "ID:${s}_${hap}" -r "SM:${s}" -r "LB:${hap}" - \
            | samtools sort -@ "$THREADS" -o "$sub_bam" -

        hap_bams+=("$sub_bam")
    done

    if [[ "${#hap_bams[@]}" -eq 0 ]]; then
        echo "  WARNING: no haplotype BAMs for $s, skipping"
        continue
    fi

    combined="$RESULTS/$s/${s}_combined.bam"
    echo "  Merging ${#hap_bams[@]} hap(s) -> $combined"
    samtools merge -f -@ "$THREADS" -o "$combined" "${hap_bams[@]}"
    combined_bams+=("$combined")
done

echo ""
echo "=== Merging ${#combined_bams[@]} samples -> $FINAL ==="
samtools merge -f -@ "$THREADS" -o "$FINAL" "${combined_bams[@]}"

echo "=== Indexing ==="
samtools index "$FINAL"

echo ""
echo "Done: $FINAL ($(du -sh "$FINAL" | cut -f1))"
