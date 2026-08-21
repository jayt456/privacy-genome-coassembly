#!/bin/bash
# SLURM worker: MAPQ-filter + merge all haplotypes for one gene × dataset.
# Final output: $BASE/cohort_bams/<DATASET>_<GENE>_complete.bam + .bai
# Intermediates (*_subsample.bam, *_combined.bam) are deleted on success.
# Required env vars: GENE (e.g. CYP2D6), DATASET (pedigree or hprc)
#
#SBATCH --job-name=merge_cohort
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --output=logs/merge_%x_%j.out

set -euo pipefail

module load awscli
module load SAMtools/1.21-GCC-13.3.0

MAPQ=21
THREADS=8
BASE="${COASM_BASE:-/nfs/roberts/project/pi_hc878/jst72}"
REF_NAME="hg38${GENE}"
OUT_DIR="$BASE/cohort_bams"
FINAL="$OUT_DIR/${DATASET}_${GENE}_complete.bam"

mkdir -p "$OUT_DIR"

if [[ -f "$FINAL" && -f "${FINAL}.bai" ]]; then
    echo "Output already exists: $FINAL — skipping"
    exit 0
fi

if [[ "$DATASET" == "pedigree" ]]; then
    RESULTS="$BASE/pedigree_results/${REF_NAME}"
    samples=(NA12877 NA12878 NA12879 NA12881 NA12882 NA12885 NA12886 NA12889 NA12890 NA12891 NA12892)
    haps=(hap1 hap2)
elif [[ "$DATASET" == "hprc" ]]; then
    RESULTS="$BASE/hprc_results/${REF_NAME}"
    samples=(
        HG002    HG005    HG00438  HG00621  HG00673  HG00733  HG00735  HG00741
        HG01071  HG01106  HG01109  HG01123  HG01175  HG01243  HG01258  HG01358
        HG01361  HG01891  HG01928  HG01952  HG01978  HG02055  HG02080  HG02109
        HG02145  HG02148  HG02257  HG02486  HG02559  HG02572  HG02622  HG02630
        HG02717  HG02723  HG02818  HG02886  HG03098  HG03453  HG03486  HG03492
        HG03516  HG03540  HG03579  NA18906  NA19240  NA20129  NA21309
    )
    haps=(paternal maternal)
else
    echo "ERROR: DATASET must be 'pedigree' or 'hprc', got: '$DATASET'"
    exit 1
fi

echo "=== Merging $DATASET ${GENE} (MAPQ >= $MAPQ) ==="
echo "Source: $RESULTS"
echo "Output: $FINAL"
echo ""

combined_bams=()
intermediate_files=()

for s in "${samples[@]}"; do
    echo "--- $s ---"
    combined="$RESULTS/$s/${s}_combined.bam"

    if [[ -f "$combined" && -s "$combined" ]]; then
        echo "  Reusing existing combined.bam"
        # Add any leftover subsample BAMs from a prior run to cleanup list
        for hap in "${haps[@]}"; do
            sub_bam="$RESULTS/$s/${s}_${hap}_subsample.bam"
            [[ -f "$sub_bam" ]] && intermediate_files+=("$sub_bam")
        done
        combined_bams+=("$combined")
        intermediate_files+=("$combined")
        continue
    fi

    hap_bams=()

    for hap in "${haps[@]}"; do
        in_bam="$RESULTS/$s/${s}_${hap}_to_${REF_NAME}.bam"
        if [[ ! -f "$in_bam" ]]; then
            echo "  WARNING: $in_bam not found, skipping"
            continue
        fi

        sub_bam="$RESULTS/$s/${s}_${hap}_subsample.bam"
        echo "  Filtering $hap (MAPQ >= $MAPQ) -> $(basename "$sub_bam")"

        samtools view -@ "$THREADS" -bhq "$MAPQ" "$in_bam" \
            | samtools addreplacerg -O bam \
                -r "ID:${s}_${hap}" -r "SM:${s}" -r "LB:${hap}" - \
            | samtools sort -@ "$THREADS" -o "$sub_bam" -

        hap_bams+=("$sub_bam")
        intermediate_files+=("$sub_bam")
    done

    if [[ "${#hap_bams[@]}" -eq 0 ]]; then
        echo "  WARNING: no haplotype BAMs for $s, skipping"
        continue
    fi

    echo "  Merging ${#hap_bams[@]} hap(s) -> $(basename "$combined")"
    samtools merge -f -@ "$THREADS" -o "$combined" "${hap_bams[@]}"
    combined_bams+=("$combined")
    intermediate_files+=("$combined")
done

echo ""
echo "=== Merging ${#combined_bams[@]} samples -> $(basename "$FINAL") ==="
samtools merge -f -@ "$THREADS" -o "$FINAL" "${combined_bams[@]}"

echo "=== Indexing ==="
samtools index "$FINAL"

echo "=== Cleaning up ${#intermediate_files[@]} intermediate files ==="
rm -f "${intermediate_files[@]}"

echo ""
echo "Done: $FINAL ($(du -sh "$FINAL" | cut -f1))"
