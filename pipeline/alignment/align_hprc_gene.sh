#!/bin/bash
# SLURM worker: align one HPRC sample (paternal + maternal) to one gene reference.
# Launched as a job array by submit_align_hprc.sh.
# Required env var: GENE (e.g. CYP2D6)
#
#SBATCH --job-name=hprc_gene
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=logs/hprc_%x_%A_%a.out

set -euo pipefail

module load awscli
module load SAMtools/1.21-GCC-13.3.0
module load minimap2/2.29-GCCcore-13.3.0

BASE="${COASM_BASE:-/nfs/roberts/project/pi_hc878/jst72}"
REF_NAME="hg38${GENE}"
REF="$BASE/reference/${REF_NAME}.fa"

if [[ ! -f "$REF" ]]; then
    echo "ERROR: reference not found: $REF"
    exit 1
fi

samples=(
    HG002    HG005    HG00438  HG00621  HG00673  HG00733  HG00735  HG00741
    HG01071  HG01106  HG01109  HG01123  HG01175  HG01243  HG01258  HG01358
    HG01361  HG01891  HG01928  HG01952  HG01978  HG02055  HG02080  HG02109
    HG02145  HG02148  HG02257  HG02486  HG02559  HG02572  HG02622  HG02630
    HG02717  HG02723  HG02818  HG02886  HG03098  HG03453  HG03486  HG03492
    HG03516  HG03540  HG03579  NA18906  NA19240  NA20129  NA21309
)
sample="${samples[$SLURM_ARRAY_TASK_ID]}"

out_dir="$BASE/hprc_results/${REF_NAME}/${sample}"
mkdir -p "$out_dir"

asm_dir="$BASE/hprc_assemblies/${sample}"

echo "=== $sample -> $REF_NAME (SLURM task $SLURM_ARRAY_TASK_ID) ==="

for hap in paternal maternal; do
    fasta="$asm_dir/${sample}.${hap}.f1_assembly_v2_genbank.fa.gz"
    if [[ ! -f "$fasta" ]]; then
        echo "[$hap] assembly not found: $fasta, skipping"
        continue
    fi

    bam="$out_dir/${sample}_${hap}_to_${REF_NAME}.bam"
    if [[ -f "$bam" && -f "${bam}.bai" ]]; then
        echo "[$hap] already done, skipping"
        continue
    fi

    echo "[$hap] aligning $(basename "$fasta") -> $bam"
    minimap2 -t 4 -ax asm5 "$REF" "$fasta" \
        | samtools sort -@ 4 -o "$bam" -

    if [[ ! -s "$bam" ]]; then
        echo "[$hap] ERROR: output BAM is empty"
        exit 1
    fi

    samtools index "$bam"
    echo "[$hap] done ($(du -sh "$bam" | cut -f1))"
done
