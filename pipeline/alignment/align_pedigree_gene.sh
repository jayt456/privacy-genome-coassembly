#!/bin/bash
# SLURM worker: align one pedigree sample (both haps) to one gene reference.
# Launched as a job array by submit_align_pedigree.sh.
# Required env var: GENE (e.g. CYP2D6)
#
#SBATCH --job-name=ped_gene
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=logs/ped_%x_%A_%a.out

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

samples=(NA12877 NA12878 NA12879 NA12881 NA12882 NA12885 NA12886 NA12889 NA12890 NA12891 NA12892)
sample="${samples[$SLURM_ARRAY_TASK_ID]}"

out_dir="$BASE/pedigree_results/${REF_NAME}/${sample}"
mkdir -p "$out_dir"

asm_dir="$BASE/pedigree_assemblies/${sample}"

echo "=== $sample -> $REF_NAME (SLURM task $SLURM_ARRAY_TASK_ID) ==="

for hap in hap1 hap2; do
    fasta=$(find "$asm_dir" -maxdepth 1 -type f \
            -name "*.${hap}.p_ctg.gfa.fasta" ! -name "*.fai" | head -1)
    if [[ -z "$fasta" ]]; then
        echo "[$hap] no FASTA found in $asm_dir, skipping"
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
