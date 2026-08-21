#!/bin/bash
# Stream and subsample HG002 HiFi reads to 12x from S3, storing the FASTQ in scratch.
# Output: $SCRATCH/hprc_results/hifiasm_12x/HG002/reads_12x/HG002_12x.fastq.gz
# Does NOT delete the FASTQ — stored for targeted assembly experiments.
#
#SBATCH --job-name=HG002_12x_store
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --output=logs/HG002_12x_store_%j.out

set -euo pipefail
module load awscli
module load SAMtools/1.21-GCC-13.3.0

SCRATCH="${COASM_SCRATCH:-/nfs/roberts/scratch/pi_hc878/jst72}"
OUT="$SCRATCH/hprc_results/hifiasm_12x/HG002/reads_12x"
FASTQ="$OUT/HG002_12x.fastq.gz"

BAM_URL="https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC_PLUS/HG002/analysis/aligned_reads/hifi/GRCh38/HG002_aligned_GRCh38_winnowmap.sorted.bam"

TARGET_COV=12
THREADS=16
SEED=42
EST_REGIONS="chr1:50000000-51000000 chr4:60000000-61000000 chr7:60000000-61000000 chr11:60000000-61000000 chr15:50000000-51000000 chr19:30000000-31000000"

mkdir -p "$OUT"

echo "=== HG002 12x FASTQ storage run ==="
echo "BAM: $BAM_URL"
echo "Output: $FASTQ"
echo "Started: $(date)"
echo ""

if [[ -f "$FASTQ" && -s "$FASTQ" ]]; then
    sz=$(du -sh "$FASTQ" | cut -f1)
    echo "FASTQ already exists ($sz) — nothing to do."
    exit 0
fi

# ---- 1. Estimate coverage ----
echo "--- Estimating coverage from 6 windows ---"
cov=$(for r in $EST_REGIONS; do
    samtools depth -a -r "$r" "$BAM_URL" 2>/dev/null
done | awk '{sum+=$3; n++} END {if(n>0) print sum/n; else print 0}')
printf "  estimated coverage: %.1fx\n" "$cov"

if (( $(echo "$cov <= 0" | bc -l) )); then
    echo "ERROR: coverage estimate 0 — BAM not accessible"
    exit 1
fi

# ---- 2. Compute subsample fraction ----
frac=$(echo "$TARGET_COV / $cov" | bc -l)
if (( $(echo "$frac >= 1" | bc -l) )); then
    echo "  coverage already <= ${TARGET_COV}x; keeping all reads"
    seed_frac=""
else
    seed_frac=$(echo "$frac" | awk -v sd="$SEED" '{printf "%d%.4f", sd, $1}')
    printf "  fraction: %.4f  (seed.frac: %s)\n" "$frac" "$seed_frac"
fi

# ---- 3. Stream + subsample → FASTQ ----
echo ""
echo "--- Streaming from S3 ---"
if [[ -z "$seed_frac" ]]; then
    samtools view -b -@ "$THREADS" "$BAM_URL" \
        | samtools fastq -@ "$THREADS" - | gzip > "$FASTQ"
else
    samtools view -b -@ "$THREADS" -s "$seed_frac" "$BAM_URL" \
        | samtools fastq -@ "$THREADS" - | gzip > "$FASTQ"
fi

if [[ ! -s "$FASTQ" ]]; then
    echo "ERROR: output FASTQ is empty"
    exit 1
fi

sz=$(du -sh "$FASTQ" | cut -f1)
nreads=$(zcat "$FASTQ" | awk 'NR%4==1{c++} END{print c}')
echo ""
echo "=== Done ==="
echo "  file:   $FASTQ"
echo "  size:   $sz"
echo "  reads:  $nreads"
echo "Finished: $(date)"
