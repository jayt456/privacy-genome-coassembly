#!/bin/bash
# SLURM worker: stream+subsample HPRC HiFi GRCh38 BAM to 12x, run hifiasm, delete FASTQ.
# Source: s3://human-pangenomics/working/HPRC_PLUS/ (18 samples) or working/HPRC/ (29 samples)
# Output: hprc_results/hifiasm_12x/<sample>/<sample>_12x.bp.hap{1,2}.p_ctg.{gfa,fa}
# Called from submit_hprc_hifiasm.sh as a SLURM array job.
#
#SBATCH --job-name=hprc_hifiasm_12x
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --output=logs/hprc_hifiasm_12x_%a_%j.out

set -euo pipefail

S3_BASE="https://s3-us-west-2.amazonaws.com/human-pangenomics/working"
TARGET_COV=12
THREADS=32
SEED=42
EST_REGIONS="chr1:50000000-51000000 chr4:60000000-61000000 chr7:60000000-61000000 chr11:60000000-61000000 chr15:50000000-51000000 chr19:30000000-31000000"
BASE="${COASM_BASE:-/nfs/roberts/project/pi_hc878/jst72}"
SCRATCH_BASE="${COASM_SCRATCH:-/nfs/roberts/scratch/pi_hc878/jst72}"

module load awscli
module load SAMtools/1.21-GCC-13.3.0
module load hifiasm/0.25.0-GCCcore-13.3.0

samples=(
    HG002    HG005    HG00438  HG00621  HG00673  HG00733  HG00735  HG00741
    HG01071  HG01106  HG01109  HG01123  HG01175  HG01243  HG01258  HG01358
    HG01361  HG01891  HG01928  HG01952  HG01978  HG02055  HG02080  HG02109
    HG02145  HG02148  HG02257  HG02486  HG02559  HG02572  HG02622  HG02630
    HG02717  HG02723  HG02818  HG02886  HG03098  HG03453  HG03486  HG03492
    HG03516  HG03540  HG03579  NA18906  NA19240  NA20129  NA21309
)

# 18 samples with GRCh38 BAMs in HPRC_PLUS; remaining 29 are in HPRC (no PLUS)
hprc_plus=(
    HG002 HG005 HG00733 HG01109 HG01243 HG02055 HG02080 HG02109
    HG02145 HG02723 HG02818 HG03098 HG03486 HG03492 NA18906 NA19240 NA20129 NA21309
)

SAMPLE="${samples[$SLURM_ARRAY_TASK_ID]}"

# Route to correct bucket prefix
if printf '%s\n' "${hprc_plus[@]}" | grep -qx "$SAMPLE"; then
    BUCKET="$S3_BASE/HPRC_PLUS"
else
    BUCKET="$S3_BASE/HPRC"
fi

OUTDIR="$SCRATCH_BASE/hprc_results/hifiasm_12x/$SAMPLE"
PREFIX="$OUTDIR/${SAMPLE}_12x"
BAM_URL="$BUCKET/$SAMPLE/analysis/aligned_reads/hifi/GRCh38/${SAMPLE}_aligned_GRCh38_winnowmap.sorted.bam"
FASTQ_TMP="$OUTDIR/${SAMPLE}_12x_tmp.fastq.gz"

mkdir -p "$OUTDIR"

echo "=== HPRC hifiasm 12x: $SAMPLE (array task $SLURM_ARRAY_TASK_ID) ==="
echo "BAM: $BAM_URL"
echo "Output: $PREFIX.*"
echo ""

# Skip if already assembled
if [[ -f "${PREFIX}.bp.hap1.p_ctg.gfa" && -f "${PREFIX}.bp.hap2.p_ctg.gfa" ]]; then
    echo "GFAs already exist — skipping"
    exit 0
fi

# ---- 1. Estimate coverage ----
echo "--- Estimating coverage from $(echo $EST_REGIONS | wc -w) windows ---"
cov=$(for r in $EST_REGIONS; do
    samtools depth -a -r "$r" "$BAM_URL" 2>/dev/null
done | awk '{sum+=$3; n++} END {if(n>0) print sum/n; else print 0}')
printf "  estimated coverage: %.1fx\n" "$cov"

if (( $(echo "$cov <= 0" | bc -l) )); then
    echo "ERROR: coverage estimate 0 — BAM may not be accessible at $BAM_URL"
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

# ---- 3. Stream + subsample → temp FASTQ ----
echo ""
if [[ -f "$FASTQ_TMP" && -s "$FASTQ_TMP" ]]; then
    sz=$(du -sh "$FASTQ_TMP" | cut -f1)
    echo "--- Reusing existing FASTQ ($sz): $(basename $FASTQ_TMP) ---"
else
    echo "--- Streaming to $(basename $FASTQ_TMP) ---"
    if [[ -z "$seed_frac" ]]; then
        samtools view -b -@ "$THREADS" "$BAM_URL" \
            | samtools fastq -@ "$THREADS" - | gzip > "$FASTQ_TMP"
    else
        samtools view -b -@ "$THREADS" -s "$seed_frac" "$BAM_URL" \
            | samtools fastq -@ "$THREADS" - | gzip > "$FASTQ_TMP"
    fi

    if [[ ! -s "$FASTQ_TMP" ]]; then
        echo "ERROR: empty output FASTQ"
        exit 1
    fi
    sz=$(du -sh "$FASTQ_TMP" | cut -f1)
    echo "  written: $sz"
fi

# ---- 4. Run hifiasm ----
echo ""
echo "--- Running hifiasm ($THREADS threads) ---"
hifiasm -o "$PREFIX" -t "$THREADS" "$FASTQ_TMP"

# ---- 5. GFA → FASTA ----
echo ""
echo "--- Converting GFAs to FASTA ---"
for gfa in "${PREFIX}".bp.hap1.p_ctg.gfa \
           "${PREFIX}".bp.hap2.p_ctg.gfa \
           "${PREFIX}".bp.p_ctg.gfa; do
    [[ ! -f "$gfa" ]] && { echo "  WARNING: $gfa missing"; continue; }
    fa="${gfa%.gfa}.fa"
    if [[ -f "$fa" ]]; then
        echo "  $(basename $fa) already exists"
    else
        awk '/^S/{print ">"$2; print $3}' "$gfa" > "$fa"
        echo "  $(basename $fa): $(grep -c '^>' "$fa") contigs"
    fi
done

# ---- 6. Cleanup: delete FASTQ + .bin files ----
echo ""
echo "--- Cleaning up ---"
rm -f "$FASTQ_TMP"
find "$OUTDIR" -name "*.bin" -delete
echo "  temp FASTQ and .bin files deleted"

# ---- Summary ----
echo ""
echo "=== Done: $SAMPLE ==="
for fa in "${PREFIX}".bp.hap1.p_ctg.fa "${PREFIX}".bp.hap2.p_ctg.fa; do
    [[ ! -f "$fa" ]] && continue
    ctgs=$(grep -c '^>' "$fa")
    total_bp=$(grep -v '^>' "$fa" | tr -d '\n' | wc -c)
    printf "  %-40s  %d contigs  %d bp\n" "$(basename $fa)" "$ctgs" "$total_bp"
done
du -sh "$OUTDIR"
