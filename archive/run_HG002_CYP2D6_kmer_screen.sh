#!/bin/bash
# run_HG002_CYP2D6_kmer_screen.sh
#
# Build a 61-mer reference set from the 4 CYP2D6 regions in the HG002 high-cov
# boundary subgraph, screen all HG002 12x p_utg segments for k-mer matches,
# then align the found segments (as contigs) to:
#   (a) Q100 HG002 CYP2D6 pat — ground truth paternal gene body (4.3 kb)
#   (b) Q100 HG002 CYP2D6 mat — ground truth maternal gene body (4.3 kb)
#   (c) hg38 CYP2D6 flank window — same coordinate space as cohort BAMs (IGV-ready)
#
# All alignments produce sorted, indexed BAMs for IGV loading.
#
#SBATCH --job-name=HG002_kmer_screen
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=logs/HG002_kmer_screen_%j.out

set -euo pipefail
export PATH="$HOME/.conda/envs/biotools/bin:$PATH"
module load awscli
module load minimap2/2.29-GCCcore-13.3.0
module load SAMtools/1.21-GCC-13.3.0

HIGHCOV=/nfs/roberts/scratch/pi_hc878/jst72/hprc_results/highcov/HG002
ASM12X=/nfs/roberts/scratch/pi_hc878/jst72/hprc_results/hifiasm_12x/HG002
REF=/nfs/roberts/project/pi_hc878/jst72/reference
SCRIPTS=/nfs/roberts/project/pi_hc878/jst72/scripts
OUT="$ASM12X/HG002_12x.CYP2D6.kmer61"
TMP="$ASM12X/sort_tmp"

echo "=== HG002 CYP2D6 k-mer screen (61-mer, high-cov boundary subgraph) ==="
echo "Started: $(date)"
echo ""

# --- Step 1: k-mer screen (skip if FASTA already exists) ---
if [[ -f "${OUT}.fa" && -s "${OUT}.fa" ]]; then
    echo "SKIP k-mer screen — ${OUT}.fa already exists"
else
    python3 "$SCRIPTS/kmer_screen_12x.py" \
        --highcov-gfa "$HIGHCOV/HG002_highcov.CYP2D6.boundary.subset.gfa" \
        --alines      "$HIGHCOV/HG002_highcov.CYP2D6.hap1.boundary.lines" \
                      "$HIGHCOV/HG002_highcov.CYP2D6.hap2.boundary.lines" \
        --gfa-12x     "$ASM12X/HG002_12x.bp.p_utg.gfa" \
        -k 61 \
        --min-hits 2 \
        --out-names   "${OUT}.names" \
        --out-fasta   "${OUT}.fa"
fi

echo "Kept segments in FASTA: $(grep -c '^>' "${OUT}.fa")"
echo ""

# Helper: align FASTA to REF, emit sorted indexed BAM of mapped reads only
# Usage: align_bam <ref.fa> <query.fa> <out.bam>
align_bam() {
    local ref=$1 query=$2 bam=$3
    if [[ -f "$bam" && -s "$bam" ]]; then
        echo "  SKIP — $bam already exists"
        return
    fi
    minimap2 -a -x asm5 -t 4 "$ref" "$query" \
        | samtools view -F 4 -b \
        | samtools sort -T "$TMP" -o "$bam"
    samtools index "$bam"
}

# --- Step 2a: vs Q100 CYP2D6 paternal (4.3 kb gene body) ---
echo "--- Aligning to Q100 HG002 CYP2D6 (pat) ---"
align_bam "$HIGHCOV/hg002.cyp2d6.pat.fa" \
          "${OUT}.fa" \
          "${OUT}.vs_Q100_pat.bam"
samtools flagstat "${OUT}.vs_Q100_pat.bam"

echo ""

# --- Step 2b: vs Q100 CYP2D6 maternal (4.3 kb gene body) ---
echo "--- Aligning to Q100 HG002 CYP2D6 (mat) ---"
align_bam "$HIGHCOV/hg002.cyp2d6.mat.fa" \
          "${OUT}.fa" \
          "${OUT}.vs_Q100_mat.bam"
samtools flagstat "${OUT}.vs_Q100_mat.bam"

echo ""

# --- Step 2c: vs hg38 CYP2D6 flank window (286 kb, IGV-ready) ---
echo "--- Aligning to hg38 CYP2D6 flank window ---"
align_bam "$REF/hg38CYP2D6.fa" \
          "${OUT}.fa" \
          "${OUT}.vs_hg38CYP2D6.bam"
samtools flagstat "${OUT}.vs_hg38CYP2D6.bam"

echo ""
echo "Done: $(date)"
echo "BAMs:"
ls -lh "${OUT}".vs_*.bam
