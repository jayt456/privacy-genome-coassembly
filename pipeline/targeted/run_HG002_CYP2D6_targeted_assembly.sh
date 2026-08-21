#!/bin/bash
# Targeted CYP2D6 assembly for HG002 (Baris suggestion).
#
# Pre-filters stored 12x FASTQ by genome-unique CYP2D6 21-mers, assembles the
# filtered reads with hifiasm, then compares the k=21 subgraph from the targeted
# assembly against the k=21 subgraph from the whole-genome 12x assembly.
#
# Input:    $SCRATCH/hprc_results/hifiasm_12x/HG002/reads_12x/HG002_12x.fastq.gz
#           (produced by run_HG002_12x_store_fastq.sh — must exist before submitting)
# Output:   $SCRATCH/hprc_results/hifiasm_12x/HG002/targeted_CYP2D6/
#
#SBATCH --job-name=HG002_CYP2D6_targeted
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --output=logs/HG002_CYP2D6_targeted_%j.out

set -euo pipefail
export PATH="$HOME/.conda/envs/biotools/bin:$PATH"
module load hifiasm/0.25.0-GCCcore-13.3.0

BASE="${COASM_BASE:-/nfs/roberts/project/pi_hc878/jst72}"
SCRATCH="${COASM_SCRATCH:-/nfs/roberts/scratch/pi_hc878/jst72}"
SCRIPTS="$BASE/scripts"
KMERS="$BASE/reference/hg38CYP2D6.unique_k21.txt"

FASTQ_IN="$SCRATCH/hprc_results/hifiasm_12x/HG002/reads_12x/HG002_12x.fastq.gz"
WHOLE_PUTG="$SCRATCH/hprc_results/hifiasm_12x/HG002/assembly/HG002_12x.bp.p_utg.gfa"

OUT="$SCRATCH/hprc_results/hifiasm_12x/HG002/targeted_CYP2D6"
mkdir -p "$OUT"

FILTERED="$OUT/HG002_CYP2D6.fastq.gz"
ASM_PREFIX="$OUT/HG002_CYP2D6_targeted"
PUTG="$ASM_PREFIX.bp.p_utg.gfa"

echo "=== HG002 CYP2D6 targeted assembly ==="
echo "Started: $(date)"
echo ""

# Require input FASTQ (produced by run_HG002_12x_store_fastq.sh)
if [[ ! -s "$FASTQ_IN" ]]; then
    echo "ERROR: input FASTQ not found: $FASTQ_IN"
    echo "Run run_HG002_12x_store_fastq.sh first."
    exit 1
fi

# ---- Step 1: filter reads ----
echo "--- Step 1: filter reads by CYP2D6 k=21 k-mers ---"
if [[ -f "$FILTERED" && -s "$FILTERED" ]]; then
    echo "  SKIP — $FILTERED exists"
else
    python3 "$SCRIPTS/filter_reads_by_kmers.py" \
        --fastq "$FASTQ_IN" \
        --kmers "$KMERS" \
        --out   "$FILTERED"
fi
nreads=$(zcat "$FILTERED" | awk 'NR%4==1{c++} END{print c}')
echo "  filtered reads: $nreads"
echo ""

# ---- Step 2: hifiasm on filtered reads ----
echo "--- Step 2: hifiasm on CYP2D6-filtered reads (8 threads) ---"
if [[ -f "$PUTG" && -s "$PUTG" ]]; then
    echo "  SKIP — $PUTG exists"
else
    hifiasm -o "$ASM_PREFIX" -t 8 "$FILTERED" 2>&1
fi
echo "  p_utg segments: $(grep -c '^S' "$PUTG")"
echo "  p_utg total bp: $(awk '/^S/{sum+=length($3)} END{print sum}' "$PUTG")"
echo ""

# GFA → FASTA for hap1/hap2/p_ctg/p_utg
echo "--- GFA → FASTA ---"
for gfa in "${ASM_PREFIX}".bp.hap1.p_ctg.gfa \
           "${ASM_PREFIX}".bp.hap2.p_ctg.gfa \
           "${ASM_PREFIX}".bp.p_ctg.gfa \
           "${ASM_PREFIX}".bp.p_utg.gfa; do
    [[ ! -f "$gfa" ]] && { echo "  WARNING: $gfa missing"; continue; }
    fa="${gfa%.gfa}.fa"
    if [[ -f "$fa" ]]; then
        echo "  $(basename $fa) already exists"
    else
        awk '/^S/{print ">"$2; print $3}' "$gfa" > "$fa"
        echo "  $(basename $fa): $(grep -c '^>' "$fa") contigs"
    fi
done
echo ""

# ---- Step 3: k=21 subgraph — targeted assembly ----
echo "--- Step 3: k=21 subgraph (targeted) ---"
SEGS_T="$OUT/HG002_CYP2D6_targeted.k21.segments"
GFA_T="$OUT/HG002_CYP2D6_targeted.k21.subset.gfa"

if [[ -f "$GFA_T" && -s "$GFA_T" ]]; then
    echo "  SKIP — $GFA_T exists"
else
    python3 "$SCRIPTS/gfa_k21.py" --kmers "$KMERS" "$PUTG" > "$SEGS_T"
    gfatools view -l "@$SEGS_T" -r 0 "$PUTG" > "$GFA_T"
fi
echo "  segments: $(grep -c '^S' "$GFA_T")"
echo "  links:    $(grep -c '^L' "$GFA_T" || true)"
echo "  total bp: $(awk '/^S/{sum+=length($3)} END{print sum}' "$GFA_T")"
echo ""

# ---- Step 4: k=21 subgraph — whole-genome 12x baseline ----
echo "--- Step 4: k=21 subgraph (whole-genome 12x baseline) ---"
SEGS_W="$OUT/HG002_12x_whole.k21.segments"
GFA_W="$OUT/HG002_12x_whole.k21.subset.gfa"

if [[ -f "$GFA_W" && -s "$GFA_W" ]]; then
    echo "  SKIP — $GFA_W exists"
else
    python3 "$SCRIPTS/gfa_k21.py" --kmers "$KMERS" "$WHOLE_PUTG" > "$SEGS_W"
    gfatools view -l "@$SEGS_W" -r 0 "$WHOLE_PUTG" > "$GFA_W"
fi
echo "  segments: $(grep -c '^S' "$GFA_W")"
echo "  links:    $(grep -c '^L' "$GFA_W" || true)"
echo "  total bp: $(awk '/^S/{sum+=length($3)} END{print sum}' "$GFA_W")"
echo ""

# ---- Summary ----
echo "=== Comparison: targeted vs whole-genome 12x ==="
printf "  %-28s  segs=%-4s  bp=%s\n" \
    "targeted (pre-filtered)" \
    "$(grep -c '^S' "$GFA_T")" \
    "$(awk '/^S/{sum+=length($3)} END{print sum}' "$GFA_T")"
printf "  %-28s  segs=%-4s  bp=%s\n" \
    "whole-genome 12x" \
    "$(grep -c '^S' "$GFA_W")" \
    "$(awk '/^S/{sum+=length($3)} END{print sum}' "$GFA_W")"
echo ""
echo "Done: $(date)"
echo "Results: $OUT"
