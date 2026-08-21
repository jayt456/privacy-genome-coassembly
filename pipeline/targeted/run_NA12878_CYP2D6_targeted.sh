#!/bin/bash
# run_NA12878_CYP2D6_targeted.sh
#
# Targeted CYP2D6 assembly for NA12878 12x reads (Baris suggestion).
#
# Hypothesis: pre-filtering reads before assembly prevents hifiasm from
# discarding low-coverage CYP2D6 data at genome scale, so both haplotypes
# are more likely to survive into the final assembly.
#
# Pipeline:
#   1. Filter NA12878 12x FASTQ — keep reads with ≥1 genome-unique CYP2D6 21-mer
#   2. Assemble filtered reads with hifiasm
#   3. Extract CYP2D6 subgraph from targeted assembly (k=21)
#   4. Extract CYP2D6 subgraph from whole-genome 12x assembly (k=21) for comparison
#
# Output: pedigree_results/hifiasm_targeted_CYP2D6/NA12878/
#
#SBATCH --job-name=NA12878_CYP2D6_targeted
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --output=logs/NA12878_CYP2D6_targeted_%j.out

set -euo pipefail
export PATH="$HOME/.conda/envs/biotools/bin:$PATH"
module load awscli
module load hifiasm/0.25.0-GCCcore-13.3.0

BASE="${COASM_BASE:-/nfs/roberts/project/pi_hc878/jst72}"
SCRIPTS="$BASE/scripts"
KMERS="$BASE/reference/hg38CYP2D6.unique_k21.txt"

FASTQ="$BASE/pedigree_reads/wholegenome_12x/NA12878/NA12878_12x.fastq.gz"
WHOLE_PUTG="$BASE/pedigree_results/hifiasm_12x/NA12878/NA12878_12x.bp.p_utg.gfa"

OUT="$BASE/pedigree_results/hifiasm_targeted_CYP2D6/NA12878"
mkdir -p "$OUT"

FILTERED="$OUT/NA12878_CYP2D6.fastq.gz"
ASM="$OUT/NA12878_CYP2D6_targeted"
PUTG="$ASM.bp.p_utg.gfa"

echo "=== NA12878 CYP2D6 targeted assembly ==="
echo "Started: $(date)"
echo ""

# --- Step 1: filter reads ---
echo "--- Step 1: filter reads ---"
if [[ -f "$FILTERED" && -s "$FILTERED" ]]; then
    echo "  SKIP — $FILTERED exists"
else
    python3 "$SCRIPTS/filter_reads_by_kmers.py" \
        --fastq "$FASTQ" \
        --kmers "$KMERS" \
        --out   "$FILTERED"
fi
echo "  filtered reads: $(zcat "$FILTERED" | awk 'NR%4==1{c++} END{print c}')"
echo ""

# --- Step 2: assemble filtered reads ---
echo "--- Step 2: hifiasm on filtered reads ---"
if [[ -f "$PUTG" && -s "$PUTG" ]]; then
    echo "  SKIP — $PUTG exists"
else
    hifiasm -o "$ASM" -t 8 "$FILTERED" 2>&1
fi
echo "  p_utg segments: $(grep -c '^S' "$PUTG")"
echo "  p_utg total bp: $(awk '/^S/{sum+=length($3)} END{print sum}' "$PUTG")"
echo ""

# --- Step 3: k=21 subgraph — targeted assembly ---
echo "--- Step 3: k=21 subgraph (targeted) ---"
SEGS_T="$OUT/NA12878_CYP2D6_targeted.k21.segments"
GFA_T="$OUT/NA12878_CYP2D6_targeted.k21.subset.gfa"

python3 "$SCRIPTS/gfa_k21.py" --kmers "$KMERS" "$PUTG" > "$SEGS_T"
gfatools view -l "@$SEGS_T" -r 0 "$PUTG" > "$GFA_T"

echo "  segments: $(grep -c '^S' "$GFA_T")"
echo "  links:    $(grep -c '^L' "$GFA_T" || true)"
echo "  total bp: $(awk '/^S/{sum+=length($3)} END{print sum}' "$GFA_T")"
echo ""

# --- Step 4: k=21 subgraph — whole-genome 12x (baseline) ---
echo "--- Step 4: k=21 subgraph (whole-genome 12x baseline) ---"
SEGS_W="$OUT/NA12878_12x_whole.k21.segments"
GFA_W="$OUT/NA12878_12x_whole.k21.subset.gfa"

python3 "$SCRIPTS/gfa_k21.py" --kmers "$KMERS" "$WHOLE_PUTG" > "$SEGS_W"
gfatools view -l "@$SEGS_W" -r 0 "$WHOLE_PUTG" > "$GFA_W"

echo "  segments: $(grep -c '^S' "$GFA_W")"
echo "  links:    $(grep -c '^L' "$GFA_W" || true)"
echo "  total bp: $(awk '/^S/{sum+=length($3)} END{print sum}' "$GFA_W")"
echo ""

# --- Summary ---
echo "=== Comparison: targeted vs whole-genome 12x ==="
printf "  %-30s  segs=%s  bp=%s\n" \
    "targeted" \
    "$(grep -c '^S' "$GFA_T")" \
    "$(awk '/^S/{sum+=length($3)} END{print sum}' "$GFA_T")"
printf "  %-30s  segs=%s  bp=%s\n" \
    "whole-genome 12x" \
    "$(grep -c '^S' "$GFA_W")" \
    "$(awk '/^S/{sum+=length($3)} END{print sum}' "$GFA_W")"
echo ""
echo "Done: $(date)"
echo "Results: $OUT"
