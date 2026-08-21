#!/bin/bash
# run_HG002_CYP2D6_trace.sh
# A-line read tracing for HG002 highcov × CYP2D6 subgraph extraction.
#
# Implements Baris's approach:
#   1. Use MAPQ60 contig-to-Q100 BAMs to find contig offsets of CYP2D6 locus
#   2. Trace reads through p_ctg.gfa A-lines (position-filtered)
#   3. Find those reads in p_utg.gfa A-lines → collect unitig names
#   4. gfatools view → subgraph.gfa
#
# Runs hap1 + hap2 sequentially, merges utg lists, extracts combined subgraph.
# Output: HG002_highcov.CYP2D6.trace.subset.gfa
#
#SBATCH --job-name=HG002_trace
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --output=logs/HG002_trace_%j.out

set -euo pipefail
export PATH="$HOME/.conda/envs/biotools/bin:$PATH"

module load awscli
module load SAMtools/1.21-GCC-13.3.0

BASE="${COASM_BASE:-/nfs/roberts/project/pi_hc878/jst72}"
SCRATCH="${COASM_SCRATCH:-/nfs/roberts/scratch/pi_hc878/jst72}"
WORKDIR="$SCRATCH/hprc_results/highcov/HG002"
SCRIPT="$BASE/scripts/trace_reads.py"
UTG="$WORKDIR/HG002.asm.bp.p_utg.gfa"
FLANK=150000  # ±150 kb around the 4.3 kb CYP2D6 gene body

# Output files
UTG_LIST="$WORKDIR/HG002_highcov.CYP2D6.trace.utgs"
GFA_OUT="$WORKDIR/HG002_highcov.CYP2D6.trace.subset.gfa"

echo "=== HG002 CYP2D6 read-tracing subgraph ==="
echo "UTG GFA: $UTG"
echo "Flank:   ±${FLANK} bp"
echo "Started: $(date)"

if [[ -f "$GFA_OUT" && -s "$GFA_OUT" ]]; then
    echo "SKIP — output already exists: $GFA_OUT"
    exit 0
fi

# --- hap1 ---
echo ""
echo "--- hap1 ---"
python3 "$SCRIPT" \
    --bam  "$WORKDIR/hg002.hap1_to_pat.CYP2D6.MAPQ60.bam" \
    --ctg-gfa "$WORKDIR/HG002.asm.bp.hap1.p_ctg.gfa" \
    --utg-gfa "$UTG" \
    --flank "$FLANK" \
    > "$WORKDIR/HG002_highcov.CYP2D6.hap1.trace.utgs"

echo "  hap1 utgs: $(wc -l < "$WORKDIR/HG002_highcov.CYP2D6.hap1.trace.utgs")"

# --- hap2 ---
echo ""
echo "--- hap2 ---"
python3 "$SCRIPT" \
    --bam  "$WORKDIR/hg002.hap2_to_mat.CYP2D6.MAPQ60.bam" \
    --ctg-gfa "$WORKDIR/HG002.asm.bp.hap2.p_ctg.gfa" \
    --utg-gfa "$UTG" \
    --flank "$FLANK" \
    > "$WORKDIR/HG002_highcov.CYP2D6.hap2.trace.utgs"

echo "  hap2 utgs: $(wc -l < "$WORKDIR/HG002_highcov.CYP2D6.hap2.trace.utgs")"

# --- merge utg lists and extract subgraph ---
echo ""
echo "--- merging utg lists ---"
sort -u \
    "$WORKDIR/HG002_highcov.CYP2D6.hap1.trace.utgs" \
    "$WORKDIR/HG002_highcov.CYP2D6.hap2.trace.utgs" \
    > "$UTG_LIST"

echo "  combined utgs (deduplicated): $(wc -l < "$UTG_LIST")"

echo ""
echo "--- extracting subgraph with gfatools ---"
gfatools view -l "@$UTG_LIST" -r 0 "$UTG" > "$GFA_OUT"

echo "  subgraph segments: $(grep -c '^S' "$GFA_OUT")"
echo "  subgraph links:    $(grep -c '^L' "$GFA_OUT" || true)"
echo ""
echo "Done: $(date)"
echo "Output: $GFA_OUT"
