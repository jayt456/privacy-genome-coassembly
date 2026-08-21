#!/bin/bash
# k=21 subgraph extraction for HG002 x all 8 genes.
# Runs both the 12x and highcov assemblies for each gene.
# Each gene gets its own subdirectory under subgraph_k21/.
#
# Output layout:
#   hifiasm_12x/HG002/subgraph_k21/<GENE>/HG002_12x.<GENE>.k21.subset.gfa
#   highcov/HG002/subgraph_k21/<GENE>/HG002_highcov.<GENE>.k21.subset.gfa
#
# Depends on build_k21_all_genes.sh having completed for the 7 new genes.
# CYP2D6 k-mer file already exists; task 0 will produce a fresh result in
# the per-gene subdirectory (consistent structure alongside the new genes).
#
# Submit after build_k21_all_genes.sh completes:
#   JOBA=$(sbatch --parsable build_k21_all_genes.sh)
#   sbatch --dependency=afterok:$JOBA run_HG002_all_genes_k21.sh
#
#SBATCH --job-name=HG002_all_k21
#SBATCH --array=0-7
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --output=logs/HG002_all_k21_%A_%a.out

set -euo pipefail
export PATH="$HOME/.conda/envs/biotools/bin:$PATH"

BASE="${COASM_BASE:-/nfs/roberts/project/pi_hc878/jst72}"
SCRATCH="${COASM_SCRATCH:-/nfs/roberts/scratch/pi_hc878/jst72}"
GFA_SCRIPT="${COASM_TOOLS:-$BASE/src}/gfa_k21.py"

genes=(CYP2D6 TERT TCF3 SLC6A3 LMF1 IGH KIR HLA)
GENE="${genes[$SLURM_ARRAY_TASK_ID]}"

KMERS="$BASE/reference/hg38${GENE}.unique_k21.txt"

# Assembly GFA paths (post-reorganization, both in assembly/ subdirs)
PUTG_12X="$SCRATCH/hprc_results/hifiasm_12x/HG002/assembly/HG002_12x.bp.p_utg.gfa"
PUTG_HC="$SCRATCH/hprc_results/highcov/HG002/assembly/HG002.asm.bp.p_utg.gfa"

# Per-gene output directories
OUT_12X="$SCRATCH/hprc_results/hifiasm_12x/HG002/subgraph_k21/$GENE"
OUT_HC="$SCRATCH/hprc_results/highcov/HG002/subgraph_k21/$GENE"

mkdir -p "$OUT_12X" "$OUT_HC"

echo "=== HG002 k=21 subgraph: $GENE ==="
echo "K-mers: $KMERS"
echo "Started: $(date)"
echo ""

if [[ ! -f "$KMERS" || ! -s "$KMERS" ]]; then
    echo "ERROR: k-mer file not found: $KMERS"
    echo "Run build_k21_all_genes.sh first."
    exit 1
fi
echo "K-mers loaded: $(wc -l < "$KMERS")"
echo ""

# ---------------------------------------------------------------
# run_subgraph LABEL PUTG OUTDIR
# ---------------------------------------------------------------
run_subgraph() {
    local label="$1"
    local putg="$2"
    local outdir="$3"

    local segs="$outdir/HG002_${label}.${GENE}.k21.segments"
    local gfa="$outdir/HG002_${label}.${GENE}.k21.subset.gfa"
    local log="$outdir/HG002_${label}.${GENE}.k21.log"

    echo "--- $label ---"

    if [[ -f "$gfa" && -s "$gfa" ]]; then
        echo "  SKIP — $gfa exists"
        echo "  segments: $(grep -c '^S' "$gfa")  links: $(grep -c '^L' "$gfa" || true)  bp: $(awk '/^S/{sum+=length($3)} END{print sum}' "$gfa")"
        echo ""
        return
    fi

    python3 "$GFA_SCRIPT" --kmers "$KMERS" "$putg" > "$segs" 2>"$log"
    cat "$log"
    echo "  segments passing threshold: $(wc -l < "$segs")"

    gfatools view -l "@$segs" -r 0 "$putg" > "$gfa" 2>>"$log"

    echo "  segments: $(grep -c '^S' "$gfa")  links: $(grep -c '^L' "$gfa" || true)  bp: $(awk '/^S/{sum+=length($3)} END{print sum}' "$gfa")"
    echo "  output:   $gfa"
    echo ""
}

run_subgraph "12x"     "$PUTG_12X" "$OUT_12X"
run_subgraph "highcov" "$PUTG_HC"  "$OUT_HC"

echo "=== $GENE done: $(date) ==="
