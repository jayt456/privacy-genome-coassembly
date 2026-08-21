#!/bin/bash
# run_HG002_CYP2D6_12x_trace.sh
# Q100-anchored read-tracing subgraph extraction for HG002 12x low-coverage assembly.
#
# Mirrors run_HG002_CYP2D6_trace.sh (highcov) exactly, but uses 12x GFAs.
# The 12x hap1/hap2 p_ctg FASTAs already exist; this script maps them to the
# Q100 CYP2D6 reference, filters to MAPQ60, then runs trace_reads.py to collect
# reads from the ctg A-lines and lifts them to unitigs in the 12x p_utg.gfa.
#
# Step 1: minimap2 asm5 hap1/hap2 p_ctg.fa → Q100 CYP2D6 → MAPQ60 sorted BAM
# Step 2: trace_reads.py (per hap) → UTG name list
# Step 3: merge UTG lists → gfatools subgraph
#
# Output: HG002_12x.CYP2D6.trace.subset.gfa
#
#SBATCH --job-name=HG002_12x_trace
#SBATCH --time=01:30:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --output=logs/HG002_12x_trace_%j.out

set -euo pipefail
export PATH="$HOME/.conda/envs/biotools/bin:$PATH"
module load awscli
module load minimap2/2.29-GCCcore-13.3.0
module load SAMtools/1.21-GCC-13.3.0

HIGHCOV=/nfs/roberts/scratch/pi_hc878/jst72/hprc_results/highcov/HG002
ASM12X=/nfs/roberts/scratch/pi_hc878/jst72/hprc_results/hifiasm_12x/HG002
SCRIPTS=/nfs/roberts/project/pi_hc878/jst72/scripts

UTG_GFA="$ASM12X/HG002_12x.bp.p_utg.gfa"
FLANK=150000   # ±150 kb around the 4.3 kb CYP2D6 gene body

OUT_PREFIX="$ASM12X/HG002_12x.CYP2D6.trace"
UTG_LIST="${OUT_PREFIX}.utgs"
GFA_OUT="${OUT_PREFIX}.subset.gfa"
TMP="$ASM12X/sort_tmp_trace"

echo "=== HG002 12x CYP2D6 read-tracing subgraph ==="
echo "Started: $(date)"
echo ""

mkdir -p "$TMP"

if [[ -f "$GFA_OUT" && -s "$GFA_OUT" ]]; then
    echo "SKIP — output already exists: $GFA_OUT"
    exit 0
fi

# Q100 references (same files used in highcov pipeline)
declare -A Q100_REF
Q100_REF[hap1]="$HIGHCOV/hg002.cyp2d6.pat.fa"
Q100_REF[hap2]="$HIGHCOV/hg002.cyp2d6.mat.fa"

# --- Step 1: align each hap's p_ctg FA to its Q100 reference ---
for hap in hap1 hap2; do
    fa="$ASM12X/HG002_12x.bp.${hap}.p_ctg.fa"
    bam="${OUT_PREFIX}.${hap}.MAPQ60.bam"
    ref="${Q100_REF[$hap]}"

    if [[ -f "$bam" && -s "$bam" ]]; then
        echo "SKIP alignment $hap — $bam already exists"
    else
        echo "--- Aligning $hap to Q100 CYP2D6 ($(basename "$ref")) ---"
        minimap2 -a -x asm5 -t 4 "$ref" "$fa" \
            | samtools view -F 4 -q 60 -b \
            | samtools sort -T "$TMP/$hap" -o "$bam"
        samtools index "$bam"
        echo "  flagstat:"
        samtools flagstat "$bam"
        echo ""
    fi
done

# --- Step 2: trace_reads.py (per hap) → UTG name lists ---
for hap in hap1 hap2; do
    bam="${OUT_PREFIX}.${hap}.MAPQ60.bam"
    ctg_gfa="$ASM12X/HG002_12x.bp.${hap}.p_ctg.gfa"
    utg_list_hap="${OUT_PREFIX}.${hap}.utgs"

    if [[ -f "$utg_list_hap" && -s "$utg_list_hap" ]]; then
        echo "SKIP trace $hap — $utg_list_hap already exists"
        echo "  utgs: $(wc -l < "$utg_list_hap")"
    else
        echo "--- Tracing $hap reads → UTGs ---"
        python3 "$SCRIPTS/trace_reads.py" \
            --bam     "$bam" \
            --ctg-gfa "$ctg_gfa" \
            --utg-gfa "$UTG_GFA" \
            --flank   "$FLANK" \
            > "$utg_list_hap"
        echo "  $hap utgs found: $(wc -l < "$utg_list_hap")"
    fi
    echo ""
done

# --- Step 3: merge UTG lists and extract subgraph ---
echo "--- Merging UTG lists ---"
sort -u \
    "${OUT_PREFIX}.hap1.utgs" \
    "${OUT_PREFIX}.hap2.utgs" \
    > "$UTG_LIST"
echo "  combined UTGs (deduplicated): $(wc -l < "$UTG_LIST")"
echo ""

echo "--- Extracting subgraph with gfatools ---"
gfatools view -l "@$UTG_LIST" -r 0 "$UTG_GFA" > "$GFA_OUT"

echo ""
echo "=== Subgraph summary ==="
echo "  Segments: $(grep -c '^S' "$GFA_OUT")"
echo "  Links:    $(grep -c '^L' "$GFA_OUT" || echo 0)"
echo "  Total bp: $(awk '/^S/{sum+=length($3)} END{print sum}' "$GFA_OUT")"
echo ""
echo "Done: $(date)"
echo "Output: $GFA_OUT"
ls -lh "${OUT_PREFIX}"*
