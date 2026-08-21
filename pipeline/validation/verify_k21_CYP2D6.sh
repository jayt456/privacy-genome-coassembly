#!/bin/bash
# verify_k21_CYP2D6.sh
# Verify HG002 k=21 CYP2D6 subgraphs by aligning extracted segment sequences
# to the HG002 Q100 v1.1 paternal + maternal CYP2D6 references.
#
# For each assembly (12x, highcov):
#   1. GFA → FASTA (awk)
#   2. minimap2 asm5 → Q100 pat   → PAF
#   3. minimap2 asm5 → Q100 mat   → PAF
#   4. Print alignment summary (query, qlen, ref, ref_start, ref_end, mapq, matches)
#
# Good result: each segment aligns cleanly to pat or mat; together they cover
# the full ~4.3 kb CYP2D6 gene body.
#
# Outputs (written alongside each subgraph GFA):
#   HG002_{12x,highcov}.CYP2D6.k21.fa
#   HG002_{12x,highcov}.CYP2D6.k21.vs_Q100_pat.paf
#   HG002_{12x,highcov}.CYP2D6.k21.vs_Q100_mat.paf
#
#SBATCH --job-name=verify_k21
#SBATCH --time=00:20:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --output=logs/verify_k21_CYP2D6_%j.out

set -euo pipefail
export PATH="$HOME/.conda/envs/biotools/bin:$PATH"
module load awscli
module load minimap2/2.29-GCCcore-13.3.0

SCRATCH="${COASM_SCRATCH:-/nfs/roberts/scratch/pi_hc878/jst72}"
Q100_PAT="$SCRATCH/hprc_results/highcov/HG002/hg002.cyp2d6.pat.fa"
Q100_MAT="$SCRATCH/hprc_results/highcov/HG002/hg002.cyp2d6.mat.fa"

declare -a LABELS=("12x" "highcov")
declare -a GFAS=(
    "$SCRATCH/hprc_results/hifiasm_12x/HG002/HG002_12x.CYP2D6.k21.subset.gfa"
    "$SCRATCH/hprc_results/highcov/HG002/HG002_highcov.CYP2D6.k21.subset.gfa"
)

echo "=== k=21 CYP2D6 subgraph verification vs HG002 Q100 ==="
echo "Q100 pat: $Q100_PAT"
echo "Q100 mat: $Q100_MAT"
echo "Started: $(date)"
echo ""

for i in 0 1; do
    label="${LABELS[$i]}"
    gfa="${GFAS[$i]}"
    outdir="$(dirname "$gfa")"
    prefix="$outdir/HG002_${label}.CYP2D6.k21"
    fa="${prefix}.fa"
    paf_pat="${prefix}.vs_Q100_pat.paf"
    paf_mat="${prefix}.vs_Q100_mat.paf"

    echo "=== $label ==="
    echo "Input GFA: $gfa"
    echo "  Segments: $(grep -c '^S' "$gfa")   Total bp: $(awk '/^S/{sum+=length($3)} END{print sum}' "$gfa")"
    echo ""

    # Step 1: GFA → FASTA
    echo "--- Extracting FASTA ---"
    awk '/^S/{print ">"$2; print $3}' "$gfa" > "$fa"
    grep -c '^>' "$fa" | xargs echo "  sequences:"
    echo ""

    # Step 2: align vs Q100 pat
    echo "--- Aligning vs Q100 pat ---"
    minimap2 -c -x asm5 --cs -t 2 "$Q100_PAT" "$fa" > "$paf_pat" 2>/dev/null
    if [[ -s "$paf_pat" ]]; then
        echo "  Alignments found:"
        awk 'BEGIN{OFS="\t"; print "  query","qlen","qstart","qend","strand","ref","rlen","rstart","rend","matches","mapq"}
             {printf "  %-30s %7d %7d %7d  %s  %-30s %5d %5d %5d  %7d  %3d\n",
              $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$12}' "$paf_pat"
    else
        echo "  No alignments."
    fi
    echo ""

    # Step 3: align vs Q100 mat
    echo "--- Aligning vs Q100 mat ---"
    minimap2 -c -x asm5 --cs -t 2 "$Q100_MAT" "$fa" > "$paf_mat" 2>/dev/null
    if [[ -s "$paf_mat" ]]; then
        echo "  Alignments found:"
        awk 'BEGIN{OFS="\t"; print "  query","qlen","qstart","qend","strand","ref","rlen","rstart","rend","matches","mapq"}
             {printf "  %-30s %7d %7d %7d  %s  %-30s %5d %5d %5d  %7d  %3d\n",
              $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$12}' "$paf_mat"
    else
        echo "  No alignments."
    fi
    echo ""
done

echo "Done: $(date)"
