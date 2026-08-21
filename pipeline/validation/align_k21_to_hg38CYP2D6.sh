#!/bin/bash
# align_k21_to_hg38CYP2D6.sh
# Align HG002 k=21 CYP2D6 subgraph segments (12x + highcov) to the hg38CYP2D6
# reference slice, producing sorted + indexed BAMs ready for IGV.
#
# Reference: reference/hg38CYP2D6.fa (chr22:41,983,509-42,269,959)
# Inputs:    k=21 subgraph GFAs (S-lines extracted to FASTA, then minimap2 asm5)
#
# Outputs in $SCRATCH/hprc_results/k21_igv/:
#   HG002_12x.CYP2D6.k21.vs_hg38CYP2D6.bam  + .bai
#   HG002_highcov.CYP2D6.k21.vs_hg38CYP2D6.bam  + .bai
#   HG002_12x.CYP2D6.k21.fa                   (intermediate FASTA)
#   HG002_highcov.CYP2D6.k21.fa               (intermediate FASTA)
#
# IGV: load hg38CYP2D6.fa as genome, then open BAMs from k21_igv/
#
#SBATCH --job-name=k21_igv
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --output=logs/k21_igv_%j.out

set -euo pipefail
export PATH="$HOME/.conda/envs/biotools/bin:$PATH"
module load awscli
module load SAMtools/1.21-GCC-13.3.0
module load minimap2/2.29-GCCcore-13.3.0

BASE="${COASM_BASE:-/nfs/roberts/project/pi_hc878/jst72}"
SCRATCH="${COASM_SCRATCH:-/nfs/roberts/scratch/pi_hc878/jst72}"
REF="$BASE/reference/hg38CYP2D6.fa"
OUTDIR="$SCRATCH/hprc_results/k21_igv"
mkdir -p "$OUTDIR"

declare -a LABELS=("12x" "highcov")
declare -a GFAS=(
    "$SCRATCH/hprc_results/hifiasm_12x/HG002/HG002_12x.CYP2D6.k21.subset.gfa"
    "$SCRATCH/hprc_results/highcov/HG002/HG002_highcov.CYP2D6.k21.subset.gfa"
)

echo "=== k=21 CYP2D6 subgraphs → hg38CYP2D6 BAMs ==="
echo "Reference: $REF"
echo "Output:    $OUTDIR"
echo "Started:   $(date)"
echo ""

for i in 0 1; do
    label="${LABELS[$i]}"
    gfa="${GFAS[$i]}"
    fa="$OUTDIR/HG002_${label}.CYP2D6.k21.fa"
    bam="$OUTDIR/HG002_${label}.CYP2D6.k21.vs_hg38CYP2D6.bam"

    echo "=== $label ==="
    echo "  GFA:  $gfa"
    echo "  BAM:  $bam"

    if [[ -f "$bam" && -s "$bam" ]]; then
        echo "  SKIP — BAM already exists"
        echo ""
        continue
    fi

    echo "  Segments in GFA: $(grep -c '^S' "$gfa")  Total bp: $(awk '/^S/{sum+=length($3)} END{print sum}' "$gfa")"

    # Step 1: GFA → FASTA
    echo "  Extracting FASTA..."
    awk '/^S/{print ">"$2; print $3}' "$gfa" > "$fa"

    # Step 2: align to hg38CYP2D6 reference slice
    echo "  Aligning with minimap2 asm5..."
    minimap2 -a -x asm5 --cs -t 4 "$REF" "$fa" \
        | samtools sort -@ 4 -o "$bam"
    samtools index "$bam"

    echo "  Done. flagstat:"
    samtools flagstat "$bam" | sed 's/^/    /'
    echo ""
done

echo "Done: $(date)"
echo ""
echo "To view in IGV:"
echo "  1. Genome → Load Genome from File → reference/hg38CYP2D6.fa"
echo "  2. File → Load from File → k21_igv/HG002_12x.CYP2D6.k21.vs_hg38CYP2D6.bam"
echo "  3. File → Load from File → k21_igv/HG002_highcov.CYP2D6.k21.vs_hg38CYP2D6.bam"
