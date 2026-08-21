#!/bin/bash
# run_gfa_subset.sh: extract gene-region subgraph from each sample's unitig GFA.
#
# Uses k=21 genome-unique k-mers (confirmed best method; eliminates paralog contamination
# that affected k=15). Scores each p_utg segment by ref_pct = hits/len(ref_kmers)*100,
# threshold >1%. Extracts passing segments + their edges with gfatools view -r 0.
#
# Pipeline per sample × gene:
#   1. python gfa_k21.py --kmers reference/hg38<GENE>.unique_k21.txt <utg.gfa>
#      → segment IDs with ref_pct > 1%
#   2. gfatools view -l @segments -r 0 <utg.gfa> → subset.gfa
#
# Requires: reference/hg38<GENE>.unique_k21.txt for each gene (all 8 complete)
#   → regenerate with: sbatch scripts/build_k21_all_genes.sh
#
# Output: {pedigree,hprc}_results/hifiasm_12x/<sample>/<sample>_12x.<GENE>.k21.{segments,subset.gfa}
# HPRC assemblies (46 non-HG002 samples): /nfs/roberts/scratch/pi_hc878/jst72/hprc_results/hifiasm_12x/
# HG002 assembly: /nfs/roberts/project/pi_hc878/jst72/hprc_results/hifiasm_12x/HG002/
#
# Usage:
#   sbatch run_gfa_subset.sh              # both datasets
#   sbatch run_gfa_subset.sh pedigree     # pedigree only
#   sbatch run_gfa_subset.sh hprc         # hprc only
#
# Skip-if-present: won't overwrite existing .k21.subset.gfa files.
#
#SBATCH --job-name=gfa_subset
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --output=logs/gfa_subset_%j.out

set -euo pipefail

BASE="${COASM_BASE:-/nfs/roberts/project/pi_hc878/jst72}"
SCRATCH_BASE="${COASM_SCRATCH:-/nfs/roberts/scratch/pi_hc878/jst72}"
SCRIPT="${COASM_TOOLS:-$BASE/src}/gfa_k21.py"
GENES_FILE="$BASE/genes.txt"

# Read gene names from column 1 of genes.txt (skip comment lines)
mapfile -t GENES < <(grep -v '^#' "$GENES_FILE" | awk '{print $1}')

pedigree_samples=(
    NA12877 NA12878 NA12879 NA12881 NA12882
    NA12885 NA12886 NA12889 NA12890 NA12891 NA12892
)

hprc_samples=(
    HG002    HG005    HG00438  HG00621  HG00673  HG00733  HG00735  HG00741
    HG01071  HG01106  HG01109  HG01123  HG01175  HG01243  HG01258  HG01358
    HG01361  HG01891  HG01928  HG01952  HG01978  HG02055  HG02080  HG02109
    HG02145  HG02148  HG02257  HG02486  HG02559  HG02572  HG02622  HG02630
    HG02717  HG02723  HG02818  HG02886  HG03098  HG03453  HG03486  HG03492
    HG03516  HG03540  HG03579  NA18906  NA19240  NA20129  NA21309
)

run_dataset () {
    local dataset="$1"   # "pedigree" or "hprc"
    local -n samples_ref=$2
    local asm_base="$3"  # directory containing hifiasm_12x/<sample>/

    echo "=== Dataset: $dataset (${#samples_ref[@]} samples × ${#GENES[@]} genes) ==="

    local done=0 skipped=0 missing=0

    for sample in "${samples_ref[@]}"; do
        local outdir="$asm_base/hifiasm_12x/$sample"
        local utg="$outdir/${sample}_12x.bp.p_utg.gfa"

        if [[ ! -f "$utg" ]]; then
            echo "  SKIP $sample — unitig GFA not found: $utg"
            (( missing++ )) || true
            continue
        fi

        for gene in "${GENES[@]}"; do
            local kmers_file="$BASE/reference/hg38${gene}.unique_k21.txt"
            local seg_out="$outdir/${sample}_12x.${gene}.k21.segments"
            local gfa_out="$outdir/${sample}_12x.${gene}.k21.subset.gfa"
            local log="$outdir/${sample}_12x.gfa_subset.log"

            if [[ ! -f "$kmers_file" ]]; then
                echo "  WARN: unique k-mers not found: $kmers_file  (run build_unique_kmers.sh first)"
                continue
            fi

            if [[ -f "$gfa_out" && -s "$gfa_out" ]]; then
                (( skipped++ )) || true
                continue
            fi

            python3 "$SCRIPT" --kmers "$kmers_file" "$utg" \
                > "$seg_out" 2>>"$log"

            if [[ ! -s "$seg_out" ]]; then
                echo "  WARN: no segments found — $sample × $gene"
                continue
            fi

            # -r 0: no hop expansion; k=21 genome-unique k-mers already identify the right segments
            gfatools view -l "@$seg_out" -r 0 "$utg" > "$gfa_out" \
                2>>"$log"

            (( done++ )) || true
            segs=$(grep -c '^S' "$gfa_out" 2>/dev/null || echo 0)
            bp=$(awk '/^S/{sum+=length($3)} END{print sum+0}' "$gfa_out")
            echo "  OK  $sample × $gene  ($segs segs, ${bp} bp)"
        done
    done

    echo "  $dataset: $done extracted, $skipped skipped (existed), $missing samples missing GFA"
}

mode="${1:-both}"

case "$mode" in
    pedigree)
        run_dataset pedigree pedigree_samples "$BASE/pedigree_results" ;;
    hprc)
        run_dataset hprc hprc_samples "$SCRATCH_BASE/hprc_results" ;;
    both)
        run_dataset pedigree pedigree_samples "$BASE/pedigree_results"
        run_dataset hprc hprc_samples "$SCRATCH_BASE/hprc_results" ;;
    *)
        echo "Usage: $0 [pedigree|hprc|both]"; exit 1 ;;
esac

echo ""
echo "Done."
