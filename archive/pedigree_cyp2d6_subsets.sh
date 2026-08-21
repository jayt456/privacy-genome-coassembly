#!/bin/bash
# pedigree_cyp2d6_subsets.sh
# Extract CYP2D6 subgraph from each pedigree sample, then combine into one GFA per coverage level.
#
# Per sample:
#   12x:     gfa.py on p_utg.gfa  → segments file → gfatools subset GFA
#   highcov: gfa.py on hap1+hap2 FASTA → minimal GFA (S-lines only)
#
# Combined (ID-prefixed to avoid collisions):
#   pedigree_results/CYP2D6.pedigree.12x.combined.gfa
#   pedigree_results/CYP2D6.pedigree.highcov.combined.gfa
#
# Usage:
#   sbatch pedigree_cyp2d6_subsets.sh          # all 11 samples
#   sbatch pedigree_cyp2d6_subsets.sh NA12878  # single sample (no combine step)
#
#SBATCH --job-name=cyp2d6_subsets
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --output=logs/cyp2d6_subsets_%j.out

set -euo pipefail

export PATH="/home/jst72/.conda/envs/biotools/bin:$PATH"

BASE="${COASM_BASE:-/nfs/roberts/project/pi_hc878/jst72}"
GENE=CYP2D6
REF_FA="$BASE/reference/hg38${GENE}.fa"
SCRIPT="$BASE/scripts/gfa.py"

all_samples=(NA12877 NA12878 NA12879 NA12881 NA12882
             NA12885 NA12886 NA12889 NA12890 NA12891 NA12892)

if [[ $# -eq 1 ]]; then
    samples=("$1")
    single_mode=1
else
    samples=("${all_samples[@]}")
    single_mode=0
fi

# ---- prefix_gfa: rewrite segment IDs in a GFA with a sample-name prefix ----
# Handles S-lines and L-lines. Output goes to stdout.
prefix_gfa () {
    local sample="$1"
    local gfa="$2"
    awk -v s="$sample" '
        /^H/ { print; next }
        /^S/ { $2 = s "_" $2; print; next }
        /^L/ { $2 = s "_" $2; $4 = s "_" $4; print; next }
        { print }
    ' OFS='\t' "$gfa"
}

echo "=== 12x low-coverage subsets ==="
for s in "${samples[@]}"; do
    utg="$BASE/pedigree_results/hifiasm_12x/$s/${s}_12x.bp.p_utg.gfa"
    seg="$BASE/pedigree_results/hifiasm_12x/$s/${s}_12x.${GENE}.segments"
    out="$BASE/pedigree_results/hifiasm_12x/$s/${s}_12x.${GENE}.subset.gfa"

    if [[ ! -f "$utg" ]]; then
        echo "  SKIP $s — unitig GFA not found"; continue
    fi
    if [[ -f "$out" && -s "$out" ]]; then
        echo "  SKIP $s — 12x subset GFA already exists"
        grep -c '^S' "$out" | xargs -I{} echo "         ({} segments)"
        continue
    fi

    printf "  %-12s " "$s"
    if [[ -f "$seg" && -s "$seg" ]]; then
        n_seg=$(wc -l < "$seg")
        printf "scoring skipped (%d IDs cached) → " "$n_seg"
    else
        printf "scoring ... "
        python3 "$SCRIPT" "$REF_FA" "$utg" > "$seg" 2>/dev/null
        n_seg=$(wc -l < "$seg")
        printf "%d matching IDs → " "$n_seg"
    fi
    gfatools view -l "@$seg" -r 1 "$utg" > "$out"
    n_s=$(grep -c '^S' "$out" || true)
    echo "${n_s} segments"
done

echo ""
echo "=== High-coverage subsets (FASTA → minimal GFA) ==="
for s in "${samples[@]}"; do
    out="$BASE/pedigree_assemblies/$s/${s}_highcov.${GENE}.subset.gfa"

    if [[ -f "$out" && -s "$out" ]]; then
        echo "  SKIP $s — highcov subset GFA already exists"
        continue
    fi

    hap1=$(ls "$BASE/pedigree_assemblies/$s/"*.hap1.p_ctg.gfa.fasta 2>/dev/null | head -1)
    hap2=$(ls "$BASE/pedigree_assemblies/$s/"*.hap2.p_ctg.gfa.fasta 2>/dev/null | head -1)

    if [[ -z "$hap1" || -z "$hap2" ]]; then
        echo "  SKIP $s — hap FASTAs not found"; continue
    fi

    printf "  %-12s scoring ... " "$s"
    {
        python3 "$SCRIPT" "$REF_FA" "$hap1" 2>/dev/null
        python3 "$SCRIPT" "$REF_FA" "$hap2" 2>/dev/null | grep -v '^H'
    } > "$out"
    n_s=$(grep -c '^S' "$out" || true)
    echo "${n_s} contigs"
done

# ---- combine step: only in full mode ----
if [[ $single_mode -eq 0 ]]; then
    echo ""
    echo "=== Combining per-sample GFAs ==="

    combined_12x="$BASE/pedigree_results/CYP2D6.pedigree.12x.combined.gfa"
    combined_hc="$BASE/pedigree_results/CYP2D6.pedigree.highcov.combined.gfa"

    echo "H\tVN:Z:1.0" > "$combined_12x"
    echo "H\tVN:Z:1.0" > "$combined_hc"

    for s in "${all_samples[@]}"; do
        gfa12x="$BASE/pedigree_results/hifiasm_12x/$s/${s}_12x.${GENE}.subset.gfa"
        gfahc="$BASE/pedigree_assemblies/$s/${s}_highcov.${GENE}.subset.gfa"

        if [[ -f "$gfa12x" && -s "$gfa12x" ]]; then
            prefix_gfa "$s" "$gfa12x" | grep -v '^H' >> "$combined_12x"
            echo "  $s 12x: $(grep -c '^S' "$gfa12x") segments added (prefixed)"
        else
            echo "  SKIP $s 12x — subset GFA missing"
        fi

        if [[ -f "$gfahc" && -s "$gfahc" ]]; then
            prefix_gfa "$s" "$gfahc" | grep -v '^H' >> "$combined_hc"
            echo "  $s highcov: $(grep -c '^S' "$gfahc") segments added (prefixed)"
        else
            echo "  SKIP $s highcov — subset GFA missing"
        fi
    done

    echo ""
    echo "Combined 12x:     $combined_12x"
    echo "  $(grep -c '^S' "$combined_12x") total segments"
    echo "Combined highcov: $combined_hc"
    echo "  $(grep -c '^S' "$combined_hc") total segments"
fi

echo ""
echo "Done."
