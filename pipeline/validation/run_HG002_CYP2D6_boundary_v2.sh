#!/bin/bash
# run_HG002_CYP2D6_boundary_v2.sh
#
# CYP2D6 subgraph extraction for HG002 highcov assembly.
#
# Refinement of boundary-crossing method: separates MAPQ60 BAM reads by
# orientation (FLAG 0 = forward, FLAG 16 = reverse) and applies the correct
# CIGAR offset formula for each before scanning p_ctg A-lines for boundary-
# crossing reads. Those reads are then traced into p_utg, and per-haplotype
# subgraphs are extracted with gfatools.
#
# CIGAR parsing (after: cut -f1,4,6 | sed 's/[SM]/\t/g'):
#   fields: $1=contig_id  $2=POS  $3=leading_S  $4=M  $5=trailing_S
#
#   Forward (FLAG 0):  gene region on contig = [$2+$3, $2+$3+$4]
#                      (POS + leading soft-clip = contig start of alignment)
#   Reverse (FLAG 16): contig is RC'd; trailing_S in CIGAR = 5' of original contig
#                      gene region on contig = [$5, $5+$4]
#
# Output directory: $SCRATCH/hprc_results/highcov/HG002/subgraph_boundary_v2/
#
#SBATCH --job-name=HG002_boundary_v2
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --output=logs/HG002_boundary_v2_%j.out

set -euo pipefail
export PATH="$HOME/.conda/envs/biotools/bin:$PATH"
module load awscli
module load SAMtools/1.21-GCC-13.3.0

SCRATCH="${COASM_SCRATCH:-/nfs/roberts/scratch/pi_hc878/jst72}"
HG002="$SCRATCH/hprc_results/highcov/HG002"
ASM="$HG002/assembly"
Q100="$HG002/Q100"
OUT="$HG002/subgraph_boundary_v2"
PUTG="$ASM/HG002.asm.bp.p_utg.gfa"

mkdir -p "$OUT"

echo "=== HG002 CYP2D6 subgraph — boundary v2 (forward/reverse CIGAR) ==="
echo "Started: $(date)"
echo ""

# ---------------------------------------------------------------
# process_hap HAP PCTG_GFA BAM PREFIX
#   HAP       - label (hap1 / hap2)
#   PCTG_GFA  - phased contig GFA for this haplotype
#   BAM       - MAPQ60 contig-to-Q100 BAM for this haplotype
#   PREFIX    - output file prefix (full path)
# ---------------------------------------------------------------
process_hap() {
    local hap="$1"
    local pctg_gfa="$2"
    local bam="$3"
    local prefix="$4"

    echo "--- $hap ---"

    # Forward reads (FLAG 0): exclude reverse + secondary + supplementary
    # Gene region on contig: [POS + leading_S, POS + leading_S + M]
    awk 'NR==FNR { id=$1; s=$2+$3; e=$2+$3+$4; next }
         $1=="A" && $2==id && (($3<=e && ($3+$7)>=e) || ($3<=s && ($3+$7)>=s))' \
      <(samtools view -F 2320 "$bam" | cut -f1,4,6 | sed 's/[SM]/\t/g') \
      "$pctg_gfa" \
      > "${prefix}.fwd.lines"

    echo "  forward  boundary A-lines: $(wc -l < "${prefix}.fwd.lines")"

    # Reverse reads (FLAG 16): exclude secondary + supplementary
    # Contig is RC'd — trailing_S in CIGAR = 5' end of original contig
    # Gene region on contig: [trailing_S, trailing_S + M]
    awk 'NR==FNR { id=$1; s=$5; e=$5+$4; next }
         $1=="A" && $2==id && (($3<=e && ($3+$7)>=e) || ($3<=s && ($3+$7)>=s))' \
      <(samtools view -f 16 -F 2304 "$bam" | cut -f1,4,6 | sed 's/[SM]/\t/g') \
      "$pctg_gfa" \
      > "${prefix}.rev.lines"

    echo "  reverse  boundary A-lines: $(wc -l < "${prefix}.rev.lines")"

    # Collect all boundary-crossing read names, trace into p_utg A-lines
    cat "${prefix}.fwd.lines" "${prefix}.rev.lines" \
      | cut -f5 | sort -u \
      > "${prefix}.read_names"

    grep -f "${prefix}.read_names" "$PUTG" > "${prefix}.utg.lines"
    echo "  p_utg A-lines matched:     $(wc -l < "${prefix}.utg.lines")"

    # Extract unique UTG names
    awk '{print $2}' "${prefix}.utg.lines" | sort -u > "${prefix}.utgs"
    echo "  UTGs found:"
    sed 's/^/    /' "${prefix}.utgs"

    # Extract subgraph
    gfatools view -l "@${prefix}.utgs" -r 0 "$PUTG" > "${prefix}.subset.gfa"
    echo "  segments: $(grep -c '^S' "${prefix}.subset.gfa")"
    echo "  links:    $(grep -c '^L' "${prefix}.subset.gfa" || true)"
    echo "  output:   ${prefix}.subset.gfa"
    echo ""
}

process_hap "hap1" \
    "$ASM/HG002.asm.bp.hap1.p_ctg.gfa" \
    "$Q100/hg002.hap1_to_pat.CYP2D6.MAPQ60.bam" \
    "$OUT/HG002_highcov.CYP2D6.hap1"

process_hap "hap2" \
    "$ASM/HG002.asm.bp.hap2.p_ctg.gfa" \
    "$Q100/hg002.hap2_to_mat.CYP2D6.MAPQ60.bam" \
    "$OUT/HG002_highcov.CYP2D6.hap2"

# ---------------------------------------------------------------
# Combined subgraph — merge UTGs from both haplotypes
# ---------------------------------------------------------------
echo "--- combined subgraph ---"

cat "$OUT/HG002_highcov.CYP2D6.hap1.utgs" \
    "$OUT/HG002_highcov.CYP2D6.hap2.utgs" \
  | sort -u \
  > "$OUT/HG002_highcov.CYP2D6.combined.utgs"

echo "  total unique UTGs: $(wc -l < "$OUT/HG002_highcov.CYP2D6.combined.utgs")"

gfatools view \
    -l "@$OUT/HG002_highcov.CYP2D6.combined.utgs" \
    -r 0 \
    "$PUTG" \
    > "$OUT/HG002_highcov.CYP2D6.boundary_v2.subset.gfa"

echo "  segments: $(grep -c '^S' "$OUT/HG002_highcov.CYP2D6.boundary_v2.subset.gfa")"
echo "  links:    $(grep -c '^L' "$OUT/HG002_highcov.CYP2D6.boundary_v2.subset.gfa" || true)"
echo "  output:   $OUT/HG002_highcov.CYP2D6.boundary_v2.subset.gfa"
echo ""
echo "Done: $(date)"
echo "Results in: $OUT"
