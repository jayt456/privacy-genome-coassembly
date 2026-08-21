#!/bin/bash
# Verify all 11 pedigree samples (22 BAMs) are aligned to a gene reference.
# Usage: bash verify_pedigree_alignments.sh <GENE>   (default: CYP2D6)
# Writes pedigree_results/hg38<GENE>/alignment_verification.txt

set -euo pipefail

GENE="${1:-CYP2D6}"
REF_NAME="hg38${GENE}"

module load SAMtools/1.21-GCC-13.3.0

PROJ=/nfs/roberts/project/pi_hc878/jst72
BASE="$PROJ/pedigree_results/$REF_NAME"
OUT="$BASE/alignment_verification.txt"

if [[ ! -d "$BASE" ]]; then
  echo "ERROR: results directory not found: $BASE"
  exit 1
fi

samples=(NA12877 NA12878 NA12879 NA12881 NA12882 NA12885 NA12886 NA12889 NA12890 NA12891 NA12892)
total=$(( ${#samples[@]} * 2 ))

{
  echo "Pedigree ${GENE} Alignment Verification"
  echo "Reference: ${REF_NAME}"
  echo "Date: $(date -u '+%Y-%m-%d')"
  echo "Aligner: minimap2 -ax asm5"
  echo ""
  printf "%-12s  %-10s  %-10s  %-12s  %-10s  %s\n" \
    "SAMPLE" "HAP" "STATUS" "SIZE" "MAPPED" "PRIMARY_MAPPED"
  echo "$(printf '%.0s-' {1..80})"

  ok=0
  missing=0

  for s in "${samples[@]}"; do
    for hap in hap1 hap2; do
      bam="$BASE/$s/${s}_${hap}_to_${REF_NAME}.bam"
      bai="${bam}.bai"
      if [[ -f "$bam" && -f "$bai" ]]; then
        size=$(du -sh "$bam" | cut -f1)
        mapped=$(samtools flagstat "$bam" | awk '/mapped \(/ {print $1; exit}')
        primary=$(samtools flagstat "$bam" | awk '/primary mapped/ {print $1; exit}')
        printf "%-12s  %-10s  %-10s  %-12s  %-10s  %s\n" \
          "$s" "$hap" "OK" "$size" "$mapped" "$primary"
        ok=$((ok+1))
      else
        printf "%-12s  %-10s  %-10s\n" "$s" "$hap" "MISSING"
        missing=$((missing+1))
      fi
    done
  done

  echo ""
  echo "$(printf '%.0s-' {1..80})"
  echo "TOTAL OK:      $ok / $total"
  echo "TOTAL MISSING: $missing / $total"
} | tee "$OUT"

echo ""
echo "Written to $OUT"
