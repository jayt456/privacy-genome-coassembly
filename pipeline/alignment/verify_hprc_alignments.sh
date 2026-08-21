#!/bin/bash
# Verify all 47 HPRC samples (94 BAMs) are aligned to a gene reference.
# Usage: bash verify_hprc_alignments.sh <GENE>   (default: CYP2D6)
# Writes hprc_results/hg38<GENE>/alignment_verification.txt

set -euo pipefail

GENE="${1:-CYP2D6}"
REF_NAME="hg38${GENE}"

module load SAMtools/1.21-GCC-13.3.0

PROJ=/nfs/roberts/project/pi_hc878/jst72
BASE="$PROJ/hprc_results/$REF_NAME"
OUT="$BASE/alignment_verification.txt"

if [[ ! -d "$BASE" ]]; then
  echo "ERROR: results directory not found: $BASE"
  exit 1
fi

samples=(
  HG002 HG005 HG00438 HG00621 HG00673 HG00733 HG00735 HG00741
  HG01071 HG01106 HG01109 HG01123 HG01175 HG01243 HG01258 HG01358
  HG01361 HG01891 HG01928 HG01952 HG01978 HG02055 HG02080 HG02109
  HG02145 HG02148 HG02257 HG02486 HG02559 HG02572 HG02622 HG02630
  HG02717 HG02723 HG02818 HG02886 HG03098 HG03453 HG03486 HG03492
  HG03516 HG03540 HG03579 NA18906 NA19240 NA20129 NA21309
)
total=$(( ${#samples[@]} * 2 ))

{
  echo "HPRC ${GENE} Alignment Verification"
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
    for hap in paternal maternal; do
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
