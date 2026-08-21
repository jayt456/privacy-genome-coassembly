#!/usr/bin/env bash
# Usage: ./extract_region.sh <gene_name> <chr:start-end>
# Example: ./extract_region.sh CYP2A6 chr19:40842850-40851193
set -uo pipefail

NAME="${1:?usage: extract_region.sh <gene_name> <chr:start-end>}"
REGION="${2:?usage: extract_region.sh <gene_name> <chr:start-end>}"
SRC="${COASM_BASE:-/nfs/roberts/project/pi_hc878/jst72}/reference/hg38_orig_reference.fasta
DEST="${COASM_BASE:-/nfs/roberts/project/pi_hc878/jst72}/reference/hg38${NAME}.fa

[[ ! -f "$SRC" ]] && { echo "ERROR: source reference not found: $SRC"; exit 1; }

if [[ -s "$DEST" ]]; then
    echo "$DEST already exists, skipping"
    exit 0
fi

echo "Extracting $REGION from GRCh38 -> $DEST (contig will be named >$NAME)"
samtools faidx "$SRC" "$REGION" | sed "s/^>.*/>${NAME}/" > "$DEST"
samtools faidx "$DEST"
echo "Done."
