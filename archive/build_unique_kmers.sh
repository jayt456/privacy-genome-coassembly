#!/bin/bash
# build_unique_kmers.sh: generate hg38<GENE>.unique_k15.txt for all genes in genes.txt.
#
# For each gene, calls ffind.py --output-kmers to collect all k-mers in the full
# flank window [sf, ef] that appear at most once genome-wide (fwd + RC count <= 1)
# using the prebuilt Jellyfish k15 DB.
#
# Output: reference/hg38<GENE>.unique_k15.txt  (one k-mer string per line)
# These files are consumed by gfa.py --kmers for genome-unique subgraph filtering.
#
# Usage:
#   sbatch build_unique_kmers.sh
#   bash build_unique_kmers.sh          # runs locally (slow — jellyfish queries)
#
# Skip-if-present: skips genes that already have a .unique_k15.txt file.
#
#SBATCH --job-name=build_kmers
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --output=logs/build_unique_kmers_%j.out

set -euo pipefail

export PATH="$HOME/.conda/envs/biotools/bin:$PATH"

BASE="${COASM_BASE:-/nfs/roberts/project/pi_hc878/jst72}"
REF="$BASE/reference/hg38_orig_reference.fasta"
JF_DB="$BASE/reference/hg38_orig_reference.k15.jf"
FFIND="$BASE/scripts/ffind.py"
GENES_FILE="$BASE/genes.txt"

echo "Building genome-unique k15 files for all genes in $GENES_FILE"
echo "Jellyfish DB: $JF_DB"
echo ""

# genes.txt columns: GENE CHROM GENE_START GENE_END FLANK_START FLANK_END
grep -v '^#' "$GENES_FILE" | while read -r gene chrom gene_start gene_end flank_start flank_end _rest; do
    out="$BASE/reference/hg38${gene}.unique_k15.txt"

    if [[ -f "$out" && -s "$out" ]]; then
        echo "SKIP $gene — $out already exists ($(wc -l < "$out") k-mers)"
        continue
    fi

    echo "Processing $gene ($chrom:$gene_start-$gene_end) ..."
    python3 "$FFIND" "$REF" "$chrom" "$gene_start" "$gene_end" \
        -k 15 -t 1 -m 10000 \
        --db "$JF_DB" \
        --output-kmers "$out" || echo "  WARNING: ffind.py exited non-zero for $gene (incomplete flank?), continuing"
    if [[ -s "$out" ]]; then
        echo "  => $(wc -l < "$out") unique k-mers written to $out"
    else
        echo "  ERROR: no output written for $gene, skipping"
        rm -f "$out"
    fi
    echo ""
done

echo "Done."
