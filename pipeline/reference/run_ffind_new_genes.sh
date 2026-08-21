#!/bin/bash
# Run ffind.py for IGH, KIR, and HLA and update genes.txt in-place with results.
# HLA is ~7.7 Mb so search may need to extend far — allow 6 hours.
#
#SBATCH --job-name=ffind_new_genes
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --output=logs/ffind_new_genes_%j.out

set -euo pipefail

BASE="${COASM_BASE:-/nfs/roberts/project/pi_hc878/jst72}"
BIOTOOLS_ENV=/home/jst72/.conda/envs/biotools
PYTHON="$BIOTOOLS_ENV/bin/python"
REF="$BASE/reference/hg38_orig_reference.fasta"
DB="$BASE/reference/hg38_orig_reference.k15.jf"
SCRIPT="$BASE/scripts/ffind.py"
GENES_TXT="$BASE/genes.txt"
SUMMARY="$BASE/reference/ffind_flanks_summary.txt"

export PATH="$BIOTOOLS_ENV/bin:$PATH"

[[ ! -f "$DB" ]]     && { echo "ERROR: jellyfish DB not found: $DB"; exit 1; }
[[ ! -f "$REF" ]]    && { echo "ERROR: reference not found: $REF"; exit 1; }
[[ ! -f "$SCRIPT" ]] && { echo "ERROR: ffind.py not found: $SCRIPT"; exit 1; }

NEW_GENES=(IGH KIR HLA)

for gene in "${NEW_GENES[@]}"; do
    # Read coords from genes.txt
    line=$(grep -P "^${gene}\t" "$GENES_TXT" | head -1)
    if [[ -z "$line" ]]; then
        echo "WARNING: $gene not found in genes.txt — skipping"
        continue
    fi

    chrom=$(echo "$line"     | awk '{print $2}')
    gene_start=$(echo "$line" | awk '{print $3}')
    gene_end=$(echo "$line"   | awk '{print $4}')
    flank_start=$(echo "$line" | awk '{print $5}')
    flank_end=$(echo "$line"   | awk '{print $6}')

    if [[ "$flank_start" != "$gene_start" || "$flank_end" != "$gene_end" ]]; then
        echo "=== $gene: flank coords already set ($chrom:$flank_start-$flank_end) — skipping ==="
        continue
    fi

    echo ""
    echo "=== $gene ($chrom:$gene_start-$gene_end, len $(( gene_end - gene_start )) bp) ==="
    echo "Running ffind.py..."

    output=$("$PYTHON" "$SCRIPT" \
        "$REF" "$chrom" "$gene_start" "$gene_end" \
        --db "$DB" -k 15 -t 1 -m 10000) || true

    echo "$output"

    # Parse "=> minimum region : chrN:START-END"
    region=$(echo "$output" | grep "=> minimum region" | grep -oP '\w+:\d+-\d+' | head -1)
    if [[ -z "$region" ]]; then
        echo "ERROR: could not parse minimum region from ffind.py output for $gene"
        continue
    fi

    new_flank_start=$(echo "$region" | cut -d: -f2 | cut -d- -f1)
    new_flank_end=$(echo "$region"   | cut -d: -f2 | cut -d- -f2)

    echo ""
    echo "  => $gene flank: $chrom:$new_flank_start-$new_flank_end"

    # Update genes.txt in-place: replace cols 5 and 6 for this gene
    awk -v gene="$gene" \
        -v fs="$new_flank_start" \
        -v fe="$new_flank_end" \
        'BEGIN{OFS="\t"}
         /^#/ { print; next }
         $1 == gene { $5 = fs; $6 = fe }
         { print }' "$GENES_TXT" > "${GENES_TXT}.tmp" && mv "${GENES_TXT}.tmp" "$GENES_TXT"

    echo "  genes.txt updated for $gene"

    # Append to summary
    echo "$gene  $chrom  $new_flank_start  $new_flank_end" >> "$SUMMARY"
done

echo ""
echo "=== Updated genes.txt ==="
cat "$GENES_TXT"

echo ""
echo "=== Summary appended to $SUMMARY ==="
cat "$SUMMARY"
