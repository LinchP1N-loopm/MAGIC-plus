#!/bin/bash
# ============================================================
# 00_get_coords.sh
# Purpose : Look up gene coordinates from GTF by gene name or
#           Ensembl gene ID (ENSMUSG…), then write a TSV used
#           by subsequent steps.
# Usage   : bash 00_get_coords.sh [gene1] [gene2] ...
# Example : bash 00_get_coords.sh Penk Gfap Apoe
#           bash 00_get_coords.sh ENSMUSG00000031980 Penk
# Output  : ./gene_coords.tsv
# ============================================================

set -euo pipefail

# ── Paths ──────────────────────────────────────────────────
GTF=/histor/zhao/pangkun/database/ref_mm10_M23_release98/mm10.gtf
OUTFILE=./gene_coords.tsv

# ── Padding around each gene (bp) ──────────────────────────
# Small padding to avoid clipping the gene edges visually.
# Set to 0 for exact gene boundaries.
PADDING=1000

# ── Per-gene custom plot regions (optional, bp, no padding) ─
# Format: CUSTOM["GeneName"]="plot_start,plot_end"
# Override when a gene is very long and only a portion matters.
declare -A CUSTOM
# CUSTOM["Mbp"]="82532000,82586000"
# CUSTOM["Gapdh"]="125163000,125166000"

# ── Input genes ────────────────────────────────────────────
if [ $# -eq 0 ]; then
    echo "[ERROR] No genes supplied."
    echo "Usage : bash 00_get_coords.sh Penk Gfap Apoe"
    echo "        bash 00_get_coords.sh ENSMUSG00000031980"
    exit 1
fi
GENES=("$@")

# ── Helper: resolve Ensembl ID → gene_name from GTF ────────
resolve_id() {
    local query="$1"
    if [[ "$query" =~ ^ENSMUSG ]]; then
        local name
        name=$(grep -m1 "gene_id \"${query}\"" "$GTF" \
               | grep -oP 'gene_name "\K[^"]+' | head -1)
        echo "${name:-$query}"
    else
        echo "$query"
    fi
}

# ── Write header ────────────────────────────────────────────
echo -e "gene\tchr\tgene_start\tgene_end\tstrand\tplot_start\tplot_end" \
    > "$OUTFILE"

echo "GTF    : $GTF"
echo "Output : $OUTFILE"
echo "Padding: ${PADDING} bp"
echo ""
printf "  %-18s %-7s %-13s %-13s %-8s %-13s %-13s %s\n" \
    "gene" "chr" "gene_start" "gene_end" "strand" "plot_start" "plot_end" "note"
echo "  $(printf '%0.s─' {1..90})"

FOUND=0; NOTFOUND=0

for raw_gene in "${GENES[@]}"; do

    # Resolve Ensembl ID to gene_name if necessary
    gene=$(resolve_id "$raw_gene")

    # Fetch the "gene" feature line from GTF
    line=$(grep -w "gene_name \"${gene}\"" "$GTF" \
           | awk '$3=="gene"' | head -1)

    if [ -z "$line" ]; then
        echo "  [NOT FOUND] ${raw_gene}"
        NOTFOUND=$((NOTFOUND+1))
        continue
    fi

    chr=$(echo    "$line" | awk '{print $1}')
    gstart=$(echo "$line" | awk '{print $4}')
    gend=$(echo   "$line" | awk '{print $5}')
    strand=$(echo "$line" | awk '{print $7}')

    if [ -n "${CUSTOM[$gene]+_}" ]; then
        plot_start=$(echo "${CUSTOM[$gene]}" | cut -d, -f1)
        plot_end=$(echo   "${CUSTOM[$gene]}" | cut -d, -f2)
        note="[custom region]"
    else
        plot_start=$(( gstart - PADDING ))
        plot_end=$(( gend + PADDING ))
        [ "$plot_start" -lt 1 ] && plot_start=1
        note=""
    fi

    printf "  %-18s %-7s %-13s %-13s %-8s %-13s %-13s %s\n" \
        "$gene" "$chr" "$gstart" "$gend" "$strand" "$plot_start" "$plot_end" "$note"
    echo -e "${gene}\t${chr}\t${gstart}\t${gend}\t${strand}\t${plot_start}\t${plot_end}" \
        >> "$OUTFILE"
    FOUND=$((FOUND+1))
done

echo ""
echo "  Found: ${FOUND}   Not found: ${NOTFOUND}"
echo "  → gene_coords.tsv written"
echo ""
echo "Next: bash 01_extract_depth.sh"
