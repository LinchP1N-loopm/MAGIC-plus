#!/usr/bin/env bash
set -euo pipefail

PROJECT="/histor/zhao/projects/lupeng/zhaolab/ZJJ/Fig4_sashimi"
SCRIPTS="${PROJECT}/scripts"
OUTDIR="${PROJECT}/sashimi_junctions/screenshot_candidates_same_y_20260926"

mkdir -p "${OUTDIR}"
bash "${SCRIPTS}/00_get_coords_regions.sh" Ncam1 Bin1 U2af1 Nrxn1 Msra Rbfox3

plot_one() {
    local gene="$1" region="$2" samples="$3" comparison="$4"

    python "${SCRIPTS}/02_plot_sashimi_junctions.py" \
        --gene "${PROJECT}/gene_coords/${gene}.tsv" \
        --region "${region}" \
        --samples "${samples}" \
        --same-y \
        --outdir "${OUTDIR}/${comparison}"
}

# AD vs Young：G21-4 = Young，G21-2 = AD
plot_one Ncam1 FT    G21-4,G21-2 AD_vs_Young
plot_one Bin1  STR   G21-4,G21-2 AD_vs_Young
plot_one U2af1 TH    G21-4,G21-2 AD_vs_Young
plot_one Nrxn1 TH    G21-4,G21-2 AD_vs_Young
plot_one Msra  CTX_1 G21-4,G21-2 AD_vs_Young

# Aged vs Young：G21-4 = Young，G21-1 = Aged
plot_one Rbfox3 CA1 G21-4,G21-1 Aged_vs_Young

echo "[DONE] ${OUTDIR}"
