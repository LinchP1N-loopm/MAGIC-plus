#!/bin/bash
# ============================================================
# 01_extract_depthV2.sh
# Purpose : Convert BAM → BigWig using deepTools bamCoverage
#           (RPKM normalisation, binSize 1, minMappingQuality 10)
#           Replaces samtools depth approach in v1.
#
# Pre-req : bash 00_get_coords.sh [genes…]  (writes gene_coords.tsv)
#           deepTools >= 3.0 in PATH
# Usage   : bash 01_extract_depthV2.sh
# Outputs : ./bigwig_data/
#               <display_name>.bw    (one per sample, whole-genome)
# ============================================================

set -euo pipefail

# ── Paths ──────────────────────────────────────────────────
COORDS=./gene_coords.tsv
OUTDIR=./bigwig_data
THREADS=8          # adjust to available CPUs

# ── Sample registry ────────────────────────────────────────
SAMPLE_ORDER=("MAGIC_plus" "MAGIC_seq" "X10_Visium")

declare -A DISPLAY_NAME
DISPLAY_NAME["MAGIC_plus"]="MAGIC_plus"
DISPLAY_NAME["MAGIC_seq"]="MAGIC-seq"
DISPLAY_NAME["X10_Visium"]="10X_Visium"

declare -A BAM
BAM["MAGIC_plus"]="/histor/zhao/pangkun/Zhaolab/Lupeng/Zhujunjie/20260305_rejiace/data/result_STARsolo/G21-4/STARsolo_pe/G21-4_Aligned_AB.sortedByCoord.out.bam"
BAM["MAGIC_seq"]="/histor/zhao/pangkun/Zhaolab/Lupeng/Zhujunjie/20260307_testdata/data/result_STARsolo/ZZJ_Reg1_10/STARsolo_se/ZZJ_Reg1_10_Aligned_AB.sortedByCoord.out.bam"
BAM["X10_Visium"]="/histor/zhao/pangkun/Zhaolab/Lupeng/Zhujunjie/20260307_testdata/data/result_STARsolo/WT1_10/STARsolo_se/WT1_10_Aligned_AB.sortedByCoord.out.bam"

# ── Sanity checks ──────────────────────────────────────────
[ ! -f "$COORDS" ] && \
    echo "[ERROR] gene_coords.tsv not found. Run 00_get_coords.sh first." && exit 1

command -v bamCoverage &>/dev/null || \
    { echo "[ERROR] bamCoverage not found. Activate deepTools env first."; exit 1; }

mkdir -p "$OUTDIR"

echo "============================================================"
echo " Step 0  Diagnose BAM files"
echo "============================================================"
for key in "${SAMPLE_ORDER[@]}"; do
    sname="${DISPLAY_NAME[$key]}"
    bam="${BAM[$key]}"
    echo ""
    echo ">>> $sname"
    if [ ! -f "$bam" ]; then
        echo "    [ERROR] BAM not found: $bam"; continue
    fi
    echo "    [OK] BAM exists"
    { [ -f "${bam}.bai" ] || [ -f "${bam%.bam}.bai" ]; } \
        && echo "    [OK] Index found" \
        || echo "    [WARNING] No .bai index — run: samtools index $bam"
done

echo ""
echo "============================================================"
echo " Step 1  BAM → BigWig  (bamCoverage, RPKM, binSize=1)"
echo "============================================================"
echo " Note: each BigWig covers the whole genome — one file per"
echo "       sample, reused for all genes. This may take ~10-30"
echo "       min per sample depending on BAM size."
echo ""

for key in "${SAMPLE_ORDER[@]}"; do
    sname="${DISPLAY_NAME[$key]}"
    bam="${BAM[$key]}"
    bw="$OUTDIR/${sname}.bw"

    [ ! -f "$bam" ] && echo "[skip] $sname: BAM not found" && continue

    # Skip if BigWig already exists
    if [ -f "$bw" ]; then
        echo "  [SKIP] $sname → $bw already exists"
        continue
    fi

    echo "  → $sname"
    echo "    BAM : $bam"
    echo "    OUT : $bw"
    echo -n "    Running bamCoverage ... "

    bamCoverage \
        --bam              "$bam" \
        --outFileName      "$bw" \
        --outFileFormat    bigwig \
        --binSize          1 \
        --normalizeUsing   RPKM \
        --exactScaling \
        --minMappingQuality 10 \
        --numberOfProcessors "$THREADS" \
        2>&1 | tail -3

    [ -f "$bw" ] && echo "    [OK] $bw" || echo "    [ERROR] output not created"
    echo ""
done

echo "============================================================"
echo " Done.  BigWig files in: $OUTDIR/"
echo " Next:  python 02_plot_tracksV2.py [gene1] [gene2] ..."
echo "============================================================"
