#!/bin/bash
#PBS -N 3tech_bam2bigwig
#PBS -j oe
#PBS -d /path/to/project/20260427_AS_analysis/Fig2_trackplot/qsub
#PBS -l nodes=1:ppn=10
#PBS -l walltime=48:00:00
#PBS -q gold

eval "$(micromamba shell hook --shell bash)"
micromamba activate biotools

set -euo pipefail

# ── 路径 ───────────────────────────────────────────────────────────────────────
OUTDIR="/path/to/project/20260427_AS_analysis/Fig2_trackplot/3tech_bigwig_data"
THREADS=10

mkdir -p "${OUTDIR}"

# ── 样本配置 ───────────────────────────────────────────────────────────────────
declare -A BAM
BAM["G21-4-0504"]="/path/to/project/20260427_AS_analysis/Fig2_rmats/ASclu_results/pseudobulk_bams/G21-4-0504/G21-4-0504_ALL_pseudobulk.bam"
BAM["WT1_10"]="/path/to/project/20260427_AS_analysis/Fig2_rmats/ASclu_results/pseudobulk_bams/WT1_10/WT1_10_ALL_pseudobulk.bam"
BAM["sample_Reg1_10"]="/path/to/project/20260427_AS_analysis/Fig2_rmats/ASclu_results/pseudobulk_bams/sample_Reg1_10/sample_Reg1_10_ALL_pseudobulk.bam"

SAMPLE_ORDER=("G21-4-0504" "WT1_10" "sample_Reg1_10")

echo "============================================================"
echo " 3tech bam2bigwig  (RPKM, binSize=1)"
echo " OUTDIR: ${OUTDIR}"
echo "============================================================"

for key in "${SAMPLE_ORDER[@]}"; do
    bam="${BAM[$key]}"
    bw="${OUTDIR}/${key}_RPKM.bw"

    echo ""
    echo ">>> ${key}"

    if [[ ! -f "${bam}" ]]; then
        echo "  [ERROR] BAM not found: ${bam}"; continue
    fi

    if [[ ! -f "${bam}.bai" ]]; then
        echo "  [INFO] indexing BAM ..."
        samtools index -@ "${THREADS}" "${bam}"
    fi

    if [[ -f "${bw}" ]]; then
        echo "  [SKIP] already exists: ${bw}"; continue
    fi

    echo "  BAM : ${bam}"
    echo "  OUT : ${bw}"

    bamCoverage \
        --bam                "${bam}" \
        --outFileName        "${bw}" \
        --outFileFormat      bigwig \
        --binSize            1 \
        --normalizeUsing     RPKM \
        --exactScaling \
        --minMappingQuality  10 \
        --numberOfProcessors "${THREADS}"

    [[ -f "${bw}" ]] && echo "  [OK] ${bw}" || echo "  [ERROR] output not created"
done

echo ""
echo "============================================================"
echo " Done → ${OUTDIR}"
echo "============================================================"