#!/bin/bash
# =============================================================
# step1_extract_barcodes_regions.sh
# =============================================================
# Purpose:
#   调用 step1_extract_barcodes_regions.py，从新的 *_regions.h5ad 中
#   按 obs["Regions"] 提取各脑区 barcode 列表。
#
# 用法：
#   bash step1_extract_barcodes_regions.sh <样本名> <h5ad路径>
#
# 示例：
#   bash step1_extract_barcodes_regions.sh G21-1 \
#     /path/to/project/20260705_Fig3/anno_h5ads/G21-1_regions.h5ad
#
# 批量运行示例：
#   ADATA_DIR="/path/to/project/20260705_Fig3/anno_h5ads"
#   for s in G21-1 G21-2 G21-4 G21-5 G21-7 G21-8; do
#       bash step1_extract_barcodes_regions.sh "$s" "${ADATA_DIR}/${s}_regions.h5ad"
#   done
# =============================================================

set -euo pipefail

# ============================================================
# Paths
# ============================================================
ADATA_DIR="/path/to/project/20260705_Fig3/anno_h5ads"
WORK_DIR="/path/to/project/20260705_Fig3/rmats_works"
LOG_DIR="${WORK_DIR}/logs"

mkdir -p "${LOG_DIR}"

# ============================================================
# Parse arguments
# ============================================================
if [[ $# -ne 2 ]]; then
    echo "用法: bash $0 <样本名> <h5ad路径>"
    echo ""
    echo "批量运行示例（六个样本）:"
    for s in G21-1 G21-2 G21-4 G21-5 G21-7 G21-8; do
        echo "  bash $0 ${s} ${ADATA_DIR}/${s}_regions.h5ad"
    done
    exit 1
fi

SAMPLE="$1"
H5AD_PATH="$2"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PY_SCRIPT="${SCRIPT_DIR}/step1_extract_barcodes_regions.py"
LOG_FILE="${LOG_DIR}/step1_${SAMPLE}.log"

# ============================================================
# Check inputs
# ============================================================
if [[ ! -f "${PY_SCRIPT}" ]]; then
    echo "[ERROR] 找不到 Python 脚本: ${PY_SCRIPT}"
    echo "  请确认 step1_extract_barcodes_regions.py 与本 sh 脚本放在同一目录。"
    exit 1
fi

if [[ ! -f "${H5AD_PATH}" ]]; then
    echo "[ERROR] h5ad 文件不存在: ${H5AD_PATH}"
    exit 1
fi

# ============================================================
# Run step 1
# ============================================================
echo "[$(date '+%F %T')] Step1 开始: ${SAMPLE}"
echo "  h5ad:  ${H5AD_PATH}"
echo "  脚本:  ${PY_SCRIPT}"
echo "  日志:  ${LOG_FILE}"
echo ""

eval "$(micromamba shell hook --shell bash)"
micromamba activate py310_st

python "${PY_SCRIPT}" "${SAMPLE}" "${H5AD_PATH}" 2>&1 | tee "${LOG_FILE}"

echo ""
echo "[$(date '+%F %T')] Step1 完成: ${SAMPLE}"
echo "  日志: ${LOG_FILE}"
echo "  输出: ${WORK_DIR}/barcodes/${SAMPLE}/"
