#!/bin/bash
# =============================================================
# step1_extract_barcodes_AS.sh （shell 包装版）
# =============================================================
# Usage:
#   bash step1_extract_barcodes_AS.sh <样本名> <h5ad路径>
#
# 示例：
#   bash step1_extract_barcodes_AS.sh G21-4 /path/to/G21-4_anno.h5ad
#   bash step1_extract_barcodes_AS.sh G21-7 /path/to/G21-7_anno.h5ad
#
# 脚本放置目录：ASclu_qsub/
# 结果输出目录：ASclu_results/barcodes/<SAMPLE>/
# =============================================================

set -euo pipefail

# ============================================================
# Parse arguments
# ============================================================
if [[ $# -ne 2 ]]; then
    echo "用法: bash $0 <样本名> <h5ad路径>"
    echo ""
    echo "示例:"
    echo "  bash $0 G21-4 /path/to/G21-4_anno.h5ad"
    echo "  bash $0 G21-7 /path/to/G21-7_anno.h5ad"
    exit 1
fi

SAMPLE="${1}"
H5AD_PATH="${2}"

# ============================================================
# 固定路径（Python 脚本与本脚本放在同一目录）
# ============================================================
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PY_SCRIPT="${SCRIPT_DIR}/step1_extract_barcodes_AS.py"

if [[ ! -f "${PY_SCRIPT}" ]]; then
    echo "[ERROR] 找不到 Python 脚本: ${PY_SCRIPT}"
    exit 1
fi

if [[ ! -f "${H5AD_PATH}" ]]; then
    echo "[ERROR] h5ad 文件不存在: ${H5AD_PATH}"
    exit 1
fi

# ============================================================
# 激活环境并运行
# ============================================================
echo "[$(date '+%F %T')] Step1 开始: ${SAMPLE}"

eval "$(micromamba shell hook --shell bash)"
micromamba activate py310_st

python "${PY_SCRIPT}" "${SAMPLE}" "${H5AD_PATH}" 2>&1 | tee "step1_${SAMPLE}.log"

echo "[$(date '+%F %T')] Step1 完成: ${SAMPLE}"
echo "  日志: ${SCRIPT_DIR}/step1_${SAMPLE}.log"
