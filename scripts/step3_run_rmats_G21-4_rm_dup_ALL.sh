#!/bin/bash
# =============================================================
# step3_run_rmats_G21-4_rm_dup_ALL.sh
# =============================================================
# 对 G21-4_rm_dup 的 ALL BAM 跑 --statoff 模式
# =============================================================

set -euo pipefail

eval "$(micromamba shell hook --shell bash)"
micromamba activate rmats_env

echo "[$(date '+%F %T')] G21-4_rm_dup ALL --statoff 开始"

BASE_DIR="/histor/zhao/pangkun/Zhaolab/Lupeng/Zhujunjie/20260427_AS_analysis/Fig2_rmats"
GTF="/histor/zhao/pangkun/database/ref_mm10_M23_release98/mm10.gtf"
THREADS=10
READ_LEN=150
FORCE=${FORCE:-0}

SAMPLE="G21-4_rm_dup"
BAM_ALL="${BASE_DIR}/pseudobulk_bams/${SAMPLE}/${SAMPLE}_ALL_pseudobulk.bam"
OUT_DIR="${BASE_DIR}/rmats_prep_only/${SAMPLE}_ALL"
TMP_DIR="${BASE_DIR}/rmats_tmp/${SAMPLE}_ALL_tmp"

RMATS_PY=$(which rmats.py 2>/dev/null || echo "")
if [[ -z "${RMATS_PY}" ]]; then
    echo "[ERROR] 找不到 rmats.py"; exit 1
fi

if [[ ! -f "${BAM_ALL}" ]]; then
    echo "[ERROR] BAM 文件不存在: ${BAM_ALL}"; exit 1
fi

if ls "${OUT_DIR}"/*.rmats &>/dev/null 2>&1; then
    if [[ "${FORCE}" == "1" ]]; then
        echo "  [FORCE] 删除旧结果重跑"
        rm -rf "${OUT_DIR}" "${TMP_DIR}"
    else
        echo "  [SKIP] .rmats 文件已存在，跳过（FORCE=1 可强制重跑）"
        exit 0
    fi
fi

mkdir -p "${OUT_DIR}" "${TMP_DIR}"
echo "${BAM_ALL}" > "${OUT_DIR}/b1.txt"
echo "  [RUN] ${SAMPLE}_ALL"
echo "  BAM:  ${BAM_ALL}"
echo "  OUT:  ${OUT_DIR}"

python "${RMATS_PY}" \
    --b1 "${OUT_DIR}/b1.txt" \
    --gtf "${GTF}" \
    --od  "${OUT_DIR}" \
    --tmp "${TMP_DIR}" \
    -t paired \
    --readLength ${READ_LEN} \
    --variable-read-length \
    --allow-clipping \
    --nthread ${THREADS} \
    --task both \
    --statoff \
    2>&1 | tee "${OUT_DIR}/rmats.log"

cp "${TMP_DIR}"/*.rmats "${OUT_DIR}/" 2>/dev/null || true
rm -rf "${TMP_DIR}"

echo ""
echo "[$(date '+%F %T')] 完成 → ${OUT_DIR}"
