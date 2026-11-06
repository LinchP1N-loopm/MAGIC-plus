#!/bin/bash
# =============================================================
# submit_rmats_pe.sh
# 双端测序样本 rMATS 批量提交脚本
# 与 submit_rmats_all.sh 的唯一区别：-t paired（双端模式）
# 运行方式：bash submit_rmats_pe.sh
# =============================================================

SCRIPT_DIR="/path/to/project/20260516_AS_saturation/qsub/qsub_rmats_pe"
STAR_BASE="/path/to/project/20260516_AS_saturation/data/result_STARsolo"
OUT_BASE="/path/to/project/20260516_AS_saturation/rmats_result_pe"
GTF="/path/to/project/database/ref_mm10_M23_release98/mm10.gtf"
THREADS=10
READ_LEN=150

mkdir -p "${SCRIPT_DIR}"

# ========== 修改样本和梯度 ==========
SAMPLES=("G21-4" "G21-7")
LABELS=("1M" "3M" "5M" "10M" "20M" "30M" "50M" "100M" "150M" "200M" "250M" "300M")
# ====================================

for SAMPLE in "${SAMPLES[@]}"; do
    for LABEL in "${LABELS[@]}"; do
        TAG="${SAMPLE}-${LABEL}"
        BAM="${STAR_BASE}/${TAG}/STARsolo_pe/${TAG}_Aligned_AB.sortedByCoord.out.bam"
        OUT_DIR="${OUT_BASE}/${TAG}"
        TMP_DIR="${OUT_BASE}/tmp/${TAG}_tmp"
        PBS_FILE="${SCRIPT_DIR}/${TAG}_rmats_pe.pbs"

        cat > "${PBS_FILE}" << EOF
#!/bin/bash
#PBS -N ${TAG}_rmats_pe
#PBS -j oe
#PBS -d ${SCRIPT_DIR}
#PBS -l nodes=1:ppn=${THREADS}
#PBS -l walltime=72:00:00
#PBS -q gold

eval "\$(micromamba shell hook --shell bash)"
micromamba activate rmats_env

ulimit -n 4096
set -euo pipefail

BAM="${BAM}"
OUT_DIR="${OUT_DIR}"
TMP_DIR="${TMP_DIR}"

echo "[\$(date '+%F %T')] START ${TAG}"

if [[ ! -f "\${BAM}" ]]; then
    echo "[ERROR] BAM不存在: \${BAM}"; exit 1
fi

if ls "\${OUT_DIR}"/*.rmats &>/dev/null 2>&1; then
    echo "  [SKIP] .rmats已存在，跳过（FORCE=1可强制重跑）"
    exit 0
fi

mkdir -p "\${OUT_DIR}" "\${TMP_DIR}"
echo "\${BAM}" > "\${OUT_DIR}/b1.txt"

RMATS_PY=\$(which rmats.py 2>/dev/null || echo "")
if [[ -z "\${RMATS_PY}" ]]; then
    echo "[ERROR] 找不到rmats.py"; exit 1
fi

python "\${RMATS_PY}" \\
    --b1 "\${OUT_DIR}/b1.txt" \\
    --gtf "${GTF}" \\
    --od  "\${OUT_DIR}" \\
    --tmp "\${TMP_DIR}" \\
    -t paired \\
    --readLength ${READ_LEN} \\
    --variable-read-length \\
    --allow-clipping \\
    --nthread ${THREADS} \\
    --task both \\
    --statoff \\
    2>&1 | tee "\${OUT_DIR}/rmats.log"

cp "\${TMP_DIR}"/*.rmats "\${OUT_DIR}/" 2>/dev/null || true
rm -rf "\${TMP_DIR}"

echo "[\$(date '+%F %T')] DONE ${TAG} -> \${OUT_DIR}"
EOF

        echo "Submitting: ${TAG}"
        qsub "${PBS_FILE}"

    done
done

echo ""
echo "全部提交完成，共 $((${#SAMPLES[@]} * ${#LABELS[@]})) 个作业"
echo "PBS脚本: ${SCRIPT_DIR}"
echo "结果输出: ${OUT_BASE}"
