#!/bin/bash
# =============================================================
# submit_rmats_all.sh
# 自动为每个样本×梯度生成PBS脚本并提交
# 运行方式：bash submit_rmats_all.sh
# =============================================================

SCRIPT_DIR="/path/to/project/20260516_AS_saturation/qsub/qsub_rmats_all"
STAR_BASE="/path/to/project/20260516_AS_saturation/data/result_STARsolo"
OUT_BASE="/path/to/project/20260516_AS_saturation/rmats_result"
GTF="/path/to/project/database/ref_mm10_M23_release98/mm10.gtf"
THREADS=10
READ_LEN=150

mkdir -p "${SCRIPT_DIR}"

SAMPLES=("sample_reg1" "sample_reg2" "G21-4" "G21-7")
LABELS=("1M" "3M" "5M" "10M" "20M" "30M" "50M" "100M" "150M" "200M" "250M" "300M")

for SAMPLE in "${SAMPLES[@]}"; do
    for LABEL in "${LABELS[@]}"; do
        TAG="${SAMPLE}-${LABEL}"
        BAM="${STAR_BASE}/${TAG}/STARsolo_se/${TAG}_Aligned_AB.sortedByCoord.out.bam"
        OUT_DIR="${OUT_BASE}/${TAG}"
        TMP_DIR="${OUT_BASE}/tmp/${TAG}_tmp"
        PBS_FILE="${SCRIPT_DIR}/${TAG}_rmats.pbs"

        cat > "${PBS_FILE}" << EOF
#!/bin/bash
#PBS -N ${TAG}_rmats
#PBS -j oe
#PBS -d ${SCRIPT_DIR}
#PBS -l nodes=1:ppn=${THREADS}
#PBS -l walltime=72:00:00
#PBS -q gold
source /path/to/project/.bashrc

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
    -t single \\
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
