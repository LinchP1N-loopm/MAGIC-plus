#!/bin/bash
# =============================================================
# step2_generate_pseudobulk_bam_AS.sh （qsub 版）
# =============================================================
# Usage:
#   bash step2_generate_pseudobulk_bam_AS.sh <样本名> <BAM路径> [BC_SRC样本名] [THREADS]
#
# 示例：
#   # 样本名与 barcode 来源相同（最常见）
#   bash step2_generate_pseudobulk_bam_AS.sh G21-7 /path/to/G21-7.bam
#
#   # 不同 BAM 复用同一 barcode（G21-4 两个版本）
#   bash step2_generate_pseudobulk_bam_AS.sh G21-4-raw  /path/to/project/20260210_2/.../G21-4_Aligned_AB.sortedByCoord.out.bam  G21-4
#   bash step2_generate_pseudobulk_bam_AS.sh G21-4-0504 /path/to/project/20260504/.../G21-4_Aligned_AB.sortedByCoord.out.bam    G21-4
#
# 参数说明：
#   $1  SAMPLE   : 输出样本名（pseudobulk_bams 下的子目录名）
#   $2  BAM_PATH : 原始 BAM 文件路径
#   $3  BC_SRC   : barcode 来源样本名（默认与 SAMPLE 相同）
#   $4  THREADS  : 线程数（默认 10）
#
# 脚本放置目录：ASclu_qsub/
# 结果输出目录：ASclu_results/pseudobulk_bams/<SAMPLE>/
# =============================================================

set -euo pipefail

# ============================================================
# Parse arguments
# ============================================================
if [[ $# -lt 2 ]]; then
    echo "用法: bash $0 <样本名> <BAM路径> [BC_SRC样本名] [THREADS]"
    echo ""
    echo "示例（样本名=barcode来源）："
    echo "  bash $0 G21-4  /path/to/G21-4.bam"
    echo ""
    echo "示例（不同BAM复用同一barcode）："
    echo "  bash $0 G21-4-raw  /path/to/20260210.bam  G21-4"
    echo "  bash $0 G21-4-0504 /path/to/20260504.bam  G21-4"
    exit 1
fi

SAMPLE="${1}"
BAM_PATH="${2}"
BC_SRC="${3:-${SAMPLE}}"
THREADS="${4:-10}"

# ============================================================
# Paths
# ============================================================
BASE_DIR="/path/to/project/20260427_AS_analysis/Fig2_rmats"
BC_BASE="${BASE_DIR}/ASclu_results/barcodes"
PB_BASE="${BASE_DIR}/ASclu_results/pseudobulk_bams"
QSUB_DIR="${BASE_DIR}/ASclu_qsub"
PBS_FILE="${QSUB_DIR}/step2_${SAMPLE}.pbs"

mkdir -p "${QSUB_DIR}"

# ============================================================
# 检查输入
# ============================================================
if [[ ! -f "${BAM_PATH}" ]]; then
    echo "[ERROR] BAM 文件不存在: ${BAM_PATH}"
    exit 1
fi

BC_DIR="${BC_BASE}/${BC_SRC}"
if [[ ! -d "${BC_DIR}" ]]; then
    echo "[ERROR] barcode 目录不存在: ${BC_DIR}"
    echo "  请先运行: bash step1_extract_barcodes_AS.sh ${BC_SRC} <h5ad路径>"
    exit 1
fi

mapfile -t BC_FILES < <(ls "${BC_DIR}/${BC_SRC}_"*"_barcodes.txt" 2>/dev/null || true)
if [[ ${#BC_FILES[@]} -eq 0 ]]; then
    echo "[ERROR] 在 ${BC_DIR} 下未找到 *_barcodes.txt 文件"
    echo "  请先运行: bash step1_extract_barcodes_AS.sh ${BC_SRC} <h5ad路径>"
    exit 1
fi

echo "[Step2] 样本:      ${SAMPLE}"
echo "  BAM:          ${BAM_PATH}"
echo "  barcode来源:  ${BC_SRC}  (${#BC_FILES[@]} 个脑区)"
echo "  输出目录:     ${PB_BASE}/${SAMPLE}/"
echo "  PBS文件:      ${PBS_FILE}"
echo "  线程数:       ${THREADS}"
echo ""

# ============================================================
# 生成 PBS 文件
# ============================================================
cat > "${PBS_FILE}" << PBSEOF
#!/bin/bash
#PBS -N step2_${SAMPLE}
#PBS -j oe
#PBS -d ${QSUB_DIR}
#PBS -l nodes=1:ppn=${THREADS}
#PBS -l walltime=24:00:00
#PBS -q gold

eval "\$(micromamba shell hook --shell bash)"
micromamba activate rmats_env

echo "[\$(date '+%F %T')] Step2 开始: ${SAMPLE}"
echo "  BAM:         ${BAM_PATH}"
echo "  barcode来源: ${BC_SRC}"
echo "  输出目录:    ${PB_BASE}/${SAMPLE}/"
echo ""

BAM="${BAM_PATH}"
BC_BASE="${BC_BASE}"
PB_BASE="${PB_BASE}"
BC_SRC="${BC_SRC}"
SAMPLE="${SAMPLE}"
THREADS="${THREADS}"

mkdir -p "\${PB_BASE}/\${SAMPLE}"

# 检查并建立 BAM 索引
if [[ ! -f "\${BAM}.bai" ]]; then
    echo "  建立 BAM 索引..."
    samtools index -@ \${THREADS} "\${BAM}"
fi

# ---- 核心过滤函数 ----
filter_one_region() {
    local REGION=\$1
    local BC_FILE="\${BC_BASE}/\${BC_SRC}/\${BC_SRC}_\${REGION}_barcodes.txt"
    local OUT_BAM="\${PB_BASE}/\${SAMPLE}/\${SAMPLE}_\${REGION}_pseudobulk.bam"

    if [[ ! -f "\${BC_FILE}" ]]; then
        echo "  [SKIP]  \${REGION}: barcode 文件不存在"
        return
    fi

    local N_BC
    N_BC=\$(wc -l < "\${BC_FILE}")
    if [[ \${N_BC} -eq 0 ]]; then
        echo "  [SKIP]  \${REGION}: barcode 文件为空"
        return
    fi

    # 跳过已完成的
    if [[ -f "\${OUT_BAM}" && -f "\${OUT_BAM}.bai" ]]; then
        local N_READS
        N_READS=\$(samtools view -c -@ \${THREADS} "\${OUT_BAM}")
        echo "  [EXIST] \${REGION}: \${N_BC} barcodes, \${N_READS} reads（跳过）"
        return
    fi

    echo "  [RUN]   \${REGION}: \${N_BC} barcodes"

    samtools view \
        -@ \${THREADS} \
        -b \
        -D CB:"\${BC_FILE}" \
        "\${BAM}" \
        -o "\${OUT_BAM}"

    samtools index -@ \${THREADS} "\${OUT_BAM}"

    local N_READS
    N_READS=\$(samtools view -c -@ \${THREADS} "\${OUT_BAM}")
    echo "    → \${N_READS} reads → \$(basename \${OUT_BAM})"
    if [[ \${N_READS} -lt 400000 ]]; then
        echo "    ⚠  reads < 400k，后续 rMATS 结果可能不可靠"
    fi
}

# ---- 遍历所有脑区 ----
for BC_FILE in \${BC_BASE}/\${BC_SRC}/\${BC_SRC}_*_barcodes.txt; do
    BASENAME=\$(basename "\${BC_FILE}" _barcodes.txt)
    REGION=\${BASENAME#\${BC_SRC}_}
    filter_one_region "\${REGION}"
done

# ---- 输出汇总 ----
echo ""
echo "--- 输出汇总: \${SAMPLE} ---"
for OUT_BAM in \${PB_BASE}/\${SAMPLE}/*.bam; do
    [[ -f "\${OUT_BAM}" ]] || continue
    N=\$(samtools view -c -@ \${THREADS} "\${OUT_BAM}" 2>/dev/null || echo "?")
    echo "  \$(basename \${OUT_BAM}): \${N} reads"
done

echo ""
echo "[\$(date '+%F %T')] Step2 完成: ${SAMPLE}"
echo "  输出目录: ${PB_BASE}/${SAMPLE}/"
PBSEOF

echo "✓ PBS 文件已生成: ${PBS_FILE}"
echo ""
echo "提交命令："
echo "  qsub ${PBS_FILE}"
