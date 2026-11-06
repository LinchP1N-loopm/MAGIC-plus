#!/bin/bash
# =============================================================
# step3_run_rmats_AS.sh （ASclu_qsub 版）
# =============================================================
# Usage:
#   bash step3_run_rmats_AS.sh <样本名> [READ_LEN] [SE|PE] [THREADS]
#
# 示例：
#   bash step3_run_rmats_AS.sh G21-4-0504
#   bash step3_run_rmats_AS.sh G21-4-0504 150 single 10
#   bash step3_run_rmats_AS.sh G21-4-raw  150 single 10
#
# 参数说明：
#   $1  SAMPLE   : 样本名（与 pseudobulk_bams 下子目录名一致）
#   $2  READ_LEN : read 长度（默认 150）
#   $3  LIB_TYPE : single 或 paired（默认 single）
#   $4  THREADS  : 线程数（默认 10）
#
# 说明：
#   - 生成两个 PBS 文件：
#       step3_<SAMPLE>_pairwise.pbs   ← 脑区两两比较
#       step3_<SAMPLE>_ALL.pbs        ← 全脑 statoff
#   - PBS 文件放在 ASclu_qsub/ 目录
#   - 结果放在 ASclu_results/ 目录
# =============================================================

set -euo pipefail

# ============================================================
# Parse arguments
# ============================================================
if [[ $# -lt 1 ]]; then
    echo "用法: bash $0 <样本名> [READ_LEN] [single|paired] [THREADS]"
    echo ""
    echo "示例:"
    echo "  bash $0 G21-4-0504"
    echo "  bash $0 G21-4-0504 150 single 10"
    echo "  bash $0 G21-4-raw  150 single 10"
    exit 1
fi

SAMPLE="${1}"
READ_LEN="${2:-150}"
LIB_TYPE="${3:-single}"
THREADS="${4:-10}"

# ============================================================
# Paths
# ============================================================
BASE_DIR="/path/to/project/20260427_AS_analysis/Fig2_rmats"
PB_BASE="${BASE_DIR}/ASclu_results/pseudobulk_bams"
PAIRWISE_OUT="${BASE_DIR}/ASclu_results/rmats_pairwise"
PREP_OUT="${BASE_DIR}/ASclu_results/rmats_prep_only"
TMP_BASE="${BASE_DIR}/ASclu_results/rmats_tmp"
QSUB_DIR="${BASE_DIR}/ASclu_qsub"
GTF="/path/to/project/database/ref_mm10_M23_release98/mm10.gtf"

PBS_PAIRWISE="${QSUB_DIR}/step3_${SAMPLE}_pairwise.pbs"
PBS_ALL="${QSUB_DIR}/step3_${SAMPLE}_ALL.pbs"

mkdir -p "${QSUB_DIR}"

# ============================================================
# 检查 BAM 目录
# ============================================================
PB_DIR="${PB_BASE}/${SAMPLE}"
if [[ ! -d "${PB_DIR}" ]]; then
    echo "[ERROR] BAM 目录不存在: ${PB_DIR}"
    echo "  请先运行: bash step2_generate_pseudobulk_bam_AS.sh ${SAMPLE} <BAM路径>"
    exit 1
fi

# 收集脑区列表（排除 ALL）
mapfile -t REGIONS < <(
    ls "${PB_DIR}/${SAMPLE}_"*"_pseudobulk.bam" 2>/dev/null \
    | xargs -I{} basename {} _pseudobulk.bam \
    | sed "s/^${SAMPLE}_//" \
    | grep -v "^ALL$" \
    | sort
)

if [[ ${#REGIONS[@]} -eq 0 ]]; then
    echo "[ERROR] 未找到 pseudobulk BAM 文件: ${PB_DIR}/"
    echo "  请先运行 step2"
    exit 1
fi

N=${#REGIONS[@]}
NPAIRS=$(( N * (N-1) / 2 ))

echo "[Step3] 样本:      ${SAMPLE}"
echo "  BAM目录:      ${PB_DIR}/"
echo "  脑区数:       ${N}  → 两两对数: ${NPAIRS}"
echo "  脑区列表:     ${REGIONS[*]}"
echo "  read长度:     ${READ_LEN}"
echo "  文库类型:     ${LIB_TYPE}"
echo "  线程数:       ${THREADS}"
echo "  pairwise PBS: ${PBS_PAIRWISE}"
echo "  ALL PBS:      ${PBS_ALL}"
echo ""

# ============================================================
# 生成 PBS1：pairwise 两两比较
# ============================================================
# 把脑区列表序列化传入PBS
REGIONS_STR="${REGIONS[*]}"

cat > "${PBS_PAIRWISE}" << PBSEOF
#!/bin/bash
#PBS -N step3_${SAMPLE}_pw
#PBS -j oe
#PBS -d ${QSUB_DIR}
#PBS -l nodes=1:ppn=${THREADS}
#PBS -l walltime=120:00:00
#PBS -q batch

eval "\$(micromamba shell hook --shell bash)"
micromamba activate rmats_env

echo "[\$(date '+%F %T')] Step3 pairwise 开始: ${SAMPLE}"
echo "  脑区列表: ${REGIONS_STR}"
echo ""

BASE_DIR="${BASE_DIR}"
PB_BASE="${PB_BASE}"
PAIRWISE_OUT="${PAIRWISE_OUT}"
TMP_BASE="${TMP_BASE}"
GTF="${GTF}"
SAMPLE="${SAMPLE}"
READ_LEN="${READ_LEN}"
LIB_TYPE="${LIB_TYPE}"
THREADS="${THREADS}"
FORCE=\${FORCE:-0}

mkdir -p "\${PAIRWISE_OUT}/\${SAMPLE}" "\${TMP_BASE}"

RMATS_PY=\$(which rmats.py 2>/dev/null || echo "")
if [[ -z "\${RMATS_PY}" ]]; then
    echo "[ERROR] 找不到 rmats.py"; exit 1
fi

REGIONS=(${REGIONS_STR})
N=\${#REGIONS[@]}
DONE=0
SKIP=0
FAIL=0

for (( i=0; i<N-1; i++ )); do
    for (( j=i+1; j<N; j++ )); do
        REG_I="\${REGIONS[\$i]}"
        REG_J="\${REGIONS[\$j]}"
        LABEL="\${REG_I}_vs_\${REG_J}"
        OUT_DIR="\${PAIRWISE_OUT}/\${SAMPLE}/\${LABEL}"
        TMP_DIR="\${TMP_BASE}/\${SAMPLE}_\${LABEL}"

        BAM_I="\${PB_BASE}/\${SAMPLE}/\${SAMPLE}_\${REG_I}_pseudobulk.bam"
        BAM_J="\${PB_BASE}/\${SAMPLE}/\${SAMPLE}_\${REG_J}_pseudobulk.bam"

        if [[ ! -f "\${BAM_I}" || ! -f "\${BAM_J}" ]]; then
            echo "  [SKIP]  \${LABEL}: BAM 缺失"
            (( SKIP++ )) || true
            continue
        fi

        if [[ -f "\${OUT_DIR}/SE.MATS.JC.txt" && "\${FORCE}" != "1" ]]; then
            echo "  [SKIP]  \${LABEL}: 已完成"
            (( SKIP++ )) || true
            continue
        fi

        [[ "\${FORCE}" == "1" ]] && rm -rf "\${OUT_DIR}" "\${TMP_DIR}"
        mkdir -p "\${OUT_DIR}" "\${TMP_DIR}"
        echo "\${BAM_I}" > "\${OUT_DIR}/b1.txt"
        echo "\${BAM_J}" > "\${OUT_DIR}/b2.txt"
        echo "  [RUN]   \${LABEL}"

        python "\${RMATS_PY}" \\
            --b1 "\${OUT_DIR}/b1.txt" \\
            --b2 "\${OUT_DIR}/b2.txt" \\
            --gtf "\${GTF}" \\
            --od  "\${OUT_DIR}" \\
            --tmp "\${TMP_DIR}" \\
            -t \${LIB_TYPE} \\
            --readLength \${READ_LEN} \\
            --variable-read-length \\
            --allow-clipping \\
            --nthread \${THREADS} \\
            --task both \\
            2>&1 | tee "\${OUT_DIR}/rmats.log"

        ok=1
        for AS in SE A5SS A3SS MXE RI; do
            [[ ! -f "\${OUT_DIR}/\${AS}.MATS.JC.txt" ]] && ok=0
        done
        if [[ \${ok} -eq 1 ]]; then
            echo "    ✓ \${LABEL} 完成"
            (( DONE++ )) || true
        else
            echo "    ⚠  \${LABEL} 输出不完整，请检查 rmats.log"
            (( FAIL++ )) || true
        fi
        rm -rf "\${TMP_DIR}"
    done
done

echo ""
echo "[\$(date '+%F %T')] Step3 pairwise 完成: ${SAMPLE}"
echo "  完成: \${DONE}  跳过: \${SKIP}  失败: \${FAIL}"
echo "  结果目录: \${PAIRWISE_OUT}/\${SAMPLE}/"
PBSEOF

# ============================================================
# 生成 PBS2：ALL statoff
# ============================================================
cat > "${PBS_ALL}" << PBSEOF
#!/bin/bash
#PBS -N step3_${SAMPLE}_ALL
#PBS -j oe
#PBS -d ${QSUB_DIR}
#PBS -l nodes=1:ppn=${THREADS}
#PBS -l walltime=24:00:00
#PBS -q batch

eval "\$(micromamba shell hook --shell bash)"
micromamba activate rmats_env

echo "[\$(date '+%F %T')] Step3 ALL statoff 开始: ${SAMPLE}"

BASE_DIR="${BASE_DIR}"
PB_BASE="${PB_BASE}"
PREP_OUT="${PREP_OUT}"
TMP_BASE="${TMP_BASE}"
GTF="${GTF}"
SAMPLE="${SAMPLE}"
READ_LEN="${READ_LEN}"
LIB_TYPE="${LIB_TYPE}"
THREADS="${THREADS}"
FORCE=\${FORCE:-0}

mkdir -p "\${PREP_OUT}" "\${TMP_BASE}"

RMATS_PY=\$(which rmats.py 2>/dev/null || echo "")
if [[ -z "\${RMATS_PY}" ]]; then
    echo "[ERROR] 找不到 rmats.py"; exit 1
fi

BAM_ALL="\${PB_BASE}/\${SAMPLE}/\${SAMPLE}_ALL_pseudobulk.bam"
OUT_DIR="\${PREP_OUT}/\${SAMPLE}_ALL"
TMP_DIR="\${TMP_BASE}/\${SAMPLE}_ALL_tmp"

if [[ ! -f "\${BAM_ALL}" ]]; then
    echo "[ERROR] ALL BAM 不存在: \${BAM_ALL}"
    exit 1
fi

if ls "\${OUT_DIR}"/*.rmats &>/dev/null 2>&1 && [[ "\${FORCE}" != "1" ]]; then
    echo "  [SKIP] .rmats 已存在，跳过（FORCE=1 可强制重跑）"
    exit 0
fi

[[ "\${FORCE}" == "1" ]] && rm -rf "\${OUT_DIR}" "\${TMP_DIR}"
mkdir -p "\${OUT_DIR}" "\${TMP_DIR}"
echo "\${BAM_ALL}" > "\${OUT_DIR}/b1.txt"
echo "  [RUN] \${SAMPLE}_ALL"

python "\${RMATS_PY}" \\
    --b1 "\${OUT_DIR}/b1.txt" \\
    --gtf "\${GTF}" \\
    --od  "\${OUT_DIR}" \\
    --tmp "\${TMP_DIR}" \\
    -t \${LIB_TYPE} \\
    --readLength \${READ_LEN} \\
    --variable-read-length \\
    --allow-clipping \\
    --nthread \${THREADS} \\
    --task both \\
    --statoff \\
    2>&1 | tee "\${OUT_DIR}/rmats.log"

cp "\${TMP_DIR}"/*.rmats "\${OUT_DIR}/" 2>/dev/null || true
rm -rf "\${TMP_DIR}"

ok=1
for AS in SE A5SS A3SS MXE RI; do
    [[ ! -f "\${OUT_DIR}/\${AS}.MATS.JC.txt" ]] && ok=0
done
[[ \${ok} -eq 1 ]] && echo "  ✓ 5种AS输出均存在" || echo "  ⚠  输出不完整，请检查 rmats.log"

echo "[\$(date '+%F %T')] Step3 ALL statoff 完成: ${SAMPLE}"
echo "  结果目录: \${OUT_DIR}"
PBSEOF

# ============================================================
# 输出提交命令
# ============================================================
echo "✓ PBS 文件已生成:"
echo "  ${PBS_PAIRWISE}"
echo "  ${PBS_ALL}"
echo ""
echo "提交命令（两者可同时提交）："
echo "  qsub ${PBS_PAIRWISE}"
echo "  qsub ${PBS_ALL}"
