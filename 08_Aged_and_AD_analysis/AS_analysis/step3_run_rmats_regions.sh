#!/bin/bash
# =============================================================
# step3_run_rmats_regions.sh
# =============================================================
# Purpose:
#   基于 Step2 生成的 pseudobulk BAM，为三组比较生成 rMATS PBS 脚本。
#
# 新版本特点：
#   1. 工作目录固定为 20260705_Fig3/rmats_works
#   2. 默认 -t paired，不再默认 single
#   3. Regions 从 Step1/Step0 输出动态读取，不再写死旧 Cortex_1/Cortex_2 名称
#   4. 输出目录为 rmats_results，临时目录为 rmats_tmp
#
# 用法：
#   bash step3_run_rmats_regions.sh [READ_LEN] [paired|single] [THREADS]
#
# 示例：
#   bash step3_run_rmats_regions.sh              # 默认: 150 paired 10
#   bash step3_run_rmats_regions.sh 150 paired 16
#
# 固定三组比较：
#   control_vs_old : b1=G21-4,G21-7   b2=G21-1,G21-8
#   control_vs_AD  : b1=G21-4,G21-7   b2=G21-2,G21-5
#   old_vs_AD      : b1=G21-1,G21-8   b2=G21-2,G21-5
# =============================================================

set -euo pipefail

# ============================================================
# Parameters
# ============================================================
READ_LEN="${1:-150}"
LIB_TYPE="${2:-paired}"
THREADS="${3:-10}"

if [[ ! "${READ_LEN}" =~ ^[0-9]+$ ]]; then
    echo "[ERROR] READ_LEN 必须是整数: ${READ_LEN}" >&2
    exit 1
fi

if [[ "${LIB_TYPE}" != "paired" && "${LIB_TYPE}" != "single" ]]; then
    echo "[ERROR] LIB_TYPE 只能是 paired 或 single: ${LIB_TYPE}" >&2
    exit 1
fi

if [[ ! "${THREADS}" =~ ^[0-9]+$ ]]; then
    echo "[ERROR] THREADS 必须是整数: ${THREADS}" >&2
    exit 1
fi

# ============================================================
# Paths
# ============================================================
BASE_DIR="/path/to/project/20260705_Fig3/rmats_works"
BC_BASE="${BASE_DIR}/barcodes"
PB_BASE="${BASE_DIR}/pseudobulk_bams"
RMATS_OUT="${BASE_DIR}/rmats_results"
TMP_BASE="${BASE_DIR}/rmats_tmp"
QSUB_DIR="${BASE_DIR}/qsub"
LOG_DIR="${BASE_DIR}/logs"
META_DIR="${BASE_DIR}/metadata"
GTF="/path/to/project/database/ref_mm10_M23_release98/mm10.gtf"

mkdir -p "${QSUB_DIR}" "${LOG_DIR}" "${RMATS_OUT}" "${TMP_BASE}" "${META_DIR}"

# ============================================================
# Sample groups
# ============================================================
SAMPLES=(G21-1 G21-2 G21-4 G21-5 G21-7 G21-8)
COMPARISONS=(
    "control_vs_old:G21-4,G21-7:G21-1,G21-8"
    "control_vs_AD:G21-4,G21-7:G21-2,G21-5"
    "old_vs_AD:G21-1,G21-8:G21-2,G21-5"
)

# ============================================================
# Helper functions
# ============================================================
shell_join() {
    printf "%q " "$@"
}

load_regions_from_metadata() {
    local counts_valid="${META_DIR}/step0_regions_spot_counts_valid.tsv"
    local presence_valid="${META_DIR}/step0_regions_presence_valid.tsv"

    # 优先使用 Step0 输出的有效区域矩阵表头
    if [[ -f "${counts_valid}" ]]; then
        head -n 1 "${counts_valid}" | awk -F'\t' '{for (i=2; i<=NF; i++) print $i}'
        return 0
    fi

    if [[ -f "${presence_valid}" ]]; then
        head -n 1 "${presence_valid}" | awk -F'\t' '{for (i=2; i<=NF; i++) print $i}'
        return 0
    fi

    return 1
}

load_regions_from_step1() {
    local first_sample="G21-1"
    local region_list=""

    if [[ -f "${BC_BASE}/${first_sample}/${first_sample}_region_list.tsv" ]]; then
        region_list="${BC_BASE}/${first_sample}/${first_sample}_region_list.tsv"
    elif [[ -f "${BC_BASE}/${first_sample}/${first_sample}_cluster_list.tsv" ]]; then
        # 兼容旧名
        region_list="${BC_BASE}/${first_sample}/${first_sample}_cluster_list.tsv"
    else
        return 1
    fi

    awk -F'\t' 'NF >= 1 && $1 != "" && $1 !~ /^#/ {print $1}' "${region_list}"
}

load_regions() {
    if load_regions_from_metadata; then
        return 0
    fi
    if load_regions_from_step1; then
        return 0
    fi
    return 1
}

make_bam_list_file() {
    local out_txt="$1"
    shift
    local bams=("$@")
    (IFS=,; echo "${bams[*]}") > "${out_txt}"
}

# ============================================================
# Read regions
# ============================================================
mapfile -t REGIONS < <(load_regions)

if [[ "${#REGIONS[@]}" -eq 0 ]]; then
    echo "[ERROR] 未能读取有效 Regions。" >&2
    echo "        请确认以下文件至少存在一个：" >&2
    echo "        ${META_DIR}/step0_regions_spot_counts_valid.tsv" >&2
    echo "        ${BC_BASE}/G21-1/G21-1_region_list.tsv" >&2
    exit 1
fi

# Drop empty fields
REGIONS_CLEAN=()
for r in "${REGIONS[@]}"; do
    [[ -n "${r}" ]] && REGIONS_CLEAN+=("${r}")
done
REGIONS=("${REGIONS_CLEAN[@]}")

# ============================================================
# Run information
# ============================================================
echo "================================================================="
echo "[Step3] rMATS PBS 生成"
echo "================================================================="
echo "  BASE_DIR:  ${BASE_DIR}"
echo "  PB_BASE:   ${PB_BASE}"
echo "  OUT:       ${RMATS_OUT}"
echo "  TMP:       ${TMP_BASE}"
echo "  GTF:       ${GTF}"
echo "  READ_LEN:  ${READ_LEN}"
echo "  LIB_TYPE:  ${LIB_TYPE}"
echo "  THREADS:   ${THREADS}"
echo "  Regions:   ${#REGIONS[@]} → ${REGIONS[*]}"
echo ""

if [[ ! -f "${GTF}" ]]; then
    echo "[ERROR] GTF 文件不存在: ${GTF}" >&2
    exit 1
fi

# ============================================================
# Check for pseudobulk BAM files
# ============================================================
echo "================================================================="
echo "[Step3] 检查 pseudobulk BAM 文件..."
echo "================================================================="

MISSING=0
for sample in "${SAMPLES[@]}"; do
    pb_dir="${PB_BASE}/${sample}"
    if [[ ! -d "${pb_dir}" ]]; then
        echo "  [MISS] 目录不存在: ${pb_dir}"
        MISSING=1
        continue
    fi

    bam_all="${pb_dir}/${sample}_ALL_pseudobulk.bam"
    if [[ ! -f "${bam_all}" ]]; then
        echo "  [MISS] ${sample}_ALL_pseudobulk.bam"
        MISSING=1
    fi

    for region in "${REGIONS[@]}"; do
        bam="${pb_dir}/${sample}_${region}_pseudobulk.bam"
        if [[ ! -f "${bam}" ]]; then
            echo "  [MISS] ${sample}_${region}_pseudobulk.bam"
            MISSING=1
        fi
    done
done

if [[ "${MISSING}" -eq 1 ]]; then
    echo ""
    echo "[ERROR] 存在缺失的 pseudobulk BAM，请先完成 Step2 再运行 Step3。" >&2
    exit 1
fi

echo "  ✓ 所有 pseudobulk BAM 文件检查通过"
echo ""

# ============================================================
# 为每个比较组生成 PBS 文件
# ============================================================
PBS_FILES=()
REGIONS_STR="$(shell_join "${REGIONS[@]}")"

for comp_def in "${COMPARISONS[@]}"; do
    comp_name="${comp_def%%:*}"
    rest="${comp_def#*:}"
    b1_samples_str="${rest%%:*}"
    b2_samples_str="${rest##*:}"

    IFS=',' read -ra B1_SAMPLES <<< "${b1_samples_str}"
    IFS=',' read -ra B2_SAMPLES <<< "${b2_samples_str}"

    b1_str="$(shell_join "${B1_SAMPLES[@]}")"
    b2_str="$(shell_join "${B2_SAMPLES[@]}")"

    echo "================================================================="
    echo "[Step3] 比较组: ${comp_name}"
    echo "  b1 样本: ${B1_SAMPLES[*]}"
    echo "  b2 样本: ${B2_SAMPLES[*]}"
    echo "================================================================="

    pbs_regions="${QSUB_DIR}/step3_${comp_name}_regions.pbs"
    pbs_all="${QSUB_DIR}/step3_${comp_name}_ALL.pbs"
    PBS_FILES+=("${pbs_regions}" "${pbs_all}")

    out_base="${RMATS_OUT}/${comp_name}"
    mkdir -p "${out_base}"

    # ----------------------------------------------------------
    # PBS1: 各 Regions 比较
    # ----------------------------------------------------------
    cat > "${pbs_regions}" <<PBSEOF
#!/bin/bash
#PBS -N step3_${comp_name}_reg
#PBS -j oe
#PBS -o ${LOG_DIR}/step3_${comp_name}_regions.pbs.log
#PBS -d ${QSUB_DIR}
#PBS -l nodes=1:ppn=${THREADS}
#PBS -l walltime=120:00:00
#PBS -q gold

set -euo pipefail

eval "\$(micromamba shell hook --shell bash)"
micromamba activate rmats_env

echo "[\$(date '+%F %T')] Step3 regions 开始: ${comp_name}"
echo "  LIB_TYPE: ${LIB_TYPE}"
echo ""

PB_BASE="${PB_BASE}"
OUT_BASE="${out_base}"
TMP_BASE="${TMP_BASE}"
GTF="${GTF}"
READ_LEN="${READ_LEN}"
LIB_TYPE="${LIB_TYPE}"
THREADS="${THREADS}"
FORCE=\${FORCE:-0}

B1_SAMPLES=(${b1_str})
B2_SAMPLES=(${b2_str})
REGIONS=(${REGIONS_STR})

RMATS_PY=\$(which rmats.py 2>/dev/null || echo "")
if [[ -z "\${RMATS_PY}" ]]; then
    echo "[ERROR] 找不到 rmats.py，请检查 rmats_env 环境。" >&2
    exit 1
fi

DONE=0; SKIP=0; FAIL=0

for REGION in "\${REGIONS[@]}"; do
    OUT_DIR="\${OUT_BASE}/\${REGION}"
    TMP_DIR="\${TMP_BASE}/${comp_name}_\${REGION}"

    if [[ -f "\${OUT_DIR}/SE.MATS.JC.txt" && "\${FORCE}" != "1" ]]; then
        echo "  [SKIP]  \${REGION}: 已完成（FORCE=1 可强制重跑）"
        (( SKIP++ )) || true
        continue
    fi

    B1_BAMS=()
    ALL_EXIST=1
    for S in "\${B1_SAMPLES[@]}"; do
        BAM="\${PB_BASE}/\${S}/\${S}_\${REGION}_pseudobulk.bam"
        if [[ ! -f "\${BAM}" ]]; then
            echo "  [SKIP]  \${REGION}: 缺少 b1 BAM: \${BAM}"
            ALL_EXIST=0
            break
        fi
        B1_BAMS+=("\${BAM}")
    done
    [[ "\${ALL_EXIST}" -eq 0 ]] && { (( SKIP++ )) || true; continue; }

    B2_BAMS=()
    for S in "\${B2_SAMPLES[@]}"; do
        BAM="\${PB_BASE}/\${S}/\${S}_\${REGION}_pseudobulk.bam"
        if [[ ! -f "\${BAM}" ]]; then
            echo "  [SKIP]  \${REGION}: 缺少 b2 BAM: \${BAM}"
            ALL_EXIST=0
            break
        fi
        B2_BAMS+=("\${BAM}")
    done
    [[ "\${ALL_EXIST}" -eq 0 ]] && { (( SKIP++ )) || true; continue; }

    if [[ "\${FORCE}" == "1" ]]; then
        rm -rf "\${OUT_DIR}" "\${TMP_DIR}"
    fi
    mkdir -p "\${OUT_DIR}" "\${TMP_DIR}"

    (IFS=,; echo "\${B1_BAMS[*]}") > "\${OUT_DIR}/b1.txt"
    (IFS=,; echo "\${B2_BAMS[*]}") > "\${OUT_DIR}/b2.txt"

    echo "  [RUN]   \${REGION}"
    echo "    b1: \$(cat "\${OUT_DIR}/b1.txt")"
    echo "    b2: \$(cat "\${OUT_DIR}/b2.txt")"

    python "\${RMATS_PY}" \
        --b1 "\${OUT_DIR}/b1.txt" \
        --b2 "\${OUT_DIR}/b2.txt" \
        --gtf "\${GTF}" \
        --od  "\${OUT_DIR}" \
        --tmp "\${TMP_DIR}" \
        -t "\${LIB_TYPE}" \
        --readLength "\${READ_LEN}" \
        --variable-read-length \
        --allow-clipping \
        --nthread "\${THREADS}" \
        --task both \
        2>&1 | tee "\${OUT_DIR}/rmats.log"

    ok=1
    for AS in SE A5SS A3SS MXE RI; do
        [[ ! -f "\${OUT_DIR}/\${AS}.MATS.JC.txt" ]] && ok=0
    done

    if [[ "\${ok}" -eq 1 ]]; then
        echo "    ✓ \${REGION} 完成"
        (( DONE++ )) || true
    else
        echo "    ⚠  \${REGION} 输出不完整，请检查 rmats.log"
        (( FAIL++ )) || true
    fi

    rm -rf "\${TMP_DIR}"
done

echo ""
echo "[\$(date '+%F %T')] Step3 regions 完成: ${comp_name}"
echo "  完成: \${DONE}  跳过: \${SKIP}  失败: \${FAIL}"
echo "  结果目录: \${OUT_BASE}/"
PBSEOF

    # ----------------------------------------------------------
    # PBS2: ALL 比较
    # ----------------------------------------------------------
    cat > "${pbs_all}" <<PBSEOF
#!/bin/bash
#PBS -N step3_${comp_name}_ALL
#PBS -j oe
#PBS -o ${LOG_DIR}/step3_${comp_name}_ALL.pbs.log
#PBS -d ${QSUB_DIR}
#PBS -l nodes=1:ppn=${THREADS}
#PBS -l walltime=48:00:00
#PBS -q batch

set -euo pipefail

eval "\$(micromamba shell hook --shell bash)"
micromamba activate rmats_env

echo "[\$(date '+%F %T')] Step3 ALL 开始: ${comp_name}"
echo "  LIB_TYPE: ${LIB_TYPE}"
echo ""

PB_BASE="${PB_BASE}"
OUT_BASE="${out_base}"
TMP_BASE="${TMP_BASE}"
GTF="${GTF}"
READ_LEN="${READ_LEN}"
LIB_TYPE="${LIB_TYPE}"
THREADS="${THREADS}"
FORCE=\${FORCE:-0}

B1_SAMPLES=(${b1_str})
B2_SAMPLES=(${b2_str})

RMATS_PY=\$(which rmats.py 2>/dev/null || echo "")
if [[ -z "\${RMATS_PY}" ]]; then
    echo "[ERROR] 找不到 rmats.py，请检查 rmats_env 环境。" >&2
    exit 1
fi

OUT_DIR="\${OUT_BASE}/ALL"
TMP_DIR="\${TMP_BASE}/${comp_name}_ALL"

if [[ -f "\${OUT_DIR}/SE.MATS.JC.txt" && "\${FORCE}" != "1" ]]; then
    echo "  [SKIP] ALL: 已完成（FORCE=1 可强制重跑）"
    exit 0
fi

B1_BAMS=()
for S in "\${B1_SAMPLES[@]}"; do
    BAM="\${PB_BASE}/\${S}/\${S}_ALL_pseudobulk.bam"
    if [[ ! -f "\${BAM}" ]]; then
        echo "[ERROR] 缺少 b1 ALL BAM: \${BAM}" >&2
        exit 1
    fi
    B1_BAMS+=("\${BAM}")
done

B2_BAMS=()
for S in "\${B2_SAMPLES[@]}"; do
    BAM="\${PB_BASE}/\${S}/\${S}_ALL_pseudobulk.bam"
    if [[ ! -f "\${BAM}" ]]; then
        echo "[ERROR] 缺少 b2 ALL BAM: \${BAM}" >&2
        exit 1
    fi
    B2_BAMS+=("\${BAM}")
done

if [[ "\${FORCE}" == "1" ]]; then
    rm -rf "\${OUT_DIR}" "\${TMP_DIR}"
fi
mkdir -p "\${OUT_DIR}" "\${TMP_DIR}"

(IFS=,; echo "\${B1_BAMS[*]}") > "\${OUT_DIR}/b1.txt"
(IFS=,; echo "\${B2_BAMS[*]}") > "\${OUT_DIR}/b2.txt"

echo "  [RUN] ALL"
echo "    b1: \$(cat "\${OUT_DIR}/b1.txt")"
echo "    b2: \$(cat "\${OUT_DIR}/b2.txt")"

python "\${RMATS_PY}" \
    --b1 "\${OUT_DIR}/b1.txt" \
    --b2 "\${OUT_DIR}/b2.txt" \
    --gtf "\${GTF}" \
    --od  "\${OUT_DIR}" \
    --tmp "\${TMP_DIR}" \
    -t "\${LIB_TYPE}" \
    --readLength "\${READ_LEN}" \
    --variable-read-length \
    --allow-clipping \
    --nthread "\${THREADS}" \
    --task both \
    2>&1 | tee "\${OUT_DIR}/rmats.log"

ok=1
for AS in SE A5SS A3SS MXE RI; do
    [[ ! -f "\${OUT_DIR}/\${AS}.MATS.JC.txt" ]] && ok=0
done

if [[ "\${ok}" -eq 1 ]]; then
    echo "  ✓ ALL 完成"
else
    echo "  ⚠  ALL 输出不完整，请检查 rmats.log"
fi

rm -rf "\${TMP_DIR}"

echo "[\$(date '+%F %T')] Step3 ALL 完成: ${comp_name}"
echo "  结果目录: \${OUT_DIR}"
PBSEOF

    echo "  ✓ 生成:"
    echo "    ${pbs_regions}"
    echo "    ${pbs_all}"
    echo ""
done

# ============================================================
# Print submission commands
# ============================================================
echo "================================================================="
echo "所有 Step3 PBS 文件已生成，提交命令："
echo "================================================================="
for pbs in "${PBS_FILES[@]}"; do
    echo "  qsub ${pbs}"
done

echo ""
echo "预期输出目录结构："
echo "  ${RMATS_OUT}/"
for comp_def in "${COMPARISONS[@]}"; do
    comp_name="${comp_def%%:*}"
    echo "  ├── ${comp_name}/"
    echo "  │   ├── ALL/"
    for region in "${REGIONS[@]}"; do
        echo "  │   ├── ${region}/"
    done
done
