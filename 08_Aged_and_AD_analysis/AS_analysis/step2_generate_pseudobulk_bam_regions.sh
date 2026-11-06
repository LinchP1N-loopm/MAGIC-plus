#!/bin/bash
# =============================================================
# step2_generate_pseudobulk_bam_regions.sh
# =============================================================
# Purpose:
#   读取 Step1 生成的 barcode 文件，用 samtools -D CB: 从原始 STARsolo BAM
#   中过滤每个 Regions 的 reads，生成 pseudobulk BAM，并生成 PBS/qsub 文件。
#
# 新版本特点：
#   1. 工作目录固定为 20260705_Fig3/rmats_works
#   2. 使用 Step1 输出的 *_region_list.tsv；兼容旧名 *_cluster_list.tsv
#   3. 内置 6 个样本的原始 BAM 路径
#   4. 支持一键为全部样本生成 PBS
#
# 用法：
#   bash step2_generate_pseudobulk_bam_regions.sh
#   bash step2_generate_pseudobulk_bam_regions.sh all [THREADS]
#   bash step2_generate_pseudobulk_bam_regions.sh <SAMPLE> [THREADS]
#   bash step2_generate_pseudobulk_bam_regions.sh <SAMPLE> <BAM_PATH> [THREADS]
#
# 示例：
#   bash step2_generate_pseudobulk_bam_regions.sh all 10
#   bash step2_generate_pseudobulk_bam_regions.sh G21-1 10
#   bash step2_generate_pseudobulk_bam_regions.sh G21-1 /path/to/G21-1.bam 16
# =============================================================

set -euo pipefail

# ============================================================
# Paths
# ============================================================
BASE_DIR="/path/to/project/20260705_Fig3/rmats_works"
BC_BASE="${BASE_DIR}/barcodes"
PB_BASE="${BASE_DIR}/pseudobulk_bams"
QSUB_DIR="${BASE_DIR}/qsub"
LOG_DIR="${BASE_DIR}/logs"
META_DIR="${BASE_DIR}/metadata"

mkdir -p "${QSUB_DIR}" "${LOG_DIR}" "${PB_BASE}" "${META_DIR}"

# ============================================================
# 样本与 BAM 路径
# ============================================================
SAMPLES=(G21-1 G21-2 G21-4 G21-5 G21-7 G21-8)

declare -A BAM_MAP
BAM_MAP["G21-1"]="/path/to/project/20260305_rejiace/data/result_STARsolo/G21-1/STARsolo_pe/G21-1_Aligned_AB.sortedByCoord.out.bam"
BAM_MAP["G21-2"]="/path/to/project/20260305_rejiace/data/result_STARsolo/G21-2/STARsolo_pe/G21-2_Aligned_AB.sortedByCoord.out.bam"
BAM_MAP["G21-4"]="/path/to/project/20260305_rejiace/data/result_STARsolo/G21-4/STARsolo_pe/G21-4_Aligned_AB.sortedByCoord.out.bam"
BAM_MAP["G21-5"]="/path/to/project/20260328/data/raw_result_STARsolo/G21-5/STARsolo_pe/G21-5_Aligned_AB.sortedByCoord.out.bam"
BAM_MAP["G21-7"]="/path/to/project/20260322/data/result_STARsolo/G21-7/STARsolo_pe/G21-7_Aligned_AB.sortedByCoord.out.bam"
BAM_MAP["G21-8"]="/path/to/project/20260328/data/raw_result_STARsolo/G21-8/STARsolo_pe/G21-8_Aligned_AB.sortedByCoord.out.bam"

# ============================================================
# Helper functions
# ============================================================
usage() {
    cat <<EOF
用法:
  bash $0
  bash $0 all [THREADS]
  bash $0 <SAMPLE> [THREADS]
  bash $0 <SAMPLE> <BAM_PATH> [THREADS]

样本:
  ${SAMPLES[*]}

默认:
  TARGET=all
  THREADS=10

示例:
  bash $0 all 10
  bash $0 G21-1 10
  bash $0 G21-1 /path/to/G21-1.bam 16
EOF
}

is_integer() {
    [[ "${1:-}" =~ ^[0-9]+$ ]]
}

sample_exists() {
    local sample="$1"
    local s
    for s in "${SAMPLES[@]}"; do
        [[ "${sample}" == "${s}" ]] && return 0
    done
    return 1
}

shell_join() {
    # 把数组安全转成 shell 可解析字符串
    printf "%q " "$@"
}

load_regions_for_sample() {
    local sample="$1"
    local bc_dir="${BC_BASE}/${sample}"
    local region_list=""

    if [[ -f "${bc_dir}/${sample}_region_list.tsv" ]]; then
        region_list="${bc_dir}/${sample}_region_list.tsv"
    elif [[ -f "${bc_dir}/${sample}_cluster_list.tsv" ]]; then
        # 兼容 Step1 早期输出名
        region_list="${bc_dir}/${sample}_cluster_list.tsv"
    else
        echo "[ERROR] 找不到区域列表文件:" >&2
        echo "        ${bc_dir}/${sample}_region_list.tsv" >&2
        echo "        ${bc_dir}/${sample}_cluster_list.tsv" >&2
        echo "        请先运行 Step1。" >&2
        return 1
    fi

    # 第一列为 safe region name；跳过空行和注释行
    awk -F'\t' 'NF >= 1 && $1 != "" && $1 !~ /^#/ {print $1}' "${region_list}"
}

generate_one_pbs() {
    local sample="$1"
    local bam_path="$2"
    local threads="$3"

    local bc_dir="${BC_BASE}/${sample}"
    local pb_dir="${PB_BASE}/${sample}"
    local pbs_file="${QSUB_DIR}/step2_${sample}.pbs"
    local read_count_tsv="${META_DIR}/step2_${sample}_pseudobulk_read_counts.tsv"

    echo "================================================================="
    echo "[Step2] 准备样本: ${sample}"
    echo "================================================================="
    echo "  BAM:          ${bam_path}"
    echo "  barcode目录:  ${bc_dir}"
    echo "  输出目录:     ${pb_dir}"
    echo "  PBS文件:      ${pbs_file}"
    echo "  线程数:       ${threads}"

    if [[ ! -f "${bam_path}" ]]; then
        echo "[ERROR] BAM 文件不存在: ${bam_path}" >&2
        return 1
    fi

    if [[ ! -d "${bc_dir}" ]]; then
        echo "[ERROR] barcode 目录不存在: ${bc_dir}" >&2
        echo "        请先运行 Step1。" >&2
        return 1
    fi

    mapfile -t REGIONS < <(load_regions_for_sample "${sample}")
    if [[ "${#REGIONS[@]}" -eq 0 ]]; then
        echo "[ERROR] ${sample} 没有读取到任何有效 Regions。" >&2
        return 1
    fi

    echo "  有效区域数:   ${#REGIONS[@]}"
    echo "  有效区域:     ${REGIONS[*]}"

    local all_bc="${bc_dir}/${sample}_ALL_barcodes.txt"
    if [[ ! -s "${all_bc}" ]]; then
        echo "[ERROR] 找不到 ALL barcode 文件或文件为空: ${all_bc}" >&2
        return 1
    fi

    local missing_bc=0
    local region
    for region in "${REGIONS[@]}"; do
        local bc_file="${bc_dir}/${sample}_${region}_barcodes.txt"
        if [[ ! -s "${bc_file}" ]]; then
            echo "  [MISS] barcode 文件不存在或为空: ${bc_file}" >&2
            missing_bc=1
        fi
    done
    if [[ "${missing_bc}" -eq 1 ]]; then
        echo "[ERROR] barcode 文件不完整，请检查 Step1 输出。" >&2
        return 1
    fi

    mkdir -p "${pb_dir}" "${QSUB_DIR}" "${LOG_DIR}" "${META_DIR}"

    local regions_str
    regions_str="$(shell_join "${REGIONS[@]}")"

    cat > "${pbs_file}" <<PBSEOF
#!/bin/bash
#PBS -N step2_${sample}
#PBS -j oe
#PBS -o ${LOG_DIR}/step2_${sample}.pbs.log
#PBS -d ${QSUB_DIR}
#PBS -l nodes=1:ppn=${threads}
#PBS -l walltime=24:00:00
#PBS -q gold

set -euo pipefail

eval "\$(micromamba shell hook --shell bash)"
micromamba activate rmats_env

echo "[\$(date '+%F %T')] Step2 开始: ${sample}"
echo "  BAM:      ${bam_path}"
echo "  BC_DIR:   ${bc_dir}"
echo "  输出目录: ${pb_dir}"
echo ""

BAM="${bam_path}"
BC_DIR="${bc_dir}"
PB_DIR="${pb_dir}"
SAMPLE="${sample}"
THREADS="${threads}"
READ_COUNT_TSV="${read_count_tsv}"
REGIONS=(${regions_str})

mkdir -p "\${PB_DIR}" "$(dirname "${read_count_tsv}")"

echo -e "sample\tregion\tn_barcodes\tn_reads\tbam" > "\${READ_COUNT_TSV}"

if ! command -v samtools >/dev/null 2>&1; then
    echo "[ERROR] 找不到 samtools，请检查 rmats_env 环境。" >&2
    exit 1
fi

# ---- 建立原始 BAM 索引（若不存在）----
if [[ ! -f "\${BAM}.bai" ]]; then
    echo "  建立原始 BAM 索引..."
    samtools index -@ "\${THREADS}" "\${BAM}"
fi

filter_region() {
    local REGION="\$1"
    local BC_FILE="\${BC_DIR}/\${SAMPLE}_\${REGION}_barcodes.txt"
    local OUT_BAM="\${PB_DIR}/\${SAMPLE}_\${REGION}_pseudobulk.bam"

    if [[ ! -s "\${BC_FILE}" ]]; then
        echo "  [SKIP]  \${REGION}: barcode 文件不存在或为空 (\${BC_FILE})"
        return 0
    fi

    local N_BC
    N_BC=\$(grep -c . "\${BC_FILE}" || true)
    if [[ "\${N_BC}" -eq 0 ]]; then
        echo "  [SKIP]  \${REGION}: barcode 文件为空"
        return 0
    fi

    if [[ -f "\${OUT_BAM}" && -f "\${OUT_BAM}.bai" ]]; then
        local N_READS
        N_READS=\$(samtools view -c -@ "\${THREADS}" "\${OUT_BAM}")
        echo "  [EXIST] \${REGION}: \${N_BC} barcodes, \${N_READS} reads（已完成，跳过）"
        echo -e "\${SAMPLE}\t\${REGION}\t\${N_BC}\t\${N_READS}\t\${OUT_BAM}" >> "\${READ_COUNT_TSV}"
        return 0
    fi

    echo "  [RUN]   \${REGION}: \${N_BC} barcodes → \$(basename "\${OUT_BAM}")"

    samtools view \
        -@ "\${THREADS}" \
        -b \
        -D CB:"\${BC_FILE}" \
        "\${BAM}" \
        -o "\${OUT_BAM}"

    samtools index -@ "\${THREADS}" "\${OUT_BAM}"

    local N_READS
    N_READS=\$(samtools view -c -@ "\${THREADS}" "\${OUT_BAM}")
    echo "    → \${N_READS} reads"
    echo -e "\${SAMPLE}\t\${REGION}\t\${N_BC}\t\${N_READS}\t\${OUT_BAM}" >> "\${READ_COUNT_TSV}"
}

# ---- 处理各 Regions ----
for REGION in "\${REGIONS[@]}"; do
    filter_region "\${REGION}"
done

# ---- 处理 ALL ----
filter_region "ALL"

# ---- 输出汇总 ----
echo ""
echo "--- 输出汇总: ${sample} ---"
column -t -s $'\t' "\${READ_COUNT_TSV}" || cat "\${READ_COUNT_TSV}"

echo ""
echo "[\$(date '+%F %T')] Step2 完成: ${sample}"
echo "  输出目录: ${pb_dir}"
echo "  read count: ${read_count_tsv}"
PBSEOF

    echo "  ✓ PBS 文件已生成: ${pbs_file}"
    echo "  提交命令: qsub ${pbs_file}"
    echo ""
}

# ============================================================
# Parse arguments
# ============================================================
TARGET="all"
THREADS="10"
BAM_OVERRIDE=""

if [[ $# -eq 0 ]]; then
    TARGET="all"
    THREADS="10"
elif [[ $# -eq 1 ]]; then
    TARGET="$1"
    THREADS="10"
elif [[ $# -eq 2 ]]; then
    TARGET="$1"
    if is_integer "$2"; then
        THREADS="$2"
    else
        BAM_OVERRIDE="$2"
        THREADS="10"
    fi
elif [[ $# -eq 3 ]]; then
    TARGET="$1"
    BAM_OVERRIDE="$2"
    THREADS="$3"
else
    usage
    exit 1
fi

if ! is_integer "${THREADS}"; then
    echo "[ERROR] THREADS 必须是整数: ${THREADS}" >&2
    usage
    exit 1
fi

if [[ "${TARGET}" == "all" ]]; then
    TARGET_SAMPLES=("${SAMPLES[@]}")
else
    if ! sample_exists "${TARGET}"; then
        echo "[ERROR] 未知样本: ${TARGET}" >&2
        usage
        exit 1
    fi
    TARGET_SAMPLES=("${TARGET}")
fi

# 只有单样本模式才允许 BAM_OVERRIDE
if [[ -n "${BAM_OVERRIDE}" && "${#TARGET_SAMPLES[@]}" -ne 1 ]]; then
    echo "[ERROR] BAM_PATH 覆盖只允许用于单个样本。" >&2
    exit 1
fi

# ============================================================
# 生成 PBS
# ============================================================
PBS_FILES=()
FAILED=0

for sample in "${TARGET_SAMPLES[@]}"; do
    bam_path="${BAM_MAP[${sample}]:-}"
    if [[ -n "${BAM_OVERRIDE}" ]]; then
        bam_path="${BAM_OVERRIDE}"
    fi

    if [[ -z "${bam_path}" ]]; then
        echo "[ERROR] 没有找到 ${sample} 的 BAM 路径。" >&2
        FAILED=1
        continue
    fi

    if generate_one_pbs "${sample}" "${bam_path}" "${THREADS}"; then
        PBS_FILES+=("${QSUB_DIR}/step2_${sample}.pbs")
    else
        FAILED=1
    fi
done

if [[ "${#PBS_FILES[@]}" -gt 0 ]]; then
    echo "================================================================="
    echo "Step2 PBS 文件已生成，提交命令："
    echo "================================================================="
    for pbs in "${PBS_FILES[@]}"; do
        echo "  qsub ${pbs}"
    done
    echo ""
fi

if [[ "${FAILED}" -ne 0 ]]; then
    echo "[ERROR] 部分 Step2 PBS 生成失败，请查看上方错误。" >&2
    exit 1
fi
