#!/bin/bash

set -euo pipefail

# ============================================================
# 路径
# ============================================================
BAM_BASE="/histor/zhao/pangkun/Zhaolab/Lupeng/Zhujunjie/20260705_Fig3/rmats_works/pseudobulk_bams"
WORK_DIR="/histor/zhao/projects/lupeng/zhaolab/ZJJ/Fig4_sashimi"

BW_BASE="${WORK_DIR}/bigwig"
QSUB_DIR="${WORK_DIR}/qsub"
LOG_DIR="${WORK_DIR}/logs"

THREADS=10

mkdir -p "${BW_BASE}" "${QSUB_DIR}" "${LOG_DIR}"

# ============================================================
# 样本与脑区
# ============================================================
SAMPLES=(G21-4 G21-7 G21-1 G21-8 G21-2 G21-5)

REGIONS=(
    CA1
    CA3
    CTX_1
    CTX_2
    CTX_3
    DG
    FT
    HIP
    HY
    STR
    TH
    VS
)

# ============================================================
# 生成 PBS
# ============================================================
for sample in "${SAMPLES[@]}"; do

    PBS_FILE="${QSUB_DIR}/bam2bigwig_${sample}.pbs"

    cat > "${PBS_FILE}" <<EOF
#!/bin/bash
#PBS -N bw_${sample}
#PBS -j oe
#PBS -o ${LOG_DIR}/bam2bigwig_${sample}.log
#PBS -d ${QSUB_DIR}
#PBS -l nodes=1:ppn=${THREADS}
#PBS -l walltime=48:00:00
#PBS -q gold

eval "\$(micromamba shell hook --shell bash)"
micromamba activate biotools

set -euo pipefail

SAMPLE="${sample}"
BAM_DIR="${BAM_BASE}/${sample}"
OUT_DIR="${BW_BASE}/${sample}"
THREADS=${THREADS}

REGIONS=(
    CA1
    CA3
    CTX_1
    CTX_2
    CTX_3
    DG
    FT
    HIP
    HY
    STR
    TH
    VS
)

mkdir -p "\${OUT_DIR}"

echo "============================================================"
echo " BAM → BigWig"
echo " SAMPLE: \${SAMPLE}"
echo " BAM:    \${BAM_DIR}"
echo " OUT:    \${OUT_DIR}"
echo "============================================================"

for region in "\${REGIONS[@]}"; do

    bam="\${BAM_DIR}/\${SAMPLE}_\${region}_pseudobulk.bam"
    bw="\${OUT_DIR}/\${SAMPLE}_\${region}_RPKM.bw"

    echo ""
    echo ">>> \${SAMPLE} | \${region}"

    if [[ ! -f "\${bam}" ]]; then
        echo "  [MISS] BAM not found: \${bam}"
        continue
    fi

    if [[ ! -f "\${bam}.bai" ]]; then
        echo "  [INFO] indexing BAM ..."
        samtools index -@ "\${THREADS}" "\${bam}"
    fi

    if [[ -f "\${bw}" ]]; then
        echo "  [SKIP] already exists: \${bw}"
        continue
    fi

    echo "  BAM : \${bam}"
    echo "  OUT : \${bw}"

    bamCoverage \
        --bam "\${bam}" \
        --outFileName "\${bw}" \
        --outFileFormat bigwig \
        --binSize 1 \
        --normalizeUsing RPKM \
        --exactScaling \
        --minMappingQuality 10 \
        --numberOfProcessors "\${THREADS}"

    if [[ -f "\${bw}" ]]; then
        echo "  [OK] \${bw}"
    else
        echo "  [ERROR] output not created"
    fi

done

echo ""
echo "============================================================"
echo " Done: \${SAMPLE}"
echo "============================================================"
EOF

    echo "✓ ${PBS_FILE}"
done

echo ""
echo "PBS files generated:"
for sample in "${SAMPLES[@]}"; do
    echo "qsub ${QSUB_DIR}/bam2bigwig_${sample}.pbs"
done
