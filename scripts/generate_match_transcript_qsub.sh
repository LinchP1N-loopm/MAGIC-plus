#!/bin/bash
set -euo pipefail

# 用法：bash generate_match_transcript_qsub.sh SAMPLE TSV ADATA OUTDIR [OBS_COLS] [QSUB_DIR] [INDEX_PKL]
SAMPLE=${1:?请提供样本名}
TSV=${2:?请提供重构TSV}
ADATA=${3:?请提供参考h5ad}
OUTDIR=${4:?请提供输出目录}
OBS_COLS=${5:-AS_cluster,brain_cluster,leiden}
QSUB_DIR=${6:-/histor/zhao/projects/lupeng/zhaolab/ZJJ/Reconstruct/qsub}
INDEX_PKL=${7:-/histor/zhao/projects/lupeng/zhaolab/ZJJ/Reconstruct/transcript_junction_index.pkl}

MATCH_SCRIPT=${QSUB_DIR}/match_transcript.py
QSUB_FILE=${QSUB_DIR}/${SAMPLE}_match.sh
OUT_H5AD=${OUTDIR}/${SAMPLE}_transcript_adata.h5ad
mkdir -p "${QSUB_DIR}" "${OUTDIR}"

cat > "${QSUB_FILE}" << EOF
#!/bin/bash
#PBS -N ${SAMPLE}_match
#PBS -j oe
#PBS -d ${OUTDIR}
#PBS -l nodes=1:ppn=1
#PBS -l mem=60gb
#PBS -l walltime=12:00:00
#PBS -q gold

set -euo pipefail
eval "\$(micromamba shell hook --shell bash)"
micromamba activate py310_st

python "${MATCH_SCRIPT}" \\
    --tsv "${TSV}" \\
    --index "${INDEX_PKL}" \\
    --adata "${ADATA}" \\
    --output "${OUT_H5AD}" \\
    --obs-cols "${OBS_COLS}" \\
    --min-spots 2
EOF

echo "已生成: ${QSUB_FILE}"
echo "输出矩阵: ${OUT_H5AD}"
echo "提交命令: qsub ${QSUB_FILE}"
