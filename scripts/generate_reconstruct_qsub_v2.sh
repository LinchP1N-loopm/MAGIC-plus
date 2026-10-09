#!/bin/bash
# 使用方法：bash generate_reconstruct_qsub_v2.sh <样本名> <BAM路径> [se|pe]
# 示例（双端）：bash generate_reconstruct_qsub_v2.sh G21-4 /path/to/G21-4.bam pe
# 示例（单端）：bash generate_reconstruct_qsub_v2.sh ZJJ_reg1 /path/to/ZJJ_reg1.bam se
# 第三个参数不填默认双端（pe）

SAMPLE=$1
BAM=$2
MODE=${3:-pe}   # 默认双端

if [ -z "$SAMPLE" ] || [ -z "$BAM" ]; then
    echo "用法: bash generate_reconstruct_qsub.sh <样本名> <BAM路径> [se|pe]"
    exit 1
fi

if [ "$MODE" != "se" ] && [ "$MODE" != "pe" ]; then
    echo "错误: 第三个参数只能是 se 或 pe"
    exit 1
fi

# 单端加参数，双端不加
if [ "$MODE" == "se" ]; then
    SE_FLAG="--single_end"
    MODE_LABEL="单端(SE)"
else
    SE_FLAG=""
    MODE_LABEL="双端(PE)"
fi

SCRIPT_DIR=/histor/zhao/projects/lupeng/zhaolab/ZJJ/Reconstruct/qsub
OUTDIR=/histor/zhao/projects/lupeng/zhaolab/ZJJ/Reconstruct/${SAMPLE}
mkdir -p "${SCRIPT_DIR}" "${OUTDIR}"

QSUB_FILE=${SCRIPT_DIR}/${SAMPLE}_reconstruct.sh

cat > ${QSUB_FILE} << EOF
#!/bin/bash
#PBS -N ${SAMPLE}_reconstruct
#PBS -j oe
#PBS -d ${OUTDIR}
#PBS -l nodes=1:ppn=1
#PBS -l walltime=72:00:00
#PBS -q gold

eval "\$(micromamba shell hook --shell bash)"
micromamba activate rmats_long

SCRIPT_DIR=${SCRIPT_DIR}
BAM=${BAM}
OUTDIR=${OUTDIR}
mkdir -p \${OUTDIR}
OUT_PREFIX=\${OUTDIR}/${SAMPLE}

echo "[\$(date '+%Y-%m-%d %H:%M:%S')] Step1: umi_reconstructV5 start"
echo "模式: ${MODE_LABEL}"

python \${SCRIPT_DIR}/umi_reconstructV5.py \\
    -i  \${BAM} \\
    -o  \${OUT_PREFIX} \\
    --threads       1 \\
    --max_umi_reads 500 \\
    ${SE_FLAG}

echo "[\$(date '+%Y-%m-%d %H:%M:%S')] ALL DONE"
echo "结果目录: \${OUTDIR}"
EOF

echo "========================================"
echo "已生成: ${QSUB_FILE}"
echo "  样本名  : ${SAMPLE}"
echo "  BAM路径 : ${BAM}"
echo "  比对模式: ${MODE_LABEL}"
echo "  结果目录: ${OUTDIR}"
echo "提交命令: qsub ${QSUB_FILE}"
echo "========================================"
