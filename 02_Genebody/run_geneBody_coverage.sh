#!/bin/bash
# =============================================================================
# geneBody_coverage - 通用化改造版本
# 使用方式: bash run_geneBody_coverage.sh bam文件路径 样本名称
# 脚本存放路径: /path/to/project/20260328_genebody/supple_genebody/qsub
# 输出根路径: /path/to/project/20260328_genebody/supple_genebody/results
# =============================================================================

# -------------------------- 1. 参数检查与定义 --------------------------
# 检查输入参数数量
if [ $# -ne 2 ]; then
    echo "错误：参数数量不正确！"
    echo "使用方式：bash $0 <bam文件路径> <样本名称>"
    exit 1
fi

# 固定参数定义
BED=/path/to/project/database/ref_mm10_M23_release98/mm10.bed
#BED=/path/to/project/database/ref_hg38_v32_release98/hg38.bed
ROOT_OUT_DIR=/path/to/project/20260328_genebody/supple_genebody/results
t_num=10
node=gold

# 输入参数赋值
BAM_FILE=$1          # 传入的BAM文件路径
SAMPLE_NAME=$2       # 传入的样本名称

# 衍生路径定义
SAMPLE_OUT_DIR=${ROOT_OUT_DIR}/${SAMPLE_NAME}  # 样本专属输出目录
QSUB_DIR=$(cd $(dirname $0) && pwd)             # 脚本所在的qsub目录（绝对路径）
QSUB_FILE=${QSUB_DIR}/${SAMPLE_NAME}.qsub       # 生成的qsub脚本文件

# -------------------------- 2. 目录创建 --------------------------
mkdir -p ${SAMPLE_OUT_DIR} ${QSUB_DIR}

# -------------------------- 3. 生成qsub脚本 --------------------------
cat > ${QSUB_FILE} << EOF
#PBS -N ${SAMPLE_NAME}_genebody
#PBS -l nodes=1:ppn=${t_num}
#PBS -l walltime=24:00:00
#PBS -j oe
#PBS -q ${node}

set -e

echo "[\$(date)] \${PBS_JOBID}@\$(hostname) \${PBS_JOBNAME} started"
start=\$(date +%s)

# 激活环境
eval "\$(micromamba shell hook --shell bash)"
micromamba activate biotools

# 调整文件描述符限制
ulimit -n 4096

# 确保样本输出目录存在（防止意外）
mkdir -p ${SAMPLE_OUT_DIR}

# 执行geneBody_coverage分析
echo '----- geneBody_coverage: ${SAMPLE_NAME} -----'
echo \$(date +%T)

geneBody_coverage.py \\
    -r ${BED} \\
    -i ${BAM_FILE} \\
    -f pdf \\
    -o ${SAMPLE_OUT_DIR}/${SAMPLE_NAME}_RSEQC

echo 'End -----------'
echo \$(date +%T)

# 设置目录权限
chmod -R 775 ${SAMPLE_OUT_DIR}

# 计算运行时间
echo "[\$(date)] \${PBS_JOBID}@\$(hostname) \${PBS_JOBNAME} finished"
end=\$(date +%s)
secs=\$((\$end-\$start))
printf 'Elapsed Time: %02dh:%02dm:%02ds\n' \$((secs/3600)) \$((secs%3600/60)) \$((secs%60))
EOF

# -------------------------- 4. 提交qsub任务 --------------------------
echo "Qsub脚本已生成：${QSUB_FILE}"
qsub ${QSUB_FILE}
echo "任务已提交，样本${SAMPLE_NAME}的输出结果将存储在：${SAMPLE_OUT_DIR}"