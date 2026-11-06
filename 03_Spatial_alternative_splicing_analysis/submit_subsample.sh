#!/bin/bash

# ============================================================
# 批量抽取 reads 脚本
# 样本：sample_reg1、sample_reg2、G21-4、G21-7
# 梯度：1M 3M 5M 10M 20M 30M 50M 100M 150M 200M 300M
# 每个样本提交一个PBS作业，11个梯度在作业内串行运行
# ============================================================

SCRIPT_DIR="/path/to/project/20260516_AS_saturation/qsub/qsub_fq"
OUT_BASE="/path/to/project/20260516_AS_saturation/data/rawdata"

mkdir -p ${SCRIPT_DIR}

# ---------- 样本信息定义 ----------
declare -A FASTQ1
declare -A FASTQ2

FASTQ1["sample_reg1"]="/path/to/project/Project/MAGIC_seq/Mouse_Adult_Brain_T3_70_50um/data/rawdata/sample_reg1/sample_reg1_R1.fastq.gz"
FASTQ2["sample_reg1"]="/path/to/project/Project/MAGIC_seq/Mouse_Adult_Brain_T3_70_50um/data/rawdata/sample_reg1/sample_reg1_R2.fastq.gz"

FASTQ1["sample_reg2"]="/path/to/project/Project/MAGIC_seq/Mouse_Adult_Brain_T3_70_50um/data/rawdata/sample_reg2/sample_reg2_R1.fastq.gz"
FASTQ2["sample_reg2"]="/path/to/project/Project/MAGIC_seq/Mouse_Adult_Brain_T3_70_50um/data/rawdata/sample_reg2/sample_reg2_R2.fastq.gz"

FASTQ1["G21-4"]="/path/to/project/20260210_2/data/rawdata/G21-4/G21-4_R1.fastq.gz"
FASTQ2["G21-4"]="/path/to/project/20260210_2/data/rawdata/G21-4/G21-4_R2.fastq.gz"

FASTQ1["G21-7"]="/path/to/project/20260328/data/rawdata/G21-7/G21-7_R1.fastq.gz"
FASTQ2["G21-7"]="/path/to/project/20260328/data/rawdata/G21-7/G21-7_R2.fastq.gz"

SAMPLES=("sample_reg1" "sample_reg2" "G21-4" "G21-7")

# ---------- 梯度定义 ----------
GRADIENTS=(1000000 3000000 5000000 10000000 20000000 30000000 50000000 100000000 150000000 200000000 300000000)
LABELS=("1M" "3M" "5M" "10M" "20M" "30M" "50M" "100M" "150M" "200M" "300M")

# ---------- 每个样本生成并提交一个PBS作业 ----------
for SAMPLE in "${SAMPLES[@]}"; do
    fq1=${FASTQ1[$SAMPLE]}
    fq2=${FASTQ2[$SAMPLE]}

    QSUB_FILE="${SCRIPT_DIR}/${SAMPLE}_subsample.pbs"

    # 写PBS头部
    cat > ${QSUB_FILE} << EOF
#!/bin/bash
#PBS -N ${SAMPLE}_subsample
#PBS -j oe
#PBS -d ${SCRIPT_DIR}
#PBS -l nodes=1:ppn=10
#PBS -l walltime=72:00:00
#PBS -q gold
source /path/to/project/.bashrc
source activate st_pipeline

ulimit -n 4096

echo "[START] Sample: ${SAMPLE}"

EOF

    # 逐梯度追加命令
    for i in "${!GRADIENTS[@]}"; do
        NUM_READ=${GRADIENTS[$i]}
        LABEL=${LABELS[$i]}

        OUT_DIR="${OUT_BASE}/${SAMPLE}-${LABEL}"
        OUT_R1="${OUT_DIR}/${SAMPLE}-${LABEL}_R1.fastq.gz"
        OUT_R2="${OUT_DIR}/${SAMPLE}-${LABEL}_R2.fastq.gz"

        cat >> ${QSUB_FILE} << EOF
echo "[\$(date '+%Y-%m-%d %H:%M:%S')] Running ${LABEL} (${NUM_READ} reads)..."
mkdir -p ${OUT_DIR}
seqkit head -n ${NUM_READ} -j 10 ${fq1} -o ${OUT_R1}
seqkit head -n ${NUM_READ} -j 10 ${fq2} -o ${OUT_R2}
echo "[\$(date '+%Y-%m-%d %H:%M:%S')] Done: ${LABEL}"

EOF
    done

    echo "echo '[DONE] All gradients finished for ${SAMPLE}'" >> ${QSUB_FILE}

    echo "Submitting: ${SAMPLE}"
    qsub ${QSUB_FILE}

done

echo ""
echo "全部作业已提交，共 ${#SAMPLES[@]} 个任务（每个样本一个PBS作业）"
echo "PBS脚本保存在: ${SCRIPT_DIR}"
echo "输出数据保存在: ${OUT_BASE}"
