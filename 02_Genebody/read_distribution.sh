#!/bin/bash
# Usage: bash read_distribution.sh <sample_name> <bam_path>
# Example: bash read_distribution.sh Sample1 /path/to/Sample1_Aligned.sortedByCoord.out.bam

# =============================================================================
# 参数接收
# =============================================================================
sample=$1
bam=$2

if [[ -z "${sample}" || -z "${bam}" ]]; then
    echo "Usage: bash $0 <sample_name> <bam_path>"
    exit 1
fi

# =============================================================================
# 路径配置
# =============================================================================
project_path=/path/to/project/20260511_readdistribution
qsub_dir=${project_path}/qsub
result_path=${project_path}/results/${sample}

# =============================================================================
# ref-info（选择 mm10 或 hg38，注释另一组）
# =============================================================================
BED=/path/to/project/database/ref_mm10_M23_release98/mm10.bed
# BED=/path/to/project/database/ref_hg38_v32_release98/hg38.bed
# BED=/path/to/project/database/ref_scrofa_release115/Sus_scrofa.Sscrofa11.1.115.bed

# =============================================================================
# 其他配置
# =============================================================================
t_num=10
node=batch

# =============================================================================
# 创建目录
# =============================================================================
mkdir -p ${qsub_dir}
mkdir -p ${result_path}

echo "Sample  : ${sample}"
echo "BAM     : ${bam}"
echo "BED     : ${BED}"
echo "Results : ${result_path}"
echo "Qsub    : ${qsub_dir}"

# =============================================================================
# 生成 .sh 脚本
# =============================================================================
cat > ${qsub_dir}/${sample}.RD.sh << EOF
#!/bin/bash
#-------Starting-Log-----#
set -e
echo "[\$(date)] \${PBS_JOBID}@\$(hostname) \${PBS_JOBNAME} started"
start=\$(date +%s)
#-------Starting-Log-----#

source /path/to/project/.bashrc
source activate st_pipeline

ulimit -n 4096

mkdir -p ${result_path}

cd ${result_path}

echo "[\$(date)] Running read_distribution.py for ${sample}"

read_distribution.py \\
    -i ${bam} \\
    -r ${BED} \\
    > ${result_path}/${sample}_read_distribution.txt 2>&1

echo "[\$(date)] Finished. Results saved to ${result_path}/${sample}_read_distribution.txt"

chmod -R 775 ${result_path}

#-------Ending-Log-----#
echo "[\$(date)] \${PBS_JOBID}@\$(hostname) \${PBS_JOBNAME} finished"
end=\$(date +%s)
secs=\$((\$end-\$start))
printf 'Elapsed Time: %02dh:%02dm:%02ds\n' \$((secs/3600)) \$((secs%3600/60)) \$((secs%60))
#-------Ending-Log-----#
EOF

# =============================================================================
# 生成 .qsub 文件
# =============================================================================
cat > ${qsub_dir}/${sample}.RD.qsub << EOF
#PBS -N ${sample}_RD
#PBS -l nodes=1:ppn=${t_num}
#PBS -l walltime=24:00:00
#PBS -j oe
#PBS -o ${result_path}/${sample}_RD.log
#PBS -q ${node}

set -e

bash ${qsub_dir}/${sample}.RD.sh

EOF

echo ""
echo "Done! To submit:"
echo "  qsub ${qsub_dir}/${sample}.RD.qsub"
