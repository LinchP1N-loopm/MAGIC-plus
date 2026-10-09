#!/bin/bash
#PBS -N G21-4-0515-pe_reconstruct
#PBS -j oe
#PBS -d /histor/zhao/projects/lupeng/zhaolab/ZJJ/Reconstruct/G21-4-0515-pe
#PBS -l nodes=1:ppn=1
#PBS -l walltime=72:00:00
#PBS -q gold

eval "$(micromamba shell hook --shell bash)"
micromamba activate rmats_long

SCRIPT_DIR=/histor/zhao/projects/lupeng/zhaolab/ZJJ/Reconstruct/qsub
BAM=/histor/zhao/pangkun/Zhaolab/Lupeng/Zhujunjie/20260515/data/raw_result_STARsolo/G21-4/STARsolo_pe/G21-4_Aligned_AB.sortedByCoord.out.bam
OUTDIR=/histor/zhao/projects/lupeng/zhaolab/ZJJ/Reconstruct/G21-4-0515-pe
mkdir -p ${OUTDIR}
OUT_PREFIX=${OUTDIR}/G21-4-0515-pe

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Step1: umi_reconstructV5 start"
echo "模式: 双端(PE)"

python ${SCRIPT_DIR}/umi_reconstructV5.py \
    -i  ${BAM} \
    -o  ${OUT_PREFIX} \
    --threads       1 \
    --max_umi_reads 500 \
    

echo "[$(date '+%Y-%m-%d %H:%M:%S')] ALL DONE"
echo "结果目录: ${OUTDIR}"
