#!/bin/bash
#PBS -N ZZJ_Reg1_10_500_1500
#PBS -l nodes=1:ppn=10
#PBS -l walltime=72:00:00
#PBS -j oe
#PBS -q gold
#PBS -o /histor/zhao/pangkun/Zhaolab/Lupeng/Zhujunjie/20260328_genebody/qsub/ZZJ_Reg1_10_500_1500.log

# 加载环境
source /histor/zhao/pangkun/.bashrc
source activate st_pipeline
ulimit -n 4096

echo "开始：ZZJ_Reg1_10 | 500_1500 | $(date)"
geneBody_coverage.py \
    -r /histor/zhao/pangkun/database/ref_mm10_M23_release98/mm10_500_1500.bed \
    -i /histor/zhao/pangkun/Zhaolab/Lupeng/Zhujunjie/20260307_testdata/data/result_STARsolo/ZZJ_Reg1_10/STARsolo_se/ZZJ_Reg1_10_Aligned_AB.sortedByCoord.out.bam \
    -o /histor/zhao/pangkun/Zhaolab/Lupeng/Zhujunjie/20260328_genebody/result/ZZJ_Reg1_10_RSeQC/500_1500/ZZJ_Reg1_10_500_1500_genebody
echo "完成：ZZJ_Reg1_10 | 500_1500 | $(date)"
