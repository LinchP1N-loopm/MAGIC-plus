#!/bin/bash
#-------Starting-Log-----#
set -e
echo "[$(date)] ${PBS_JOBID}@$(hostname) ${PBS_JOBNAME} started"
start=$(date +%s)
#-------Starting-Log-----#

source /histor/zhao/pangkun/.bashrc
source activate st_pipeline

ulimit -n 4096

mkdir -p /histor/zhao/pangkun/Zhaolab/Lupeng/Zhujunjie/20260511_readdistribution/results/MAGIC-seq1

cd /histor/zhao/pangkun/Zhaolab/Lupeng/Zhujunjie/20260511_readdistribution/results/MAGIC-seq1

echo "[$(date)] Running read_distribution.py for MAGIC-seq1"

read_distribution.py \
    -i /histor/zhao/pangkun/Zhaolab/Lupeng/Zhujunjie/20260307_testdata/data/result_STARsolo/ZJJ_reg1/STARsolo_se/ZJJ_reg1_Aligned_AB.sortedByCoord.out.bam \
    -r /histor/zhao/pangkun/database/ref_mm10_M23_release98/mm10.bed \
    > /histor/zhao/pangkun/Zhaolab/Lupeng/Zhujunjie/20260511_readdistribution/results/MAGIC-seq1/MAGIC-seq1_read_distribution.txt 2>&1

echo "[$(date)] Finished. Results saved to /histor/zhao/pangkun/Zhaolab/Lupeng/Zhujunjie/20260511_readdistribution/results/MAGIC-seq1/MAGIC-seq1_read_distribution.txt"

chmod -R 775 /histor/zhao/pangkun/Zhaolab/Lupeng/Zhujunjie/20260511_readdistribution/results/MAGIC-seq1

#-------Ending-Log-----#
echo "[$(date)] ${PBS_JOBID}@$(hostname) ${PBS_JOBNAME} finished"
end=$(date +%s)
secs=$(($end-$start))
printf 'Elapsed Time: %02dh:%02dm:%02ds\n' $((secs/3600)) $((secs%3600/60)) $((secs%60))
#-------Ending-Log-----#
