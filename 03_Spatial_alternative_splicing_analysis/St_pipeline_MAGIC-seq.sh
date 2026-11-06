#!/bin/bash
# =============================================================================
# MAGIC-seq Pipeline Script (mawk streaming version)
# Replace: seqkit split2 + subseq + concat (slow, disk-heavy)
# With:    zcat | mawk | pigz  (fast, streaming, no temp files)
# =============================================================================

sample=$1

sampath=`dirname ${sample}`
ifile1=`basename ${sample}`
array2=(${ifile1//_R1/ })
sample=${array2[0]}

# =============================================================================
# >>>>>>>>>> 需要修改 <<<<<<<<<<
# =============================================================================
project_path=/path/to/project/20260516_AS_saturation
# Barcode-info
whitelist=/path/to/project/20260307_testdata/data/Barcode/Barcode-T3-70_MAGIC-seq/whitelist.txt
ID=/path/to/project/20260307_testdata/data/Barcode/Barcode-T3-70_MAGIC-seq/T3-ids-barcode.txt
# result_path
st_path=${project_path}/data/result_STARsolo/${sample}

# ref-info (选择 mm10 或 hg38，注释另一组)
MAP=/path/to/project/database/ref_mm10_M23_release98/ref_mm10_st_pipeline_150_index
ANN=/path/to/project/database/ref_mm10_M23_release98/mm10.gtf
BED=/path/to/project/database/ref_mm10_M23_release98/mm10.bed
# MAP=/path/to/project/database/ref_hg38_v32_release98/ref_hg38_v32_release98_st_pipeline_150_index
# ANN=/path/to/project/database/ref_hg38_v32_release98/hg38.gtf
# BED=/path/to/project/database/ref_hg38_v32_release98/hg38.bed

# Tool paths
MAWK=/path/to/software/mawk-1.3.4/bin/mawk
PIGZ=/path/to/software/pigz-2.8/pigz
FASTP=/path/to/software/fastp/fastp
STAR=/path/to/project/software/STAR/source/STAR

# Other-info
log_file=${st_path}/${sample}_st_log.txt
t_num=10
node=batch
# =============================================================================

echo $ifile1
echo $sample

fastq1=${sampath}/${sample}_R1.fastq.gz
fastq2=${sampath}/${sample}_R2.fastq.gz

echo $fastq1

echo "
#-------Starting-Log-----#
set -e
echo \"[\$(date)] ${PBS_JOBID}@\$(hostname) \${PBS_JOBNAME} started\"
start=\$(date +%s)
#-------Starting-Log-----#

source /path/to/project/.bashrc
source activate st_pipeline

ulimit -n 4096

mkdir -p ${st_path}

# =============================================================================
echo 'Step0 ----- QC reports (raw fastq) -----'
echo \$(date +%T)
# =============================================================================

${FASTP} \\
    -i ${fastq1} -I ${fastq2} \\
    -o ${st_path}/${sample}_R1.cleaned.fq.gz \\
    -O ${st_path}/${sample}_R2.cleaned.fq.gz \\
    -Q -L -G -A \\
    -m --merged_out ${st_path}/${sample}_se.cleaned.fq.gz \\
    -w ${t_num} --detect_adapter_for_pe \\
    -j ${st_path}/${sample}.raw.fastp.json \\
    -h ${st_path}/${sample}.raw.fastp.html
rm ${st_path}/*fq.gz

# =============================================================================
echo 'Step1 ----- MAGIC-seq Barcode Reformat (mawk streaming) -----'
echo \$(date +%T)
# =============================================================================
#
# MAGIC-seq R1 barcode 结构：
#   位置  1- 8 : Barcode 前段 (8 bp)
#   位置  9-26 : Linker 序列  (18 bp，跳过)
#   位置 27-46 : Barcode 后段 + UMI (20 bp)
#   拼接后 R1  : 28 bp  (对应 STARsolo --soloCBlen 16 + --soloUMIlen 12)
#
# 原方案 (慢): seqkit split2 → loop subseq → concat → cat → rm
# 新方案 (快): zcat | mawk 流式处理 | pigz 多线程压缩，零中间文件
# =============================================================================

zcat ${fastq1} | \\
${MAWK} '{
    if (NR % 4 == 1) {
        print
    }
    else if (NR % 4 == 2) {
        print substr(\$0, 1, 8) substr(\$0, 27, 20)
    }
    else if (NR % 4 == 3) {
        print
    }
    else if (NR % 4 == 0) {
        print substr(\$0, 1, 8) substr(\$0, 27, 20)
    }
}' | ${PIGZ} -p ${t_num} -f > ${st_path}/${sample}_reformat_R1.fastq.gz

# R2 无需改动，直接软链节省磁盘空间
ln -sf ${fastq2} ${st_path}/${sample}_reformat_R2.fastq.gz

# =============================================================================
echo 'Step2 ----- Cutadapt (TA 建库) -----'
echo \$(date +%T)
# =============================================================================

# TA 建库去除 TSO adapter
cutadapt \\
    -g AAGCAGTGGTATCAACGCAGAGTGAATGGG \\
    -e 0.01 \\
    -j ${t_num} \\
    -o ${st_path}/${sample}_reformat_cutadapt_R2.fastq.gz \\
    ${st_path}/${sample}_reformat_R2.fastq.gz

# =============================================================================
echo 'Step3 ----- fastp QC (barcode reads) -----'
echo \$(date +%T)
# =============================================================================

${FASTP} \\
    -i ${st_path}/${sample}_reformat_R1.fastq.gz \\
    -I ${st_path}/${sample}_reformat_cutadapt_R2.fastq.gz \\
    -o ${st_path}/${sample}_R1.cleaned.fq.gz \\
    -O ${st_path}/${sample}_R2.cleaned.fq.gz \\
    -l 28 -x -g \\
    -w ${t_num} --detect_adapter_for_pe \\
    -j ${st_path}/${sample}.barcode.fastp.json \\
    -h ${st_path}/${sample}.barcode.fastp.html

# 清理 reformat 中间文件
rm -f ${st_path}/${sample}_reformat_R1.fastq.gz \\
      ${st_path}/${sample}_reformat_R2.fastq.gz \\
      ${st_path}/${sample}_reformat_cutadapt_R2.fastq.gz

# =============================================================================
echo 'Step4 ----- STARsolo -----'
echo \$(date +%T)
# =============================================================================

${STAR} --genomeDir ${MAP} \\
    --outFileNamePrefix ${st_path}/STARsolo_se/${sample}_ \\
    --readFilesCommand zcat \\
    --readFilesIn ${st_path}/${sample}_R2.cleaned.fq.gz \\
                  ${st_path}/${sample}_R1.cleaned.fq.gz \\
    --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM sF \\
    --outSAMtype BAM SortedByCoordinate \\
    --limitBAMsortRAM 121539607552 \\
    --soloType CB_UMI_Simple \\
    --soloCBwhitelist ${whitelist} \\
    --soloCBstart 1 \\
    --soloCBlen 16 \\
    --soloUMIstart 17 \\
    --soloUMIlen 12 \\
    --soloBarcodeReadLength 0 \\
    --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \\
    --soloFeatures Gene GeneFull SJ Velocyto \\
    --soloMultiMappers EM \\
    --soloUMIdedup 1MM_CR \\
    --soloUMIfiltering MultiGeneUMI_CR \\
    --soloCellFilter EmptyDrops_CR \\
    --soloCellReadStats Standard \\
    --clipAdapterType CellRanger4 \\
    --outFilterScoreMin 30 \\
    --outReadsUnmapped Fastx \\
    --runThreadN ${t_num}
    #--soloStrand Reverse

samtools index -@ ${t_num} ${st_path}/STARsolo_se/${sample}_Aligned.sortedByCoord.out.bam

time ${PIGZ} -p ${t_num} ${st_path}/STARsol*/*_Solo.out/*/*/*mtx

# =============================================================================
echo 'Step5 ----- Cleanup & rename -----'
echo \$(date +%T)
# =============================================================================

rm -f ${st_path}/STARsolo_se/*_Unmapped.out.mate*

# =============================================================================
echo 'Step6 ----- QC stats (RSeQC) -----'
echo \$(date +%T)
# =============================================================================

conda deactivate

cd ${st_path}/STARsolo_se

samtools view -h ${sample}_Aligned.sortedByCoord.out.bam | \\
    grep -v 'CB:Z:-' > ${sample}_Aligned_AB.sortedByCoord.out.sam
samtools view -@ ${t_num} -b ${sample}_Aligned_AB.sortedByCoord.out.sam \\
    > ${sample}_Aligned_AB.sortedByCoord.out.bam
rm ${sample}_Aligned_AB.sortedByCoord.out.sam
samtools index -@ ${t_num} ${sample}_Aligned_AB.sortedByCoord.out.bam
samtools flagstat -@ ${t_num} ${sample}_Aligned_AB.sortedByCoord.out.bam \\
    > ${sample}_flagstat.txt

# read_distribution.py  -i ${sample}_Aligned_AB.sortedByCoord.out.bam -r ${BED}
# geneBody_coverage.py  -r ${BED} -i ${sample}_Aligned_AB.sortedByCoord.out.bam \\
#     -o ${sample}_RSEQC_

echo 'End -----------'
echo \$(date +%T)

chmod -R 775 ${st_path}

#-------Ending-Log-----#
echo \"[\$(date)] \${PBS_JOBID}@\$(hostname) \$PBS_JOBNAME finished\"
end=\$(date +%s)
secs=\$((\$end-\$start))
printf 'Elapsed Time: %02dh:%02dm:%02ds\n' \$((secs/3600)) \$((secs%3600/60)) \$((secs%60))
#-------Ending-Log-----#

" > ${project_path}/qsub/qsub_solo/${sample}.new_St.sh

echo "#PBS -N ${sample}
#PBS -l nodes=1:ppn=${t_num}
#PBS -l walltime=72:00:00
#PBS -j oe
#PBS -q ${node}

set -e

cd ${project_path}/qsub/qsub_solo/

bash ${sample}.new_St.sh

" > ${project_path}/qsub/qsub_solo/${sample}.new_St.qsub
