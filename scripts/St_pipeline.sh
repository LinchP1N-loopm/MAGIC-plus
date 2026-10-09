#!/bin/bash
# =============================================================================
# GGG-seq Pipeline Script (mawk streaming version)
# Replace: seqkit split2 + subseq loop + concat  (slow, disk-heavy)
# With:    zcat | mawk | pigz                     (fast, streaming, zero temp files)
#
# GGG-seq R1 Barcode 结构：
#   位置  1- 8  : Barcode 前段         ( 8 bp)
#   位置  9-32  : Linker-1             (24 bp, 跳过)
#   位置 33-44  : Barcode 后段B        (12 bp)
#   位置 45-52  : Barcode 后段A        ( 8 bp)
#   位置 53-59  : Linker-2             ( 7 bp, 跳过)
#   位置 60-150 : RNA 序列             (91 bp)
#
# 拼接顺序（保持与原 seqkit concat 一致）：
#   test1(1-8) + test2(45-52) + test3(33-44) + test4(60-150)
#   = 8 + 8 + 12 + 91 = 119 bp
#   前 28 bp = CB(16 bp) + UMI(12 bp)，后 91 bp = RNA
# =============================================================================

sample=$1

sampath=`dirname ${sample}`
ifile1=`basename ${sample}`
array2=(${ifile1//_R1/ })
sample=${array2[0]}

# =============================================================================
# >>>>>>>>>> 需要修改 <<<<<<<<<<
# =============================================================================
project_path=/histor/zhao/pangkun/Zhaolab/Lupeng/Zhujunjie/20260515

# Barcode-info
whitelist=/histor/zhao/pangkun/Zhaolab/Lupeng/Zhujunjie/20251203/data/Barcode/whitelist.txt
ID=/histor/zhao/pangkun/Zhaolab/Lupeng/Zhujunjie/20251203/data/Barcode/T3-ids-barcode.txt

# result_path
st_path=${project_path}/data/result_STARsolo/${sample}

# ref-info（选择 mm10 或 hg38，注释另一组）
MAP=/histor/zhao/pangkun/database/ref_mm10_M23_release98/ref_mm10_st_pipeline_150_index
ANN=/histor/zhao/pangkun/database/ref_mm10_M23_release98/mm10.gtf
BED=/histor/zhao/pangkun/database/ref_mm10_M23_release98/mm10.bed
# MAP=/histor/zhao/pangkun/database/ref_hg38_v32_release98/ref_hg38_v32_release98_st_pipeline_150_index
# ANN=/histor/zhao/pangkun/database/ref_hg38_v32_release98/hg38.gtf
# BED=/histor/zhao/pangkun/database/ref_hg38_v32_release98/hg38.bed
# MAP=/histor/zhao/pangkun/Zhaolab/Lupeng/database/ref_scrofa_release115/STAR150
# ANN=/histor/zhao/pangkun/Zhaolab/Lupeng/database/ref_scrofa_release115/Sus_scrofa.Sscrofa11.1.115.gtf
# BED=/histor/zhao/pangkun/Zhaolab/Lupeng/database/ref_scrofa_release115/Sus_scrofa.Sscrofa11.1.115.bed


# Tool paths
MAWK=/histor/public/software/mawk-1.3.4/bin/mawk
PIGZ=/histor/public/software/pigz-2.8/pigz
FASTP=/histor/public/software/fastp/fastp
STAR=/histor/zhao/pangkun/software/STAR/source/STAR

# Other-info
log_file=${st_path}/${sample}_st_log.txt
t_num=20
node=gold
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

source /histor/zhao/pangkun/.bashrc
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
echo 'Step1 ----- GGG-seq Barcode Reformat (mawk streaming) -----'
echo \$(date +%T)
# =============================================================================
#
# 原方案 (慢):
#   seqkit split2 → 循环10次 subseq(×4) + concat(×3) → cat → rm
#   共约 10×7 = 70 次读写磁盘，大量临时文件
#
# 新方案 (快):
#   zcat | mawk 流式处理 | pigz 多线程压缩，零中间文件，一次扫描完成
#
# 拼接逻辑（与原 seqkit 完全一致）：
#   substr(1,8)  +  substr(45,8)  +  substr(33,12)  +  substr(60,91)
# =============================================================================

zcat ${fastq1} | \\
${MAWK} '{
    if (NR % 4 == 1) {
        print
    }
    else if (NR % 4 == 2) {
        print substr(\$0, 1, 8) substr(\$0, 45, 8) substr(\$0, 33, 12) substr(\$0, 60, 91)
    }
    else if (NR % 4 == 3) {
        print
    }
    else if (NR % 4 == 0) {
        print substr(\$0, 1, 8) substr(\$0, 45, 8) substr(\$0, 33, 12) substr(\$0, 60, 91)
    }
}' | ${PIGZ} -p ${t_num} -f > ${st_path}/${sample}_reformat_R1.fastq.gz

# R2 无需改动，软链节省空间
ln -sf ${fastq2} ${st_path}/${sample}_reformat_R2.fastq.gz

# =============================================================================
echo 'Step2 ----- Cutadapt (TA 建库去除 TSO adapter) -----'
echo \$(date +%T)
# =============================================================================

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
    -l 30 -x -g \\
    -w ${t_num} --detect_adapter_for_pe \\
    -j ${st_path}/${sample}.barcode.fastp.json \\
    -h ${st_path}/${sample}.barcode.fastp.html

# 清理 reformat 中间文件
rm -f ${st_path}/${sample}_reformat_R1.fastq.gz \\
      ${st_path}/${sample}_reformat_R2.fastq.gz \\
      ${st_path}/${sample}_reformat_cutadapt_R2.fastq.gz

# =============================================================================
echo 'Step4 ----- Split barcode / RNA from cleaned R1 (mawk streaming) -----'
echo \$(date +%T)
# =============================================================================
#
# 原方案: seqkit subseq 两次（1:28 和 29:150），各自读一遍文件
# 新方案: mawk 两路输出，只读一遍，各自 pigz 压缩
#   Barcode+UMI : 位置 1-28  (CB 16bp + UMI 12bp)
#   RNA         : 位置 29-end
# =============================================================================

zcat ${st_path}/${sample}_R1.cleaned.fq.gz | \\
${MAWK} -v bc_out="${st_path}/${sample}_R1.barcode.fq" \\
        -v rna_out="${st_path}/${sample}_R1.rna.fq" '{
    if (NR % 4 == 1) {
        header = \$0
    }
    else if (NR % 4 == 2) {
        bc_seq  = substr(\$0, 1, 28)
        rna_seq = substr(\$0, 29)
    }
    else if (NR % 4 == 3) {
        plus = \$0
    }
    else {
        print header  > bc_out
        print bc_seq  > bc_out
        print plus    > bc_out
        print substr(\$0, 1, 28) > bc_out

        print header  > rna_out
        print rna_seq > rna_out
        print plus    > rna_out
        print substr(\$0, 29) > rna_out
    }
}'

${PIGZ} -p ${t_num} -f ${st_path}/${sample}_R1.barcode.fq
${PIGZ} -p ${t_num} -f ${st_path}/${sample}_R1.rna.fq
# 生成：${sample}_R1.barcode.fq.gz  和  ${sample}_R1.rna.fq.gz

# =============================================================================
# echo 'Step5 ----- STARsolo (se: R2 cDNA + reformat R1 barcode) -----'
# echo \$(date +%T)
# # =============================================================================

# mkdir -p ${st_path}/STARsolo_se
# ${STAR} --genomeDir ${MAP} \\
#     --outFileNamePrefix ${st_path}/STARsolo_se/${sample}_ \\
#     --readFilesCommand zcat \\
#     --readFilesIn ${st_path}/${sample}_R2.cleaned.fq.gz \\
#                   ${st_path}/${sample}_R1.cleaned.fq.gz \\
#     --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM sF \\
#     --outSAMtype BAM SortedByCoordinate \\
#     --limitBAMsortRAM 121539607552 \\
#     --soloType CB_UMI_Simple \\
#     --soloCBwhitelist ${whitelist} \\
#     --soloCBstart 1 \\
#     --soloCBlen 16 \\
#     --soloUMIstart 17 \\
#     --soloUMIlen 12 \\
#     --soloBarcodeReadLength 0 \\
#     --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \\
#     --soloFeatures Gene GeneFull SJ Velocyto \\
#     --soloMultiMappers EM \\
#     --soloUMIdedup 1MM_CR \\
#     --soloStrand Reverse \\
#     --soloUMIfiltering MultiGeneUMI_CR \\
#     --soloCellFilter EmptyDrops_CR \\
#     --soloCellReadStats Standard \\
#     --clipAdapterType CellRanger4 \\
#     --outFilterScoreMin 30 \\
#     --outReadsUnmapped Fastx \\
#     --runThreadN ${t_num}

# # # =============================================================================
# # echo 'Step6 ----- STARsolo (r1: R1 RNA + R1 barcode) -----'
# # echo \$(date +%T)
# # # =============================================================================

# mkdir -p ${st_path}/STARsolo_r1
# ${STAR} --genomeDir ${MAP} \\
#     --outFileNamePrefix ${st_path}/STARsolo_r1/${sample}_ \\
#     --readFilesCommand zcat \\
#     --readFilesIn ${st_path}/${sample}_R1.rna.fq.gz \\
#                   ${st_path}/${sample}_R1.barcode.fq.gz \\
#     --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM sF \\
#     --outSAMtype BAM SortedByCoordinate \\
#     --limitBAMsortRAM 121539607552 \\
#     --soloType CB_UMI_Simple \\
#     --soloCBwhitelist ${whitelist} \\
#     --soloCBstart 1 \\
#     --soloCBlen 16 \\
#     --soloUMIstart 17 \\
#     --soloUMIlen 12 \\
#     --soloBarcodeReadLength 0 \\
#     --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \\
#     --soloFeatures Gene GeneFull SJ Velocyto \\
#     --soloMultiMappers EM \\
#     --soloUMIdedup 1MM_CR \\
#     --soloUMIfiltering MultiGeneUMI_CR \\
#     --soloCellFilter EmptyDrops_CR \\
#     --soloCellReadStats Standard \\
#     --clipAdapterType CellRanger4 \\
#     --outFilterScoreMin 30 \\
#     --outReadsUnmapped Fastx \\
#     --runThreadN ${t_num}

# =============================================================================
echo 'Step7 ----- STARsolo (pe: paired-end R1+R2) -----'
echo \$(date +%T)
# =============================================================================

mkdir -p ${st_path}/STARsolo_pe
${STAR} --genomeDir ${MAP} \\
    --outFileNamePrefix ${st_path}/STARsolo_pe/${sample}_ \\
    --readFilesCommand zcat \\
    --readFilesIn ${st_path}/${sample}_R1.cleaned.fq.gz \\
                  ${st_path}/${sample}_R2.cleaned.fq.gz \\
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
    --soloBarcodeMate 1 \\
    --clip5pNbases 28 0 \\
    --soloUMIfiltering MultiGeneUMI_CR \\
    --soloCellFilter EmptyDrops_CR \\
    --soloCellReadStats Standard \\
    --outFilterScoreMin 30 \\
    --outReadsUnmapped Fastx \\
    --runThreadN ${t_num}

# =============================================================================
echo 'Step8 ----- Index BAM & compress mtx -----'
echo \$(date +%T)
# =============================================================================

samtools index -@ ${t_num} ${st_path}/STARsolo_se/${sample}_Aligned.sortedByCoord.out.bam
samtools index -@ ${t_num} ${st_path}/STARsolo_r1/${sample}_Aligned.sortedByCoord.out.bam
samtools index -@ ${t_num} ${st_path}/STARsolo_pe/${sample}_Aligned.sortedByCoord.out.bam

time ${PIGZ} -p ${t_num} ${st_path}/STARsol*/*_Solo.out/*/*/*mtx

# =============================================================================
echo 'Step9 ----- Cleanup Unmapped -----'
echo \$(date +%T)
# =============================================================================

rm -f ${st_path}/STARsolo_se/*_Unmapped.out.mate* \\
      ${st_path}/STARsolo_r1/*_Unmapped.out.mate* \\
      ${st_path}/STARsolo_pe/*_Unmapped.out.mate*

# =============================================================================
echo 'Step10 ----- QC stats (RSeQC) - STARsolo_se -----'
echo \$(date +%T)
# =============================================================================

conda deactivate

# # # # --- STARsolo_se ---
# cd ${st_path}/STARsolo_se
# samtools view -h ${sample}_Aligned.sortedByCoord.out.bam | \\
#     grep -v 'CB:Z:-' > ${sample}_Aligned_AB.sortedByCoord.out.sam
# samtools view -@ ${t_num} -b ${sample}_Aligned_AB.sortedByCoord.out.sam \\
#     > ${sample}_Aligned_AB.sortedByCoord.out.bam
# rm ${sample}_Aligned_AB.sortedByCoord.out.sam
# samtools index   -@ ${t_num} ${sample}_Aligned_AB.sortedByCoord.out.bam
# samtools flagstat -@ ${t_num} ${sample}_Aligned_AB.sortedByCoord.out.bam > ${sample}_flagstat.txt
# # read_distribution.py  -i ${sample}_Aligned_AB.sortedByCoord.out.bam -r ${BED}
# # geneBody_coverage.py  -r ${BED} -i ${sample}_Aligned_AB.sortedByCoord.out.bam -o ${sample}_RSEQC

# # # # --- STARsolo_r1 ---
# cd ${st_path}/STARsolo_r1
# samtools view -h ${sample}_Aligned.sortedByCoord.out.bam | \\
#     grep -v 'CB:Z:-' > ${sample}_Aligned_AB.sortedByCoord.out.sam
# samtools view -@ ${t_num} -b ${sample}_Aligned_AB.sortedByCoord.out.sam \\
#     > ${sample}_Aligned_AB.sortedByCoord.out.bam
# rm ${sample}_Aligned_AB.sortedByCoord.out.sam
# samtools index   -@ ${t_num} ${sample}_Aligned_AB.sortedByCoord.out.bam
# samtools flagstat -@ ${t_num} ${sample}_Aligned_AB.sortedByCoord.out.bam > ${sample}_flagstat.txt
# # read_distribution.py  -i ${sample}_Aligned_AB.sortedByCoord.out.bam -r ${BED}
# # geneBody_coverage.py  -r ${BED} -i ${sample}_Aligned_AB.sortedByCoord.out.bam -o ${sample}_RSEQC

# --- STARsolo_pe ---
cd ${st_path}/STARsolo_pe
samtools view -h ${sample}_Aligned.sortedByCoord.out.bam | \\
    grep -v 'CB:Z:-' > ${sample}_Aligned_AB.sortedByCoord.out.sam
samtools view -@ ${t_num} -b ${sample}_Aligned_AB.sortedByCoord.out.sam \\
    > ${sample}_Aligned_AB.sortedByCoord.out.bam
rm ${sample}_Aligned_AB.sortedByCoord.out.sam
samtools index   -@ ${t_num} ${sample}_Aligned_AB.sortedByCoord.out.bam
samtools flagstat -@ ${t_num} ${sample}_Aligned_AB.sortedByCoord.out.bam > ${sample}_flagstat.txt
read_distribution.py  -i ${sample}_Aligned_AB.sortedByCoord.out.bam -r ${BED}
geneBody_coverage.py  -r ${BED} -i ${sample}_Aligned_AB.sortedByCoord.out.bam -o ${sample}_RSEQC

echo 'End -----------'
echo \$(date +%T)

chmod -R 775 ${st_path}

#-------Ending-Log-----#
echo \"[\$(date)] \${PBS_JOBID}@\$(hostname) \$PBS_JOBNAME finished\"
end=\$(date +%s)
secs=\$((\$end-\$start))
printf 'Elapsed Time: %02dh:%02dm:%02ds\n' \$((secs/3600)) \$((secs%3600/60)) \$((secs%60))
#-------Ending-Log-----#

" > ${project_path}/qsub/qsub_solo/${sample}.St.sh

echo "#PBS -N ${sample}
#PBS -l nodes=1:ppn=${t_num}
#PBS -l walltime=72:00:00
#PBS -j oe
#PBS -q ${node}

set -e

cd ${project_path}/qsub/qsub_solo/

bash ${sample}.St.sh

" > ${project_path}/qsub/qsub_solo/${sample}.St.qsub
