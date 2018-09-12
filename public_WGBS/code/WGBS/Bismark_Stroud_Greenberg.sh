#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=10G,h_vmem=13G,h_fsize=500G
#$ -N Bismark_Stroud_Greenberg_PE
#$ -pe local 5
#$ -t 5-21
#$ -tc 10
#$ -m a


FILELIST=/dcl02/lieber/WGBS/PublicData/Stroud_Greenberg_Cell_2017/FASTQ/WGBS/SRR_Acc_List.txt
ID=$(awk "NR==$SGE_TASK_ID" $FILELIST)
read1=/dcl02/lieber/WGBS/PublicData/Stroud_Greenberg_Cell_2017/FASTQ/WGBS/${ID}_1.fastq.gz
read2=/dcl02/lieber/WGBS/PublicData/Stroud_Greenberg_Cell_2017/FASTQ/WGBS/${ID}_2.fastq.gz

module load bowtie2/2.2.5
module load samtools/1.1


/dcl01/lieber/WGBS/Software/Bismark_v0.19.0/bismark --multicore 4 -o /dcl02/lieber/WGBS/PublicData/Stroud_Greenberg_Cell_2017/BAMS/WGBS /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/GENCODE/GRCm38_mm10 -1 $read1 -2 $read2
