#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=10G,h_vmem=13G,h_fsize=500G
#$ -N Bismark_Gabel_Greenberg_PE
#$ -pe local 5
#$ -t 1-4
#$ -tc 10
#$ -m a

FILELIST=/dcl02/lieber/WGBS/PublicData/Gabel_Greenberg_Nature_2015/FASTQ/WGBS/SRR_PE.txt
ID=$(awk "NR==$SGE_TASK_ID" $FILELIST)
read1=/dcl02/lieber/WGBS/PublicData/Gabel_Greenberg_Nature_2015/FASTQ/WGBS/${ID}_1.fastq.gz
read2=/dcl02/lieber/WGBS/PublicData/Gabel_Greenberg_Nature_2015/FASTQ/WGBS/${ID}_2.fastq.gz

module load bowtie2/2.2.5
module load samtools/1.1

cd /dcl02/lieber/WGBS/PublicData/Gabel_Greenberg_Nature_2015/BAMS/WGBS/paired

/dcl01/lieber/WGBS/Software/Bismark_v0.19.0/bismark --multicore 4 -o /dcl02/lieber/WGBS/PublicData/Gabel_Greenberg_Nature_2015/BAMS/WGBS/paired /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/GENCODE/GRCm38_mm10 -1 $read1 -2 $read2
