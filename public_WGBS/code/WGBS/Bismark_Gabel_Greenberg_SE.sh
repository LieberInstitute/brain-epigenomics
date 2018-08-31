#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=15G,h_vmem=15G,h_fsize=500G
#$ -N Bismark_Gabel_Greenberg_SE
#$ -pe local 5
#$ -t 1-8
#$ -tc 10
#$ -m a

FILELIST=/dcl02/lieber/WGBS/PublicData/Gabel_Greenberg_Nature_2015/FASTQ/WGBS/SRR_SE.txt
ID=$(awk "NR==$SGE_TASK_ID" $FILELIST)

module load bowtie2/2.2.5
module load samtools/1.1

cd /dcl02/lieber/WGBS/PublicData/Gabel_Greenberg_Nature_2015/BAMS/WGBS/single/

/dcl01/lieber/WGBS/Software/Bismark_v0.19.0/bismark --multicore 4 -o /dcl02/lieber/WGBS/PublicData/Gabel_Greenberg_Nature_2015/BAMS/WGBS/single /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/GENCODE/GRCm38_mm10 /dcl02/lieber/WGBS/PublicData/Gabel_Greenberg_Nature_2015/FASTQ/WGBS/${ID}.fastq.gz
