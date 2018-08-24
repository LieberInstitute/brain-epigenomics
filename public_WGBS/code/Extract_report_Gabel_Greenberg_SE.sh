#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=16G,h_vmem=20G,h_fsize=100G
#$ -N Extract_report_Gabel_Greenberg_SE
#$ -pe local 4
#$ -t 1-8
#$ -tc 10
#$ -m a

FILELIST=/dcl02/lieber/WGBS/PublicData/Gabel_Greenberg_Nature_2015/FASTQ/WGBS/SRR_SE.txt
ID=$(awk "NR==$SGE_TASK_ID" $FILELIST)

module load bowtie2/2.2.5
module load samtools/1.1

cd /dcl02/lieber/WGBS/PublicData/Gabel_Greenberg_Nature_2015/BAMS/WGBS/single/

/dcl01/lieber/WGBS/Software/Bismark_v0.19.0/bismark_methylation_extractor  --single-end \
	--cytosine_report --genome_folder /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/GENCODE/GRCm38_mm10 \
	--gzip --multicore 4 --bedGraph  \
	--output /dcl02/lieber/WGBS/PublicData/Gabel_Greenberg_Nature_2015/BAMS/WGBS/reports \
	--CX_context --split_by_chromosome ${ID}_bismark_bt2.bam
