#!/bin/bash
#$ -cwd
#$ -l mem_free=40G,h_vmem=40G,h_fsize=100G
#$ -N nonCG_highCov_array
#$ -o ./logs/nonCG_highCov.$TASK_ID.txt
#$ -e ./logs/nonCG_highCov.$TASK_ID.txt
#$ -m e
#$ -t 1-25

echo "**** Job starts ****"
date

Rscript filter_highCov_nonCG.R

echo "**** Job ends ****"
date
