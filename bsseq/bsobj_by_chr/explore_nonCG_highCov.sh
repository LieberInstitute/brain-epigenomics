#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=150G,h_vmem=150G,h_fsize=100G
#$ -N explore_nonCG_highCov
#$ -o ./logs/explore_nonCG_highCov.txt
#$ -e ./logs/explore_nonCG_highCov.txt
#$ -m e
#$ -hold_jid nonCG_highCov

echo "**** Job starts ****"
date

Rscript explore_nonCG_highCov.R

echo "**** Job ends ****"
date
