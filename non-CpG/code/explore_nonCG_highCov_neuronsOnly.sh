#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=150G,h_vmem=150G,h_fsize=100G
#$ -N explore_nonCG_highCov_neuronsOnly
#$ -o ./logs/explore_nonCG_highCov_neuronsOnly.txt
#$ -e ./logs/explore_nonCG_highCov_neuronsOnly.txt
#$ -m e
#$ -M amanda.joy.price@gmail.com
#$ -hold_jid nonCG_highCov

echo "**** Job starts ****"
date

Rscript explore_nonCG_highCov_neuronsOnly.R

echo "**** Job ends ****"
date
