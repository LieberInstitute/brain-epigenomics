#!/bin/bash
#$ -cwd
#$ -l mem_free=130G,h_vmem=130G,h_fsize=100G
#$ -N nonCG_highCov
#$ -o ./logs/nonCG_highCov.txt
#$ -e ./logs/nonCG_highCov.txt
#$ -m e
#$ -hold_jid nonCG_highCov_array

echo "**** Job starts ****"
date

Rscript filter_highCov_nonCG.R

echo "**** Job ends ****"
date
