#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=40G,h_vmem=40G,h_fsize=100G
#$ -N create_bsobj_merge_and_cpgs
#$ -o ./logs/create_bsobj_merge_and_cpgs.txt
#$ -e ./logs/create_bsobj_merge_and_cpgs.txt
#$ -m e

echo "**** Job starts ****"
date

Rscript create_bsobj_merge_and_cpgs.R

echo "**** Job ends ****"
date
