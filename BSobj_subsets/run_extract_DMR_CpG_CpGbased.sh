#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=80G,h_vmem=80G,h_fsize=200G
#$ -N extract_dmr_CpG
#$ -o ./logs/extract_dmr_CpG.txt
#$ -e ./logs/extract_dmr_CpG.txt
#$ -m e

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${TASK_ID}"

module load conda_R/3.4.x
Rscript extract_DMR_CpG_CpGbased.R

echo "**** Job ends ****"
date
