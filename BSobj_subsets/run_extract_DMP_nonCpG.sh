#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=80G,h_vmem=80G,h_fsize=200G
#$ -N extract_dmp_nonCpG
#$ -o ./logs/extract_dmp_nonCpG.txt
#$ -e ./logs/extract_dmp_nonCpG.txt
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
Rscript extract_DMP_nonCpG.R

echo "**** Job ends ****"
date
