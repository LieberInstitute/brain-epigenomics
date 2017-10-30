#!/bin/bash
#$ -cwd
#$ -N atacFullCov
#$ -o ./logs/atacFullCov.txt
#$ -e ./logs/atacFullCov.txt
#$ -m e
#$ -l mem_free=55G,h_vmem=55G,h_fsize=100G
#$ -pe local 5

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${TASK_ID}"

Rscript atacFullCov.R

echo "**** Job ends ****"
date
