#!/bin/bash
#$ -cwd
#$ -l mem_free=350G,h_vmem=350G,h_fsize=300G
#$ -N create_matched_versions
#$ -o ./logs/create_matched_versions.txt
#$ -e ./logs/create_matched_versions.txt
#$ -m e

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${TASK_ID}"

Rscript create_matched_versions.R

echo "**** Job ends ****"
date
