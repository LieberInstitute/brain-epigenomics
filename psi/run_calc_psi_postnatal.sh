#!/bin/bash
#$ -cwd
#$ -l mem_free=21G,h_vmem=21G,h_fsize=200G
#$ -N calc_psi_postnatal
#$ -pe local 10
#$ -o ./logs/calc_psi_postnatal.txt
#$ -e ./logs/calc_psi_postnatal.txt
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
Rscript calc_psi_postnatal.R

echo "**** Job ends ****"
date
