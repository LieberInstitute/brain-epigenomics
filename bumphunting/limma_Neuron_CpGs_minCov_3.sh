#!/bin/bash
#$ -cwd
#$ -N limma_Neuron_CpGs_minCov_3
#$ -o ./limma_Neuron_CpGs_minCov_3_log.txt
#$ -e ./limma_Neuron_CpGs_minCov_3_log.txt
#$ -m e
#$ -l bluejay,mem_free=150G,h_vmem=150G,h_fsize=100G
echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${TASK_ID}"

module load conda_R/3.4.x
Rscript limma_Neuron_CpGs_minCov_3.R

echo "**** Job ends ****"
date
