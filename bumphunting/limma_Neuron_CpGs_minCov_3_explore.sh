#!/bin/bash
#$ -cwd
#$ -N limma_Neuron_CpGs_minCov_3_explore
#$ -o ./limma_Neuron_CpGs_minCov_3_explore_log.txt
#$ -e ./limma_Neuron_CpGs_minCov_3_explore_log.txt
#$ -m e
#$ -l mem_free=250G,h_vmem=250G,h_fsize=100G
echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${TASK_ID}"

module load conda_R/3.4.x
Rscript limma_Neuron_CpGs_minCov_3_explore.R

echo "**** Job ends ****"
date
