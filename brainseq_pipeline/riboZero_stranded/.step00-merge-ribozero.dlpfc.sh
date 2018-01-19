#!/bin/bash
#$ -cwd
#$ -l mem_free=6G,h_vmem=10G,h_fsize=150G
#$ -N step00-merge-ribozero.dlpfc
#$ -pe local 8
#$ -o ./logs/merge-ribozero.txt
#$ -e ./logs/merge-ribozero.txt
#$ -hold_jid pipeline_setup
#$ -m a

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "****"

Rscript /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/sh/step00-merge.R -s /dcl01/lieber/ajaffe/lab/brain-epigenomics/brainseq_pipeline/riboZero_stranded/.samples_unmerged.manifest -o /dcl01/lieber/ajaffe/lab/brain-epigenomics/brainseq_pipeline/riboZero_stranded/merged_fastq -c 8

echo "**** Job ends ****"
date
