#!/bin/bash
#$ -cwd
#$ -pe local 5
#$ -l mem_free=28G,h_vmem=30G,h_fsize=200G
#$ -N step7-Rcounts-polyA.dlpfc
#$ -o ./logs/Rcounts-polyA.txt
#$ -e ./logs/Rcounts-polyA.txt
#$ -hold_jid pipeline_setup,step4-featCounts-polyA.dlpfc,step6-txQuant-polyA.dlpfc
#$ -m a
echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "****"

Rscript /dcl01/lieber/ajaffe/lab/brain-epigenomics/brainseq_pipeline/polyA_unstranded/.step7-create_count_objects-human.R -o hg19 -m /dcl01/lieber/ajaffe/lab/brain-epigenomics/brainseq_pipeline/polyA_unstranded -e polyA -p dlpfc -l TRUE -c TRUE -t 5 -s FALSE

echo "**** Job ends ****"
date
