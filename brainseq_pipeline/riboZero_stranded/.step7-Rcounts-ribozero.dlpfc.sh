#!/bin/bash
#$ -cwd
#$ -pe local 5
#$ -l mem_free=28G,h_vmem=30G,h_fsize=200G
#$ -N step7-Rcounts-ribozero.dlpfc
#$ -o ./logs/Rcounts-ribozero.txt
#$ -e ./logs/Rcounts-ribozero.txt
#$ -hold_jid pipeline_setup,step4-featCounts-ribozero.dlpfc,step6-txQuant-ribozero.dlpfc
#$ -m a
echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "****"

Rscript /dcl01/lieber/ajaffe/lab/brain-epigenomics/brainseq_pipeline/riboZero_stranded/.step7-create_count_objects-human.R -o hg19 -m /dcl01/lieber/ajaffe/lab/brain-epigenomics/brainseq_pipeline/riboZero_stranded -e ribozero -p dlpfc -l TRUE -c TRUE -t 5 -s reverse

echo "**** Job ends ****"
date
