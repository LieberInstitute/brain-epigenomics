#!/bin/bash
#$ -cwd
#$ -N step4-featCounts-ribozero.dlpfc_clean
#$ -o ./logs/featCounts-ribozero_clean.txt
#$ -e ./logs/featCounts-ribozero_clean.txt
#$ -hold_jid pipeline_setup,step4-featCounts-ribozero.dlpfc
#$ -m a
echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "****"

## Delete temporary files after they have been used
rm -rf /dcl01/lieber/ajaffe/lab/brain-epigenomics/brainseq_pipeline/riboZero_stranded/Counts/junction/tmpdir

echo "**** Job ends ****"
date
