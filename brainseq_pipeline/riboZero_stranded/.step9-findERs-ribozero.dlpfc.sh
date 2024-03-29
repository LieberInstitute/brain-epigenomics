#!/bin/bash
#$ -cwd
#$ -pe local 10
#$ -l mem_free=25G,h_vmem=25G,h_fsize=200G
#$ -N step9-findERs-ribozero.dlpfc
#$ -o ./logs/findERs-ribozero.txt
#$ -e ./logs/findERs-ribozero.txt
#$ -hold_jid pipeline_setup,step5b-meanCoverage-ribozero.dlpfc
#$ -m a
echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "****"

for meanFile in Coverage/mean*.bw
do
    echo "************************************"
    date
    echo "Initializing script for ${meanFile}"
    echo "************************************"
    Rscript /dcl01/lieber/ajaffe/lab/brain-epigenomics/brainseq_pipeline/riboZero_stranded/.step9-find_expressed_regions.R -m ${meanFile} -o /dcl01/lieber/ajaffe/lab/brain-epigenomics/brainseq_pipeline/riboZero_stranded/ERs -i /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/hg19.chrom.sizes.gencode -c 10
done

echo "**** Job ends ****"
date
