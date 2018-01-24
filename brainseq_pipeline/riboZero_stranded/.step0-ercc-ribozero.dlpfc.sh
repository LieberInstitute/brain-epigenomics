#!/bin/bash
#$ -cwd
#$ -l mem_free=6G,h_vmem=10G,h_fsize=100G
#$ -N step0-ercc-ribozero.dlpfc
#$ -pe local 8
#$ -o ./logs/ercc-ribozero.$TASK_ID.txt
#$ -e ./logs/ercc-ribozero.$TASK_ID.txt
#$ -t 1-31
#$ -tc 20
#$ -hold_jid pipeline_setup,step00-merge-ribozero.dlpfc
#$ -m a
echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"
echo "****"
echo "Sample id: $(cat /dcl01/lieber/ajaffe/lab/brain-epigenomics/brainseq_pipeline/riboZero_stranded/samples.manifest | awk '{print $NF}' | awk "NR==${SGE_TASK_ID}")"
echo "****"

## Locate file and ids
FILE1=$(awk 'BEGIN {FS="\t"} {print $1}' /dcl01/lieber/ajaffe/lab/brain-epigenomics/brainseq_pipeline/riboZero_stranded/samples.manifest | awk "NR==${SGE_TASK_ID}")
if [ TRUE == "TRUE" ] 
then
    FILE2=$(awk 'BEGIN {FS="\t"} {print $3}' /dcl01/lieber/ajaffe/lab/brain-epigenomics/brainseq_pipeline/riboZero_stranded/samples.manifest | awk "NR==${SGE_TASK_ID}")
fi
ID=$(cat /dcl01/lieber/ajaffe/lab/brain-epigenomics/brainseq_pipeline/riboZero_stranded/samples.manifest | awk '{print $NF}' | awk "NR==${SGE_TASK_ID}")
mkdir -p /dcl01/lieber/ajaffe/lab/brain-epigenomics/brainseq_pipeline/riboZero_stranded/Ercc/${ID}

if [ TRUE == "TRUE" ]
then 
    /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Software/kallisto quant     -i /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/ERCC/ERCC92.idx     -o /dcl01/lieber/ajaffe/lab/brain-epigenomics/brainseq_pipeline/riboZero_stranded/Ercc/${ID} -t 8 --rf-stranded     ${FILE1} ${FILE2}
else
    /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Software/kallisto quant     -i /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/ERCC/ERCC92.idx     -o /dcl01/lieber/ajaffe/lab/brain-epigenomics/brainseq_pipeline/riboZero_stranded/Ercc/${ID} -t 8 --single --rf-stranded ${FILE1}
fi

echo "**** Job ends ****"
date
