#!/bin/bash
#$ -cwd
#$ -l mem_free=10G,h_vmem=14G,h_fsize=100G
#$ -N step1-fastqc-Stroud_Greenberg_ChIPseq.Aug212018
#$ -o ./logs/fastqc-Stroud_Greenberg_ChIPseq.$TASK_ID.txt
#$ -e ./logs/fastqc-Stroud_Greenberg_ChIPseq.$TASK_ID.txt
#$ -t 1-0
#$ -tc 100
#$ -hold_jid pipeline_setup,step00-merge-Stroud_Greenberg_ChIPseq.Aug212018
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
echo "Sample id: $(cat /dcl01/lieber/ajaffe/lab/brain-epigenomics/public_WGBS/code/samples.manifest | awk '{print $NF}' | awk "NR==${SGE_TASK_ID}")"
echo "****"

## Locate file and ids
FILE1=$(awk 'BEGIN {FS="\t"} {print $1}' /dcl01/lieber/ajaffe/lab/brain-epigenomics/public_WGBS/code/samples.manifest | awk "NR==${SGE_TASK_ID}")
if [ FALSE == "TRUE" ] 
then
    FILE2=$(awk 'BEGIN {FS="\t"} {print $3}' /dcl01/lieber/ajaffe/lab/brain-epigenomics/public_WGBS/code/samples.manifest | awk "NR==${SGE_TASK_ID}")
fi
ID=$(cat /dcl01/lieber/ajaffe/lab/brain-epigenomics/public_WGBS/code/samples.manifest | awk '{print $NF}' | awk "NR==${SGE_TASK_ID}")

mkdir -p /dcl01/lieber/ajaffe/lab/brain-epigenomics/public_WGBS/code/FastQC/Untrimmed/${ID}

if [ FALSE == "TRUE" ]
then 
    /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Software/FastQC_v0.11.5/fastqc ${FILE1} ${FILE2} --outdir=/dcl01/lieber/ajaffe/lab/brain-epigenomics/public_WGBS/code/FastQC/Untrimmed/${ID} --extract
else
    /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Software/FastQC_v0.11.5/fastqc ${FILE1} --outdir=/dcl01/lieber/ajaffe/lab/brain-epigenomics/public_WGBS/code/FastQC/Untrimmed/${ID} --extract
fi

echo "**** Job ends ****"
date
