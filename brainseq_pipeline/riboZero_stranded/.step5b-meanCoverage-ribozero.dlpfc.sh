#!/bin/bash
#$ -cwd
#$ -l mem_free=200G,h_vmem=240G,h_fsize=100G
#$ -N step5b-meanCoverage-ribozero.dlpfc
#$ -o ./logs/meanCoverage-ribozero.txt
#$ -e ./logs/meanCoverage-ribozero.txt
#$ -hold_jid pipeline_setup,step5-coverage-ribozero.dlpfc
#$ -m a
echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "****"

if [ ! -f "inferred_strandness_pattern.txt" ]
then
    echo "Missing the file inferred_strandness_pattern.txt"
    exit 1
fi

## Load required software
module load wiggletools/default
module load ucsctools

## Read the strandness information
STRANDRULE=$(cat inferred_strandness_pattern.txt)

if [ ${STRANDRULE} == "none" ]
then
    echo "*****************************"
    date
    echo "Processing unstranded data"
    echo "*****************************"
    ## Locate normalized BigWig files and concatenate them in a space separated list
    BIGWIGS=$(cat /dcl01/lieber/ajaffe/lab/brain-epigenomics/brainseq_pipeline/riboZero_stranded/samples.manifest | awk '{print "/dcl01/lieber/ajaffe/lab/brain-epigenomics/brainseq_pipeline/riboZero_stranded/Coverage/"$NF".bw"}' | paste -sd " ")
    
    ## Create mean of normalized bigwigs
    wiggletools write /dcl01/lieber/ajaffe/lab/brain-epigenomics/brainseq_pipeline/riboZero_stranded/Coverage/mean.wig mean ${BIGWIGS}
    wigToBigWig /dcl01/lieber/ajaffe/lab/brain-epigenomics/brainseq_pipeline/riboZero_stranded/Coverage/mean.wig /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/hg19.chrom.sizes.gencode /dcl01/lieber/ajaffe/lab/brain-epigenomics/brainseq_pipeline/riboZero_stranded/Coverage/mean.bw
else
    for strand in Forward Reverse
    do
        echo "*****************************"
        date
        echo "Processing strand ${strand}"
        echo "*****************************"
        ## Locate normalized BigWig files and concatenate them in a space separated list
        BIGWIGS=$(cat /dcl01/lieber/ajaffe/lab/brain-epigenomics/brainseq_pipeline/riboZero_stranded/samples.manifest | awk -v strand="${strand}" '{print "/dcl01/lieber/ajaffe/lab/brain-epigenomics/brainseq_pipeline/riboZero_stranded/Coverage/"$NF"."strand".bw"}' | paste -sd " ")
    
        ## Create mean of normalized bigwigs
        wiggletools write /dcl01/lieber/ajaffe/lab/brain-epigenomics/brainseq_pipeline/riboZero_stranded/Coverage/mean.${strand}.wig mean ${BIGWIGS}
        wigToBigWig /dcl01/lieber/ajaffe/lab/brain-epigenomics/brainseq_pipeline/riboZero_stranded/Coverage/mean.${strand}.wig /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/hg19.chrom.sizes.gencode /dcl01/lieber/ajaffe/lab/brain-epigenomics/brainseq_pipeline/riboZero_stranded/Coverage/mean.${strand}.bw
    done
fi

## Remove temp files
rm /dcl01/lieber/ajaffe/lab/brain-epigenomics/brainseq_pipeline/riboZero_stranded/Coverage/mean*.wig

echo "**** Job ends ****"
date
