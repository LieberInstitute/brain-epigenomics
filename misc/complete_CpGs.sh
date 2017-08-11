#!/bin/bash
#$ -cwd
#$ -N complete_CpGs
#$ -o ./complete_CpGs.$TASK_ID.txt
#$ -e ./complete_CpGs.$TASK_ID.txt
#$ -m e
#$ -l bluejay,mem_free=20G,h_vmem=20G,h_fsize=100G
#$ -t 2,15
echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${TASK_ID}"

FILELIST=/dcl01/lieber/ajaffe/Amanda/ATAC/FinalPeaks/scripts/finalizedList-pooled.txt
ID=$(awk "NR==$SGE_TASK_ID" $FILELIST)

CPGs=/dcl01/lieber/ajaffe/Amanda/ATAC/Coverage/coverage_at_CpGs/BSobj_bsseqSmooth_Neuron_minCov_3_resized_sorted.bed
sortedbed=/dcl01/lieber/ajaffe/Amanda/ATAC/BAMs/combined_set/${ID}.final.sorted.bed

module load bedtools/2.24.0

echo ${ID}

coverageBed -sorted -g ../../genome_files/hg19.genome.reordered2 -counts -a $CPGs -b $sortedbed > ${ID}_CpGs_coverageBed_versionByLeo.txt

echo "**** Job ends ****"
date
