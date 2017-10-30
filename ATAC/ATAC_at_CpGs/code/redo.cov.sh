# qsub -tc 245 -l mf=13G,h_vmem=13G,h_stack=256M -m e -M amanda.joy.price@gmail.com -t 1-11 redo.cov.sh

FILELIST=/dcl01/lieber/ajaffe/lab/brain-epigenomics/ATAC/ATAC_at_CpGs/redocov.txt
ID=$(awk "NR==$SGE_TASK_ID" $FILELIST)

CPGs=/dcl01/lieber/ajaffe/Amanda/ATAC/Coverage/coverage_at_CpGs/BSobj_bsseqSmooth_Neuron_minCov_3_resized_sorted.bed
sortedbed=/dcl01/lieber/ajaffe/Amanda/ATAC/BAMs/combined_set/${ID}.final.sorted.bed
col6=/dcl01/lieber/ajaffe/Amanda/ATAC/BAMs/combined_set/${ID}.final.sorted.6col.bed

cd /dcl01/lieber/ajaffe/Amanda/ATAC/Coverage/coverage_at_CpGs
module load bedtools/2.24.0

echo ${ID}

awk -F'\t' '{print NF}' $sortedbed | grep -v 6 | wc -l
awk 'NF==6{print}{}' $sortedbed > $col6

coverageBed -sorted -g ../../genome_files/hg19.genome.reordered2 -counts -a $CPGs -b $col6 > ${ID}_CpGs_coverageBed.txt
