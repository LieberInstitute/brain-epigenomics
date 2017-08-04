# qsub -tc 245 -l mf=10G,h_vmem=10G,h_fsize=15G,h_stack=256M -m e -M amanda.joy.price@gmail.com -t 1-34 calculate.ATAC.coverage.at.CpGs.sh

FILELIST=/dcl01/lieber/ajaffe/Amanda/ATAC/FinalPeaks/scripts/finalizedList-pooled.txt
ID=$(awk "NR==$SGE_TASK_ID" $FILELIST)

CPGs=/dcl01/lieber/ajaffe/Amanda/ATAC/Coverage/coverage_at_CpGs/BSobj_bsseqSmooth_Neuron_minCov_3_resized_sorted.bed
sortedbed=/dcl01/lieber/ajaffe/Amanda/ATAC/BAMs/combined_set/${ID}.final.sorted.bed

cd /dcl01/lieber/ajaffe/Amanda/ATAC/Coverage/coverage_at_CpGs
module load bedtools/2.24.0

coverageBed -sorted -g ../../genome_files/hg19.genome.reordered2 -counts -a $CPGs -b $sortedbed > ${ID}_CpGs_coverageBed.txt
