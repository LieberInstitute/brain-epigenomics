# qsub -tc 245 -l mf=10G,h_vmem=10G,h_fsize=20G,h_stack=256M -m e -M amanda.joy.price@gmail.com -t 1-34 sort.pooled.bed.files.sh

FILELIST=/dcl01/lieber/ajaffe/Amanda/ATAC/FinalPeaks/scripts/finalizedList-pooled.txt
ID=$(awk "NR==$SGE_TASK_ID" $FILELIST)

cd /dcl01/lieber/ajaffe/Amanda/ATAC/BAMs/combined_set/

sort -k1,1 -k2,2n ${ID}.final.10bp-added2.bed > ${ID}.final.sorted.bed
