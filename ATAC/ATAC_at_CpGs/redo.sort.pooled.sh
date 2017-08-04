# qsub -tc 245 -l mf=20G,h_vmem=20G,h_fsize=30G,h_stack=256M -m e -M amanda.joy.price@gmail.com -t 1-7 redo.sort.pooled.sh

FILELIST=/dcl01/lieber/ajaffe/lab/brain-epigenomics/ATAC/ATAC_at_CpGs/redo.sort.txt
ID=$(awk "NR==$SGE_TASK_ID" $FILELIST)

cd /dcl01/lieber/ajaffe/Amanda/ATAC/BAMs/combined_set/

sort -k1,1 -k2,2n ${ID}.final.10bp-added2.bed > ${ID}.final.sorted.bed
