# qsub -tc 245 -l mf=5G,h_vmem=5G,h_stack=256M -m e -M amanda.joy.price@gmail.com -t 1-38 merge_hypoDMRs.sh

FILELIST=/dcl01/lieber/ajaffe/Amanda/sorted_nuclear_RNA/merged/hypoDMR_IDs.txt
ID=$(awk "NR==$SGE_TASK_ID" $FILELIST)
bed=/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/PWMEnrich/hypoDMRs_${ID}.tab

echo ${ID}

module load bedtools/2.24.0

bedtools merge -d 1000 -i $bed > /dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/PWMEnrich/${ID}_merged_hypoDMRs.tab
