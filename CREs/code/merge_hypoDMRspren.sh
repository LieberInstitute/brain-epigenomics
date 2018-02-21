# qsub -tc 245 -l mf=5G,h_vmem=5G,h_stack=256M -m e -M amanda.joy.price@gmail.com -t 1-20 merge_hypoDMRspren.sh

FILELIST=/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/hypoDMR_prenatal_IDs.txt
ID=$(awk "NR==$SGE_TASK_ID" $FILELIST)
bed=/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/hypoDMRspren_${ID}.tab

echo ${ID}

module load bedtools/2.24.0

bedtools merge -d 1000 -i $bed > /dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/${ID}_merged_hypoDMRspren.tab
