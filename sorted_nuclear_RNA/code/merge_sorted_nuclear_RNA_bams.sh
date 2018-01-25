# qsub -tc 245 -l mf=5G,h_vmem=5G,h_stack=256M -m e -M amanda.joy.price@gmail.com -t 1-6 merge_sorted_nuclear_RNA_bams.sh

FILELIST=/dcl01/lieber/ajaffe/Amanda/sorted_nuclear_RNA/merged/sorted_nuclear_RNA_IDs.txt
ID=$(awk "NR==$SGE_TASK_ID" $FILELIST)
polyA=/dcl01/lieber/ajaffe/CellSorting/RNAseq_pipeline/HISAT2_out/${ID}_PolyA_accepted_hits.sorted.bam
Ribo=/dcl01/lieber/ajaffe/CellSorting/RNAseq_pipeline/HISAT2_out/${ID}_Ribo_accepted_hits.sorted.bam

echo ${ID}


cd /dcl01/lieber/ajaffe/Amanda/sorted_nuclear_RNA/merged/
module load samtools

samtools merge -f ${ID}_combinedLibraries.bam $polyA $Ribo
