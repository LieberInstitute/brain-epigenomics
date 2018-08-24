# qsub -cwd -l bluejay,mf=10G,h_vmem=10G,h_fsize=500G,h_stack=256M -t 1-12 -tc 30 Bowtie2_Gabel_Greenberg_ChIPseq.sh

FILELIST=/dcl02/lieber/WGBS/PublicData/Gabel_Greenberg_Nature_2015/FASTQ/ChIPseq/SRR_Acc_List.txt
ID=$(awk "NR==$SGE_TASK_ID" $FILELIST)
fastq=/dcl02/lieber/WGBS/PublicData/Gabel_Greenberg_Nature_2015/FASTQ/ChIPseq/${ID}.fastq.gz

module load bowtie2/2.2.5
module load samtools/1.1

cd /dcl02/lieber/WGBS/PublicData/Gabel_Greenberg_Nature_2015/BAMS/ChIPseq/

bowtie2 -x /dcl01/lieber/ajaffe/Amanda/Mus_musculus/Ensembl/GRCm38/Sequence/Bowtie2Index/genome -U $fastq -S /dcl02/lieber/WGBS/PublicData/Gabel_Greenberg_Nature_2015/BAMS/ChIPseq/${ID}.sam

samtools view -b ${ID}.sam > ${ID}.bam

rm ${ID}.sam
