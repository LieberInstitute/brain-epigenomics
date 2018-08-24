# qsub -cwd -l bluejay,mf=10G,h_vmem=10G,h_fsize=500G,h_stack=256M -t 1-11 -tc 30 Bowtie2_Stroud_Greenberg_ChIPseq.sh

FILELIST=/dcl02/lieber/WGBS/PublicData/Stroud_Greenberg_Cell_2017/FASTQ/ChIPseq/SRR_Acc_List.txt
ID=$(awk "NR==$SGE_TASK_ID" $FILELIST)
fastq=/dcl02/lieber/WGBS/PublicData/Stroud_Greenberg_Cell_2017/FASTQ/ChIPseq/${ID}.fastq.gz

module load bowtie2/2.2.5
module load samtools/1.1

cd /dcl02/lieber/WGBS/PublicData/Stroud_Greenberg_Cell_2017/BAMS/ChIPseq

bowtie2 -x /dcl01/lieber/ajaffe/Amanda/Mus_musculus/Ensembl/GRCm38/Sequence/Bowtie2Index/genome -U $fastq -S ${ID}.sam

samtools view -b ${ID}.sam > ${ID}.bam

rm ${ID}.sam
