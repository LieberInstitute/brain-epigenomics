# qsub -cwd -l bluejay,mf=5G,h_vmem=5G,h_fsize=100G,h_stack=256M -t 1-3 -tc 30 MACS2_Gabel_Greenberg_ChIPseq.sh

CHIP=/dcl02/lieber/WGBS/PublicData/Gabel_Greenberg_Nature_2015/FASTQ/ChIPseq/ChIP_ids.txt
CID=$(awk "NR==$SGE_TASK_ID" $CHIP)

INPUT=/dcl02/lieber/WGBS/PublicData/Gabel_Greenberg_Nature_2015/FASTQ/ChIPseq/Input_ids.txt
IID=$(awk "NR==$SGE_TASK_ID" $INPUT)

NAME=/dcl02/lieber/WGBS/PublicData/Gabel_Greenberg_Nature_2015/FASTQ/ChIPseq/Test_ids.txt
NID=$(awk "NR==$SGE_TASK_ID" $NAME)


module load macs/2.1.0

cd /dcl02/lieber/WGBS/PublicData/Gabel_Greenberg_Nature_2015/BAMS/ChIPseq

macs2 callpeak -t ${CID} -c ${IID} -g mm -n ${NID} --broad --outdir /dcl02/lieber/WGBS/PublicData/Gabel_Greenberg_Nature_2015/BAMS/ChIPseq/Peaks
