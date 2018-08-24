# qsub -cwd -l bluejay,mf=10G,h_vmem=10G,h_fsize=500G -t 1-11 -tc 30 download_Stroud_Greenberg_Cell_2017_ChIPseq.sh

FILELIST=/dcl02/lieber/WGBS/PublicData/Stroud_Greenberg_Cell_2017/FASTQ/ChIPseq/SRR_Acc_List.txt
ID=$(awk "NR==$SGE_TASK_ID" $FILELIST)

echo ${ID}

/users/ajaffe/software/sratoolkit.2.8.1-3-centos_linux64/bin/fastq-dump -O /dcl02/lieber/WGBS/PublicData/Stroud_Greenberg_Cell_2017/FASTQ/ChIPseq/ --split-3 --gzip $ID
