# qsub -cwd -l bluejay,mf=10G,h_vmem=10G,h_fsize=100G -t 1-2 -tc 30 download_Stroud_Greenberg_Cell_2017_RNAseq.sh

FILELIST=/dcl02/lieber/WGBS/PublicData/Stroud_Greenberg_Cell_2017/FASTQ/RNAseq/SRR_Acc_List2.txt
ID=$(awk "NR==$SGE_TASK_ID" $FILELIST)

/users/ajaffe/software/sratoolkit.2.8.1-3-centos_linux64/bin/fastq-dump -O /dcl02/lieber/WGBS/PublicData/Stroud_Greenberg_Cell_2017/FASTQ/RNAseq/ --split-3 --gzip $ID
