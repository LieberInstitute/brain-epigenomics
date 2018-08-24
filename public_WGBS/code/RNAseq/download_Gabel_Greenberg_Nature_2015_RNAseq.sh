# qsub -cwd -l bluejay,mf=30G,h_vmem=30G,h_fsize=100G,h_stack=256M -M amanda.joy.price@gmail.com -t 1 -tc 30 download_Gabel_Greenberg_Nature_2015_RNAseq.sh

FILELIST=/dcl02/lieber/WGBS/PublicData/Gabel_Greenberg_Nature_2015/FASTQ/RNAseq/SRR_Acc_List2.txt
ID=$(awk "NR==$SGE_TASK_ID" $FILELIST)

/users/ajaffe/software/sratoolkit.2.8.1-3-centos_linux64/bin/fastq-dump -O /dcl01/lieber/ajaffe/Amanda/WGBS/ --split-3 --gzip SRR1930034
