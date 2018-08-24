# qsub -cwd -l bluejay,mf=10G,h_vmem=10G,h_fsize=100G,h_stack=256M -M amanda.joy.price@gmail.com -t 1-12 -tc 30 download_Gabel_Greenberg_Nature_2015_WGBS.sh

FILELIST=/dcl02/lieber/WGBS/PublicData/Gabel_Greenberg_Nature_2015/FASTQ/WGBS/SRR_Acc_List.txt
ID=$(awk "NR==$SGE_TASK_ID" $FILELIST)

/users/ajaffe/software/sratoolkit.2.8.1-3-centos_linux64/bin/fastq-dump -O /dcl02/lieber/WGBS/PublicData/Gabel_Greenberg_Nature_2015/FASTQ/WGBS/ --split-3 --gzip $ID
