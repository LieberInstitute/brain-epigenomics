# qsub -cwd -l bluejay,mf=20G,h_vmem=20G,h_fsize=500G,h_stack=256M -t 1-21 -tc 30 download_Stroud_Greenberg_Cell_2017_WGBS.sh

FILELIST=/dcl02/lieber/WGBS/PublicData/Stroud_Greenberg_Cell_2017/FASTQ/WGBS/SRR_Acc_List.txt
ID=$(awk "NR==$SGE_TASK_ID" $FILELIST)

echo $ID

/users/ajaffe/software/sratoolkit.2.8.1-3-centos_linux64/bin/fastq-dump -O /dcl01/lieber/ajaffe/Amanda/WGBS/wgbs/ --split-3 --gzip $ID
