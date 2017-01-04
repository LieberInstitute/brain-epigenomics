#!/bin/sh

## Usage
# sh BismarkMethylationExtractor_report.sh

## Based on /users/ajaffe/Lieber/Projects/WGBS/Analysis/BismarkMethylationExtractor_report.sh

SHORT="bismark-non-CpG"

mkdir -p logs

# Construct shell files
FILELIST=/dcl01/lieber/WGBS/LIBD_Data/WGC_IDs.txt
NUM=$(cat $FILELIST | awk '{print $NF}' | uniq | wc -l)
echo "Creating script ${SHORT}"

cat > .${SHORT}.sh <<EOF
#!/bin/bash
#$ -cwd
#$ -l mem_free=16G,h_vmem=20G,h_fsize=100G
#$ -N ${SHORT}
#$ -pe local 4
#$ -o ./logs/${SHORT}.o.\$TASK_ID.txt
#$ -e ./logs/${SHORT}.e.\$TASK_ID.txt
#$ -t 1-${NUM}
#$ -tc 20
#$ -m e
echo "**** Job starts ****"
date

ID=\$(awk "NR==$SGE_TASK_ID" $FILELIST )

## Create ouput directory
mkdir -p Reports
mkdir -p Reports/\${ID}

/users/ajaffe/software/bismark_v0.16.3/bismark_methylation_extractor --single-end \
	--cytosine_report --genome_folder /dcl01/lieber/ajaffe/Annotation/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/ \
	--gzip --multicore 4 --bedGraph  \
	-o Reports/\${ID} --CX_context --split_by_chromosome \
	/dcl01/lieber/WGBS/LIBD_Data/BAM/\${ID}.concatenated.sorted.duplicatesRemoved.bam

echo "**** Job ends ****"
date
EOF

call="qsub .${SHORT}.sh"
echo $call
$call
