#!/bin/bash

## Usage:
# sh compute_aucs.sh

mkdir -p logs
mkdir -p auc_files

CORES=4
FILELIST=/dcl01/lieber/ajaffe/lab/brain-epigenomics/misc/cov/WGC_IDs_subset.txt
NUM=$(cat $FILELIST | awk '{print $NF}' | uniq | wc -l)

for bamtype in duplicatesRemoved Marked_duplicates
do
    
    SHORT="compute_aucs_${bamtype}"

    # Construct shell file
    echo "Creating script for BAM files ${bamtype}"
    cat > .${SHORT}.sh <<EOF
#!/bin/bash
#$ -cwd
#$ -l mem_free=3G,h_vmem=3G,h_fsize=100G
#$ -pe local ${CORES}
#$ -N ${SHORT}
#$ -o ./logs/${SHORT}.\$TASK_ID.txt
#$ -e ./logs/${SHORT}.\$TASK_ID.txt
#$ -m e
#$ -t 1-${NUM}

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: \${USER}"
echo "Job id: \${JOB_ID}"
echo "Job name: \${JOB_NAME}"
echo "Hostname: \${HOSTNAME}"
echo "Task id: \${TASK_ID}"


ID=\$(awk "NR==\${SGE_TASK_ID}" $FILELIST )

## Load dependencies
module load bamcount/0.2.6

## List current modules
module list

## compute aucs
ls -lh /dcl02/lieber/WGBS/LIBD_Data/BAM/\${ID}*${bamtype}.bam
bamcount /dcl02/lieber/WGBS/LIBD_Data/BAM/\${ID}*${bamtype}.bam --threads 4 --no_head --auc auc_files/\${ID}_${bamtype}

echo "**** Job ends ****"
date
EOF
    call="qsub .${SHORT}.sh"
    echo $call
    $call
done
