#!/bin/bash

## Usage:
# sh create_bsobj_by_chr.sh

mkdir -p logs

#for chrnum in {1..22} X Y M
#for chrnum in M Y 21 22
#for chrnum in {3..20} X
for chrnum in 1
do 

chr="chr${chrnum}"
SHORT="bsseq_${chr}"

# Construct shell file
echo "Creating script for chromosome ${chr}"
cat > .${SHORT}.sh <<EOF
#!/bin/bash
#$ -cwd
#$ -l mem_free=35G,h_vmem=40G,h_fsize=100G
#$ -N ${SHORT}
#$ -pe local 4
#$ -o ./logs/${SHORT}.txt
#$ -e ./logs/${SHORT}.txt
#$ -m e

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: \${USER}"
echo "Job id: \${JOB_ID}"
echo "Job name: \${JOB_NAME}"
echo "Hostname: \${HOSTNAME}"
echo "Task id: \${TASK_ID}"

Rscript create_bsobj_by_chr.R -c ${chr} -t 4

echo "**** Job ends ****"
date
EOF
    call="qsub .${SHORT}.sh"
    echo $call
    $call
done
