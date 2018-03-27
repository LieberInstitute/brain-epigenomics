#!/bin/bash

## Usage:
# sh create_matched_versions.sh

mkdir -p logs

for chrnum in {1..22} X Y M
do 

chr="chr${chrnum}"
SHORT="filter_${chr}"

# Construct shell file
echo "Creating script for chromosome ${chr}"
cat > .${SHORT}.sh <<EOF
#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=150G,h_vmem=150G,h_fsize=200G
#$ -N ${SHORT}
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

Rscript create_matched_versions.R -c ${chr}

echo "**** Job ends ****"
date
EOF
    call="qsub .${SHORT}.sh"
    echo $call
    $call
done


for type in CpG CpH
do 
SHORT="create_matched_versions_${type}"

# Construct shell file
echo "Creating script for type ${type}"
cat > .${SHORT}.sh <<EOF
#!/bin/bash
#$ -cwd
#$ -l mem_free=300G,h_vmem=300G,h_fsize=300G
#$ -N ${SHORT}
#$ -o ./logs/${SHORT}.txt
#$ -e ./logs/${SHORT}.txt
#$ -m e
#$ -hold_jid filter_chr*

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: \${USER}"
echo "Job id: \${JOB_ID}"
echo "Job name: \${JOB_NAME}"
echo "Hostname: \${HOSTNAME}"
echo "Task id: \${TASK_ID}"

Rscript create_matched_versions.R -c "all" -t ${type}

echo "**** Job ends ****"
date
EOF
    call="qsub .${SHORT}.sh"
    echo $call
    $call
done
