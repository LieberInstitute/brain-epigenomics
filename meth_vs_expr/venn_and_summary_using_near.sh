#!/bin/bash

## Usage:
# mkdir -p logs
# sh venn_and_summary_using_near.sh

mkdir -p logs

for type in jxleft jxright
    do

SHORT="venn_and_summary_using_near_${type}"

# Construct shell file
echo "Creating scriptfeature type ${type} side"

cat > .${SHORT}.sh <<EOF
#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=210G,h_vmem=210G,h_fsize=200G
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

module load conda_R/3.4.x
Rscript venn_and_summary_using_near.R -f ${type}

echo "**** Job ends ****"
date
EOF

call="qsub .${SHORT}.sh"
echo $call
$call

done
