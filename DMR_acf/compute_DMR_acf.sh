#!/bin/bash

## Usage:
# mkdir -p logs
# sh compute_DMR_acf.sh

mkdir -p logs

for model in age cell interaction
do
    #for context in nonCG CG CHG CHH
    for context in all
    do

SHORT="compute_DMR_acf_${model}_${context}"

# Construct shell file
echo "Creating script for model ${model} under the context of ${context}"

cat > .${SHORT}.sh <<EOF
#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=5G,h_vmem=5G,h_fsize=100G
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
Rscript compute_DMR_acf.R -t 1 -m ${model} -x ${context}

echo "**** Job ends ****"
date
EOF

call="qsub .${SHORT}.sh"
echo $call
$call
    done
done
