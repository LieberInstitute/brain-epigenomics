#!/bin/bash

## Usage:
# mkdir -p logs
# sh compute_DMR_lag_dist.sh

mkdir -p logs

for model in age cell interaction
do
    for context in all nonCG CG CHG CHH
    do

SHORT="compute_DMR_lag_dist_${model}_${context}"

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

Rscript compute_DMR_lag_dist.R -m ${model} -x ${context}

echo "**** Job ends ****"
date
EOF

    call="qsub .${SHORT}.sh"
    echo $call
    $call
    done
done
