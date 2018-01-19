#!/bin/bash

## Usage:
# mkdir -p logs
# sh meth_vs_expr_QTL.sh

mkdir -p logs

for cpg in FALSE
#for cpg in TRUE FALSE
do
    for type in gene exon jx #psi
    do

SHORT="meth_vs_expr_QTL_cpg${cpg}_${type}"

# Construct shell file
echo "Creating script for cpg ${cpg} for feature type ${type}"

cat > .${SHORT}.sh <<EOF
#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=100G,h_vmem=100G,h_fsize=100G
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
Rscript meth_vs_expr_QTL.R -c ${cpg} -f ${type}

echo "**** Job ends ****"
date
EOF

call="qsub .${SHORT}.sh"
echo $call
$call
    done
done
