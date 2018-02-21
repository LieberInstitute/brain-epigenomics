#!/bin/bash

## Usage:
# mkdir -p logs
# sh meth_vs_expr_QTL_near_nonCpG_meQTLs.sh

mkdir -p logs

for cpg in TRUE
do
    for type in jx #gene exon jx psi
    do

SHORT="meth_vs_expr_QTL_cpg${cpg}_${type}_near_nonCpG_meQTLs"

# Construct shell file
echo "Creating script for cpg ${cpg} for feature type ${type}"

cat > .${SHORT}.sh <<EOF
#!/bin/bash
#$ -cwd
#$ -l mem_free=350G,h_vmem=350G,h_fsize=200G
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
Rscript meth_vs_expr_QTL_near_nonCpG_meQTLs.R -c ${cpg} -f ${type}

echo "**** Job ends ****"
date
EOF

call="qsub .${SHORT}.sh"
echo $call
$call
    done
done
