#!/bin/bash

## Usage:
# mkdir -p logs
# sh meth_vs_expr_QTL_jxonly.sh

mkdir -p logs

for cpg in TRUE
do
    for type in jx
    do
        for jxside in left right
        do

SHORT="meth_vs_expr_QTL_cpg${cpg}_${type}_near_nonCpG_meQTLs_${jxside}"

# Construct shell file
echo "Creating script for cpg ${cpg} for feature type ${type} side ${jxside} (near nonCpGs)"

cat > .${SHORT}.sh <<EOF
#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=140G,h_vmem=140G,h_fsize=200G
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
Rscript meth_vs_expr_QTL_near_nonCpG_meQTLs.R -c ${cpg} -f ${type} -j ${jxside}

echo "**** Job ends ****"
date
EOF

call="qsub .${SHORT}.sh"
echo $call
$call

SHORT="meth_vs_expr_QTL_cpg${cpg}_${type}_${jxside}"

# Construct shell file
echo "Creating script for cpg ${cpg} for feature type ${type} side ${jxside} (all)"

cat > .${SHORT}.sh <<EOF
#!/bin/bash
#$ -cwd
#$ -l mem_free=180G,h_vmem=180G,h_fsize=200G
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
Rscript meth_vs_expr_QTL.R -c ${cpg} -f ${type} -j ${jxside}

echo "**** Job ends ****"
date
EOF

call="qsub .${SHORT}.sh"
echo $call
$call
        done
    done
done
