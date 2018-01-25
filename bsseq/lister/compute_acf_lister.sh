#!/bin/bash

## Usage:
# sh compute_acf_lister.sh

mkdir -p logs

for chrnum in {1..22} X Y M
do 

    chr="chr${chrnum}"
    
    for context in nonCG CG CHG CHH
    do

    
    SHORT="lister_acf_${chr}_context_${context}"

    # Construct shell file
    echo "Creating script for chromosome ${chr}"
    cat > .${SHORT}.sh <<EOF
#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=30G,h_vmem=30G,h_fsize=100G
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
Rscript compute_acf_lister.R -c ${chr} -t 1 -x ${context}

echo "**** Job ends ****"
date
EOF
        call="qsub .${SHORT}.sh"
        echo $call
        $call
    done
done
