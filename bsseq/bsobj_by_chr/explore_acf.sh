#!/bin/bash

## Usage:
# sh explore_acf.sh

mkdir -p logs

#for chrnum in {1..22} X Y M
#for chrnum in M Y 21 22
#for chrnum in {2..20} X
#for chrnum in {2..22} X Y M
for chrnum in {1..22} X Y M
do 

    chr="chr${chrnum}"
    
    if [[ "$chrnum" == "1" ]]
    then
        MEM=80
    elif [[ "$chrnum" == "2" ]] || [[ "$chrnum" == "3" ]]
    then
        MEM=75
    else
        MEM=55
    fi
    
    for context in nonCG CG CHG CHH
    do

    
    SHORT="acf_${chr}_context_${context}"

    # Construct shell file
    echo "Creating script for chromosome ${chr} with ${MEM}G of RAM"
    cat > .${SHORT}.sh <<EOF
#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=${MEM}G,h_vmem=${MEM}G,h_fsize=100G
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
Rscript explore_acf.R -c ${chr} -t 1 -x ${context}

echo "**** Job ends ****"
date
EOF
        call="qsub .${SHORT}.sh"
        echo $call
        $call
    done
done
