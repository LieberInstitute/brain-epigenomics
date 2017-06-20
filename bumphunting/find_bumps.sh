#!/bin/bash

## Usage:
# sh find_bumps.sh

CORES=8
mkdir -p logs

for cell in Neuron
do
    for model in cell age interaction
    do
        for permutations in 5
        do

SHORT="finding_bumps_${cell}_${model}_${permutations}"

# Construct shell file
echo "Creating script for chromosome Glia + ${cell} using model ${model} with ${permutations} permutations"

cat > .${SHORT}.sh <<EOF
#!/bin/bash
#$ -cwd
#$ -l mem_free=10G,h_vmem=12G,h_fize=100G
#$ -N ${SHORT}
#$ -pe local ${CORES}
#$ -o ./logs/${SHORT}.o.txt
#$ -e ./logs/${SHORT}.e.txt
#$ -m e

echo "**** Job starts ****"
date

Rscript find_bumps.R -m ${model} -s ${cell} -p ${permutations} -t ${CORES}

echo "**** Job ends ****"
date
EOF

call="qsub .${SHORT}.sh"
echo $call
$call
        done
    done
done
