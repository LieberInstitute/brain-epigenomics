#!/bin/bash

## Usage:
# sh find_bumps_bsseqSmooth.sh

CORES=4
mkdir -p logs

for cell in Neuron
do
    for model in cell age interaction
    do
        for permutations in 0
        do

SHORT="finding_bumps_bsseqSmooth_${cell}_${model}_${permutations}"

# Construct shell file
echo "Creating script for chromosome Glia + ${cell} using model ${model} with ${permutations} permutations"

cat > .${SHORT}.sh <<EOF
#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=45G,h_vmem=45G,h_fsize=100G
#$ -N ${SHORT}
#$ -pe local ${CORES}
#$ -o ./logs/${SHORT}.txt
#$ -e ./logs/${SHORT}.txt
#$ -m e

echo "**** Job starts ****"
date

Rscript find_bumps_bsseqSmooth.R -m ${model} -s ${cell} -p ${permutations} -t ${CORES}

echo "**** Job ends ****"
date
EOF

call="qsub .${SHORT}.sh"
echo $call
$call
        done
    done
done