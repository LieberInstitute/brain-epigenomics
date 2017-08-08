#!/bin/bash

## Usage:
# sh plot_bumps_bsseqSmooth.sh

mkdir -p logs

for cell in Neuron
do
    for model in cell age interaction
    do
        for permutations in 250
        do

SHORT="plotting_bumps_bsseqSmooth_${cell}_${model}_${permutations}"

# Construct shell file
echo "Creating script for chromosome Glia + ${cell} using model ${model} with ${permutations} permutations"

cat > .${SHORT}.sh <<EOF
#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=55G,h_vmem=55G,h_fsize=100G
#$ -N ${SHORT}
#$ -o ./logs/${SHORT}.txt
#$ -e ./logs/${SHORT}.txt
#$ -m e

echo "**** Job starts ****"
date

Rscript plot_bumps_bsseqSmooth.R -m ${model} -s ${cell} -p ${permutations} -b FALSE -i FALSE

echo "**** Job ends ****"
date
EOF

call="qsub .${SHORT}.sh"
echo $call
$call
        done
    done
done
