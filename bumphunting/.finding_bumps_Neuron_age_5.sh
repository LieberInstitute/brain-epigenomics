#!/bin/bash
#$ -cwd
#$ -l mem_free=10G,h_vmem=12G,h_fsize=100G
#$ -N finding_bumps_Neuron_age_5
#$ -pe local 8
#$ -o ./logs/finding_bumps_Neuron_age_5.txt
#$ -e ./logs/finding_bumps_Neuron_age_5.txt
#$ -m e

echo "**** Job starts ****"
date

Rscript find_bumps.R -m age -s Neuron -p 5 -t 8

echo "**** Job ends ****"
date
