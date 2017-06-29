#!/bin/bash
#$ -cwd
#$ -l mem_free=35G,h_vmem=40G,h_fsize=100G
#$ -N finding_bumps_Neuron_cell_0
#$ -pe local 4
#$ -o ./logs/finding_bumps_Neuron_cell_0.txt
#$ -e ./logs/finding_bumps_Neuron_cell_0.txt
#$ -m e

echo "**** Job starts ****"
date

Rscript find_bumps.R -m cell -s Neuron -p 0 -t 4

echo "**** Job ends ****"
date
