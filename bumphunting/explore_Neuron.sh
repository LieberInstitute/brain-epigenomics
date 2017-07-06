#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=150G,h_vmem=150G,h_fsize=100G
#$ -N explore_Neuron
#$ -o ./logs/explore_Neuron.txt
#$ -e ./logs/explore_Neuron.txt
#$ -m e

echo "**** Job starts ****"
date

Rscript explore_Neuron.R

echo "**** Job ends ****"
date
