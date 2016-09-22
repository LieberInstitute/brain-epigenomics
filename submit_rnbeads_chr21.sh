#!/bin/bash
#$ -cwd
#$ -m e
#$ -M amanda.price@libd.org
#$ -l mem_free=400G,h_vmem=400G,h_fsize=125G
#$ -N RnBeads-chr21-job

echo "**** Job starts ****"
date

# Make logs directory
mkdir -p logs

# run RNBeads.chr21.test.R
Rscript RNBeads.chr21.test.R

# Move log files into the logs directory
mv RnBeads-chr21-job.* logs/

echo "**** Job ends ****"
date
