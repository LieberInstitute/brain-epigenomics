# Run LDSC
# Peter Hickey
# 2017-11-21

### =========================================================================
### SGE variables
### -------------------------------------------------------------------------
###

#$ -l mem_free=4G
#$ -l h_vmem=4.1G
#$ -m n
#$ -l h_fsize=500G
#$ -l cegs
#$ -pe local 10

### =========================================================================
### Load modules
### -------------------------------------------------------------------------
###

module load python/2.7.6

### =========================================================================
### Run permutation script
### -------------------------------------------------------------------------
###

Rscript run_LDScore_estimation.R ${SGE_TASK_ID}