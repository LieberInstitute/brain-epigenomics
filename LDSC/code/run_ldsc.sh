# Run LDSC


### =========================================================================
### SGE variables
### -------------------------------------------------------------------------
###

#$ -l mem_free=4G
#$ -l h_vmem=4.1G
#$ -l h_fsize=500G
#$ -l bluejay
#$ -m n
#$ -o ./logs/
#$ -e ./logs/ 
#$ -cwd
#$ -M amanda.joy.price@gmail.com
#$ -t 1-36


### =========================================================================
### Load modules
### -------------------------------------------------------------------------
###

source /users/aprice26/biotools/ENTER/bin/activate ldsc

FILELIST=/dcl01/lieber/ajaffe/lab/brain-epigenomics/LDSC/categories.txt
cn=$(awk "NR==$SGE_TASK_ID" $FILELIST)


### =========================================================================
### Adjusting for baseline
### -------------------------------------------------------------------------
###

for gwas in /dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/ldsc/munge_sumstats/gwas.basenames.txt
do

if [ ${gwas} == "CNS" ]
then
    python /users/aprice26/biotools/ldsc/ldsc.py \
        --h2 /dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/ldsc/munge_sumstats/${gwas}.sumstats.gz \
        --w-ld-chr /dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/ldsc/Phase1/weights_hm3_no_hla/weights. \
        --ref-ld-chr /dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/ldsc/Phase1/cell_type_groups/CNS., \
                     /dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/ldsc/Phase1/baseline/baseline. \
        --overlap-annot \
        --frqfile-chr /dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/ldsc/Phase1/1000G_frq/1000G.mac5eur. \
        --out /dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/ldsc/output/CNS.${gwas}.Phase1 \
        --print-coefficients

else
    python /users/aprice26/biotools/ldsc/ldsc.py \
        --h2 /dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/ldsc/munge_sumstats/${gwas}.sumstats.gz \
        --w-ld-chr /dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/ldsc/Phase1/weights_hm3_no_hla/weights. \
        --ref-ld-chr /dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/ldsc/${cn}.Phase1., \
                     /dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/ldsc/Phase1/baseline/baseline. \
        --overlap-annot \
        --frqfile-chr /dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/ldsc/Phase1/1000G_frq/1000G.mac5eur. \
        --out /dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/ldsc/output/${cn}.${gwas}.Phase1 \
        --print-coefficients
fi

done