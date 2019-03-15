### =========================================================================
### SGE variables
### -------------------------------------------------------------------------
###

#$ -l mem_free=4G
#$ -l h_vmem=4.1G
#$ -l bluejay
#$ -m n
#$ -l h_fsize=500G
#$ -pe local 10
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
### Begin code
### -------------------------------------------------------------------------
###

for sl in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
do

python /users/aprice26/biotools/ldsc/ldsc.py \
	--l2 \
	--bfile /dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/ldsc/Phase1/1000G_plinkfiles/1000G.mac5eur.${sl} \
	--ld-wind-cm 1 \
	--annot /dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/ldsc/annotation/${cn}.Phase1.${sl}.annot.gz \
	--out /dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/ldsc/${cn}.Phase1.${sl} \
	--print-snps /dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/ldsc/Phase1/hapmap3_snps/hm.${sl}.snp \

done
