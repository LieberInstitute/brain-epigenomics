# Run LDSC


### =========================================================================
### SGE variables
### -------------------------------------------------------------------------
###

#$ -l mem_free=6G
#$ -l h_vmem=6.1G
#$ -l h_fsize=500G
#$ -l bluejay
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


### ============================================================================
### Adjusting for baseline + CNS + chromHMM_union + Prenatal_LMRs + Neuronal_LMRs + Glial_LMRs
###

for gwas in ADHD Agreeableness Alzheimers_disease Anorexia_nervosa Anxiety_disorder Autism_spectrum_disorder Bipolar_disorder BMI Cardioembolic_stroke Childhood_cognitive_performance Cigarettes_per_day College_attainment Conscientiousness Coronary_artery_disease Crohns_disease Depressive_symptoms Epilepsy Ever_smoked Extraversion Focal_epilepsy Generalized_epilepsy Height Intracarebral_hemorrhage IQ Ischemic_stroke Large-vessel_disease Major_depressive_disorder Neuroticism Openness PTSD Schizophrenia Small-vessel_disease Subjective_well-being Years_of_education
do

if [ ! -e ${cn}.${gwas}.adjusting_for_CNS_LMRs_chromHMM.Phase1.results ]
then

if [[ ! "${cn}" =~ ^(CNS|Prenatal_LMRs|Neuronal_LMRs|Glial_LMRs|chromHMM_union)$ ]]
then


python /users/aprice26/biotools/ldsc/ldsc.py \
	--h2 /dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/ldsc/munge_sumstats/${gwas}.sumstats.gz \
        --w-ld-chr /dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/ldsc/Phase1/weights_hm3_no_hla/weights. \
	--ref-ld-chr /dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/ldsc/Phase1/baseline/baseline.,/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/ldsc/Phase1/cell_type_groups/CNS.,/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/ldsc/chromHMM_union.Phase1.,/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/ldsc/Prenatal_LMRs.Phase1.,/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/ldsc/Neuronal_LMRs.Phase1.,/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/ldsc/Glial_LMRs.Phase1.,/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/ldsc/${cn}.Phase1. \
        --overlap-annot \
        --frqfile-chr /dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/ldsc/Phase1/1000G_frq/1000G.mac5eur. \
        --out /dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/ldsc/output/${cn}.${gwas}.adjusting_for_CNS_LMRs_chromHMM.Phase1 \
        --print-coefficients

else
	echo Nothing to do for ${cn}

fi

fi

done
