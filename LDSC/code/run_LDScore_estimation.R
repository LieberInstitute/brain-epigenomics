# Run LD score estimation with functional categories (made in make-annot.R)
# Code adapted from Peter Hickey (https://github.com/hansenlab/BrainEpigenomeNN/blob/master/SLDSR/scripts/run_LDScore_estimation.R)


categories = readRDS("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/ldsc/categories.rds")

seqlevels = 1:22

#------------------------------------------------------------------------------
# 'base' category
#
# Ran interactively

mclapply(seqlevels, function(sl) {
 cmd <- paste0("python ",
               "/users/aprice26/biotools/ldsc/ldsc.py ",
               "--l2 ",
               "--bfile /dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/ldsc/Phase1/1000G_plinkfiles/1000G.mac5eur.", sl, " ",
               "--ld-wind-cm 1 ",
               "--annot /dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/ldsc/annotation/base.Phase1.", sl, ".annot.gz ",
               "--out /dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/ldsc/base.Phase1.", sl, " ",
               "--print-snps /dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/ldsc/Phase1/hapmap3_snps/hm.", sl, ".snp")
 print(cmd)
 system(cmd)
}, mc.cores = 10)
