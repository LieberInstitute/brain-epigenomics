####
####

library("limma")
library("edgeR")
library(bsseq)
library(recount)

## load expression data
load("/dcl01/lieber/ajaffe/CellSorting/RNAseq_pipeline/rse_jx_CellSorting_July5_n12.Rdata")

## load CpG and non-CpG data
load('/dcl01/lieber/ajaffe/lab/brain-epigenomics/bumphunting/BSobj_bsseqSmooth_Neuron_minCov_3.Rdata')
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/non-CpG/nonCG_highCov_neuronsOnly_pca_pd_methTable.Rdata")
meth_ch = meth
methMap_ch = gr

meth_cg =getMeth(BSobj, type = 'raw')
methMap_cg =granges(BSobj)

# get expression
rowData(rse_jx)$bp_length = 100
jRp10m = getRPKM(rse_jx)
jMap = rowRanges(rse_jx)

## filter to the same people
keepIndex = which(pd$Working.Num %in% 12:14)
meth_ch = meth_ch[,keepIndex]
meth_cg = meth_cg[,keepIndex]
pd = pd[keepIndex,]
jRp10m = jRp10m[,c(2,6,4)] #polyA, neun+

## get change over age
mod = model.matrix(~pd$Age)

fit_ch = lmFit(meth_ch, mod)