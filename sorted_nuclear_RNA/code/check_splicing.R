####

library("limma")
library("edgeR")
library(bsseq)
library(recount)

## load expression data
load("/dcl01/lieber/ajaffe/CellSorting/RNAseq_pipeline/rse_gene_CellSorting_July5_n12.Rdata")
load("/dcl01/lieber/ajaffe/CellSorting/RNAseq_pipeline/rse_exon_CellSorting_July5_n12.Rdata")
load("/dcl01/lieber/ajaffe/CellSorting/RNAseq_pipeline/rse_jx_CellSorting_July5_n12.Rdata")

## load CpG data
load('/dcl01/lieber/ajaffe/lab/brain-epigenomics/bumphunting/BSobj_bsseqSmooth_Neuron_minCov_3.Rdata')
keepIndex = which(pData(BSobj)$Working.Num %in% 12:14)
BSobj_cpg = BSobj[,keepIndex]
pd_cpg <- pData(BSobj_cpg)
gr_cpg <- granges(BSobj_cpg)
meth_cpg <- getMeth(BSobj_cpg, type = 'smooth')
rm(BSobj)

## load non-CpG data
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/bsseq/bsobj_by_chr/allChrs_postNatal_cleaned_nonCG_noHomogenate_highCov.Rdata")
keepIndex = which(pData(BSobj)$Working.Num %in% 12:14)
BSobj_non = BSobj[,keepIndex]
pd_non <- pData(BSobj_non)
gr_non <- granges(BSobj_non)
meth_non <- getMeth(BSobj_non, type = 'raw')


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