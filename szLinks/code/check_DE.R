###
library(GenomicRanges)

## read in case-control results
load("/users/ajaffe/Lieber/Projects/RNAseq/firstRnaSeqPaper/caseControl/rdas/expressed_de_features.rda")

## read in DMRs
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/bumphunting/rda/limma_Neuron_CpGs_minCov_3_ageInfo_dmrs.Rdata")
dmrList = split(dmrs, dmrs$k6cluster_label)

#### gene stats
geneStats = outStatsExprs$Gene

dmrs$overlapExprs = countOverlaps(dmrs, geneStats, maxgap = 10000) > 0
geneStats$overlapDMRs = countOverlaps(geneStats, dmrs, maxgap = 10000) > 0
geneStats$sigAndRepl = geneStats$fdr_qsva < 0.1 & 
		sign(geneStats$CMC_log2FC_qsva) == sign(geneStats$log2FC_qsva) & 
		geneStats$CMC_pval_qsva < 0.05

## test overlap overall
chisq.test(table(geneStats$overlapDMRs, geneStats$fdr_qsva < 0.05)) # no overlap
chisq.test(table(geneStats$overlapDMRs, geneStats$sigAndRepl)) # no overlap

## by type
oMat = sapply(dmrList, function(x) countOverlaps(geneStats, x, maxgap=10000) > 0)
apply(oMat, 2, function(x) {
	chisq.test(table(x, geneStats$fdr_qsva < 0.05))
})
apply(oMat, 2, function(x) {
	chisq.test(table(x, geneStats$sigAndRepl))
})

#############
### exon
exonStats = outStatsExprs$Exon

dmrs$overlapExprs = countOverlaps(dmrs, exonStats, maxgap = 50000) > 0
exonStats$overlapDMRs = countOverlaps(exonStats, dmrs, maxgap = 10000) > 0
exonStats$sigAndRepl = exonStats$fdr_qsva < 0.1 & 
		sign(exonStats$CMC_log2FC_qsva) == sign(exonStats$log2FC_qsva) & 
		exonStats$CMC_pval_qsva < 0.05
		
#############################
#### check isoform shifts ###
#############################
library(jaffelab)

isoSwitch = read.csv("/users/ajaffe/Lieber/Projects/RNAseq/firstRnaSeqPaper/revision/suppTable5_isoformSwitches_withCoord.csv",as.is=TRUE)
isoList = split(isoSwitch, isoSwitch$Type)

hasSwitch = sapply(isoList, function(x) names(geneStats) %in% x$EnsemblID)
rownames(hasSwitch) = names(geneStats)

lapply(as.data.frame(hasSwitch), function(x) table(geneStats$overlapDMRs,x))
lapply(as.data.frame(hasSwitch), function(x) chisq.test(table(geneStats$overlapDMRs,x)))
lapply(as.data.frame(hasSwitch), function(x) getOR(table(geneStats$overlapDMRs,x)))
cdDMRs_pval = sapply(as.data.frame(hasSwitch), function(x)
	chisq.test(table(geneStats$overlapDMRs,x))$p.value)
cdDMRs_or = sapply(as.data.frame(hasSwitch), function(x)
	getOR(table(geneStats$overlapDMRs,x)))

## by dmr type

orMat = pvalMat = matrix(NA, nrow = ncol(oMat), ncol = ncol(hasSwitch),
	dimnames= list(colnames(oMat), colnames(hasSwitch)))
for(i in 1:ncol(oMat)) {
	for(j in 1:ncol(hasSwitch)) {
		tt = table(oMat[,i], hasSwitch[,j])
		orMat[i,j] = getOR(tt)
		pvalMat[i,j] = chisq.test(tt)$p.value
	}
}
signif(rbind(cdDMRs_pval, pvalMat),3)
signif(rbind(cdDMRs_or, orMat),3)
out = cbind(rbind(cdDMRs_pval, pvalMat), rbind(cdDMRs_or, orMat))
colnames(out) = paste0(colnames(out), "_", rep(c("pval", "or"),each=4))
write.csv(out ,file="isoformShiftEnrichment.csv")

lapply(as.data.frame(oMat),function(x) {
	sapply(as.data.frame(hasSwitch), function(y) chisq.test(table(x, y)))
})