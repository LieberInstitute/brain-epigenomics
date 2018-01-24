## 
library(GenomicRanges)
library('bsseq')
library('devtools')
library('limma')

## Load the raw data, need the map/annotation/coordinates
load('../bumphunting/BSobj_bsseqSmooth_Neuron_minCov_3.Rdata')
map = rowRanges(BSobj)

## load WGBS model fits
load('../bumphunting/limma_Neuron_CpGs_minCov_3.Rdata')

## interaction model:
wgbsStats = cbind(fits$interaction$coefficients[,2:4],
	fits$interaction$t[,2:4],fits$interaction$p.value[,2:4])
colnames(wgbsStats) = paste0(rep(c("dm", "t","p"), each =3), "_",
	rep(c("Age","Cell","Int"), times=3))
wgbsStats = as.data.frame(wgbsStats)

## array stats	
load("LIBD_DNAm_age_stats_postnatal_GRanges.rda")
load("brainspan_DNAm_age_stats_postnatal.rda")
outSpan = outSpan[names(outLibd),]
colnames(mcols(outLibd))[10:12] = paste0("LIBD_",colnames(mcols(outLibd))[10:12])

colnames(outSpan) = paste0("Span_", colnames(outSpan))
mcols(outLibd) = cbind(mcols(outLibd), outSpan)

outLibd$ProbeSeqA = outLibd$ProbeSeqB = NULL
outLibd$LIBD_fdr= p.adjust(outLibd$LIBD_pval, "fdr")

## match up ###
outLibd$chrpos = paste0(seqnames(outLibd), ":", start(outLibd))
map$chrpos = paste0(seqnames(map), ":", start(map))

table(outLibd$chrpos %in% map$chrpos) 
# 407880 w/ coverage

mm = match(outLibd$chrpos, map$chrpos)
outLibd = outLibd[!is.na(mm),]
wgbsStats = wgbsStats[mm[!is.na(mm)],]

## plots?
plot(outLibd$LIBD_tstat, wgbsStats$t_Age, pch = 21, bg="grey")
cor(outLibd$LIBD_tstat, wgbsStats$t_Age)

plot(outLibd$LIBD_ageChange, wgbsStats$dm_Age, pch = 21, bg="grey")
cor(outLibd$LIBD_ageChange, wgbsStats$dm_Age)
plot(outLibd$Span_ageChange_dfc, wgbsStats$dm_Age, pch = 21, bg="grey")

## break into age and cell type
sigInd = which(outLibd$LIBD_fdr < 0.05)
mean(wgbsStats$p_Age[sigInd] < 0.05)
mean(wgbsStats$p_Age[-sigInd] < 0.05)

colSums(wgbsStats[sigInd,7:9] < 0.05)
table(rowSums(wgbsStats[sigInd,7:9] < 0.05))
table(rowSums(wgbsStats[sigInd,8:9] < 0.05))

## break into age and cell type
sigInd2 = which(p.adjust(outLibd$LIBD_pval,"bonf", 
	n=length(mm)) < 0.05)
mean(wgbsStats$p_Age[sigInd2] < 0.05)
mean(wgbsStats$p_Age[-sigInd2] < 0.05)

colSums(wgbsStats[sigInd2,7:9] < 0.05)
colSums(wgbsStats[sigInd2,7:9] < 0.01)
table(rowSums(wgbsStats[sigInd2,7:9] < 0.05))
table(rowSums(wgbsStats[sigInd2,8:9] < 0.05))