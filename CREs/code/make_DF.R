##

library(GenomicRanges)
library(jaffelab)

# load data
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/PMDs_methylSeekR.rda")

# make list
PMDsegments.CG = GRangesList(PMDsegments.CG)

# only keep autosomes
PMDsegments.CG = keepSeqlevels(PMDsegments.CG, 
	paste0("chr", c(1:22,"X","Y","M")))

PMDsegments.CG = endoapply(PMDsegments.CG , function(x) {
	mcols(x)$PMD = ifelse(mcols(x)$type =="PMD", TRUE, FALSE)
	return(x)
})

## merge together
pmdCovDf = lapply(seqlevels(PMDsegments.CG), function(chr) {
	cat(".")
	byChr = lapply(PMDsegments.CG, function(x) x[seqnames(x)==chr])
	stateByChr = lapply(byChr, function(x) Rle(x$PMD, width(x)))
	DataFrame(stateByChr)
})
names(pmdCovDf) = seqlevels(PMDsegments.CG)


######
## total lengths per state
stateStatsList = lapply(PMDsegments.CG, function(x) {
	xSplit = split(x, x$PMD)
	data.frame(numState = elementNROWS(xSplit),
		widthState = sapply(xSplit, function(y) 
			sum(as.numeric(width(y)))))
})
stateStats = do.call("cbind", stateStatsList)
theWidths = stateStats[,grep("width",names(stateStats))]/1e6
colnames(theWidths) = ss(colnames(theWidths), "\\.")
totWidth = t(round(100*theWidths / colSums(theWidths),2))

## get age
library(bsseq)
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/BSobj_subsets/rda/DMP_nonCpG_interaction.Rdata")
pd = pData(DMP_nonCpG)[rownames(totWidth),]
pd$AgeGroup = cut(pd$Age, c(0,1,10,30))
boxplot(totWidth[,"TRUE"] ~ pd$Cell.Type)
boxplot(totWidth[,"TRUE"] ~ pd$Cell.Type*pd$AgeGroup)

plot(totWidth[,"TRUE"] ~ pd$Age, pch = 21, bg =factor(pd$Cell.Type))

###############
# read DMRs ###
###############
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/DMR/DMR_objects.rda")
grDMR = GRangesList(lapply(DMR, makeGRangesFromDataFrame,keep=TRUE))
sigDMR = endoapply(grDMR, function(x) x[x$fwer < 0.05,])

## main chr
mainChr = paste0("chr", c(1:22,"X","Y"))

## cut PMDs around DMRs
# cell type
cellTypeByChr = split(sigDMR$CellType, seqnames(sigDMR$CellType))
bpOverlapByChr_cellType = mapply(function(regByChr, stateByChr) {
	cat(".")
	enr = stateByChr[ranges(regByChr),]
	colSums(as.data.frame(enr))
}, cellTypeByChr[mainChr], pmdCovDf[mainChr],SIMPLIFY=TRUE)
pmdMember_cellType = rowSums(bpOverlapByChr_cellType) / sum(width(sigDMR$CellType))
pmdMember_cellType
boxplot(pmdMember_cellType~ pd$Cell.Type)

# age
ageByChr = split(sigDMR$Age, seqnames(sigDMR$Age))
bpOverlapByChr_age = mapply(function(regByChr, stateByChr) {
	cat(".")
	enr = stateByChr[ranges(regByChr),]
	colSums(as.data.frame(enr))
}, ageByChr[mainChr], pmdCovDf[mainChr],SIMPLIFY=TRUE)
pmdMember_age = rowSums(bpOverlapByChr_age) / sum(width(sigDMR$Age))
pmdMember_age 
boxplot(pmdMember_age~ pd$Cell.Type)

# interaction
intByChr = split(sigDMR$Int, seqnames(sigDMR$Int))
bpOverlapByChr_int= mapply(function(regByChr, stateByChr) {
	cat(".")
	enr = stateByChr[ranges(regByChr),]
	colSums(as.data.frame(enr))
}, intByChr[mainChr], pmdCovDf[mainChr],SIMPLIFY=TRUE)
pmdMember_int = rowSums(bpOverlapByChr_int) / sum(width(sigDMR$Int))
pmdMember_int
boxplot(pmdMember_int~ pd$Cell.Type)
plot(pmdMember_int~ pd$Age)

### combine
pmdDf = data.frame(genome = totWidth[,"TRUE"],
	cellDMR = 100*pmdMember_cellType,
	ageDMR = 100*pmdMember_age,
	intDMR = 100*pmdMember_int)
pmdDf = round(pmdDf,2)

gIndexes = splitit(paste0(pd$Cell.Type,":", pd$AgeGroup))
round(t(sapply(gIndexes, function(ii) colMeans(pmdDf[ii,]))),2)