########
library('bsseq')
library('bumphunter')
library('doParallel')
library(limma)
library(rtracklayer)
library(jaffelab)

## load BSobj
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/processed_beta_values_plusMap.rda")

## atac-seq 
peaks = import("/dcl01/lieber/ajaffe/Amanda/ATAC/FinalPeaks/pooled_sample_peaks/mergedPeaks-shifted.bed")
peaks$name=NULL

## mean meth within each peak
oo = findOverlaps(peaks, gr)

rIndexes = split(subjectHits(oo), queryHits(oo))
meanMeth_inPeak= sapply(rIndexes, function(ii) 
	colMeans(t(t(meth[ii,]))))
meanMeth_inPeak = do.call("rbind", meanMeth_inPeak)

## get flanking regions 
bpShift = 500 
leftShift = shift(peaks,-1*bpShift)
end(leftShift) = start(peaks)-1
rightShift = shift(peaks, bpShift)
start(rightShift) = end(peaks)+1

## match up
ooLeft = findOverlaps(leftShift, gr)
ooRight = findOverlaps(rightShift, gr)
ooOut = rbind(as.matrix(ooLeft), as.matrix(ooRight))
rIndexesOut = split(ooOut[,"subjectHits"], ooOut[,"queryHits"])

## get mean meth
meanMeth_outPeak= sapply(rIndexesOut, function(ii) 
	colMeans(t(t(meth[ii,]))))
meanMeth_outPeak = do.call("rbind", meanMeth_outPeak)

## put in order
peakMeth_inPeak = peakMeth_outPeak = matrix(
	NA, nrow = length(peaks),ncol = ncol(meanMeth_inPeak))
colnames(peakMeth_inPeak) = colnames(peakMeth_outPeak) = colnames(meanMeth_inPeak)
peakMeth_inPeak[as.numeric(rownames(meanMeth_inPeak)),] = meanMeth_inPeak
peakMeth_outPeak[as.numeric(rownames(meanMeth_outPeak)),] = meanMeth_outPeak

## differences
peakDiff = rowMeans(peakMeth_inPeak - peakMeth_outPeak)

pdf("plots/atac_DNAm_inVsOut.pdf",w=12)
par(mar=c(5,6,2,2), cex.axis=2,cex.lab=2,cex.main=2)
hist(-1*peakDiff, xlab = "(Outside Peak) - (Inside Peak)",	
	main = "DNAm Differences",col="grey")
dev.off()

peaks_methDiffOrdered = peaks[order(peakDiff)]
save(peaks_methDiffOrdered, file="rdas/ATAC_peaks_methDiffOrdered.rda")

## by cell type
cIndexes = splitit(pd$Cell.Type)
peakDiff_cellType = sapply(cIndexes, function(ii) 
	rowMeans(peakMeth_inPeak[,ii] - peakMeth_outPeak[,ii]))