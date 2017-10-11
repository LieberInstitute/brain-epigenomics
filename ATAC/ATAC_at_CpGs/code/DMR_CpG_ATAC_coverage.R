library('bsseq')
library('bumphunter')
library('doParallel')
library('limma')
library('rtracklayer')
library('jaffelab')

load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/non-CpG/CH_object.rda")
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/DMR/DMR_objects.rda")
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/ATAC/CpG.single.site.ATAC.coverageMatrix.rda")

#For all DMRs, compare ATAC-seq coverage within DMR to 1kbp up- and down-stream

DMRgr = lapply(DMR, function(x) makeGRangesFromDataFrame(x, keep.extra.columns = T))
ATAC = as.matrix(ATACcovCpGs/colSums(ATACcovCpGs)/(10^6))
ATACgr = GRanges(rownames(ATAC))

# Calculate mean coverage within each DMR
oo = lapply(DMRgr, function(x) findOverlaps(x, ATACgr))
rIndexes = lapply(oo, function(x) split(subjectHits(x), queryHits(x)))
meanCov_inDMR= list(CellType = matrix(nrow = length(rIndexes$CellType), ncol = length(colnames(ATAC)), dimnames = list(DMR$CellType$regionID, colnames(ATAC))),
                    Age = matrix(nrow = length(rIndexes$Age), ncol = length(colnames(ATAC)), dimnames = list(DMR$Age$regionID, colnames(ATAC))),
                    Interaction = matrix(nrow = length(rIndexes$Interaction), ncol = length(colnames(ATAC)), dimnames = list(DMR$Interaction$regionID, colnames(ATAC))))
for (i in 1:length(meanCov_inDMR)){
  tmp = meanCov_inDMR[[i]]
  for (j in 1:nrow(tmp)){
    row = rIndexes[[i]][[j]]
    for (k in 1:ncol(tmp)){
      tmp[j,k] = mean(ATAC[row,k])
      meanCov_inDMR[[i]] = tmp
    }
  }
}  

## get flanking regions 
bpShift = 1000 
leftShift = lapply(DMRgr, function(x) shift(x,-1*bpShift))
end(leftShift) = lapply(DMRgr, function(x) start(x)-1)
rightShift = lapply(DMRgr, function(x) shift(x, bpShift))
start(rightShift) = lapply(DMRgr, function(x) end(x)+1)

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

#Subset DNAm to set in ATAC-seq, then plot mean coverage in ATAC-seq vs mean methylation for a given DMR at a time. Take into account library size. 
#`peakMeth_inPeak` vs `peakMeth_outPeak` from the previous code (something like that)

## load BSobj
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/processed_beta_values_plusMap.rda")

