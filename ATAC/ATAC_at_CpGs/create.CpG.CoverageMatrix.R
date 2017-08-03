# load matrix of peaks as rows and columns as samples
samps = scan("/dcl01/lieber/ajaffe/Amanda/ATAC/FinalPeaks/scripts/finalizedList-all.txt", what = "character")

covMat = list()
for (i in 1:length(samps)){
  tmp = samps[i]
  covMat[[i]] = read.table(paste0("/dcl01/lieber/ajaffe/Amanda/ATAC/Coverage/coverage_at_CpGs/",
                                  samps[i], "_CpGs_coverageBed.txt"),
                           col.names = c("chr","start","end", "ReadName", "ReadCount"))
}
coord = lapply(covMat, function(x) paste0(x$chr, ":", x$start, "-", x$end))
covMat = Map(cbind, covMat, coord = coord)
names(covMat) = samps
ReadCount = lapply(covMat, function(x) x$ReadCount)
CpGcov = as.data.frame(ReadCount, row.names = coord[[1]])
CpGcov = CpGcov[,order(colnames(CpGcov))]
pd = read.table("/dcl01/lieber/ajaffe/Amanda/ATAC/pd.txt", header=TRUE)
rownames(pd) = pd$ID

save(CpGcov, pd, file="/dcl01/lieber/ajaffe/Amanda/ATAC/RDAs/CpG.single.site.coverageMatrix.rda")
