library('bsseq')
library('bumphunter')
library('doParallel')
library('limma')
library('rtracklayer')
library('jaffelab')
library("ggplot2")
library(reshape2)
library(data.table)

load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/non-CpG/CH_object.rda") # includes CH
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/DMR/DMR_objects.rda") # includes DMR, geneMap
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/ATAC/CpG.single.site.ATAC.coverageMatrix.rda") # includes ATACcovCpGs, ATACpd
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/processed_beta_values_plusMap.rda") # contains gr, meth, and pd

#For all DMRs, compare ATAC-seq coverage within DMR to 1kbp up- and down-stream

DMRgr = lapply(DMR, function(x) makeGRangesFromDataFrame(x, keep.extra.columns = T))
ATAC = as.matrix(ATACcovCpGs/colSums(ATACcovCpGs)*(10^6))
ATACgr = GRanges(rownames(ATAC))

# Calculate mean coverage within each DMR
oo = lapply(DMRgr, function(x) findOverlaps(x, ATACgr))
rIndexes = lapply(oo, function(x) split(subjectHits(x), queryHits(x)))

meanCov_inDMR= lapply(rIndexes, function(x) sapply(x, function(ii) colMeans(t(t(ATAC[ii,])))))
meanCov_inDMR = lapply(meanCov_inDMR, function(x) do.call("rbind", x))

## get flanking regions 
leftShift = lapply(DMR, function(x) data.frame(seqnames = x$seqnames, start = (x$start-1000), end = x$start-1))
rightShift = lapply(DMR, function(x) data.frame(seqnames = x$seqnames, start = (x$end+1), end = x$end+1000))
leftShift = lapply(leftShift, makeGRangesFromDataFrame)
rightShift = lapply(rightShift, makeGRangesFromDataFrame)

## match up
ooLeft = lapply(leftShift, function(x) findOverlaps(x, ATACgr))
ooRight = lapply(rightShift, function(x) findOverlaps(x, ATACgr))
ooOut = mapply(function(x,y) rbind(as.matrix(x), as.matrix(y)), ooLeft, ooRight)
rIndexesOut = lapply(ooOut, function(x) split(x[,"subjectHits"], x[,"queryHits"]))

## get mean coverage
meanCov_outDMR= lapply(rIndexesOut, function(x) lapply(x, function(ii) colMeans(t(t(ATAC[ii,])))))
meanCov_outDMR = lapply(meanCov_outDMR, function(x) do.call("rbind", x))
for (i in 1:length(meanCov_outDMR)){rownames(meanCov_inDMR[[i]]) = rownames(meanCov_outDMR[[i]]) = DMR[[i]][,"regionID"]}

## differences in all regions
covDiff = mapply(function(x,y) rowMeans(x - y), meanCov_inDMR, meanCov_outDMR)

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/ATAC/ATAC_at_CpGs/figures/DMR_CpG_ATAC_coverage_vs_outDMR_allregions.pdf",w=12)
par(mar=c(5,6,2,2), cex.axis=2,cex.lab=2,cex.main=2)
for (i in 1:3){
  p = hist(-1*covDiff[[i]], xlab = "(Outside Region) - (Inside Region)",	
	main = paste0("ATAC-seq Coverage Differences: ", names(covDiff)[i], ", All Regions"), col="grey")
  print(p)
}
dev.off()

## Only DMR fwer <= 0.05
covDiff.05 = mapply(function(x,y,z) rowMeans(x[which(rownames(x) %in% y[which(y$sig=="FWER < 0.05"),"regionID"]),] - 
                                            z[which(rownames(x) %in% y[which(y$sig=="FWER < 0.05"),"regionID"]),]), meanCov_inDMR, DMR, meanCov_outDMR)

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/ATAC/ATAC_at_CpGs/figures/DMR_CpG_ATAC_coverage_vs_outDMR_fwer_0.05.pdf",w=12)
par(mar=c(5,6,2,2), cex.axis=2,cex.lab=2,cex.main=2)
for (i in 1:3){
  p = hist(-1*covDiff[[i]], xlab = "(Outside Region) - (Inside Region)",	
           main = paste0("ATAC-seq Coverage Differences: ", names(covDiff)[i], ", FWER < 0.05"), col="grey")
  print(p)
}
dev.off()


### Break this down by cell type, age bins, significance and direction of methylation change
cIndexes = splitit(ATACpd$CellType)
covDiff_cellType = list(Neuron = mapply(function(x,y) rowMeans(x[,cIndexes$Neuron] - y[,cIndexes$Neuron]), meanCov_inDMR, meanCov_outDMR),
                        Glia = mapply(function(x,y) rowMeans(x[,cIndexes$Glia] - y[,cIndexes$Glia]), meanCov_inDMR, meanCov_outDMR))
covDiff_cellType = lapply(covDiff_cellType, function(x) lapply(x, function(y) data.frame(Difference = y)))
covDiff_cellType$Neuron = Map(cbind, covDiff_cellType$Neuron, regionID = lapply(covDiff_cellType$Neuron, rownames), CellType = "Neuron")
covDiff_cellType$Glia = Map(cbind, covDiff_cellType$Glia, regionID = lapply(covDiff_cellType$Glia, rownames), CellType = "Glia")
covDiff_cellType = list("CellType" = rbind(covDiff_cellType$Neuron$CellType, covDiff_cellType$Glia$CellType), 
                        "Age" = rbind(covDiff_cellType$Neuron$Age, covDiff_cellType$Glia$Age), 
                        "Interaction" = rbind(covDiff_cellType$Neuron$Interaction, covDiff_cellType$Glia$Interaction))
covDiff_cellType = Map(cbind, covDiff_cellType, 
                       Sig = mapply(function(x,y) ifelse(x$regionID %in% y[which(y$sig=="FWER < 0.05"),"regionID"], "FWER < 0.05", "FWER > 0.05"), covDiff_cellType, DMR),
                       Direction = mapply(function(x,y) ifelse(x$regionID %in% y[which(y$Dir=="pos"),"regionID"], "More Methylated\nIn Neurons", "More Methylated\nIn Glia"), covDiff_cellType, DMR))

aIndexes = splitit(ATACpd$Age.Bins)
covDiff_age = list(Neonate = mapply(function(x,y) rowMeans(x[,aIndexes$Neonate] - y[,aIndexes$Neonate]), meanCov_inDMR, meanCov_outDMR),
                   Toddler = mapply(function(x,y) rowMeans(x[,aIndexes$Toddler] - y[,aIndexes$Toddler]), meanCov_inDMR, meanCov_outDMR),
                   Child = mapply(function(x,y) rowMeans(x[,aIndexes$Child] - y[,aIndexes$Child]), meanCov_inDMR, meanCov_outDMR),
                   Early.Teen = mapply(function(x,y) rowMeans(x[,aIndexes$Early.Teen] - y[,aIndexes$Early.Teen]), meanCov_inDMR, meanCov_outDMR),
                   Teen = mapply(function(x,y) rowMeans(x[,aIndexes$Teen] - y[,aIndexes$Teen]), meanCov_inDMR, meanCov_outDMR),     
                   Young.Adult = mapply(function(x,y) rowMeans(x[,aIndexes$Young.Adult] - y[,aIndexes$Young.Adult]), meanCov_inDMR, meanCov_outDMR))
covDiff_age = lapply(covDiff_age, function(x) lapply(x, function(y) data.frame(Difference = y)))

covDiff_age$Neonate = Map(cbind, covDiff_age$Neonate, regionID = lapply(covDiff_age$Neonate, rownames), AgeBin = "Neonate")
covDiff_age$Toddler = Map(cbind, covDiff_age$Toddler, regionID = lapply(covDiff_age$Toddler, rownames), AgeBin = "Toddler")
covDiff_age$Child = Map(cbind, covDiff_age$Child, regionID = lapply(covDiff_age$Child, rownames), AgeBin = "Child")
covDiff_age$Early.Teen = Map(cbind, covDiff_age$Early.Teen, regionID = lapply(covDiff_age$Early.Teen, rownames), AgeBin = "Early Teen")
covDiff_age$Teen = Map(cbind, covDiff_age$Teen, regionID = lapply(covDiff_age$Teen, rownames), AgeBin = "Teen")
covDiff_age$Young.Adult = Map(cbind, covDiff_age$Young.Adult, regionID = lapply(covDiff_age$Young.Adult, rownames), AgeBin = "Young Adult")
covDiff_age = list("CellType" = rbind(covDiff_age$Neonate$CellType, covDiff_age$Toddler$CellType, covDiff_age$Child$CellType, covDiff_age$Early.Teen$CellType, covDiff_age$Teen$CellType, covDiff_age$Young.Adult$CellType), 
                   "Age" = rbind(covDiff_age$Neonate$Age, covDiff_age$Toddler$Age, covDiff_age$Child$Age, covDiff_age$Early.Teen$Age, covDiff_age$Teen$Age, covDiff_age$Young.Adult$Age), 
                   "Interaction" = rbind(covDiff_age$Neonate$Interaction, covDiff_age$Toddler$Interaction, covDiff_age$Child$Interaction, covDiff_age$Early.Teen$Interaction, covDiff_age$Teen$Interaction, covDiff_age$Young.Adult$Interaction))
covDiff_age = Map(cbind, covDiff_age, 
                  Sig = mapply(function(x,y) ifelse(x$regionID %in% y[which(y$sig=="FWER < 0.05"),"regionID"], "FWER < 0.05", "FWER > 0.05"), covDiff_age, DMR),
                  Direction = mapply(function(x,y) ifelse(x$regionID %in% y[which(y$Dir=="pos"),"regionID"], "Positive", "Negative"), covDiff_age, DMR))

ATACpd = ATACpd[,1:15]
ATACpd$Group = paste(ATACpd$CellType, ATACpd$Age.Bins, sep = "_")
iIndexes = splitit(ATACpd$Group)
covDiff_int = list("Neuron_Neonate" = mapply(function(x,y) rowMeans(x[,iIndexes$Neuron_Neonate] - y[,iIndexes$Neuron_Neonate]), meanCov_inDMR, meanCov_outDMR),
                   "Neuron_Toddler" = mapply(function(x,y) rowMeans(x[,iIndexes$Neuron_Toddler] - y[,iIndexes$Neuron_Toddler]), meanCov_inDMR, meanCov_outDMR),
                   "Neuron_Child" = mapply(function(x,y) rowMeans(x[,iIndexes$Neuron_Child] - y[,iIndexes$Neuron_Child]), meanCov_inDMR, meanCov_outDMR),
                   "Neuron_Early.Teen" = mapply(function(x,y) rowMeans(x[,iIndexes$Neuron_Early.Teen] - y[,iIndexes$Neuron_Early.Teen]), meanCov_inDMR, meanCov_outDMR),
                   "Neuron_Teen" = mapply(function(x,y) rowMeans(x[,iIndexes$Neuron_Teen] - y[,iIndexes$Neuron_Teen]), meanCov_inDMR, meanCov_outDMR),
                   "Neuron_Young.Adult" = mapply(function(x,y) rowMeans(x[,iIndexes$Neuron_Young.Adult] - y[,iIndexes$Neuron_Young.Adult]), meanCov_inDMR, meanCov_outDMR),
                   "Glia_Neonate" = mapply(function(x,y) rowMeans(x[,iIndexes$Glia_Neonate] - y[,iIndexes$Glia_Neonate]), meanCov_inDMR, meanCov_outDMR),
                   "Glia_Toddler" = mapply(function(x,y) (x[,iIndexes$Glia_Toddler] - y[,iIndexes$Glia_Toddler]), meanCov_inDMR, meanCov_outDMR),
                   "Glia_Child" = mapply(function(x,y) (x[,iIndexes$Glia_Child] - y[,iIndexes$Glia_Child]), meanCov_inDMR, meanCov_outDMR),
                   "Glia_Early.Teen" = mapply(function(x,y) rowMeans(x[,iIndexes$Glia_Early.Teen] - y[,iIndexes$Glia_Early.Teen]), meanCov_inDMR, meanCov_outDMR),
                   "Glia_Teen" = mapply(function(x,y) (x[,iIndexes$Glia_Teen] - y[,iIndexes$Glia_Teen]), meanCov_inDMR, meanCov_outDMR),
                   "Glia_Young.Adult" = mapply(function(x,y) rowMeans(x[,iIndexes$Glia_Young.Adult] - y[,iIndexes$Glia_Young.Adult]), meanCov_inDMR, meanCov_outDMR))
covDiff_int = lapply(covDiff_int, function(x) lapply(x, function(y) data.frame(Difference = y)))

covDiff_int$Neuron_Neonate = Map(cbind, covDiff_int$Neuron_Neonate, regionID = lapply(covDiff_int$Neuron_Neonate, rownames), "Group" = "Neuron:Neonate")
covDiff_int$Neuron_Toddler = Map(cbind, covDiff_int$Neuron_Toddler, regionID = lapply(covDiff_int$Neuron_Toddler, rownames), "Group" = "Neuron:Toddler")
covDiff_int$Neuron_Child = Map(cbind, covDiff_int$Neuron_Child, regionID = lapply(covDiff_int$Neuron_Child, rownames), "Group" = "Neuron:Child")
covDiff_int$Neuron_Early.Teen = Map(cbind, covDiff_int$Neuron_Early.Teen, regionID = lapply(covDiff_int$Neuron_Early.Teen, rownames), "Group" = "Neuron:Early Teen")
covDiff_int$Neuron_Teen = Map(cbind, covDiff_int$Neuron_Teen, regionID = lapply(covDiff_int$Neuron_Teen, rownames), "Group" = "Neuron:Teen")
covDiff_int$Neuron_Young.Adult = Map(cbind, covDiff_int$Neuron_Young.Adult, regionID = lapply(covDiff_int$Neuron_Young.Adult, rownames), "Group" = "Neuron:Young Adult")
covDiff_int$Glia_Neonate = Map(cbind, covDiff_int$Glia_Neonate, regionID = lapply(covDiff_int$Glia_Neonate, rownames), "Group" = "Glia:Neonate")
covDiff_int$Glia_Toddler = Map(cbind, covDiff_int$Glia_Toddler, regionID = lapply(covDiff_int$Glia_Toddler, rownames), "Group" = "Glia:Toddler")
covDiff_int$Glia_Child = Map(cbind, covDiff_int$Glia_Child, regionID = lapply(covDiff_int$Glia_Child, rownames), "Group" = "Glia:Child")
covDiff_int$Glia_Early.Teen = Map(cbind, covDiff_int$Glia_Early.Teen, regionID = lapply(covDiff_int$Glia_Early.Teen, rownames), "Group" = "Glia:Early Teen")
covDiff_int$Glia_Teen = Map(cbind, covDiff_int$Glia_Teen, regionID = lapply(covDiff_int$Glia_Teen, rownames), "Group" = "Glia:Teen")
covDiff_int$Glia_Young.Adult = Map(cbind, covDiff_int$Glia_Young.Adult, regionID = lapply(covDiff_int$Glia_Young.Adult, rownames), "Group" = "Glia:Young Adult")
covDiff_int = list("CellType" = rbind(covDiff_int$Neuron_Neonate$CellType, covDiff_int$Neuron_Toddler$CellType, covDiff_int$Neuron_Child$CellType, covDiff_int$Neuron_Early.Teen$CellType, covDiff_int$Neuron_Teen$CellType, covDiff_int$Neuron_Young.Adult$CellType,
                                      covDiff_int$Glia_Neonate$CellType, covDiff_int$Glia_Toddler$CellType, covDiff_int$Glia_Child$CellType, covDiff_int$Glia_Early.Teen$CellType, covDiff_int$Glia_Teen$CellType, covDiff_int$Glia_Young.Adult$CellType), 
                   "Age" = rbind(covDiff_int$Neuron_Neonate$Age, covDiff_int$Neuron_Toddler$Age, covDiff_int$Neuron_Child$Age, covDiff_int$Neuron_Early.Teen$Age, covDiff_int$Neuron_Teen$Age, covDiff_int$Neuron_Young.Adult$Age,
                                 covDiff_int$Glia_Neonate$Age, covDiff_int$Glia_Toddler$Age, covDiff_int$Glia_Child$Age, covDiff_int$Glia_Early.Teen$Age, covDiff_int$Glia_Teen$Age, covDiff_int$Glia_Young.Adult$Age), 
                   "Interaction" = rbind(covDiff_int$Neuron_Neonate$Interaction, covDiff_int$Neuron_Toddler$Interaction, covDiff_int$Neuron_Child$Interaction, covDiff_int$Neuron_Early.Teen$Interaction, covDiff_int$Neuron_Teen$Interaction, covDiff_int$Neuron_Young.Adult$Interaction,
                                         covDiff_int$Glia_Neonate$Interaction, covDiff_int$Glia_Toddler$Interaction, covDiff_int$Glia_Child$Interaction, covDiff_int$Glia_Early.Teen$Interaction, covDiff_int$Glia_Teen$Interaction, covDiff_int$Glia_Young.Adult$Interaction))
covDiff_int = Map(cbind, covDiff_int, Sig = mapply(function(x,y) ifelse(rownames(x) %in% y[which(y$sig=="FWER < 0.05"),"regionID"], "FWER < 0.05", "FWER > 0.05"), covDiff_int, DMR),
                  Direction = mapply(function(x,y) ifelse(rownames(x) %in% y[which(y$Dir=="pos"),"regionID"], "Positive", "Negative"), covDiff_int, DMR))

# Save objects
save(covDiff, covDiff_cellType, covDiff_age, covDiff_int, file="/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/ATAC/DMR_ATACcovDiff_objects.rda")

covDiff_cellType = Map(cbind, covDiff_cellType, direction = lapply(covDiff_cellType, function(x) ifelse(x$Direction=="Positive", "More Methylated\nIn Neurons", "More Methylated\nIn Glia")))


## Plot density of distance by cell type, age, direction of methylation change and whether the DMR is significant or not

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/ATAC/ATAC_at_CpGs/figures/DMR_CpG_ATAC_coverage_vs_outDMR_byCellType_redo.pdf",w=12)
par(mar=c(5,6,2,2), cex.axis=2,cex.lab=2,cex.main=2)
for (i in 1:3){
  p = ggplot(covDiff_cellType[[i]], aes(x=Difference)) + geom_density(aes(group=CellType, colour=CellType)) + 
    facet_grid(Direction ~ Sig) +
    ylab("Density") + ylim(0,30) +
    xlab("(Inside Region) - (Outside Region)") + 
    ggtitle(paste0("Mean ATAC-seq Coverage Difference at CpGs\nBetween DMR and Flanking Regions: " , names(covDiff_cellType)[i], " Regions")) + 
    theme(title = element_text(size = 20)) +
    theme(text = element_text(size = 20))
  print(p)
}
  dev.off()

######## THIS ONE ########

######## HERE ######

x = Map(cbind, covDiff_cellType, Sig = mapply(function(x,y) ifelse(rownames(x) %in% y[which(y$sig=="FWER < 0.05"),"regionID"], "FWER < 0.05", "FWER > 0.05"), covDiff_cellType, DMR), covDiff_cellType, DMR)
x = Map(cbind, x, CT = lapply(x, function(y) factor(x$CellType, levels = c("Neuron","Glia"))))


x = covDiff_cellType$CellType
x$regionID = rownames(x)
x[x$CellType=="Glia","regionID"] = substr(x[x$CellType=="Glia","regionID"],1,nchar(x[x$CellType=="Glia","regionID"])-1)
cell = DMR$CellType
x$sig = ifelse(x$regionID %in% cell[which(cell$sig=="FWER < 0.05"),"regionID"], "FWER < 0.05","FWER > 0.05")
x$CellType = factor(x$CellType, levels = c("Glia", "Neuron"))
x$dir = ifelse(x$regionID %in% cell[which(cell$Dir=="pos"),"regionID"], "More Methylated\nIn Neurons", "More Methylated\nIn Glia")


pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/ATAC/ATAC_at_CpGs/figures/DMR_CpG_ATAC_coverage_vs_outDMR_byCellType_fwer.0.05_redo.pdf",w=12)
ggplot(x[which(x$sig=="FWER < 0.05"),], aes(x=Difference)) + geom_density(aes(group=CellType, colour=CellType), size = 2) + 
    facet_grid(. ~ dir) + scale_color_brewer(palette = "Dark2") +
    theme_classic() +
    ylab("Density") + ylim(0,30) + 
    xlab("(Inside Region) - (Outside Region)") + xlim(-0.5,0.5) +
    ggtitle("Mean ATAC-seq Coverage Difference at CpGs\nBetween DMR and Flanking Regions: Cell Type Regions") + 
    theme(title = element_text(size = 20)) +
    theme(text = element_text(size = 20), legend.title=element_blank()) + theme(legend.position="bottom")
dev.off()
  
  

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/ATAC/ATAC_at_CpGs/figures/DMR_CpG_ATAC_coverage_vs_outDMR_byAge.pdf",w=12)
par(mar=c(5,6,2,2), cex.axis=2,cex.lab=2,cex.main=2)
for (i in 1:3){
    p = ggplot(covDiff_age[[i]], aes(x=Difference)) + geom_density(aes(group=AgeBin, colour=AgeBin)) + 
      facet_grid(Direction ~ Sig) +
      ylab("Density") + ylim(0,30) +
      xlab("(Inside Region) - (Outside Region)") +
      ggtitle(paste0("Mean ATAC-seq Coverage Difference at CpGs\nBetween DMR and Flanking Regions: " , names(covDiff_age)[i], " Regions")) + 
      theme(title = element_text(size = 20)) +
      theme(text = element_text(size = 20))
    print(p)
  }
dev.off()
  
pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/ATAC/ATAC_at_CpGs/figures/DMR_CpG_ATAC_coverage_vs_outDMR_byInteraction.pdf",w=12)
par(mar=c(5,6,2,2), cex.axis=2,cex.lab=2,cex.main=2)
for (i in 1:3){
  p = ggplot(covDiff_int[[i]], aes(x=Difference)) + geom_density(aes(group=Group, colour=Group)) + 
      facet_grid(Direction ~ Sig) +
      ylab("Density") + ylim(0,30) +
      xlab("(Inside Region) - (Outside Region)") +
      ggtitle(paste0("Mean ATAC-seq Coverage Difference at CpGs\nBetween DMR and Flanking Regions: " , names(covDiff)[i], " Regions")) + 
      theme(title = element_text(size = 20)) +
      theme(text = element_text(size = 20))
  print(p)
}
for (i in 1:3){
  p = ggplot(covDiff_int[[i]], aes(x=Difference)) + geom_density(aes(group=Group, colour=Group)) + 
    facet_grid(. ~ Sig) +
    ylab("Density") + ylim(0,30) +
    xlab("(Inside Region) - (Outside Region)") +
    ggtitle(paste0("Mean ATAC-seq Coverage Difference at CpGs\nBetween DMR and Flanking Regions: " , names(covDiff)[i], " Regions")) + 
    theme(title = element_text(size = 20)) +
    theme(text = element_text(size = 20))
  print(p)
}
dev.off()


### Plot each DMR as a scatterplot, with each point being a sample, one axis being mean DNAm, and the other axis being mean ATAC coverage
## Make list of data frames containing methylation and library normalized ATAC coverage for each sample at each DMR
save(ATAC,meth,gr,ATACpd,pd,DMR, file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/ATAC/DMR_ATAC_plot_objects.rda")

ATAC = as.matrix(ATACcovCpGs/colSums(ATACcovCpGs)*(10^6))

sig = lapply(DMR, function(x) x[which(x$fwer<=0.05),])
sig = lapply(sig, function(x) x[order(x$fwer),])

sigDMRgr = lapply(sig, function(x) makeGRangesFromDataFrame(x, keep.extra.columns = T))
ATACgr = GRanges(rownames(ATAC))

oo = findOverlaps(gr, ATACgr)
subsetATACgr = ATACgr[subjectHits(oo)]
subsetgr = gr[queryHits(oo)]
subsetATAC = ATAC[subjectHits(oo),]
subsetmeth = meth[queryHits(oo),]

cov = cbind(subsetATAC, subsetmeth)
oo.DMRs = lapply(sigDMRgr, function(x) findOverlaps(x, subsetgr))
rIndexes.DMRs = lapply(oo.DMRs, function(x) split(subjectHits(x), queryHits(x)))
covDMRs = lapply(rIndexes.DMRs, function(x) lapply(x, function(ii) cov[ii,]))
covDMRs = lapply(covDMRs, function(x) lapply(x, melt))
CellType = mapply(function(x,y) c(x,y), split(as.character(ATACpd$SampleID),ATACpd$CellType),
                 split(pd$Data.ID, pd$Cell.Type))
AgeBin = mapply(function(x,y) c(x,y), split(as.character(ATACpd$SampleID), ATACpd$Age.Bin),
                split(pd$Data.ID, pd$Age.Bin))
keys = unique(c(ATACpd$BrainID,pd$Brain.ID))
ID = setNames(mapply(c, split(as.character(ATACpd$SampleID), ATACpd$BrainID)[keys], split(pd$Data.ID, pd$Brain.ID)[keys]), keys)

for (i in 1:length(sig)){ 
  names(covDMRs[[i]]) = sig[[i]][,"regionID"]
  for (j in 1:nrow(sig[[i]])){
    tmp = covDMRs[[i]][[j]]
    tmp[grep("WGC", tmp$Var2),"Experiment"] = "WGBS"
    tmp[-grep("WGC", tmp$Var2),"Experiment"] = "ATAC"
    tmp[which(tmp$Var2 %in% CellType$Neuron),"CellType"] = "Neuron"
    tmp[which(tmp$Var2 %in% CellType$Glia),"CellType"] = "Glia"
    tmp[which(tmp$Var2 %in% AgeBin$Child),"AgeBin"] = "Child"
    tmp[which(tmp$Var2 %in% AgeBin$Early.Teen),"AgeBin"] = "Early.Teen"
    tmp[which(tmp$Var2 %in% AgeBin$Neonate),"AgeBin"] = "Neonate"
    tmp[which(tmp$Var2 %in% AgeBin$Teen),"AgeBin"] = "Teen"
    tmp[which(tmp$Var2 %in% AgeBin$Toddler),"AgeBin"] = "Toddler"
    tmp[which(tmp$Var2 %in% AgeBin$Young.Adult),"AgeBin"] = "Young.Adult"
    tmp$Group = paste(tmp$Experiment, tmp$CellType, tmp$AgeBin, sep=":")
    tmp$CTGroup = paste(tmp$Experiment, tmp$CellType, sep=":")
    tmp$AgeGroup = paste(tmp$Experiment, tmp$AgeBin, sep=":")
    covDMRs[[i]][[j]] = tmp
  }}
covDMRs = mapply(function(x,y) Map(cbind, x, DMR.Direction = as.list(y$Dir)), covDMRs, sig)

# take means of groups before calculating correlation
lapply(covDMRs, function(x) head(lapply(x,class)))
covDMRs_dt = lapply(covDMRs, function(x) lapply(x,as.data.table))
covDMRs_dt = lapply(covDMRs_dt, function(x) Map(cbind, x, "CpG" = lapply(x, function(z) z[,1])))

means = list(allMean = lapply(covDMRs_dt, function(x) lapply(x, function(y) y[,mean(value), by=c("Var1","Experiment")])),
             byCT = lapply(covDMRs_dt, function(x) lapply(x, function(y) y[,mean(value), by=c("Var1","Experiment","CellType", "CTGroup")])),
             byAge = lapply(covDMRs_dt, function(x) lapply(x, function(y) y[,mean(value), by=c("Var1","Experiment","AgeBin", "AgeGroup")])))

celltype = covDMRs_dt[["CellType"]]
a = sapply(celltype, function(x) "Var1" %in% colnames(x))



head(lapply(celltype, function(y) y[,mean(value), by=c("Var1","Experiment")]))




celltype = Map(cbind, celltype, "CpG" = lapply(celltype, function(x) x[,1]))

celltypebyAll = lapply(celltypebyAll, function(x) split(x, x$Experiment))
celltypebyAll = lapply(celltypebyAll, function(y) cbind(y$ATAC[match(y$ATAC$CpG, y$WGBS$CpG),,],
                                                      y$WGBS[match(y$ATAC$CpG, y$WGBS$CpG),,]))
celltypebyAll = lapply(celltypebyAll, as.data.frame)
celltypebyCT = lapply(celltypebyCT, function(x) split(x, x$Experiment))
celltypebyCT = lapply(celltypebyCT, function(y) cbind(y$ATAC[match(y$ATAC$CpG, y$WGBS$CpG),,],
                                                      y$WGBS[match(y$ATAC$CpG, y$WGBS$CpG),,]))
celltypebyCT = lapply(celltypebyCT, as.data.frame)
for (i in 1:length(celltypebyCT)) {colnames(celltypebyCT[[i]]) = c("CpG","Experiment","CellType","CTGroup","ATACmean","CpG","Experiment","CellType","CTGroup","WGBSmean")}

celltypebyAge = lapply(celltypebyAge, function(x) split(x, x$Experiment))
celltypebyAge = lapply(celltypebyAge, function(y) cbind(y$ATAC[match(y$ATAC$CpG, y$WGBS$CpG),,],
                                                        y$WGBS[match(y$ATAC$CpG, y$WGBS$CpG),,]))
celltypebyAge = lapply(celltypebyAge, as.data.frame)
for (i in 1:length(celltypebyAge)) {colnames(celltypebyAge[[i]]) = c("CpG","Experiment","CellType","CTGroup","ATACmean","CpG","Experiment","CellType","CTGroup","WGBSmean")}


toy = lapply(celltypebyCT[1:3],head)
corr.CellType = list(byCT = lapply(toy, function(x) cor(x = x$WGBSmean, y = x$ATACmean)),
                     byAge = lapply(celltypebyAge, function(x) cor(x = x$WGBSmean, y = x$ATACmean, use="complete.obs")))


byCT = lapply(celltypebyCT, function(x) cor(x = x$WGBSmean, y = x$ATACmean))


## Plot the top 100 regions for each comparison

top = lapply(sig, function(x) x[1:100,"regionID"])
top.celltypebyCT = celltypebyCT[which(names(celltypebyCT) %in% top$CellType)]
top.celltypebyAge = celltypebyAge[which(names(celltypebyAge) %in% top$CellType)]

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/ATAC/ATAC_at_CpGs/figures/meth_vs.ATAC_inDMRs_CellType.pdf")
for (i in 1:2){
  g = ggplot(top.celltypebyCT[[i]],  aes(x = ATACmean, y =WGBSmean)) + geom_point(aes(colour = CellType)) +
    geom_smooth(method = "lm", fullrange=TRUE) +
    ylab("Mean Beta") + xlab("Mean ATAC Coverage/Million Reads") +
    ggtitle(paste0(names(celltype)[1],": Methylation vs Accessibility")) +
    theme(title = element_text(size = 20)) +
    theme(text = element_text(size = 20))
  print(g)
  g = ggplot(top.celltypebyAge[[i]],  aes(x = ATACmean, y =WGBSmean)) + geom_point(aes(colour = AgeBin)) +
    geom_smooth(method = "lm", fullrange=TRUE) +
    ylab("Mean Beta") + xlab("Mean ATAC Coverage/Million Reads") +
    ggtitle(paste0(names(celltype)[i],": Methylation vs Accessibility")) +
    theme(title = element_text(size = 20)) +
    theme(text = element_text(size = 20))
  print(g)
}














