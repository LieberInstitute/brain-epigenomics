library(data.table)
library(ggplot2)
library(plyr)
library('bumphunter')
library(GenomicRanges)
library('jaffelab')

load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/non-CpG/CH_object.rda")
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/non-CpG/CHneurons_object.rda")
load('/dcl01/lieber/ajaffe/lab/brain-epigenomics/bsseq/bsobj_by_chr/allChrs_postNatal_cleaned_nonCG_noHomogenate_highCov.Rdata')


table(CH$CT.dir=="pos")
# FALSE     TRUE 
# 7630202 33188540 
# 18.7% of CH is greater methylated in glia

table(CH[which(CH$padj.CellType<=0.05),"CT.dir"]=="pos")
#  FALSE    TRUE 
#  79239 7602836
# 99% of all significantly regulated sites by cell type are more methylated in neurons

table(CH$Age.dir=="pos")
#FALSE     TRUE 
#9843075 30975667
#75.9% are increasing over age

table(CH[which(CH$padj.Age<=0.05),"CT.dir"]=="pos")
#FALSE    TRUE 
#24967 3169651 
# 99.2% are increasing over age in significantly dmCH sites.


## Check coverage by CH context

# load data

rm(list= ls()[!(ls() %in% c('CGlist','CGPrenlist'))])
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/UMRs_LMRs_annotated.rda")
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/dmvs_annotated.rda")
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/PMDs_annotated.rda")

cov = getCoverage(BSobj, type = "Cov")
dim(cov)
gr = granges(BSobj)

## collate regions

all = list(LMR.All = Reduce(intersect, lapply(lDMR, makeGRangesFromDataFrame)),
           UMR.All = Reduce(intersect, lapply(uDMR, makeGRangesFromDataFrame)),
           DMV.All = Reduce(intersect, lapply(dmvDMR, makeGRangesFromDataFrame)),
           PMD.All = Reduce(intersect, lapply(pmdDMR, makeGRangesFromDataFrame)))
tab = list(LMR.Prenatal = lDMR[-which(names(lDMR) %in% postnatalpd$Data.ID)],
           LMR.Neuron = lDMR[which(names(lDMR) %in% postnatalpd[postnatalpd$Cell.Type=="Neuron","Data.ID"])],
           LMR.Glia = lDMR[which(names(lDMR) %in% postnatalpd[postnatalpd$Cell.Type=="Glia","Data.ID"])],
           UMR.Prenatal = uDMR[-which(names(uDMR) %in% postnatalpd$Data.ID)],
           UMR.Neuron = uDMR[which(names(uDMR) %in% postnatalpd[postnatalpd$Cell.Type=="Neuron","Data.ID"])],
           UMR.Glia = uDMR[which(names(uDMR) %in% postnatalpd[postnatalpd$Cell.Type=="Glia","Data.ID"])],
           DMV.Prenatal = dmvDMR[-which(names(dmvDMR) %in% postnatalpd$Data.ID)],
           DMV.Neuron = dmvDMR[which(names(dmvDMR) %in% postnatalpd[postnatalpd$Cell.Type=="Neuron","Data.ID"])],
           DMV.Glia = dmvDMR[which(names(dmvDMR) %in% postnatalpd[postnatalpd$Cell.Type=="Glia","Data.ID"])],
           PMD.Prenatal = pmdDMR[-which(names(pmdDMR) %in% postnatalpd$Data.ID)],
           PMD.Neuron = pmdDMR[which(names(pmdDMR) %in% postnatalpd[postnatalpd$Cell.Type=="Neuron","Data.ID"])],
           PMD.Glia = pmdDMR[which(names(pmdDMR) %in% postnatalpd[postnatalpd$Cell.Type=="Glia","Data.ID"])])
tab = lapply(tab, function(y) lapply(y, makeGRangesFromDataFrame, keep=TRUE))
red = lapply(tab, function(x) reduce(do.call(getMethod(c, "GenomicRanges"), GRangesList(x))))
tab = c(all, list(LMR.Prenatal = setdiff(red$LMR.Prenatal,all$LMR.All), 
                  LMR.Neuron = setdiff(red$LMR.Neuron,all$LMR.All),
                  LMR.Glia = setdiff(red$LMR.Glia,all$LMR.All), 
                  UMR.Prenatal = setdiff(red$UMR.Prenatal,all$UMR.All),
                  UMR.Neuron = setdiff(red$UMR.Neuron,all$UMR.All), 
                  UMR.Glia = setdiff(red$UMR.Glia,all$UMR.All),
                  DMV.Prenatal = setdiff(red$DMV.Prenatal,all$DMV.All), 
                  DMV.Neuron = setdiff(red$DMV.Neuron,all$DMV.All),
                  DMV.Glia = setdiff(red$DMV.Glia,all$DMV.All), 
                  PMD.Prenatal = setdiff(red$PMD.Prenatal,all$PMD.All),
                  PMD.Neuron = setdiff(red$PMD.Neuron,all$PMD.All), 
                  PMD.Glia = setdiff(red$PMD.Glia,all$PMD.All)))

## get mean coverage in regions

ov = lapply(tab, function(x) findOverlaps(gr, x))
rIndexesIn = lapply(ov, function(x) split(queryHits(x), subjectHits(x)))
tabIDs = lapply(tab, function(x) paste0(data.frame(x)$seqnames, ":", data.frame(x)$start, "-", data.frame(x)$end))
for (i in 1:length(rIndexesIn)) { names(rIndexesIn[[i]]) = tabIDs[[i]][as.numeric(names(rIndexesIn[[i]]))] }

meanCov = lapply(rIndexesIn, function(x) sapply(x, function(ii) colMeans(t(t(cov[ii,])))))
meanCov = lapply(meanCov, function(x) if (class(x)=="list") { do.call("rbind", x) } else { t(x) } )
meanCov = lapply(meanCov, function(x) as.data.frame(t(data.frame(x))))
meanCov = Map(cbind, lapply(meanCov, function(x) cbind(x, id = rownames(x))), Group = as.list(names(meanCov)))
for (i in 1:length(meanCov)) {
  meanCov[[i]][which(meanCov[[i]][,"id"] %in% pd[pd$Cell.Type=="Neuron","Data.ID"]),"celltype"] = "Neuron"
  meanCov[[i]][which(meanCov[[i]][,"id"] %in% pd[pd$Cell.Type=="Glia","Data.ID"]),"celltype"] = "Glia"
  meanCov[[i]][,"Age"] = pd[match(meanCov[[i]][,"id"], pd$Data.ID), "Age"]
}

trespostnatal = lapply(meanCov, function(x) t.test(x = rowMeans(x[which(x$celltype=="Neuron"),c(1:(ncol(x)-4))]),
                                                   y = rowMeans(x[which(x$celltype=="Glia"),c(1:(ncol(x)-4))])))
write.csv(data.frame(Regions = names(trespostnatal), 
                     Neuron.Glia.Tstat = unlist(lapply(trespostnatal, function(x) x$statistic)),
                     Neuron.Glia.pval = unlist(lapply(trespostnatal, function(x) x$p.value)), 
                     Neuron.mean = unlist(lapply(trespostnatal, function(x) x$estimate[1])), 
                     Glia.mean = unlist(lapply(trespostnatal, function(x) x$estimate[2]))), quote = F,
          file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/cellType_CRE_coverage_tstats_nonCpG.csv")

df = do.call(rbind, Map(cbind, lapply(meanCov, function(x) 
  rbind(data.frame(rowmeans = rowMeans(x[which(x$celltype=="Neuron"),c(1:(ncol(x)-4))]), celltype = "Neuron"),
        data.frame(rowmeans = rowMeans(x[which(x$celltype=="Glia"),c(1:(ncol(x)-4))]), celltype = "Glia"))), 
  Regions = as.list(names(meanCov))))

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/CREs/figures/nonCpG_postnatal_coverage_inCREs.pdf", width = 26)
ggplot(df, aes(x = celltype, y = rowmeans, colour = celltype)) + 
  theme_classic() + geom_boxplot() +
  labs(fill="") + facet_grid(. ~ Regions) +
  ylab("Mean Coverage") + 
  xlab("") + theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  ggtitle("Mean nonCpG Coverage by Cell Type in CRE sequence") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20), legend.title=element_blank()) + theme(legend.position="bottom")
dev.off()

for (i in 1:length(meanCov)) {
  meanCov[[i]] = cbind(meanCov[[i]], Age.bin = pd[match(meanCov[[i]][,"id"], pd$Data.ID),"Age.Bin"])
  meanCov[[i]][,"Age.bin"] = as.character(meanCov[[i]][,"Age.bin"])
}

df = do.call(rbind, Map(cbind, lapply(meanCov, function(x) 
  rbind(data.frame(rowmeans = rowMeans(x[which(x$Age.bin=="Neonate"),c(1:(ncol(x)-5))]), Age.bin = "Neonate"),
        data.frame(rowmeans = rowMeans(x[which(x$Age.bin=="Toddler"),c(1:(ncol(x)-5))]), Age.bin = "Toddler"),
        data.frame(rowmeans = rowMeans(x[which(x$Age.bin=="Child"),c(1:(ncol(x)-5))]), Age.bin = "Child"),
        data.frame(rowmeans = rowMeans(x[which(x$Age.bin=="Early.Teen"),c(1:(ncol(x)-5))]), Age.bin = "Early.Teen"),
        data.frame(rowmeans = rowMeans(x[which(x$Age.bin=="Teen"),c(1:(ncol(x)-5))]), Age.bin = "Teen"),
        data.frame(rowmeans = rowMeans(x[which(x$Age.bin=="Young.Adult"),c(1:(ncol(x)-5))]), Age.bin = "Young.Adult"))), 
  Regions = as.list(names(meanCov))))


pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/CREs/figures/nonCpG_acrossAge_coverage_inCREs.pdf", width = 52)
ggplot(df, aes(x = Age.bin, y = rowmeans, colour = Age.bin)) + 
  theme_classic() + geom_boxplot() +
  labs(fill="") + facet_grid(. ~ Regions) +
  ylab("Mean Coverage") + 
  xlab("") + theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  ggtitle("Mean Coverage by Age Bin in CRE sequence") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20), legend.title=element_blank()) + theme(legend.position="bottom")
dev.off()


#### Are CREs depleted for mCH compared to surrounding regions, and is this affected by cell type or age? 

## Get methylation
meth = getMeth(BSobj, type = "raw")
dim(meth)
gr = granges(BSobj)

## get mean methylation
meanMeth = lapply(rIndexesIn, function(x) sapply(x, function(ii) colMeans(t(t(meth[ii,])))))
meanMeth = lapply(meanMeth, function(x) if (class(x)=="list") { do.call("rbind", x) } else { t(x) } )

## get flanking regions 
dftab = lapply(tab, as.data.frame)
leftShift = lapply(dftab, function(x) makeGRangesFromDataFrame(data.frame(seqnames = x$seqnames, start = (x$start-1000), end = x$start-1)))
rightShift = lapply(dftab, function(x) makeGRangesFromDataFrame(data.frame(seqnames = x$seqnames, start = (x$end+1), end = x$end+1000)))

## match up
ooLeft = lapply(leftShift, function(x) findOverlaps(x, gr))
ooRight = lapply(rightShift, function(x) findOverlaps(x, gr))
ooOut = mapply(function(x,y) rbind(as.matrix(x), as.matrix(y)), ooLeft, ooRight)
rIndexesOut = lapply(ooOut, function(x) split(x[,"subjectHits"], x[,"queryHits"]))
for (i in 1:length(rIndexesOut)) { names(rIndexesOut[[i]]) = tabIDs[[i]][as.numeric(names(rIndexesOut[[i]]))] }

## get mean coverage
outCRE = lapply(rIndexesOut, function(x) lapply(x, function(ii) colMeans(t(t(meth[ii,])))))
outCRE = lapply(outCRE, function(x) do.call("rbind", x))

## differences in all regions

outCRE = mapply(function(x,y) x[which(rownames(x) %in% rownames(y)),], outCRE, meanMeth)
meanMeth = mapply(function(x,y) y[which(rownames(y) %in% rownames(x)),], outCRE, meanMeth)
df = mapply(function(x,y) data.frame(methDiff = -1*rowMeans(x - y)), meanMeth, outCRE, SIMPLIFY = F)
tres = mapply(function(inn, out) t.test(rowMeans(inn), rowMeans(out)), meanMeth, outCRE, SIMPLIFY = F)
write.csv(data.frame(Group = names(tres), Tstat = unlist(lapply(tres, function(x) x$statistic)),
                     pval = unlist(lapply(tres, function(x) x$p.value)), 
                     In.mean = unlist(lapply(tres, function(x) x$estimate[1])), 
                     Out.mean = unlist(lapply(tres, function(x) x$estimate[2]))), quote = F,
          file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/nonCpG_meth_in.VS.out_CREs_tstats.csv")

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/non-CpG/figures/nonCpG_meth_in.VS.out_CREs.pdf",w=12)
for (i in 1:length(df)) { 
  g = ggplot(df, aes(methDiff)) + theme_classic() + geom_histogram() +
  labs(fill="") +
  ylab("") + 
  xlab("(Outside CRE) - (Inside CRE)") + xlim(-0.4,0.4) +
  ggtitle(paste0("mCH Within VS. Outside Region: ", names(df)[i])) +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20), legend.title=element_blank()) + theme(legend.position="bottom")
  print(g)
}
dev.off()


### Break this down by cell type and age bins
cIndexes = splitit(postnatalpd$Cell.Type)
cmethDiff = list(Neuron = mapply(function(x,y) rowMeans(x[,cIndexes$Neuron] - y[,cIndexes$Neuron]), meanMeth, outCRE),
                        Glia = mapply(function(x,y) rowMeans(x[,cIndexes$Glia] - y[,cIndexes$Glia]), meanMeth, outCRE))
cmethDiff = lapply(cmethDiff, function(x) lapply(x, function(y) data.frame(Difference = y)))
cmethDiff = do.call(rbind, Map(cbind, lapply(cmethDiff, function(x) 
  do.call(rbind, Map(cbind, x, Group = as.list(names(x))))), CellType = as.list(names(cmethDiff))))

aIndexes = splitit(postnatalpd$Age.Bin)
amethDiff = list(Neonate = mapply(function(x,y) rowMeans(x[,aIndexes$Neonate] - y[,aIndexes$Neonate]), meanMeth, outCRE),
                   Toddler = mapply(function(x,y) rowMeans(x[,aIndexes$Toddler] - y[,aIndexes$Toddler]), meanMeth, outCRE),
                   Child = mapply(function(x,y) rowMeans(x[,aIndexes$Child] - y[,aIndexes$Child]), meanMeth, outCRE),
                   Early.Teen = mapply(function(x,y) rowMeans(x[,aIndexes$Early.Teen] - y[,aIndexes$Early.Teen]), meanMeth, outCRE),
                   Teen = mapply(function(x,y) rowMeans(x[,aIndexes$Teen] - y[,aIndexes$Teen]), meanMeth, outCRE),     
                   Young.Adult = mapply(function(x,y) rowMeans(x[,aIndexes$Young.Adult] - y[,aIndexes$Young.Adult]), meanMeth, outCRE))
amethDiff = lapply(amethDiff, function(x) lapply(x, function(y) data.frame(Difference = y)))
amethDiff = do.call(rbind, Map(cbind, lapply(amethDiff, function(x) 
  do.call(rbind, Map(cbind, x, Group = as.list(names(x))))), Age.Bin = as.list(names(amethDiff))))

postnatalpd$Group = paste(postnatalpd$Cell.Type, postnatalpd$Age.Bin, sep = "_")
iIndexes = splitit(postnatalpd$Group)
imethDiff = list()
for (i in 1:length(unique(postnatalpd$Group))) {
imethDiff[[i]] = mapply(function(x,y) rowMeans(data.frame(x[,iIndexes[[unique(postnatalpd$Group)[i]]]]) - data.frame(y[,iIndexes[[unique(postnatalpd$Group)[i]]]])), 
                        meanMeth, outCRE, SIMPLIFY = F)
}
names(imethDiff) = unique(postnatalpd$Group)
imethDiff = lapply(imethDiff, function(x) lapply(x, function(y) data.frame(Difference = y)))
imethDiff = do.call(rbind, Map(cbind, lapply(imethDiff, function(x) 
  do.call(rbind, Map(cbind, x, Group = as.list(names(x))))), Int.Group = as.list(names(imethDiff))))


## Plot density of distance by cell type, age, direction of methylation change and whether the DMR is significant or not

cmethDiff$CellType = factor(cmethDiff$CellType, levels = c("Glia","Neuron"))
amethDiff$Age.Bin = factor(amethDiff$Age.Bin, levels = c("Neonate","Toddler","Child","Early.Teen","Teen","Young.Adult"))
imethDiff$Int.Group = factor(imethDiff$Int.Group, levels = c("Glia_Neonate","Glia_Toddler","Glia_Child","Glia_Early.Teen","Glia_Teen","Glia_Young.Adult",
                                                             "Neuron_Neonate","Neuron_Toddler","Neuron_Child","Neuron_Early.Teen","Neuron_Teen","Neuron_Young.Adult"))

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/non-CpG/figures/nonCpG_meth_in.VS.out_CREs.pdf")
for (i in 1:length(unique(cmethDiff$Group))) {
  g = ggplot(cmethDiff[which(cmethDiff$Group==unique(cmethDiff$Group)[i]),], aes(x=Difference)) + geom_density(aes(group=CellType, colour=CellType)) + 
    scale_colour_brewer(8, palette="Dark2") + 
    ylab("Density") + xlim(-0.4,0.4) +
    xlab("(Inside Region) - (Outside Region)") + 
    ggtitle(paste0("Mean mCH Difference by Cell Type\nBetween CREs and Flanking Regions:\n", unique(imethDiff$Group)[i])) + 
    theme(title = element_text(size = 20)) +
    theme(text = element_text(size = 20), legend.title=element_blank()) + theme(legend.position="bottom")
  print(g)
palette(c(brewer.pal(8, "Paired")[c(5:6, 7:8, 3:4, 1:2)], 'grey50'))
g = ggplot(amethDiff[which(amethDiff$Group==unique(amethDiff$Group)[i]),], aes(x=Difference)) + geom_density(aes(group=Age.Bin, colour=Age.Bin)) + 
  ylab("Density") + xlim(-0.4,0.4) +
  xlab("(Inside Region) - (Outside Region)") + 
  ggtitle(paste0("Mean mCH Difference by Age\nBetween CREs and Flanking Regions:\n", unique(imethDiff$Group)[i])) + 
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20), legend.title=element_blank()) + theme(legend.position="bottom")
print(g)
g = ggplot(imethDiff[which(imethDiff$Group==unique(imethDiff$Group)[i]),], aes(x=Difference)) + geom_density(aes(group=Int.Group, colour=Int.Group)) + 
  scale_colour_brewer(8, palette="Dark2") + 
  ylab("Density") + xlim(-0.4,0.4) +
  xlab("(Inside Region) - (Outside Region)") + 
  ggtitle(paste0("Mean mCH Difference by Age\nand Cell Type Between CREs\nand Flanking Regions: ", unique(imethDiff$Group)[i])) + 
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20), legend.title=element_blank()) + theme(legend.position="bottom")
print(g)
}
dev.off()


ctres = list(Neuron = mapply(function(inn, out) t.test(rowMeans(inn[,cIndexes$Neuron]), rowMeans(out[,cIndexes$Neuron])), meanMeth, outCRE, SIMPLIFY = F),
             Glia = mapply(function(inn, out) t.test(rowMeans(inn[,cIndexes$Glia]), rowMeans(out[,cIndexes$Glia])), meanMeth, outCRE, SIMPLIFY = F))
atres = list(Neonate = mapply(function(inn, out) t.test(rowMeans(inn[,aIndexes$Neonate]), rowMeans(out[,aIndexes$Neonate])), meanMeth, outCRE, SIMPLIFY = F),
             Toddler = mapply(function(inn, out) t.test(rowMeans(inn[,aIndexes$Toddler]), rowMeans(out[,aIndexes$Toddler])), meanMeth, outCRE, SIMPLIFY = F),
             Child = mapply(function(inn, out) t.test(rowMeans(inn[,aIndexes$Child]), rowMeans(out[,aIndexes$Child])), meanMeth, outCRE, SIMPLIFY = F),
             Early.Teen = mapply(function(inn, out) t.test(rowMeans(inn[,aIndexes$Early.Teen]), rowMeans(out[,aIndexes$Early.Teen])), meanMeth, outCRE, SIMPLIFY = F),
             Teen = mapply(function(inn, out) t.test(rowMeans(inn[,aIndexes$Teen]), rowMeans(out[,aIndexes$Teen])), meanMeth, outCRE, SIMPLIFY = F),     
             Young.Adult = mapply(function(inn, out) t.test(rowMeans(inn[,aIndexes$Young.Adult]), rowMeans(out[,aIndexes$Young.Adult])), meanMeth, outCRE, SIMPLIFY = F))
itres = list()
for (i in 1:length(unique(postnatalpd$Group))) {
  itres[[i]] = mapply(function(inn, out) t.test(rowMeans(data.frame(inn[,iIndexes[[unique(postnatalpd$Group)[i]]]])), 
                                                rowMeans(data.frame(out[,iIndexes[[unique(postnatalpd$Group)[i]]]]))), 
                          meanMeth, outCRE, SIMPLIFY = F)
}
names(itres) = unique(postnatalpd$Group)

ctres = do.call(rbind, Map(cbind, lapply(ctres, function(x) do.call(rbind, Map(cbind, lapply(x, function(y) 
                data.frame(Tstat = y$statistic, pval = y$p.value, In.mean = y$estimate[1], Out.mean = y$estimate[2])), 
                Group = as.list(names(x))))), Comparison = as.list(names(ctres))))
atres = do.call(rbind, Map(cbind, lapply(atres, function(x) do.call(rbind, Map(cbind, lapply(x, function(y) 
                data.frame(Tstat = y$statistic, pval = y$p.value, In.mean = y$estimate[1], Out.mean = y$estimate[2])), 
                Group = as.list(names(x))))), Comparison = as.list(names(atres))))
itres = do.call(rbind, Map(cbind, lapply(itres, function(x) do.call(rbind, Map(cbind, lapply(x, function(y) 
                data.frame(Tstat = y$statistic, pval = y$p.value, In.mean = y$estimate[1], Out.mean = y$estimate[2])), 
                Group = as.list(names(x))))), Comparison = as.list(names(itres))))
write.csv(rbind(cbind(ctres, padj = p.adjust(ctres$pval, method = "fdr")), cbind(atres, padj = p.adjust(atres$pval, method = "fdr")), 
                cbind(itres, padj = p.adjust(itres$pval, method = "fdr"))), quote = F, file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/nonCpG_meth_in.VS.out_CREs_tstats_byCellType_byAge.csv")

itres = cbind(itres, padj = p.adjust(itres$pval, method = "fdr"))
itres$CellType[grep("Neuron",itres$Comparison)] = "Neuron"
itres$CellType[grep("Glia",itres$Comparison)] = "Glia"
itres$Age = gsub("Neuron_", "", itres$Comparison)
itres$Age = gsub("Glia_", "", itres$Age)
itres$CRE = sub('\\..*', '', itres$Group)