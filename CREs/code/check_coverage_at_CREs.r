library('bsseq')
library('bumphunter')
library('doParallel')
library("bsseq")
library("GenomicFeatures")
library("MethylSeekR")
library("BSgenome")
library("BSgenome.Hsapiens.UCSC.hg19")


# load data

rm(list= ls()[!(ls() %in% c('CGlist','CGPrenlist'))])
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/UMRs_LMRs_annotated.rda")
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/dmvs_annotated.rda")
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/PMDs_annotated.rda")
load('/dcl01/lieber/ajaffe/lab/brain-epigenomics/bumphunting/BSobj_bsseqSmooth_Neuron_minCov_3.Rdata')

pd = pData(BSobj)
pd$Race[pd$Race== "CAUC "] <- 'CAUC'
pd$Sex[pd$Sex == " M"] <- 'M'
pd$RIN <- as.numeric(gsub(" ", "", pd$RIN))
pd$pg.DNA.nuclei.input <- as.numeric(pd$pg.DNA.nuclei.input)
pd$Reads <- as.numeric(pd$Reads)
pd$Percent.GreaterThan.Q30 <- as.numeric(pd$Percent.GreaterThan.Q30)


cov = getCoverage(BSobj, type = "Cov")
dim(cov)
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/bumphunting/BSobj_bsseqSmooth_Neuron_minCov_3_prenatal.Rdata")
covPren = getCoverage(BSobj, type = "Cov")
dim(covPren)
cov = cbind(cov, covPren)
cov = as.matrix(as.data.frame(cov))
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
ind = lapply(ov, function(x) split(queryHits(x), subjectHits(x)))

meanCov = lapply(ind, function(x) sapply(x, function(ii) colMeans(t(t(cov[ii,])))))
meanCov = lapply(meanCov, function(x) if (class(x)=="list") { do.call("rbind", x) } else { t(x) } )
meanCov = lapply(meanCov, function(x) as.data.frame(t(data.frame(x))))
meanCov = Map(cbind, lapply(meanCov, function(x) cbind(x, id = rownames(x))), Group = as.list(names(meanCov)))
for (i in 1:length(meanCov)) {
  meanCov[[i]][-which(meanCov[[i]][,"id"] %in% postnatalpd$Data.ID),"celltype"] = "Prenatal"
  meanCov[[i]][which(meanCov[[i]][,"id"] %in% postnatalpd[postnatalpd$Cell.Type=="Neuron","Data.ID"]),"celltype"] = "Neuron"
  meanCov[[i]][which(meanCov[[i]][,"id"] %in% postnatalpd[postnatalpd$Cell.Type=="Glia","Data.ID"]),"celltype"] = "Glia"
  meanCov[[i]][,"Age"] = postnatalpd[match(meanCov[[i]][,"id"], postnatalpd$Data.ID), "Age"]
}

corr = vector("list", length=length(meanCov))
for (i in 1:length(meanCov)) {
  for (j in 1:elementNROWS(ind)[i]) {
    corr[[i]][[j]] = cor.test(x = meanCov[[i]][which(meanCov[[i]][,"celltype"]!="Prenatal"),j], 
                              y = meanCov[[i]][which(meanCov[[i]][,"celltype"]!="Prenatal"),"Age"])
  }
}
save(corr, file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/coverage_check.r")

x = lapply(corr, function(x) do.call(rbind, lapply(x, function(y) data.frame(pval = y$p.value, estimate = y$estimate))))
x = do.call(rbind, Map(cbind, x, Regions = as.list(names(meanCov))))
x = data.table(x)
x[,list(max = max(estimate), min = min(estimate), mean = mean(estimate), median = median(estimate), sd = sd(estimate)),by="Regions"]
#        Regions       max        min        mean      median         sd
#1:      LMR.All 0.1381719 -0.4009710 -0.10216588 -0.10992318 0.14186114
#2:      UMR.All        NA         NA          NA          NA         NA
#3:      DMV.All 0.1177752 -0.3678718 -0.15409766 -0.15318706 0.07637643
#4:      PMD.All 0.4540102 -0.3966330  0.10347465  0.10911761 0.13486947
#5: LMR.Prenatal        NA         NA          NA          NA         NA
#6:   LMR.Neuron        NA         NA          NA          NA         NA
#7:     LMR.Glia        NA         NA          NA          NA         NA
#8: UMR.Prenatal        NA         NA          NA          NA         NA
#9:   UMR.Neuron        NA         NA          NA          NA         NA
#10:     UMR.Glia        NA         NA          NA          NA         NA
#11: DMV.Prenatal        NA         NA          NA          NA         NA
#12:   DMV.Neuron        NA         NA          NA          NA         NA
#13:     DMV.Glia        NA         NA          NA          NA         NA
#14: PMD.Prenatal 0.5455393 -0.2704574  0.22828721  0.25222361 0.10881298
#15:   PMD.Neuron 0.5455393 -0.2704574  0.09358677  0.09302474 0.11710940
#16:     PMD.Glia 0.5455393 -0.2704574  0.22695230  0.24657647 0.10687669

x[pval<=0.01, list(length()), by="Regions"]

tres = lapply(meanCov, function(x) t.test(x = rowMeans(x[which(x$celltype!="Prenatal"),c(1:(ncol(x)-4))]),
                                           y = rowMeans(x[which(x$celltype=="Prenatal"),c(1:(ncol(x)-4))])))
trespostnatal = lapply(meanCov, function(x) t.test(x = rowMeans(x[which(x$celltype=="Neuron"),c(1:(ncol(x)-4))]),
                                          y = rowMeans(x[which(x$celltype=="Glia"),c(1:(ncol(x)-4))])))
write.csv(data.frame(Regions = names(tres), Postnatal.Prenatal.Tstat = unlist(lapply(tres, function(x) x$statistic)),
                     Postnatal.Prenatal.pval = unlist(lapply(tres, function(x) x$p.value)), 
                     Postnatal.mean = unlist(lapply(tres, function(x) x$estimate[1])), 
                     Prenatal.mean = unlist(lapply(tres, function(x) x$estimate[2])),
                     Neuron.Glia.Tstat = unlist(lapply(trespostnatal, function(x) x$statistic)),
                     Neuron.Glia.pval = unlist(lapply(trespostnatal, function(x) x$p.value)), 
                     Neuron.mean = unlist(lapply(trespostnatal, function(x) x$estimate[1])), 
                     Glia.mean = unlist(lapply(trespostnatal, function(x) x$estimate[2]))), quote = F,
          file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/cellType_CRE_coverage_tstats.csv")

df = do.call(rbind, Map(cbind, lapply(meanCov, function(x) 
  rbind(data.frame(rowmeans = rowMeans(x[which(x$celltype=="Prenatal"),c(1:(ncol(x)-4))]), celltype = "Prenatal"),
  data.frame(rowmeans = rowMeans(x[which(x$celltype=="Neuron"),c(1:(ncol(x)-4))]), celltype = "Neuron"),
  data.frame(rowmeans = rowMeans(x[which(x$celltype=="Glia"),c(1:(ncol(x)-4))]), celltype = "Glia"))), 
  Regions = as.list(names(meanCov))))

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/CREs/figures/prenatal_postnatal_coverage_inCREs.pdf", width = 26)
ggplot(df, aes(x = celltype, y = rowmeans, colour = celltype)) + 
  theme_classic() + geom_boxplot() +
  labs(fill="") + facet_grid(. ~ Regions) +
  ylab("Mean Coverage") + 
  xlab("") + theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  ggtitle("Mean Coverage by Cell Type in CRE sequence") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20), legend.title=element_blank()) + theme(legend.position="bottom")
dev.off()

for (i in 1:length(meanCov)) {
  meanCov[[i]] = cbind(meanCov[[i]], Age.bin = postnatalpd[match(meanCov[[i]][,"id"], postnatalpd$Data.ID),"Age.Bin"])
  meanCov[[i]][,"Age.bin"] = as.character(meanCov[[i]][,"Age.bin"])
  meanCov[[i]][is.na(meanCov[[i]][,"Age.bin"]),"Age.bin"] = "Prenatal"
}

df = do.call(rbind, Map(cbind, lapply(meanCov, function(x) 
  rbind(data.frame(rowmeans = rowMeans(x[which(x$Age.bin=="Prenatal"),c(1:(ncol(x)-5))]), Age.bin = "Prenatal"),
        data.frame(rowmeans = rowMeans(x[which(x$Age.bin=="Neonate"),c(1:(ncol(x)-5))]), Age.bin = "Neonate"),
        data.frame(rowmeans = rowMeans(x[which(x$Age.bin=="Toddler"),c(1:(ncol(x)-5))]), Age.bin = "Toddler"),
        data.frame(rowmeans = rowMeans(x[which(x$Age.bin=="Child"),c(1:(ncol(x)-5))]), Age.bin = "Child"),
        data.frame(rowmeans = rowMeans(x[which(x$Age.bin=="Early.Teen"),c(1:(ncol(x)-5))]), Age.bin = "Early.Teen"),
        data.frame(rowmeans = rowMeans(x[which(x$Age.bin=="Teen"),c(1:(ncol(x)-5))]), Age.bin = "Teen"),
        data.frame(rowmeans = rowMeans(x[which(x$Age.bin=="Young.Adult"),c(1:(ncol(x)-5))]), Age.bin = "Young.Adult"))), 
  Regions = as.list(names(meanCov))))


pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/CREs/figures/acrossAge_coverage_inCREs.pdf", width = 52)
ggplot(df, aes(x = Age.bin, y = rowmeans, colour = Age.bin)) + 
  theme_classic() + geom_boxplot() +
  labs(fill="") + facet_grid(. ~ Regions) +
  ylab("Mean Coverage") + 
  xlab("") + theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  ggtitle("Mean Coverage by Age Bin in CRE sequence") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20), legend.title=element_blank()) + theme(legend.position="bottom")
dev.off()