library(GenomicRanges)
library(GenomicFeatures)
library(data.table)
library(ggplot2)
library(clusterProfiler)
require(org.Hs.eg.db)
library(RColorBrewer)
library(bumphunter)


load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/bumphunting/rda/limma_Neuron_CpGs_minCov_3_ageInfo_dmrs.Rdata")
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/DMR/DMR_objects.rda")
load('/dcl01/lieber/ajaffe/lab/brain-epigenomics/bumphunting/BSobj_bsseqSmooth_Neuron_minCov_3.Rdata')


## Separate Interaction DMRs into kmeans clusters

dtinteraction = data.table(DMR$Interaction)
dmrs = split(dmrs, dmrs$k6cluster_label)
oo = lapply(dmrs, function(x) findOverlaps(x, makeGRangesFromDataFrame(dtinteraction[sig=="FWER < 0.05",,])))
dtinteraction = dtinteraction[sig=="FWER < 0.05",,]
intclusters = lapply(oo, function(x) dtinteraction[subjectHits(x),,])

## Get CpGs in each DMR

gr = granges(BSobj)
meth = getMeth(BSobj, type = 'raw')

oo = lapply(intclusters, function(x) findOverlaps(reduce(makeGRangesFromDataFrame(x)),gr))

rIndexes = lapply(oo, function(x) split(subjectHits(x), queryHits(x)))
splmeth = lapply(rIndexes, function(x) sapply(x, function(ii) colMeans(t(t(meth[ii,])))))
meanCov = lapply(splmeth, function(x) do.call("rbind", x))
meanCov = lapply(meanCov, reshape2::melt)
meanCov = do.call(rbind, Map(cbind, meanCov, Cluster = as.list(names(meanCov))))
meanCov$Age = factor(pd[match(meanCov$Var2, pd$Data.ID),"Age"])
meanCov$Cluster = factor(meanCov$Cluster, levels = c("1:G-N+","2:G0N+","3:G0N-","4:G+N0","5:G+N-","6:G-N0"))
meanCov$CellType = factor(pd[match(meanCov$Var2, pd$Data.ID),"Cell.Type"])

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/DMR/CT_Age_Interaction/figures/interaction_kmeansClusters_mean_mCpG_byAge.pdf", width = 12)
ggplot(meanCov, aes(x = Age, y = value, fill= Cluster)) + geom_boxplot() +
  theme_classic() +  
  ylab("mCpG / CpG") + xlab("") + facet_grid(CellType ~ Cluster) +
  ggtitle("Mean mCpG Within cdDMRs by Cluster") + 
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20), legend.title=element_blank()) + theme(legend.position="bottom")
dev.off()


## Add CpH levels and plot as scatterplot 

meanCov = lapply(splmeth, function(x) colMeans(do.call("rbind", x)))
meanCov = do.call("rbind", meanCov)
meanCov = as.data.frame(meanCov)
meanCov$Cluster = rownames(meanCov)
meanCov = reshape2::melt(meanCov)
meanCov$Age = pd[match(meanCov$variable, pd$Data.ID),"Age"]
meanCov$Cluster = factor(meanCov$Cluster, levels = c("1:G-N+","2:G0N+","3:G0N-","4:G+N0","5:G+N-","6:G-N0"))
meanCov$CellType = factor(pd[match(meanCov$variable, pd$Data.ID),"Cell.Type"])


load('/dcl01/lieber/ajaffe/lab/brain-epigenomics/bsseq/bsobj_by_chr/allChrs_postNatal_cleaned_nonCG_noHomogenate_highCov.Rdata')

Hmeth <- getMeth(BSobj, type = 'raw')
Hgr = granges(BSobj)
Hmeth = getMeth(BSobj, type = 'raw')

oo = lapply(intclusters, function(x) findOverlaps(reduce(makeGRangesFromDataFrame(x)),Hgr))

rIndexes = lapply(oo, function(x) split(subjectHits(x), queryHits(x)))
splHmeth = lapply(rIndexes, function(x) sapply(x, function(ii) colMeans(t(t(Hmeth[ii,])))))
HmeanCov = lapply(splHmeth, function(x) colMeans(do.call("rbind", x)))
HmeanCov = do.call("rbind", HmeanCov)
HmeanCov = as.data.frame(HmeanCov)
HmeanCov$Cluster = rownames(HmeanCov)
HmeanCov = reshape2::melt(HmeanCov)
HmeanCov$Age = pd[match(HmeanCov$variable, pd$Data.ID),"Age"]
HmeanCov$Cluster = factor(HmeanCov$Cluster, levels = c("1:G-N+","2:G0N+","3:G0N-","4:G+N0","5:G+N-","6:G-N0"))
HmeanCov$CellType = factor(pd[match(HmeanCov$variable, pd$Data.ID),"Cell.Type"])

meanCov = rbind(cbind(meanCov, Context = "CpG"), cbind(HmeanCov, Context = "CpH"))
write.csv(meanCov, quote=F, file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/DMR/CT_Age_Interaction/interaction_kmeansClusters_mean_mC.csv")

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/DMR/CT_Age_Interaction/figures/interaction_kmeansClusters_mean_mCpG_byAge_scatterplot.pdf", 
    height = 4, width = 12)
ggplot(meanCov[meanCov$Context=="CpG",], aes(x = Age, y = value, colour = CellType)) + geom_point() + facet_grid(. ~ Cluster) +
  geom_smooth(method = "lm") +  
  theme_classic() + scale_colour_brewer(8, palette="Dark2") +
  ylab("mCpG / CpG") + xlab("Age") + ylim(0,1) +
  ggtitle("Mean mCpG Within cdDMRs by Cluster") + 
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20), legend.title=element_blank()) + theme(legend.position="bottom")
ggplot(meanCov[meanCov$Context=="CpG",], aes(x = Age, y = value, colour = CellType)) + facet_grid(. ~ Cluster) +
  geom_path() + geom_point() + 
  theme_classic() + scale_colour_brewer(8, palette="Dark2") +
  ylab("mCpG / CpG") + xlab("Age") + ylim(0,1) +
  ggtitle("Mean mCpG Within cdDMRs by Cluster") + 
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20), legend.title=element_blank()) + theme(legend.position="bottom")
ggplot(meanCov[meanCov$Context=="CpH",], aes(x = Age, y = value, colour = CellType)) + facet_grid(. ~ Cluster) +
  geom_path() + geom_point() + 
  theme_classic() + scale_colour_brewer(8, palette="Dark2") +
  ylab("mCpH / CpH") + xlab("Age") + ylim(0,.1) +
  ggtitle("Mean mCpH Within cdDMRs by Cluster") + 
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20), legend.title=element_blank()) + theme(legend.position="bottom")
dev.off()


## Correlate mean mCG and mCH methylation within cdDMRs

oo = lapply(intclusters, function(x) findOverlaps(reduce(makeGRangesFromDataFrame(x)),c(granges(gr),granges(Hgr))))
rIndexes = lapply(oo, function(x) split(subjectHits(x), queryHits(x)))
meth = as.data.frame(meth)
Hmeth = as.data.frame(Hmeth)
meth$context = "CpG"
Hmeth$context = "CpH"

spl = lapply(rIndexes, function(x) sapply(x, function(ii) rbind(meth,Hmeth)[ii,]))



corr = mapply(function(G,H) mapply(function(g,h) data.frame(CpG = g, CpH = h), G,H,SIMPLIFY = F), splmeth,splHmeth,SIMPLIFY = F)












