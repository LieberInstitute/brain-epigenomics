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

x = split(meanCov[meanCov$CellType=="Neuron",], meanCov[meanCov$CellType=="Neuron","Cluster"])
x = lapply(x, function(y) split(y, y$Context))
lapply(x, function(d) match(d$CpH$Age, d$CpG$Age))

meth.cor = lapply(x, function(d) cor.test(d$CpG$value, d$CpH$value))
meth.cor = do.call(rbind, Map(cbind, Cluster = as.list(names(meth.cor)), 
                              lapply(meth.cor, function(x) data.frame(cor = x$estimate, Tstat = x$statistic, pval = x$p.value))))
meth.cor$N = c("+","+","-","0","-","0")

t.test(meth.cor[meth.cor$N=="+","cor"],meth.cor[meth.cor$N=="-","cor"], alternative = "greater")
#t = 5.784, df = 1.0013, p-value = 0.05439
#alternative hypothesis: true difference in means is greater than 0
#95 percent confidence interval:
#  -0.1035194        Inf
#sample estimates:
#  mean of x  mean of y 
#0.9711352 -0.1886908 

t.test(meth.cor[meth.cor$N=="+","cor"],meth.cor[meth.cor$N=="0","cor"], alternative = "greater")
#t = 1.014, df = 1.0004, p-value = 0.2478
#alternative hypothesis: true difference in means is greater than 0
#95 percent confidence interval:
#  -1.802356       Inf
#sample estimates:
#  mean of x mean of y 
#0.9711352 0.6259981 

t.test(meth.cor[meth.cor$N=="0","cor"],meth.cor[meth.cor$N=="-","cor"])
#t = 2.0625, df = 1.6193, p-value = 0.2044
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -1.331503  2.960881
#sample estimates:
#  mean of x  mean of y 
#0.6259981 -0.1886908 


meth.cor.6 = lapply(x, function(d) cor.test(d$CpG[d$CpG$Age<=6,"value"], d$CpH[d$CpH$Age<=6,"value"]))
meth.cor.6 = do.call(rbind, Map(cbind, Cluster = as.list(names(meth.cor.6)),
                                lapply(meth.cor.6, function(x) data.frame(cor = x$estimate, Tstat = x$statistic, pval = x$p.value))))
meth.cor.6$N = c("+","+","-","0","-","0")

t.test(meth.cor.6[meth.cor.6$N=="+","cor"],meth.cor.6[meth.cor.6$N=="-","cor"], alternative = "greater")
# t = 3.4082, df = 1.007, p-value = 0.09016
#alternative hypothesis: true difference in means is greater than 0
#95 percent confidence interval:
#  -0.6126862        Inf
#sample estimates:
#  mean of x mean of y 
#0.9471962 0.2097767 

t.test(meth.cor.6[meth.cor.6$N=="+","cor"],meth.cor.6[meth.cor.6$N=="0","cor"], alternative = "greater")
# t = 0.97075, df = 1.0007, p-value = 0.2547
#alternative hypothesis: true difference in means is greater than 0
#95 percent confidence interval:
#  -3.662687       Inf
#sample estimates:
#  mean of x mean of y 
#0.9471962 0.2808161 

t.test(meth.cor.6[meth.cor.6$N=="0","cor"],meth.cor.6[meth.cor.6$N=="-","cor"])
# t = 0.098731, df = 1.1961, p-value = 0.9353
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -6.189754  6.331833
#sample estimates:
#  mean of x mean of y 
#0.2808161 0.2097767

slope = lapply(x, function(d) lapply(d, function(y) (y[y$Age==5.34,"value"]-y[y$Age==0.20,"value"])/(5.34-0.2)))
slope = do.call(rbind, Map(cbind, Cluster = as.list(names(slope)), 
                           lapply(slope, function(x) do.call(rbind, Map(cbind, Context = as.list(names(x)),
                                                                        lapply(x, function(y) data.frame(slope = y)))))))
slope$N = c("+","+","+","+","-","-","0","0","-","-","0","0")

t.test(slope[slope$N=="+" & slope$Context=="CpG","slope"],slope[slope$N=="+" & slope$Context=="CpH","slope"])
#t = 2.8549, df = 1.0006, p-value = 0.2144
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -0.1342821  0.2122615
#sample estimates:
#  mean of x   mean of y 
#0.043152489 0.004162783 

t.test(slope[slope$N=="-" & slope$Context=="CpG","slope"],slope[slope$N=="-" & slope$Context=="CpH","slope"])
#t = -3.7049, df = 1.0008, p-value = 0.1677
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -0.2389328  0.1308965
#sample estimates:
#  mean of x    mean of y 
#-0.055589213 -0.001571071 

t.test(slope[slope$N=="0" & slope$Context=="CpG","slope"],slope[slope$N=="0" & slope$Context=="CpH","slope"])
#t = 0.14343, df = 1.2723, p-value = 0.9053
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -0.04821896  0.05002962
#sample estimates:
#  mean of x   mean of y 
#0.002625927 0.001720598 

                                                        
slope.perc = lapply(x, function(d) lapply(d, function(y) ((y[y$Age==5.34,"value"]/max(y$value))-(y[y$Age==0.20,"value"]/max(y$value)))/(5.34-0.2)))
slope.perc = do.call(rbind, Map(cbind, Cluster = as.list(names(slope.perc)), 
                           lapply(slope.perc, function(x) do.call(rbind, Map(cbind, Context = as.list(names(x)),
                                                                        lapply(x, function(y) data.frame(slope = y)))))))
slope.perc$N = c("+","+","+","+","-","-","0","0","-","-","0","0")

t.test(slope.perc[slope.perc$N=="+" & slope.perc$Context=="CpG","slope"],slope.perc[slope.perc$N=="+" & slope.perc$Context=="CpH","slope"])
# t = -0.09166, df = 1.0247, p-value = 0.9415
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -0.2295401  0.2260613
#sample estimates:
#  mean of x  mean of y 
#0.05666395 0.05840340 

t.test(slope.perc[slope.perc$N=="-" & slope.perc$Context=="CpG","slope"],slope.perc[slope.perc$N=="-" & slope.perc$Context=="CpH","slope"])
# t = -1.3793, df = 1.9632, p-value = 0.3039
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -0.12060332  0.06284357
#sample estimates:
#  mean of x   mean of y 
#-0.07478515 -0.04590527 

t.test(slope.perc[slope.perc$N=="0" & slope.perc$Context=="CpG","slope"],slope.perc[slope.perc$N=="0" & slope.perc$Context=="CpH","slope"])
# t = -0.60062, df = 1.1624, p-value = 0.6439
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -0.3191516  0.2799994
#sample estimates:
#  mean of x   mean of y 
#0.001284608 0.020860714 



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

ggplot(meanCov[meanCov$Context=="CpG",], aes(x = Age, y = value, colour = CellType)) + facet_grid(. ~ Cluster) +
  geom_point() + geom_smooth() +
  theme_classic() + scale_colour_brewer(8, palette="Dark2") +
  ylab("mCpG / CpG") + xlab("Age") + ylim(0,1) +
  ggtitle("Mean mCpG Within cdDMRs by Cluster") + 
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20), legend.title=element_blank()) + theme(legend.position="bottom")
ggplot(meanCov[meanCov$Context=="CpH",], aes(x = Age, y = value, colour = CellType)) + facet_grid(. ~ Cluster) +
  geom_smooth() + geom_point() + 
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












