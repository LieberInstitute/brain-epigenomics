library(ggplot2)
library(GenomicRanges)

load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/DMR/DMR_objects.rda", verbose = T)
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/PMDs_100kb_noGaps.rda", verbose = T)
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/UMRs_LMRs_annotated.rda", verbose = T)
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/dmvs_annotated.rda", verbose = T)
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/bumphunting/rda/limma_Neuron_CpGs_minCov_3_ageInfo_dmrs.Rdata")


## Prepare DMRs

dmrs = split(dmrs, dmrs$k6cluster_label)
oo = lapply(dmrs, function(x) findOverlaps(x, makeGRangesFromDataFrame(DMR$Interaction)))
dmrs = c(lapply(DMR[which(names(DMR) %in% c("CellType","Age"))], function(x) x[which(x$sig=="FWER < 0.05"),]), lapply(oo, function(x) DMR$Interaction[subjectHits(x),]))
names(dmrs) = c("CellType","Age","Gr1","Gr2","Gr3","Gr4","Gr5","Gr6")
rdmrs = lapply(dmrs, function(x) reduce(makeGRangesFromDataFrame(x)))
elementNROWS(rdmrs)
#CellType      Age      Gr1      Gr2      Gr3      Gr4      Gr5      Gr6 
#   11179      129      567      423      362       96      170      560 


## Prepare methylation features

methfeatures = list(UMR = lapply(uDMR, makeGRangesFromDataFrame, keep=T), LMR = lapply(lDMR, makeGRangesFromDataFrame, keep=T), 
                    PMD = lapply(total.nogaps[-grep("GSM", names(total.nogaps))], function(x) x[which(x$type=="PMD")]), 
                    DMV = lapply(dmvDMR, makeGRangesFromDataFrame, keep=T))
methfeatures$UMR = mapply(function(u,d) u[!u %in% d],  methfeatures$UMR, methfeatures$DMV, SIMPLIFY = F) # because DMV is a special case of UMR, remove DMVs from UMR list

rmethfeatures = lapply(methfeatures, function(m) lapply(m, reduce))
do.call(cbind, lapply(rmethfeatures, elementNROWS))


## Count the hits for each type of DMR for each person, for each 

oo = lapply(rdmrs, function(d) lapply(rmethfeatures, function(m) lapply(m, function(x) findOverlaps(makeGRangesFromDataFrame(d),x))))

df = lapply(oo, function(d) lapply(d, function(m) do.call(rbind, Map(cbind, lapply(m, function(i) 
  data.frame(dmrHits = length(unique(queryHits(i))), featureHits = length(unique(subjectHits(i))))), id = as.list(names(m))))))  
df = lapply(df, function(d) Map(cbind, d, num.features = lapply(rmethfeatures, elementNROWS), 
                                num.bases = lapply(rmethfeatures, function(m) unlist(lapply(m, function(i) sum(as.numeric(width(i))))))))
df = lapply(df, function(d) do.call(rbind, Map(cbind, d, feature = as.list(names(d)))))  
df = do.call(rbind, Map(cbind, df, group = as.list(names(df)), num.dmrs = as.list(elementNROWS(rdmrs))))


df$perc.dmrHits = round(df$dmrHits / df$num.dmrs * 100, 1)
df$perc.dmrHits.transf = df$perc.dmrHits / df$num.bases * 1000000
df$perc.featureHits = round(df$featureHits / df$num.features * 100, 1)
df$group = gsub("CellType","Cell Type (11179)", df$group)
df$group = gsub("Age","Age (129)", df$group)
df$group = gsub("Gr1","1:G-N+ (567)", df$group)
df$group = gsub("Gr2","2:G0N+ (423)", df$group)
df$group = gsub("Gr3","3:G0N- (362)", df$group)
df$group = gsub("Gr4","4:G+N0 (96)", df$group)
df$group = gsub("Gr5","5:G+N- (170)", df$group)
df$group = gsub("Gr6","6:G-N0 (560)", df$group)
df$celltype = ifelse(df$id %in% pd$Data.ID, pd[match(df$id, pd$Data.ID),"Cell.Type"], "Prenatal")
df$age = pd[match(df$id, pd$Data.ID),"Age"]
df$feature = factor(df$feature, levels = c("UMR","LMR","DMV","PMD"))
df$group = factor(df$group, levels = c("Cell Type (11179)","Age (129)","1:G-N+ (567)","2:G0N+ (423)","3:G0N- (362)","4:G+N0 (96)","5:G+N- (170)","6:G-N0 (560)"))
df$celltype = factor(df$celltype, levels = c("Prenatal","Glia","Neuron"))

write.csv(df, quote=F,file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/methFeatures_Overlap_withDMRs.csv")

df = read.csv("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/methFeatures_Overlap_withDMRs.csv")

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/CREs/figures/methFeatures_Overlap_withDMRs.pdf", width = 12, height = 8)
ggplot(df[which(df$group %in% c("Cell Type (11179)","Age (129)")),], aes(x = celltype, y = perc.dmrHits)) + geom_boxplot() +
  labs(fill="") + facet_grid(group ~ feature) +
  ylab("Percent") + 
  xlab("") +
  ggtitle("Percent DMRs Overlapped by Methylation Features") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20), legend.title=element_blank()) + theme(legend.position="bottom")
ggplot(df[-which(df$group %in% c("Cell Type (11179)","Age (129)")),], aes(x = celltype, y = perc.dmrHits, fill = group)) + geom_boxplot() + 
  labs(fill="") + facet_grid(feature ~ group) +
  scale_fill_brewer(8, palette="Dark2") +
  ylab("Percent") + xlab("") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Percent cdDMRs Overlapped by Methylation Features") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20), legend.title=element_blank()) + theme(legend.position="bottom")
w = df[-which(df$group %in% c("Cell Type (11179)","Age (129)")),]
ggplot(w[which(w$celltype!="Prenatal"),], aes(x = age, y = perc.dmrHits, colour = celltype)) + geom_point() +
  geom_smooth(method="loess", se=T, fullrange=TRUE) + scale_colour_brewer(8, palette="Dark2") + 
  labs(fill="") + facet_grid(feature ~ group) +
  ylab("Percent") + xlab("Age (Years)") +
  ggtitle("Percent cdDMRs Overlapped by Methylation Features") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20), legend.title=element_blank()) + theme(legend.position="bottom")
ggplot(w[which(w$celltype!="Prenatal"),], aes(x = age, y = perc.dmrHits, colour = celltype)) + geom_point() +
  geom_smooth(method=lm, se=T, fullrange=TRUE) + scale_colour_brewer(8, palette="Dark2") + 
  labs(fill="") + facet_grid(feature ~ group) +
  ylab("Percent") + xlab("Age (Years)") +
  ggtitle("Percent cdDMRs Overlapped by Methylation Features") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20), legend.title=element_blank()) + theme(legend.position="bottom")
ggplot(df[-which(df$group %in% c("Cell Type (11179)","Age (129)")),], aes(x = celltype, y = perc.dmrHits.transf, fill = group)) + geom_boxplot() + 
  labs(fill="") + facet_grid(feature ~ group) +
  scale_fill_brewer(8, palette="Dark2") +
  ylab("Percent / Mb") + xlab("") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Percent cdDMRs Overlapped by Methylation Features") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20), legend.title=element_blank()) + theme(legend.position="bottom")
w = df[-which(df$group %in% c("Cell Type (11179)","Age (129)")),]
ggplot(w[which(w$celltype!="Prenatal"),], aes(x = age, y = perc.dmrHits.transf, colour = celltype)) + geom_point() +
  geom_smooth(method="loess", se=T, fullrange=TRUE) + scale_colour_brewer(8, palette="Dark2") + 
  labs(fill="") + facet_grid(feature ~ group) +
  ylab("Percent / Mb") + xlab("Age (Years)") +
  ggtitle("Percent cdDMRs Overlapped by Methylation Features") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20), legend.title=element_blank()) + theme(legend.position="bottom")
ggplot(w[which(w$celltype!="Prenatal"),], aes(x = age, y = perc.dmrHits.transf, colour = celltype)) + geom_point() +
  geom_smooth(method=lm, se=T, fullrange=TRUE) + scale_colour_brewer(8, palette="Dark2") + 
  labs(fill="") + facet_grid(feature ~ group) +
  ylab("Percent / Mb") + xlab("Age (Years)") +
  ggtitle("Percent cdDMRs Overlapped by Methylation Features") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20), legend.title=element_blank()) + theme(legend.position="bottom")
dev.off()


## Calculate correlation from last plot

w = w[which(w$celltype!="Prenatal"),]
f = unique(w$feature)
group = unique(w$group)
cor = list()
for (i in 1:length(group)) {
  cor[[i]] = list()
  for (j in 1:length(f)) {
    cor[[i]][[j]] = list(Neuron = cor.test(w[which(w[,"group"]==group[i] & w[,"feature"]==f[j] & w$celltype=="Neuron"),"perc.dmrHits"],
                                           w[which(w[,"group"]==group[i] & w[,"feature"]==f[j] & w$celltype=="Neuron"),"age"]),
                         Glia = cor.test(w[which(w[,"group"]==group[i] & w[,"feature"]==f[j] & w$celltype=="Glia"),"perc.dmrHits"],
                                         w[which(w[,"group"]==group[i] & w[,"feature"]==f[j] & w$celltype=="Glia"),"age"]))
  }
  names(cor[[i]]) = f
}
names(cor) = group
cor = lapply(cor, function(g) lapply(g, function(f) 
  do.call(rbind, Map(cbind, lapply(f, function(cell) data.frame(pval = cell$p.value, tstat = cell$statistic, cor = cell$estimate)), celltype = as.list(names(f))))))
cor = do.call(rbind, Map(cbind, lapply(cor, function(g) do.call(rbind, Map(cbind, g, feature = as.list(names(g))))), group = names(cor)))
cor$fdr = p.adjust(cor$pval, method = "fdr")
rownames(cor) = c()
write.csv(cor, quote=F,file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/methFeatures_intCluster_overlap_correlation_byAge.csv")

cor = cor[order(cor$celltype,cor$feature,cor$group),]
cor[which(cor$fdr<=0.05),]

means = data.table(df)[,mean(perc.dmrHits), by=c("celltype","group","feature")]
means = means[order(feature),,]
write.csv(means, quote=F,file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/methFeatures_meanOverlap_withDMRs.csv")