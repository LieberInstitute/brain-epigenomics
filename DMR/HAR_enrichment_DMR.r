library(GenomicRanges)
library(bumphunter)
library(RColorBrewer)

load('/dcl01/lieber/ajaffe/lab/brain-epigenomics/bumphunting/BSobj_bsseqSmooth_Neuron_minCov_3.Rdata')
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/DMR/DMR_objects.rda")
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/bumphunting/rda/limma_Neuron_CpGs_minCov_3_ageInfo_dmrs.Rdata")


# load HARs
HARs = openxlsx::read.xlsx('/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/HARs_hg19_Doan_Walsh_Table_S1.xlsx')
hars = makeGRangesFromDataFrame(HARs, keep.extra.columns=T)
length(hars) # 2737

# Identify all CpG clusters in the genome
gr <- granges(BSobj)
cl = clusterMaker( chr = as.character(seqnames(gr)), 
                   pos = start(gr),   maxGap = 1000)
gr.clusters = split(gr, cl)
gr.clusters = unlist(range(gr.clusters))
df.clusters = as.data.frame(gr.clusters)

# Find overlaps with DMRS in all three models

dmrs = split(dmrs, dmrs$k6cluster_label)
names(dmrs) = c("Gr1","Gr2","Gr3","Gr4","Gr5","Gr6")
oo = lapply(dmrs, function(x) findOverlaps(x, makeGRangesFromDataFrame(DMR$Interaction)))
dmrs = lapply(oo, function(x) DMR$Interaction[subjectHits(x),])

DMRgr = lapply(c(DMR, dmrs), function(x) makeGRangesFromDataFrame(x[which(x$fwer<=0.05),], keep.extra.columns=T))


oo = lapply(DMRgr, function(x) findOverlaps(gr.clusters, x))
lapply(oo, function(x) length(unique(queryHits(x))))

harOverlap = findOverlaps(hars, gr.clusters)

df.clusters$regionID = paste0(df.clusters$seqnames,":",df.clusters$start,"-",df.clusters$end)
df.clusters$rnum = 1:length(gr.clusters)
df.clusters$CellType = ifelse(df.clusters$rnum %in% queryHits(oo$CellType), "CellType","no")
df.clusters$Age = ifelse(df.clusters$rnum %in% queryHits(oo$Age), "Age","no")
df.clusters$Interaction = ifelse(df.clusters$rnum %in% queryHits(oo$Interaction), "Interaction","no")
df.clusters$Gr1 = ifelse(df.clusters$rnum %in% queryHits(oo$Gr1), "Gr1","no")
df.clusters$Gr2 = ifelse(df.clusters$rnum %in% queryHits(oo$Gr2), "Gr2","no")
df.clusters$Gr3 = ifelse(df.clusters$rnum %in% queryHits(oo$Gr3), "Gr3","no")
df.clusters$Gr4 = ifelse(df.clusters$rnum %in% queryHits(oo$Gr4), "Gr4","no")
df.clusters$Gr5 = ifelse(df.clusters$rnum %in% queryHits(oo$Gr5), "Gr5","no")
df.clusters$Gr6 = ifelse(df.clusters$rnum %in% queryHits(oo$Gr6), "Gr6","no")
df.clusters$HARs = ifelse(df.clusters$rnum %in% subjectHits(harOverlap), "HAR","no")


## make contingency tables

tables = list()
for (i in 1:length(names(DMRgr))) {
  tables[[i]] = data.frame(YesHAR = c(nrow(df.clusters[df.clusters$HARs=="HAR" & df.clusters[,colnames(df.clusters)==names(DMRgr)[i]]==names(DMRgr)[i],]),
                                      nrow(df.clusters[df.clusters$HARs=="HAR" & df.clusters[,colnames(df.clusters)==names(DMRgr)[i]]=="no",])),
                           NoHAR = c(nrow(df.clusters[df.clusters$HARs=="no" & df.clusters[,colnames(df.clusters)==names(DMRgr)[i]]==names(DMRgr)[i],]),
                                     nrow(df.clusters[df.clusters$HARs=="no" & df.clusters[,colnames(df.clusters)==names(DMRgr)[i]]=="no",])), 
                           row.names = c("YesDMR","NoDMR"))
}
names(tables) = names(DMRgr)

fisher = lapply(tables, fisher.test)
df = do.call(rbind, Map(cbind, lapply(fisher, function(x) data.frame(OR = x$estimate, pval = x$p.value)), Model = as.list(names(fisher))))
df$fdr = p.adjust(df$pval, method= "fdr")


## test for enrichment of conserved, evolutionarily dated enhancers

noonan = openxlsx::read.xlsx('/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/Emera_Noonan_supptable_1.xlsx')
x = GRanges(noonan[,1])
mcols(x) = noonan[,2:3]
noonan = x
length(noonan) # 30526
noonan = c(list(enhancers = noonan), as.list(split(noonan, noonan$Phylogenetic.age.assignment)))

noonanOv = lapply(noonan, function(x) findOverlaps(x, gr.clusters))

hits = lapply(noonanOv, function(x) ifelse(df.clusters$rnum %in% subjectHits(x), "yes","no"))

for (i in 1:length(hits)) { df.clusters[,17+i] = hits[[i]] }
colnames(df.clusters)[18:28] = names(hits)


## make contingency tables

tables = list()
for (i in 1:length(names(DMRgr))) {
  tables[[i]] = list(enhancers = data.frame(YesHit = c(nrow(df.clusters[df.clusters$enhancers=="yes" & df.clusters[,colnames(df.clusters)==names(DMRgr)[i]]==names(DMRgr)[i],]),
                                                       nrow(df.clusters[df.clusters$enhancers=="yes" & df.clusters[,colnames(df.clusters)==names(DMRgr)[i]]=="no",])),
                                            NoHit = c(nrow(df.clusters[df.clusters$enhancers=="no" & df.clusters[,colnames(df.clusters)==names(DMRgr)[i]]==names(DMRgr)[i],]),
                                                      nrow(df.clusters[df.clusters$enhancers=="no" & df.clusters[,colnames(df.clusters)==names(DMRgr)[i]]=="no",])), 
                                            row.names = c("YesDMR","NoDMR")),
                     Amniota = data.frame(YesHit = c(nrow(df.clusters[df.clusters$Amniota=="yes" & df.clusters[,colnames(df.clusters)==names(DMRgr)[i]]==names(DMRgr)[i],]),
                                                     nrow(df.clusters[df.clusters$Amniota=="yes" & df.clusters[,colnames(df.clusters)==names(DMRgr)[i]]=="no",])),
                                          NoHit = c(nrow(df.clusters[df.clusters$Amniota=="no" & df.clusters[,colnames(df.clusters)==names(DMRgr)[i]]==names(DMRgr)[i],]),
                                                    nrow(df.clusters[df.clusters$Amniota=="no" & df.clusters[,colnames(df.clusters)==names(DMRgr)[i]]=="no",])), 
                                          row.names = c("YesDMR","NoDMR")),
                     Eutheria = data.frame(YesHit = c(nrow(df.clusters[df.clusters$Eutheria=="yes" & df.clusters[,colnames(df.clusters)==names(DMRgr)[i]]==names(DMRgr)[i],]),
                                                      nrow(df.clusters[df.clusters$Eutheria=="yes" & df.clusters[,colnames(df.clusters)==names(DMRgr)[i]]=="no",])),
                                           NoHit = c(nrow(df.clusters[df.clusters$Eutheria=="no" & df.clusters[,colnames(df.clusters)==names(DMRgr)[i]]==names(DMRgr)[i],]),
                                                     nrow(df.clusters[df.clusters$Eutheria=="no" & df.clusters[,colnames(df.clusters)==names(DMRgr)[i]]=="no",])), 
                                           row.names = c("YesDMR","NoDMR")),
                     Gnathostomata = data.frame(YesHit = c(nrow(df.clusters[df.clusters$Gnathostomata=="yes" & df.clusters[,colnames(df.clusters)==names(DMRgr)[i]]==names(DMRgr)[i],]),
                                                           nrow(df.clusters[df.clusters$Gnathostomata=="yes" & df.clusters[,colnames(df.clusters)==names(DMRgr)[i]]=="no",])),
                                                NoHit = c(nrow(df.clusters[df.clusters$Gnathostomata=="no" & df.clusters[,colnames(df.clusters)==names(DMRgr)[i]]==names(DMRgr)[i],]),
                                                          nrow(df.clusters[df.clusters$Gnathostomata=="no" & df.clusters[,colnames(df.clusters)==names(DMRgr)[i]]=="no",])), 
                                                row.names = c("YesDMR","NoDMR")),
                     Human = data.frame(YesHit = c(nrow(df.clusters[df.clusters$Human=="yes" & df.clusters[,colnames(df.clusters)==names(DMRgr)[i]]==names(DMRgr)[i],]),
                                                   nrow(df.clusters[df.clusters$Human=="yes" & df.clusters[,colnames(df.clusters)==names(DMRgr)[i]]=="no",])),
                                        NoHit = c(nrow(df.clusters[df.clusters$Human=="no" & df.clusters[,colnames(df.clusters)==names(DMRgr)[i]]==names(DMRgr)[i],]),
                                                  nrow(df.clusters[df.clusters$Human=="no" & df.clusters[,colnames(df.clusters)==names(DMRgr)[i]]=="no",])), 
                                        row.names = c("YesDMR","NoDMR")),     
                     Mammalia = data.frame(YesHit = c(nrow(df.clusters[df.clusters$Mammalia=="yes" & df.clusters[,colnames(df.clusters)==names(DMRgr)[i]]==names(DMRgr)[i],]),
                                                      nrow(df.clusters[df.clusters$Mammalia=="yes" & df.clusters[,colnames(df.clusters)==names(DMRgr)[i]]=="no",])),
                                           NoHit = c(nrow(df.clusters[df.clusters$Mammalia=="no" & df.clusters[,colnames(df.clusters)==names(DMRgr)[i]]==names(DMRgr)[i],]),
                                                     nrow(df.clusters[df.clusters$Mammalia=="no" & df.clusters[,colnames(df.clusters)==names(DMRgr)[i]]=="no",])), 
                                           row.names = c("YesDMR","NoDMR")), 
                     "No age assignment" = data.frame(YesHit = c(nrow(df.clusters[df.clusters$"No age assignment"=="yes" & df.clusters[,colnames(df.clusters)==names(DMRgr)[i]]==names(DMRgr)[i],]),
                                                                 nrow(df.clusters[df.clusters$"No age assignment"=="yes" & df.clusters[,colnames(df.clusters)==names(DMRgr)[i]]=="no",])),
                                                      NoHit = c(nrow(df.clusters[df.clusters$"No age assignment"=="no" & df.clusters[,colnames(df.clusters)==names(DMRgr)[i]]==names(DMRgr)[i],]),
                                                                nrow(df.clusters[df.clusters$"No age assignment"=="no" & df.clusters[,colnames(df.clusters)==names(DMRgr)[i]]=="no",])), 
                                                      row.names = c("YesDMR","NoDMR")),
                     Primate = data.frame(YesHit = c(nrow(df.clusters[df.clusters$Primate=="yes" & df.clusters[,colnames(df.clusters)==names(DMRgr)[i]]==names(DMRgr)[i],]),
                                                     nrow(df.clusters[df.clusters$Primate=="yes" & df.clusters[,colnames(df.clusters)==names(DMRgr)[i]]=="no",])),
                                          NoHit = c(nrow(df.clusters[df.clusters$Primate=="no" & df.clusters[,colnames(df.clusters)==names(DMRgr)[i]]==names(DMRgr)[i],]),
                                                    nrow(df.clusters[df.clusters$Primate=="no" & df.clusters[,colnames(df.clusters)==names(DMRgr)[i]]=="no",])), 
                                          row.names = c("YesDMR","NoDMR")),         
                     Tetrapoda = data.frame(YesHit = c(nrow(df.clusters[df.clusters$Tetrapoda=="yes" & df.clusters[,colnames(df.clusters)==names(DMRgr)[i]]==names(DMRgr)[i],]),
                                                       nrow(df.clusters[df.clusters$Tetrapoda=="yes" & df.clusters[,colnames(df.clusters)==names(DMRgr)[i]]=="no",])),
                                            NoHit = c(nrow(df.clusters[df.clusters$Tetrapoda=="no" & df.clusters[,colnames(df.clusters)==names(DMRgr)[i]]==names(DMRgr)[i],]),
                                                      nrow(df.clusters[df.clusters$Tetrapoda=="no" & df.clusters[,colnames(df.clusters)==names(DMRgr)[i]]=="no",])), 
                                            row.names = c("YesDMR","NoDMR")),
                     Theria = data.frame(YesHit = c(nrow(df.clusters[df.clusters$Theria=="yes" & df.clusters[,colnames(df.clusters)==names(DMRgr)[i]]==names(DMRgr)[i],]),
                                                    nrow(df.clusters[df.clusters$Theria=="yes" & df.clusters[,colnames(df.clusters)==names(DMRgr)[i]]=="no",])),
                                         NoHit = c(nrow(df.clusters[df.clusters$Theria=="no" & df.clusters[,colnames(df.clusters)==names(DMRgr)[i]]==names(DMRgr)[i],]),
                                                   nrow(df.clusters[df.clusters$Theria=="no" & df.clusters[,colnames(df.clusters)==names(DMRgr)[i]]=="no",])), 
                                         row.names = c("YesDMR","NoDMR")),        
                     Vertebrata = data.frame(YesHit = c(nrow(df.clusters[df.clusters$Vertebrata=="yes" & df.clusters[,colnames(df.clusters)==names(DMRgr)[i]]==names(DMRgr)[i],]),
                                                        nrow(df.clusters[df.clusters$Vertebrata=="yes" & df.clusters[,colnames(df.clusters)==names(DMRgr)[i]]=="no",])),
                                             NoHit = c(nrow(df.clusters[df.clusters$Vertebrata=="no" & df.clusters[,colnames(df.clusters)==names(DMRgr)[i]]==names(DMRgr)[i],]),
                                                       nrow(df.clusters[df.clusters$Vertebrata=="no" & df.clusters[,colnames(df.clusters)==names(DMRgr)[i]]=="no",])), 
                                             row.names = c("YesDMR","NoDMR")))
}
names(tables) = names(DMRgr)

fisher = lapply(tables, function(x) lapply(x, fisher.test))
res = do.call(rbind, Map(cbind, lapply(fisher, function(y) 
           do.call(rbind, Map(cbind, lapply(y, function(x) data.frame(OR = x$estimate, pval = x$p.value)), Feature = as.list(names(y))))),
           Model = as.list(names(fisher))))
res$fdr = p.adjust(res$pval, method= "fdr")
res$Feature = gsub("enhancers","All Enhancers", res$Feature)
res$Feature = gsub("assignment","assigned", res$Feature)
res$Feature = factor(res$Feature, levels=c("All Enhancers","No age assigned","Vertebrata","Gnathostomata","Tetrapoda",     
                                           "Amniota","Mammalia","Theria","Eutheria","Primate","Human"))     
                      
[5]            "Mammalia"                
[9] "Tetrapoda"                                                                               
                                                    
pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/DMR/figures/HARs_DMR_oddsRatios.pdf", height = 4)
ggplot(df[which(df$Model %in% c("CellType", "Age", "Interaction")),], aes(Model, OR, fill = Model)) + 
  geom_col() + theme_classic() + 
  ylab("Odds Ratio") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("") + geom_hline(yintercept=1, linetype="dotted") +
  ggtitle("Enrichment for HARs") + 
  theme(title = element_text(size = 20)) + theme(text = element_text(size = 20), legend.title=element_blank()) + theme(legend.position="bottom")
ggplot(df[which(df$Model %in% c("Gr1", "Gr2", "Gr3","Gr4","Gr5","Gr6")),], aes(Model, OR, fill = Model)) + 
  geom_col() + scale_fill_brewer(8, palette="Dark2") +
  theme_classic() + geom_hline(yintercept=1, linetype="dotted") +
  ylab("Odds Ratio") + 
  xlab("") + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Enrichment for HARs") + 
  theme(title = element_text(size = 20)) + theme(text = element_text(size = 20), legend.title=element_blank()) + theme(legend.position="bottom")
dev.off()

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/DMR/figures/conservedEnhancers_DMR_oddsRatios.pdf", width = 16, height = 4)
ggplot(res[which(res$Model %in% c("CellType", "Age", "Interaction") & res$fdr<=0.05),], aes(Model,OR, fill = Model)) + 
  geom_col() + theme_classic() + 
  facet_grid(. ~ Feature) +
  ylab("Odds Ratio") + 
  xlab("") + geom_hline(yintercept=1, linetype="dotted") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Enrichment for Enhancers") + 
  theme(title = element_text(size = 20)) + theme(text = element_text(size = 20), legend.title=element_blank()) + theme(legend.position="bottom")
ggplot(res[which(res$Model %in% c("Gr1", "Gr2", "Gr3","Gr4","Gr5","Gr6") & res$fdr<=0.05),], aes(Model,OR, fill = Model)) + 
  geom_col() + scale_fill_brewer(8, palette="Dark2") +
  theme_classic() + geom_hline(yintercept=1, linetype="dotted") +
  facet_grid(. ~ Feature) +
  ylab("Odds Ratio") + 
  xlab("") + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Enrichment for Enhancers") + 
  theme(title = element_text(size = 20)) + theme(text = element_text(size = 20), legend.title=element_blank()) + theme(legend.position="bottom")
dev.off()


res = rbind(res, cbind(df[,1:2], Feature = "HAR", df[,3:4]))
rownames(res) = NULL
write.csv(res, quote=F, file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/DMR/HAR_conservedEnhancers_fisher.results.csv")



