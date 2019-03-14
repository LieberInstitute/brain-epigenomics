library(GenomicRanges)
library(bumphunter)
library(ggplot2)
library(openxlsx)

load('/dcl01/lieber/ajaffe/lab/brain-epigenomics/bumphunting/BSobj_bsseqSmooth_Neuron_minCov_3.Rdata')
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/DMR/DMR_objects.rda")
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/bumphunting/rda/limma_Neuron_CpGs_minCov_3_ageInfo_dmrs.Rdata")

## Read in the Dracheva DMRs (Kozlenkov, Dracheva, Science Advances 2018)

paperDMRs = list()
for (i in 1:9) {
  paperDMRs[[i]] = read.xlsx("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/Kozlenkov_ScienceAdvances2018_TableS4.xlsx", sheet = i)
}
names(paperDMRs) = c("GLU hypo GABA tmCG DMRs", "GLU hypo GABA mCG DMRs","GLU hypo GABA hmCG DMRs",
                     "GABA hypo GLU tmCG DMRs","GABA hypo GLU mCG DMRs","GABA hypo GLU hmCG DMRs",
                     "OLIG hypo neuron tmCG DMRs","OLIG hypo neuron mCG DMRs","OLIG hypo neuron hmCG DMRs")
paperDMRs = lapply(paperDMRs, makeGRangesFromDataFrame)
elementNROWS(paperDMRs)
#GLU hypo GABA tmCG DMRs     GLU hypo GABA mCG DMRs    GLU hypo GABA hmCG DMRs    GABA hypo GLU tmCG DMRs     GABA hypo GLU mCG DMRs 
#                  54178                       9998                      23495                      14310                      40644 
#GABA hypo GLU hmCG DMRs OLIG hypo neuron tmCG DMRs  OLIG hypo neuron mCG DMRs OLIG hypo neuron hmCG DMRs 
#                   1934                      64435                      66704                        969

# GLU hypo GABA means lower methylation in GLU
# OLIG hypo neuron means lower methylation in OLIG
# GABA hypo GLU means lower methylation in GABA


## Identify all CpG clusters in the genome

gr = granges(BSobj)
cl = clusterMaker(chr = as.character(seqnames(gr)), pos = start(gr),   maxGap = 1000)
gr.clusters = split(gr, cl)
gr.clusters = unlist(range(gr.clusters))
df.clusters = as.data.frame(gr.clusters)


## Find overlaps with DMRS in all three models

dmrs = split(dmrs, dmrs$k6cluster_label)
names(dmrs) = c("Group 1 (G-N+)","Group 2 (G0N+)","Group 3 (G0N-)","Group 4 (G+N0)","Group 5 (G+N-)","Group 6 (G-N0)")
dmrs = c(dmrs, list(cdDMRs = unique(makeGRangesFromDataFrame(DMR$Interaction[which(DMR$Interaction$fwer<=0.05),]))))

ooDMR = lapply(dmrs, function(x) findOverlaps(gr.clusters, x))
ooGG = lapply(paperDMRs, function(x) findOverlaps(gr.clusters,x))

df.clusters$regionID = paste0(df.clusters$seqnames,":",df.clusters$start,"-",df.clusters$end)
df.clusters$rnum = 1:length(gr.clusters)
for (i in 1:length(dmrs)) { df.clusters[,names(dmrs)[i]] = ifelse(df.clusters$rnum %in% queryHits(ooDMR[[i]]), "yes","no") }
for (i in 1:length(ooGG)) { df.clusters[,names(ooGG)[i]] = ifelse(df.clusters$rnum %in% queryHits(ooGG[[i]]), "yes","no") }


## make contingency tables

gg = names(ooGG)
ours = names(ooDMR)

tables = list()
for (i in 1:length(gg)){
  tables[[i]] = list()
  for (j in 1:length(ours)) {
    tables[[i]][[j]] = data.frame(c(nrow(df.clusters[which(df.clusters[,which(colnames(df.clusters)==gg[i])]=="yes" & 
                                                             df.clusters[,which(colnames(df.clusters)==ours[j])]=="yes"),]),
                                    nrow(df.clusters[which(df.clusters[,which(colnames(df.clusters)==gg[i])]=="yes" & 
                                                             df.clusters[,which(colnames(df.clusters)==ours[j])]=="no"),])),
                                  c(nrow(df.clusters[which(df.clusters[,which(colnames(df.clusters)==gg[i])]=="no" & 
                                                             df.clusters[,which(colnames(df.clusters)==ours[j])]=="yes"),]),
                                    nrow(df.clusters[which(df.clusters[,which(colnames(df.clusters)==gg[i])]=="no" & 
                                                             df.clusters[,which(colnames(df.clusters)==ours[j])]=="no"),])),
                                  row.names = c(paste0("Yes", ours[j]), paste0("No", ours[j])))
    colnames(tables[[i]][[j]]) = c(paste0("Yes",gg[i]),paste0("No",gg[i]))
  }
  names(tables[[i]]) = ours
}
names(tables) = gg

res = lapply(tables, function(x) lapply(x, fisher.test))
df = lapply(res, function(x) lapply(x, function(y) data.frame(pval = y$p.value, OR = y$estimate)))
df = do.call(rbind, Map(cbind, PaperGroups = as.list(names(res)), lapply(df, function(x) do.call(rbind, Map(cbind, DMRgroups = as.list(names(x)), x)))))
df$FDR = p.adjust(df$pval, method = "fdr")

write.csv(df,file="/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/DMR/DMR_Kozlenkov_Dracheva_ScienceAdvances_2018_hydroxymeth_Overlap.csv", quote=F)



## Count the number of DMRs that overlap DMRs from Stella's paper

oo = lapply(paperDMRs, function(x) lapply(dmrs, function(d) findOverlaps(d, x)))

df = list()
for (i in 1:length(oo)) {
  df[[i]] = list()
  for (j in 1:length(oo[[i]])) {
    if ( length(queryHits(oo[[i]][[j]]))==0 ) {
      df[[i]][[j]] = data.frame(DMRgroup = names(oo[[i]])[j], queryHits = 0, subjectHits = 0, queryLength = queryLength(oo[[i]][[j]]), subjectLength = subjectLength(oo[[i]][[j]])) } 
    else {
      df[[i]][[j]] = data.frame(DMRgroup = names(oo[[i]])[j], queryHits = length(unique(queryHits(oo[[i]][[j]]))), subjectHits = length(unique(subjectHits(oo[[i]][[j]]))),
                                queryLength = queryLength(oo[[i]][[j]]), subjectLength = subjectLength(oo[[i]][[j]]))
      }
  }
}

df = do.call(rbind, Map(cbind, PaperGroup = as.list(names(oo)), lapply(df, function(x) do.call(rbind, x))))
df$PaperType = df$PaperComparison = NA
df[grep("tmCG", df$PaperGroup), "PaperType"] = "tmCG"
df[grep(" mCG", df$PaperGroup), "PaperType"] = "mCG"
df[grep("hmCG", df$PaperGroup), "PaperType"] = "hmCG"
df$percentDMR = df$queryHits/df$queryLength * 100
df$percentPaper = df$subjectHits/df$subjectLength * 100
df[grep("GLU hypo GABA", df$PaperGroup), "PaperComparison"] = "GLU < GABA"
df[grep("GABA hypo GLU", df$PaperGroup), "PaperComparison"] = "GABA < GLU"
df[grep("OLIG hypo neuron", df$PaperGroup), "PaperComparison"] = "OLIG < neuron"


write.csv(df,file="/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/DMR/DMR_Kozlenkov_Dracheva_ScienceAdvances_2018_hydroxymeth_Overlap_Counts.csv", quote=F)


## Are the number of DMRs between the three types in the Science Advances paper similar?

elementNROWS(paperDMRs)
res = list(tmCG.v.hmCG = t.test(unique(df[which(df$PaperType=="tmCG"),"subjectLength"]), unique(df[which(df$PaperType=="hmCG"),"subjectLength"])),
           tmCG.v.mCG = t.test(unique(df[which(df$PaperType=="tmCG"),"subjectLength"]), unique(df[which(df$PaperType=="mCG"),"subjectLength"])),
           mCG.v.hmCG = t.test(unique(df[which(df$PaperType=="mCG"),"subjectLength"]), unique(df[which(df$PaperType=="hmCG"),"subjectLength"])))

res = do.call(rbind, Map(cbind, Group = as.list(names(res)), 
                         lapply(res, function(x) data.frame(Tstat = x$statistic, mean1 = x$estimate[1], mean2 = x$estimate[2], pval = x$p.value))))
res$FDR = p.adjust(res$pval, method="fdr")
res
#         Group     Tstat    mean1     mean2      pval       FDR
#t  tmCG.v.hmCG 2.0930739 44307.67  8799.333 0.1312315 0.2961011
#t1  tmCG.v.mCG 0.2316804 44307.67 39115.333 0.8282077 0.8282077
#t2  mCG.v.hmCG 1.6878292 39115.33  8799.333 0.1974007 0.2961011

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/DMR/figures/DMR_Kozlenkov_Dracheva_ScienceAdvances_2018_numberDMRs.pdf",width=5.5,height=3.5)
ggplot(unique(df[,colnames(df) %in% c("PaperGroup","PaperType","subjectLength")]), aes(x = PaperType, y = subjectLength)) + 
  geom_boxplot(outlier.shape = NA) +
  labs(fill="") + ylab("Count") + xlab("") +
  ggtitle("Number of DMRs Identified\n(Kozlenkov, 2018)") +
  theme(title = element_text(size = 20), text = element_text(size = 20))
dev.off()


## Plot

df$PaperType = factor(df$PaperType, levels = c("hmCG", "mCG", "tmCG"))  
df$DMRgroup = gsub(" (", "\n(", df$DMRgroup, fixed=T)
df$DMRgroup = gsub("Group ", "Gr", df$DMRgroup, fixed=T)
df$DMRgroup = factor(df$DMRgroup, levels = c("Gr1\n(G-N+)","Gr2\n(G0N+)","Gr3\n(G0N-)",
                                             "Gr4\n(G+N0)","Gr5\n(G+N-)","Gr6\n(G-N0)","cdDMRs"))

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/DMR/figures/DMR_Kozlenkov_Dracheva_ScienceAdvances_2018_Overlap.pdf",width=18,height=6)
ggplot(df, aes(x = DMRgroup, y = percentDMR, fill=DMRgroup)) + geom_bar(stat = "identity") +
  geom_text(aes( label = queryHits), vjust = -.5) +
  scale_fill_brewer(8, palette="Dark2") + facet_grid(PaperType ~ PaperComparison) + theme_classic() +
  labs(fill="") + ylab("Percent") + ylim(0,100) + xlab("") +
  ggtitle("Overlap with DMRs defined in Kozlenkov (2018)") +
  theme(title = element_text(size = 20), text = element_text(size = 20), legend.position = "none")
dev.off()

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/DMR/figures/DMR_Kozlenkov_Dracheva_ScienceAdvances_2018_Overlap2.pdf",width=8,height=8)
ggplot(df, aes(x = DMRgroup, y = percentDMR, fill=PaperType)) + geom_bar(stat = "identity", position="dodge") +
  scale_fill_brewer(8, palette="Set1") + facet_grid(PaperComparison ~ .) +
  labs(fill="") + ylab("Percent") + xlab("") + ylim(0,82) +
  ggtitle("Overlap with DMRs defined in Kozlenkov (2018)") +
  geom_text(aes(label=queryHits), position=position_dodge(width=0.9), vjust=-0.25) +
  theme(title = element_text(size = 20), text = element_text(size = 20), legend.position = "bottom",
        legend.title = element_blank())
dev.off()


## find enrichment for interaction clusters

oo = lapply(dmrs, function(x) lapply(paperDMRs, function(d) findOverlaps(d, x)))

DMRs = list()
for (i in 1:length(dmrs)) {
  DMRs[[i]] = as.data.frame(dmrs[[i]])
  DMRs[[i]][,"rnum"] = 1:nrow(DMRs[[i]])
  for (j in 1:length(oo[[i]])) {
    DMRs[[i]][,names(oo[[i]])[j]] = ifelse(DMRs[[i]][,"rnum"] %in% subjectHits(oo[[i]][[j]]), "yes","no")
  }
}

set = list("GLU hypo GABA" = c("GLU hypo GABA tmCG DMRs", "GLU hypo GABA hmCG DMRs"),
           "GABA hypo GLU" = c("GABA hypo GLU tmCG DMRs", "GABA hypo GLU hmCG DMRs"),
           "OLIG hypo neuron" = c("OLIG hypo neuron tmCG DMRs", "OLIG hypo neuron hmCG DMRs"),
           "GLU hypo GABA" = c("GLU hypo GABA mCG DMRs", "GLU hypo GABA hmCG DMRs"),
           "GABA hypo GLU" = c("GABA hypo GLU mCG DMRs", "GABA hypo GLU hmCG DMRs"),
           "OLIG hypo neuron" = c("OLIG hypo neuron mCG DMRs", "OLIG hypo neuron hmCG DMRs"))


table = list()
for (i in 1:length(dmrs)) {
  table[[i]] = list()
  for (j in 1:length(set)) {
    table[[i]][[j]] = data.frame(c(nrow(DMRs[[i]][which(DMRs[[i]][,which(colnames(DMRs[[i]])==set[[j]][1])]=="yes" & 
                                                           DMRs[[i]][,which(colnames(DMRs[[i]])==set[[j]][2])]=="yes"),]),
                                    nrow(DMRs[[i]][which(DMRs[[i]][,which(colnames(DMRs[[i]])==set[[j]][1])]=="yes" & 
                                                           DMRs[[i]][,which(colnames(DMRs[[i]])==set[[j]][2])]=="no"),])),
                                  c(nrow(DMRs[[i]][which(DMRs[[i]][,which(colnames(DMRs[[i]])==set[[j]][1])]=="no" & 
                                                           DMRs[[i]][,which(colnames(DMRs[[i]])==set[[j]][2])]=="yes"),]),
                                    nrow(DMRs[[i]][which(DMRs[[i]][,which(colnames(DMRs[[i]])==set[[j]][1])]=="no" & 
                                                           DMRs[[i]][,which(colnames(DMRs[[i]])==set[[j]][2])]=="no"),])),
                                  row.names = c("Yes.hmCG", "No.hmCG"))
    if (j %in% c(1:3)) { colnames(table[[i]][[j]]) = c("Yes.tmCG", "No.tmCG") }
    else { colnames(table[[i]][[j]]) = c("Yes.mCG", "No.mCG") }
    }
  names(table[[i]]) = names(set)
}
names(table) = names(dmrs)


res = lapply(table, function(x) lapply(x, fisher.test))
df = do.call(rbind, Map(cbind, DMRgroup = as.list(names(res)), 
                        lapply(lapply(res, function(x) Map(cbind, Comparison = as.list(c(rep.int("hmCG.vs.tmCG", 3), rep.int("hmCG.vs.mCG", 3))),
                                                            lapply(x, function(y) data.frame(pval = y$p.value, OR = y$estimate)))), 
                               function(x) do.call(rbind, Map(cbind, PaperComparison = as.list(names(set)), x)))))
df$FDR = p.adjust(df$pval, method = "fdr")
rownames(df) = NULL

write.csv(df,file="/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/DMR/Kozlenkov_ScienceAdvances_2018_hmCG_vs_mCG_or_tmCG_fishertest.csv", quote=F)

df[df$FDR<=0.05,]
#     Comparison       DMRgroup PaperComparison          pval        OR
#1  hmCG.vs.tmCG Group 1 (G-N+)   GLU hypo GABA  1.173035e-11 17.977691
#4  hmCG.vs.tmCG Group 2 (G0N+)   GLU hypo GABA  4.551132e-03 12.947219
#7  hmCG.vs.tmCG Group 3 (G0N-)   GLU hypo GABA  1.657188e-10  4.235920
#10 hmCG.vs.tmCG Group 4 (G+N0)   GLU hypo GABA  9.329989e-04  9.234284
#13 hmCG.vs.tmCG Group 5 (G+N-)   GLU hypo GABA  2.505506e-06  4.985189
#16 hmCG.vs.tmCG Group 6 (G-N0)   GLU hypo GABA  1.254871e-12 21.025757
#17 hmCG.vs.tmCG Group 6 (G-N0)   GABA hypo GLU  9.775109e-04       Inf
#19 hmCG.vs.tmCG         cdDMRs   GLU hypo GABA 7.781311e-103 17.812102
#20 hmCG.vs.tmCG         cdDMRs   GABA hypo GLU  5.860258e-03 30.753890
#             FDR
#1   1.642249e-10
#4   2.389344e-02
#7   1.740047e-09
#10  5.865065e-03
#13  2.104625e-05
#16  2.635229e-11
#17  5.865065e-03
#19 3.268151e-101
#20  2.734787e-02

table$cdDMRs
#$`GLU hypo GABA`
#         Yes.tmCG No.tmCG
#Yes.hmCG      220     102
#No.hmCG       200    1656

#$`GABA hypo GLU`
#         Yes.tmCG No.tmCG
#Yes.hmCG        2       2
#No.hmCG        68    2106

#$`OLIG hypo neuron`
#         Yes.tmCG No.tmCG
#Yes.hmCG        5       1
#No.hmCG       879    1293