library('RColorBrewer')
library(ggplot2)
library(GenomicRanges)

load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/DMR/DMR_objects.rda")
load('/dcl01/lieber/ajaffe/lab/brain-epigenomics/bumphunting/BSobj_bsseqSmooth_Neuron_minCov_3.Rdata')

# Identify all CpG clusters in the genome
gr <- granges(BSobj)
cl = clusterMaker( chr = as.character(seqnames(gr)), 
                   pos = start(gr),   maxGap = 1000)
gr.clusters = split(gr, cl)
gr.clusters = unlist(range(gr.clusters))
df.clusters = as.data.frame(gr.clusters)

# Find overlaps with DMRS in all three models
DMRgr = lapply(DMR, function(x) makeGRangesFromDataFrame(x[which(x$fwer<=0.05),], keep.extra.columns=T))

oo = lapply(DMRgr, function(x) findOverlaps(gr.clusters, x))
lapply(oo, function(x) length(unique(queryHits(x))))

df.clusters$regionID = paste0(df.clusters$seqnames,":",df.clusters$start,"-",df.clusters$end)
df.clusters$rnum = 1:length(gr.clusters)
df.clusters$CellType = ifelse(df.clusters$rnum %in% queryHits(oo$CellType), "CellType","no")
df.clusters$Age = ifelse(df.clusters$rnum %in% queryHits(oo$Age), "Age","no")
df.clusters$Interaction = ifelse(df.clusters$rnum %in% queryHits(oo$Interaction), "Interaction","no")
df.clusters$DMRs = paste(df.clusters$CellType, df.clusters$Age, df.clusters$Interaction, sep=":")

CellTypexAge = findOverlaps(DMRgr$CellType, DMRgr$Age)
CellTypexInt = findOverlaps(DMRgr$CellType, DMRgr$Interaction)
AgexInt = findOverlaps(DMRgr$Age, DMRgr$Interaction)

length(unique(queryHits(CellTypexAge))) # 55 of 11179 Cell type overlap age
length(unique(subjectHits(CellTypexAge))) # 52 of 129 age overlap cell type
length(unique(queryHits(CellTypexInt))) # 775 of 11179 cell type overlap interaction
length(unique(subjectHits(CellTypexInt))) # 799 of 2178 interaction overlap cell type
length(unique(queryHits(AgexInt))) # 59 of 129 age overlap interaction
length(unique(subjectHits(AgexInt))) # 66 of 2178 interaction overlap age

## make contingency tables

tables = list(CellTypexAge = data.frame(YesCT = c(nrow(df.clusters[df.clusters$CellType=="CellType" & df.clusters$Age=="Age",]),
                                                  nrow(df.clusters[df.clusters$CellType=="CellType" & df.clusters$Age=="no",])),
                                        NoCT = c(nrow(df.clusters[df.clusters$CellType=="no" & df.clusters$Age=="Age",]),
                                                 nrow(df.clusters[df.clusters$CellType=="no" & df.clusters$Age=="no",])), row.names = c("YesAge","NoAge")),
              CellTypexInt = data.frame(YesCT = c(nrow(df.clusters[df.clusters$CellType=="CellType" & df.clusters$Interaction=="Interaction",]),
                                                  nrow(df.clusters[df.clusters$CellType=="CellType" & df.clusters$Interaction=="no",])),
                                        NoCT = c(nrow(df.clusters[df.clusters$CellType=="no" & df.clusters$Interaction=="Interaction",]),
                                                 nrow(df.clusters[df.clusters$CellType=="no" & df.clusters$Interaction=="no",])), row.names = c("YesInt","NoInt")),
              AgexInt = data.frame(YesAge = c(nrow(df.clusters[df.clusters$Age=="Age" & df.clusters$Interaction=="Interaction",]),
                                              nrow(df.clusters[df.clusters$Age=="Age" & df.clusters$Interaction=="no",])),
                                   NoAge = c(nrow(df.clusters[df.clusters$Age=="no" & df.clusters$Interaction=="Interaction",]),
                                             nrow(df.clusters[df.clusters$Age=="no" & df.clusters$Interaction=="no",])), row.names = c("YesInt","NoInt")))

fisher = lapply(tables, fisher.test)
df = t(data.frame(pval = unlist(lapply(fisher, function(x) x$p.value)),
                OR = unlist(lapply(fisher, function(x) x$estimate))))

#     CellTypexAge CellTypexInt       AgexInt
#pval 2.390711e-119      0.00000 5.161157e-152
#OR   1.334017e+02     83.48967  4.590351e+02





