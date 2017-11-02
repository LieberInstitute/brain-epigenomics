library(GenomicRanges)
library(bumphunter)

load('/dcl01/lieber/ajaffe/lab/brain-epigenomics/bumphunting/BSobj_bsseqSmooth_Neuron_minCov_3.Rdata')
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/DMR/DMR_objects.rda")

# load HARs
HARs = openxlsx::read.xlsx('/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/HARs_hg19_Doan_Walsh_Table_S1.xlsx')
hars = makeGRangesFromDataFrame(HARS, keep.extra.columns=T)
length(hars) # 2737

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

harOverlap = findOverlaps(hars, gr.clusters)

df.clusters$regionID = paste0(df.clusters$seqnames,":",df.clusters$start,"-",df.clusters$end)
df.clusters$rnum = 1:length(gr.clusters)
df.clusters$CellType = ifelse(df.clusters$rnum %in% queryHits(oo$CellType), "CellType","no")
df.clusters$Age = ifelse(df.clusters$rnum %in% queryHits(oo$Age), "Age","no")
df.clusters$Interaction = ifelse(df.clusters$rnum %in% queryHits(oo$Interaction), "Interaction","no")
df.clusters$DMRs = paste(df.clusters$CellType, df.clusters$Age, df.clusters$Interaction, sep=":")
df.clusters$HARs = ifelse(df.clusters$rnum %in% subjectHits(harOverlap), "HAR","no")


## make contingency tables

tables = list(CellType = data.frame(YesHAR = c(nrow(df.clusters[df.clusters$HARs=="HAR" & df.clusters$CellType=="CellType",]),
                                                  nrow(df.clusters[df.clusters$HARs=="HAR" & df.clusters$CellType=="no",])),
                                    NoHAR = c(nrow(df.clusters[df.clusters$HARs=="no" & df.clusters$CellType=="CellType",]),
                                                 nrow(df.clusters[df.clusters$HARs=="no" & df.clusters$CellType=="no",])), row.names = c("YesCT","NoCT")),
              Age = data.frame(YesHAR = c(nrow(df.clusters[df.clusters$HARs=="HAR" & df.clusters$Age=="Age",]),
                                                  nrow(df.clusters[df.clusters$HARs=="HAR" & df.clusters$Age=="no",])),
                               NoHAR = c(nrow(df.clusters[df.clusters$HARs=="no" & df.clusters$Age=="Age",]),
                                                 nrow(df.clusters[df.clusters$HARs=="no" & df.clusters$Age=="no",])), row.names = c("YesAge","NoAge")),
              Interaction = data.frame(YesHAR = c(nrow(df.clusters[df.clusters$HARs=="HAR" & df.clusters$Interaction=="Interaction",]),
                                              nrow(df.clusters[df.clusters$HARs=="HAR" & df.clusters$Interaction=="no",])),
                                       NoHAR = c(nrow(df.clusters[df.clusters$HARs=="no" & df.clusters$Interaction=="Interaction",]),
                                             nrow(df.clusters[df.clusters$HARs=="no" & df.clusters$Interaction=="no",])), row.names = c("YesInt","NoInt")))

fisher = lapply(tables, fisher.test)
df = t(data.frame(pval = unlist(lapply(fisher, function(x) x$p.value)),
                  OR = unlist(lapply(fisher, function(x) x$estimate))))

#          CellType          Age  Interaction
#pval 5.160241e-150 1.621188e-04 1.128312e-37
#OR    9.657996e+00 1.038447e+01 9.253464e+00

df = t(data.frame(pval = round(unlist(lapply(fisher, function(x) x$p.value)),3),
                  OR = round(unlist(lapply(fisher, function(x) x$estimate)),3)))

#     CellType    Age Interaction
#pval    0.000  0.000       0.000
#OR      9.658 10.384       9.253









## Enrichment in DMRs by cell type, age and interaction

geneuniverse = lapply(DMR, function(x) na.omit(unique(x$EntrezID)))
aej_sets_expressed = lapply(geneuniverse, function(x) aej_sets[which(aej_sets$EntrezGene.ID %in% x), ]) # drop genes that are not present in the test set
splitSets = lapply(aej_sets_expressed, function(x) split(x, x$Gene.Set))
sig = lapply(DMR, function(x) na.omit(unique(x[which(x$fwer<=0.05),"EntrezID"])))
notsig = mapply(function(x,y) x[!(x %in% y)], geneuniverse, sig)

DMRenrich = mapply(function(sig,notsig,sets) sapply(sets, function(x) {
  DE_OVERLAP = c( sum( sig %in% x$EntrezGene.ID),sum(!(sig %in% x$EntrezGene.ID)))
  NOT_DE_OVERLAP= c(sum(notsig %in% x$EntrezGene.ID), sum(!(notsig %in% x$EntrezGene.ID)))
  enrich_table = cbind(DE_OVERLAP, NOT_DE_OVERLAP)
  res = fisher.test(enrich_table)
  dat=c(res$p.value, res$estimate)
  names(dat) <- c("P.Value","Odds Ratio")
  return(dat)
}), sig, notsig,splitSets,SIMPLIFY =F) 

mapply(function(x,y,z) write.csv(x,file=paste0("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/DMR/",y,"/Birnbaum_geneSet_enrichment_DMR_",z,".csv"),
                                                 quote=F),DMRenrich,c("CellType","Age","CT_Age_Interaction"),names(DMRenrich))
                   
