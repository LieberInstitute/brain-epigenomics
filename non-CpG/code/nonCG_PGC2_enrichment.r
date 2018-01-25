library(GenomicRanges)
library(jaffelab)
library(bumphunter)

load('/dcl01/lieber/ajaffe/lab/brain-epigenomics/bumphunting/BSobj_bsseqSmooth_Neuron_minCov_3.Rdata')
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/non-CpG/CH_object.rda")
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/non-CpG/CHneurons_object.rda")


## create pgc2 loci granges object

pgc = read.delim("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/tableS3_Nature2014_pgc2_loci.txt",as.is=TRUE)
pgc$chr = ss(pgc$Position..hg19.,":")
tmp = ss(pgc$Position..hg19.,":",2)
pgc$start = as.numeric(ss(tmp, "-"))
pgc$end = as.numeric(ss(tmp, "-",2))
pgcGR = GRanges(pgc$chr, IRanges(pgc$start,pgc$end))


## CH granges

CHgr = list("CellType" = makeGRangesFromDataFrame(CH[which(CH$padj.CellType<=0.01),], keep.extra.columns = T),
            "Age" = makeGRangesFromDataFrame(CH[which(CH$padj.Age<=0.01),], keep.extra.columns = T), 
            "Interaction" = makeGRangesFromDataFrame(CH[which(CH$padj.Interaction<=0.01),], keep.extra.columns = T),
            "Age (Neurons)" = makeGRangesFromDataFrame(CHneurons[which(CHneurons$padj<=0.01),], keep.extra.columns = T))


## Test for enrichment of dmCH in PGC loci by CpG clusters

# Identify all CpG clusters in the genome

gr <- granges(BSobj)
cl = clusterMaker( chr = as.character(seqnames(gr)), 
                   pos = start(gr),   maxGap = 1000)
gr.clusters = split(gr, cl)
gr.clusters = unlist(range(gr.clusters))
df.clusters = as.data.frame(gr.clusters)


# Find overlaps with CH in all models

oo = lapply(CHgr, function(x) findOverlaps(gr.clusters, x))
lapply(oo, function(x) length(unique(queryHits(x))))

pgcOverlap = findOverlaps(pgcGR, gr.clusters)

df.clusters$regionID = paste0(df.clusters$seqnames,":",df.clusters$start,"-",df.clusters$end)
df.clusters$rnum = 1:length(gr.clusters)
df.clusters$CellType = ifelse(df.clusters$rnum %in% queryHits(oo$CellType), "CellType","no")
df.clusters$Age = ifelse(df.clusters$rnum %in% queryHits(oo$Age), "Age","no")
df.clusters$Interaction = ifelse(df.clusters$rnum %in% queryHits(oo$Interaction), "Interaction","no")
df.clusters$DMRs = paste(df.clusters$CellType, df.clusters$Age, df.clusters$Interaction, sep=":")
df.clusters$AgeNeurons = ifelse(df.clusters$rnum %in% queryHits(oo$"Age (Neurons)"), "AgeNeurons","no")
df.clusters$PGC = ifelse(df.clusters$rnum %in% subjectHits(pgcOverlap), "PGC","no")


# make contingency tables

tables = list(CellType = data.frame(YesPGC = c(nrow(df.clusters[df.clusters$PGC=="PGC" & df.clusters$CellType=="CellType",]),
                                               nrow(df.clusters[df.clusters$PGC=="PGC" & df.clusters$CellType=="no",])),
                                    NoPGC = c(nrow(df.clusters[df.clusters$PGC=="no" & df.clusters$CellType=="CellType",]),
                                              nrow(df.clusters[df.clusters$PGC=="no" & df.clusters$CellType=="no",])), row.names = c("YesCT","NoCT")),
              Age = data.frame(YesPGC = c(nrow(df.clusters[df.clusters$PGC=="PGC" & df.clusters$Age=="Age",]),
                                          nrow(df.clusters[df.clusters$PGC=="PGC" & df.clusters$Age=="no",])),
                               NoPGC = c(nrow(df.clusters[df.clusters$PGC=="no" & df.clusters$Age=="Age",]),
                                         nrow(df.clusters[df.clusters$PGC=="no" & df.clusters$Age=="no",])), row.names = c("YesAge","NoAge")),
              Interaction = data.frame(YesPGC = c(nrow(df.clusters[df.clusters$PGC=="PGC" & df.clusters$Interaction=="Interaction",]),
                                                  nrow(df.clusters[df.clusters$PGC=="PGC" & df.clusters$Interaction=="no",])),
                                       NoPGC = c(nrow(df.clusters[df.clusters$PGC=="no" & df.clusters$Interaction=="Interaction",]),
                                                 nrow(df.clusters[df.clusters$PGC=="no" & df.clusters$Interaction=="no",])), row.names = c("YesInt","NoInt")),
              AgeNeurons = data.frame(YesPGC = c(nrow(df.clusters[df.clusters$PGC=="PGC" & df.clusters$AgeNeurons=="AgeNeurons",]),
                                                  nrow(df.clusters[df.clusters$PGC=="PGC" & df.clusters$AgeNeurons=="no",])),
                                       NoPGC = c(nrow(df.clusters[df.clusters$PGC=="no" & df.clusters$AgeNeurons=="AgeNeurons",]),
                                                 nrow(df.clusters[df.clusters$PGC=="no" & df.clusters$AgeNeurons=="no",])), row.names = c("YesAgeNeuron","NoAgeNeuron")))

fisher = lapply(tables, fisher.test)
df = data.frame(pval = unlist(lapply(fisher, function(x) x$p.value)),
                OR = unlist(lapply(fisher, function(x) x$estimate)), row.names = names(fisher))

#                    pval        OR
#CellType    1.264233e-11 0.7943156
#Age         3.690218e-20 0.6793588
#Interaction 1.000000e+00 0.0000000
#AgeNeurons  8.784457e-18 0.7194306


## Test for enrichment of DMRs in PGC loci genes

# Get genes in the PGC2 loci

geneMapgr = makeGRangesFromDataFrame(geneMap, keep.extra.columns = T)
oo = findOverlaps(pgcGR, geneMapgr)
PGCgenes = geneMap[subjectHits(oo),]


## Enrichment in CH by cell type, age and interaction

geneuniverse = unique(CH$nearestID)
sigCH = list("CellType" = unique(CH[which(CH$padj.CellType<=0.01),"nearestID"]),
             "Age" = unique(CH[which(CH$padj.Age<=0.01),"nearestID"]), 
             "Interaction" = unique(CH[which(CH$padj.Interaction<=0.01),"nearestID"]))
notsigCH = lapply(sigCH, function(x) geneuniverse[!(geneuniverse %in% x)])

CHenrich = mapply(function(sig,notsig) {
  DE_OVERLAP = c( sum( sig %in% PGCgenes$gencodeID),sum(!(sig %in% PGCgenes$gencodeID)))
  NOT_DE_OVERLAP= c(sum(notsig %in% PGCgenes$gencodeID), sum(!(notsig %in% PGCgenes$gencodeID)))
  enrich_table = cbind(DE_OVERLAP, NOT_DE_OVERLAP)
  res = fisher.test(enrich_table)
  dat=c(res$p.value, res$estimate)
  names(dat) <- c("P.Value","Odds Ratio")
  return(dat)
}, sigCH, notsigCH, SIMPLIFY = F) 

# Enrichment by age in neurons

geneuniverse = unique(CHneurons$nearestID)
sig = unique(CHneurons[which(CHneurons$padj<=0.01),"nearestID"])
notsig = geneuniverse[!(geneuniverse %in% sig)]

enrich_table = cbind(c( sum( sig %in% PGCgenes$gencodeID),sum(!(sig %in% PGCgenes$gencodeID))),
                     c(sum(notsig %in% PGCgenes$gencodeID), sum(!(notsig %in% PGCgenes$gencodeID))))
res = fisher.test(enrich_table)
CHneuronEnrich=c(res$p.value, res$estimate)
names(CHneuronEnrich) <- c("P.Value","Odds Ratio")


data.frame(do.call(rbind, c(CHenrich, list(CHneuronEnrich))), row.names=c(names(CHenrich), "AgeNeurons"))
#                 P.Value Odds.Ratio
#CellType    1.398847e-02  0.7111053
#Age         8.031735e-09  0.5662367
#Interaction 1.000000e+00  0.0000000
#AgeNeurons  5.142322e-06  0.6116810