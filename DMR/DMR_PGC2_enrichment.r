library(GenomicRanges)
library(jaffelab)
library(bumphunter)

load('/dcl01/lieber/ajaffe/lab/brain-epigenomics/bumphunting/BSobj_bsseqSmooth_Neuron_minCov_3.Rdata')
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/DMR/DMR_objects.rda")


## create pgc2 loci granges object

pgc = read.delim("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/tableS3_Nature2014_pgc2_loci.txt",as.is=TRUE)
pgc$chr = ss(pgc$Position..hg19.,":")
tmp = ss(pgc$Position..hg19.,":",2)
pgc$start = as.numeric(ss(tmp, "-"))
pgc$end = as.numeric(ss(tmp, "-",2))
pgcGR = GRanges(pgc$chr, IRanges(pgc$start,pgc$end))


## DMR granges

DMRgr = lapply(DMR, function(x) makeGRangesFromDataFrame(x, keep.extra.columns = T))


## Test for enrichment of DMRs in PGC loci by CpG clusters

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

pgcOverlap = findOverlaps(pgcGR, gr.clusters)

df.clusters$regionID = paste0(df.clusters$seqnames,":",df.clusters$start,"-",df.clusters$end)
df.clusters$rnum = 1:length(gr.clusters)
df.clusters$CellType = ifelse(df.clusters$rnum %in% queryHits(oo$CellType), "CellType","no")
df.clusters$Age = ifelse(df.clusters$rnum %in% queryHits(oo$Age), "Age","no")
df.clusters$Interaction = ifelse(df.clusters$rnum %in% queryHits(oo$Interaction), "Interaction","no")
df.clusters$DMRs = paste(df.clusters$CellType, df.clusters$Age, df.clusters$Interaction, sep=":")
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
                                                 nrow(df.clusters[df.clusters$PGC=="no" & df.clusters$Interaction=="no",])), row.names = c("YesInt","NoInt")))

fisher = lapply(tables, fisher.test)
df = data.frame(pval = unlist(lapply(fisher, function(x) x$p.value)),
                OR = unlist(lapply(fisher, function(x) x$estimate)), row.names = names(fisher))

#         CellType        Age  Interaction
#pval 1.594936e-08 0.01686004 2.125739e-07
#OR   1.817179e+00 4.22238034 2.684525e+00


## Test for enrichment of DMRs in PGC loci genes

# Get genes in the PGC2 loci

geneMapgr = makeGRangesFromDataFrame(geneMap, keep.extra.columns = T)
oo = findOverlaps(pgcGR, geneMapgr)
PGCgenes = geneMap[subjectHits(oo),]


## Enrichment in DMRs by cell type, age and interaction

geneuniverse = unique(geneMap$gencodeID)
sig = lapply(DMR, function(x) unique(x[which(x$fwer<=0.05),"nearestID"]))
notsig = lapply(sig, function(x) geneuniverse[!(geneuniverse %in% x)])

DMRenrich = mapply(function(sig,notsig) {
  DE_OVERLAP = c( sum( sig %in% PGCgenes$gencodeID),sum(!(sig %in% PGCgenes$gencodeID)))
  NOT_DE_OVERLAP= c(sum(notsig %in% PGCgenes$gencodeID), sum(!(notsig %in% PGCgenes$gencodeID)))
  enrich_table = cbind(DE_OVERLAP, NOT_DE_OVERLAP)
  res = fisher.test(enrich_table)
  dat=c(res$p.value, res$estimate)
  names(dat) <- c("P.Value","Odds Ratio")
  return(dat)
}, sig, notsig, SIMPLIFY = F) 

data.frame(do.call(rbind, DMRenrich), row.names=names(DMRenrich))
#                 P.Value Odds.Ratio
#CellType    2.512262e-06   1.694533
#Age         3.665844e-02   3.308604
#Interaction 7.462028e-05   2.087669