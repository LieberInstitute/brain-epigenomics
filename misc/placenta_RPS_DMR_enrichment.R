library(GenomicRanges)
library(bumphunter)
library(SNPlocs.Hsapiens.dbSNP144.GRCh37)

load('/dcl01/lieber/ajaffe/lab/brain-epigenomics/bumphunting/BSobj_bsseqSmooth_Neuron_minCov_3.Rdata')
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/DMR/DMR_objects.rda")
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/bumphunting/rda/limma_Neuron_CpGs_minCov_3_ageInfo_dmrs.Rdata")


# Load SNPs

plac = scan("/dcl01/lieber/ajaffe/lab/brain-epigenomics/misc/rps.snpZ4D.p5e-08.snp", what = "character")
nonplac = scan("/dcl01/lieber/ajaffe/lab/brain-epigenomics/misc/rps.snpZ4Dexc.p5e-08.snp", what = "character")

plac = snpsById(SNPlocs.Hsapiens.dbSNP144.GRCh37, plac)
nonplac = snpsById(SNPlocs.Hsapiens.dbSNP144.GRCh37, nonplac)
plac = as.data.frame(plac)
plac$seqnames = paste0("chr", plac$seqnames)
plac = makeGRangesFromDataFrame(plac, start.field = "pos", end.field = "pos", keep.extra.columns = T)
nonplac = as.data.frame(nonplac)
nonplac$seqnames = paste0("chr", nonplac$seqnames)
nonplac = makeGRangesFromDataFrame(nonplac, start.field = "pos", end.field = "pos", keep.extra.columns = T)


### test LD boundaries, +/- 500kb, and +/- 1Mb

plac500 = as.data.frame(plac)
plac500$start = plac500$start-500000
plac500$end = plac500$end+500000
plac500 = makeGRangesFromDataFrame(plac500, keep.extra.columns = T)
plac1000 = as.data.frame(plac)
plac1000$start = plac1000$start-1000000
plac1000$end = plac1000$end+1000000
plac1000 = makeGRangesFromDataFrame(plac1000, keep.extra.columns = T)

nonplac500 = as.data.frame(nonplac)
nonplac500$start = nonplac500$start-500000
nonplac500$end = nonplac500$end+500000
nonplac500 = makeGRangesFromDataFrame(nonplac500, keep.extra.columns = T)
nonplac1000 = as.data.frame(nonplac)
nonplac1000$start = nonplac1000$start-1000000
nonplac1000$end = nonplac1000$end+1000000
nonplac1000 = makeGRangesFromDataFrame(nonplac1000, keep.extra.columns = T)


## DMR granges

DMRgr = lapply(DMR, function(x) makeGRangesFromDataFrame(x[which(x$sig=="FWER < 0.05"),], keep.extra.columns = T))

dmrs = split(dmrs, dmrs$k6cluster_label)
names(dmrs) = c("Gr1","Gr2","Gr3","Gr4","Gr5","Gr6")
oo = lapply(dmrs, function(x) findOverlaps(x, makeGRangesFromDataFrame(DMR$Interaction)))
dmrs = lapply(oo, function(x) DMR$Interaction[subjectHits(x),])

DMRgr = c(DMRgr, list("CellType.N+" = DMRgr$CellType[which(DMRgr$CellType$Dir=="pos" & DMRgr$CellType$sig=="FWER < 0.05"),],
                      "CellType.N-" = DMRgr$CellType[which(DMRgr$CellType$Dir=="neg" & DMRgr$CellType$sig=="FWER < 0.05"),]),
          lapply(dmrs, function(x) makeGRangesFromDataFrame(x[which(x$sig=="FWER < 0.05"),], keep.extra.columns = T)))
DMRgr = lapply(DMRgr, reduce)


## Test for enrichment of DMRs in regions by CpG clusters

# Identify all CpG clusters in the genome
gr <- granges(BSobj)
cl = clusterMaker( chr = as.character(seqnames(gr)), 
                   pos = start(gr),   maxGap = 1000)
gr.clusters = split(gr, cl)
gr.clusters = unlist(range(gr.clusters))
df.clusters = as.data.frame(gr.clusters)


# Find overlaps with DMRS in all three models

oo = lapply(DMRgr, function(x) findOverlaps(gr.clusters, x))
lapply(oo, function(x) length(unique(queryHits(x))))

df.clusters$regionID = paste0(df.clusters$seqnames,":",df.clusters$start,"-",df.clusters$end)
df.clusters$rnum = 1:length(gr.clusters)
df.clusters$CellType = ifelse(df.clusters$rnum %in% queryHits(oo$CellType), 1,0)
df.clusters$"CellType.N+" = ifelse(df.clusters$rnum %in% queryHits(oo$"CellType.N+"), 1,0)
df.clusters$"CellType.N-" = ifelse(df.clusters$rnum %in% queryHits(oo$"CellType.N-"), 1,0)
df.clusters$Age = ifelse(df.clusters$rnum %in% queryHits(oo$Age), 1,0)
df.clusters$Interaction = ifelse(df.clusters$rnum %in% queryHits(oo$Interaction), 1,0)
df.clusters$Gr1 = ifelse(df.clusters$rnum %in% queryHits(oo$Gr1), 1,0)
df.clusters$Gr2 = ifelse(df.clusters$rnum %in% queryHits(oo$Gr2), 1,0)
df.clusters$Gr3 = ifelse(df.clusters$rnum %in% queryHits(oo$Gr3), 1,0)
df.clusters$Gr4 = ifelse(df.clusters$rnum %in% queryHits(oo$Gr4), 1,0)
df.clusters$Gr5 = ifelse(df.clusters$rnum %in% queryHits(oo$Gr5), 1,0)
df.clusters$Gr6 = ifelse(df.clusters$rnum %in% queryHits(oo$Gr6), 1,0)
df.clusters$DMRs = paste(df.clusters$CellType, df.clusters$Age, df.clusters$Interaction,
                         df.clusters$Gr1,df.clusters$Gr2,df.clusters$Gr3,df.clusters$Gr4,df.clusters$Gr5,df.clusters$Gr6, sep=":")
ooG = findOverlaps(gr.clusters, makeGRangesFromDataFrame(geneMap))
df.clusters$Gene = ifelse(df.clusters$rnum %in% queryHits(ooG), 1,0)
ooPC = findOverlaps(gr.clusters, makeGRangesFromDataFrame(geneMap[which(geneMap$gene_type=="protein_coding"),]))
df.clusters$ProteinCoding = ifelse(df.clusters$rnum %in% queryHits(ooPC), 1,0)

ooP = findOverlaps(plac1000, gr.clusters)
ooNP = findOverlaps(nonplac1000, gr.clusters)
df.clusters$plac1000 = ifelse(df.clusters$rnum %in% subjectHits(ooP), 1,0)
df.clusters$nonplac1000 = ifelse(df.clusters$rnum %in% subjectHits(ooNP), 1,0)

ooP = findOverlaps(plac500, gr.clusters)
ooNP = findOverlaps(nonplac500, gr.clusters)
df.clusters$plac500 = ifelse(df.clusters$rnum %in% subjectHits(ooP), 1,0)
df.clusters$nonplac500 = ifelse(df.clusters$rnum %in% subjectHits(ooNP), 1,0)


# make contingency tables

cols = c("CellType","CellType.N+","CellType.N-","Age","Interaction","Gr1","Gr2","Gr3","Gr4","Gr5","Gr6")
tables = list()
for (i in 1:length(cols)) {
  tables[[i]] = list(plac500kb = data.frame(YesPGC = c(nrow(df.clusters[df.clusters$plac500=="plac500" & df.clusters[,colnames(df.clusters)==cols[i]]!="no",]),
                                                       nrow(df.clusters[df.clusters$plac500=="plac500" & df.clusters[,colnames(df.clusters)==cols[i]]=="no",])),
                                            NoPGC = c(nrow(df.clusters[df.clusters$plac500=="no" & df.clusters[,colnames(df.clusters)==cols[i]]!="no",]),
                                                      nrow(df.clusters[df.clusters$plac500=="no" & df.clusters[,colnames(df.clusters)==cols[i]]=="no",])), 
                                            row.names = c("YesDMR","NoDMR")),
                     nonplac500kb = data.frame(YesPGC = c(nrow(df.clusters[df.clusters$nonplac500=="nonplac500" & df.clusters[,colnames(df.clusters)==cols[i]]!="no",]),
                                                          nrow(df.clusters[df.clusters$nonplac500=="nonplac500" & df.clusters[,colnames(df.clusters)==cols[i]]=="no",])),
                                               NoPGC = c(nrow(df.clusters[df.clusters$nonplac500=="no" & df.clusters[,colnames(df.clusters)==cols[i]]!="no",]),
                                                         nrow(df.clusters[df.clusters$nonplac500=="no" & df.clusters[,colnames(df.clusters)==cols[i]]=="no",])), 
                                               row.names = c("YesDMR","NoDMR")),
                     plac1Mb = data.frame(YesPGC = c(nrow(df.clusters[df.clusters$plac1000=="plac1000" & df.clusters[,colnames(df.clusters)==cols[i]]!="no",]),
                                                     nrow(df.clusters[df.clusters$plac1000=="plac1000" & df.clusters[,colnames(df.clusters)==cols[i]]=="no",])),
                                          NoPGC = c(nrow(df.clusters[df.clusters$plac1000=="no" & df.clusters[,colnames(df.clusters)==cols[i]]!="no",]),
                                                    nrow(df.clusters[df.clusters$plac1000=="no" & df.clusters[,colnames(df.clusters)==cols[i]]=="no",])), 
                                          row.names = c("YesDMR","NoDMR")),
                     nonplac1Mb = data.frame(YesPGC = c(nrow(df.clusters[df.clusters$nonplac1000=="nonplac1000" & df.clusters[,colnames(df.clusters)==cols[i]]!="no",]),
                                                        nrow(df.clusters[df.clusters$nonplac1000=="nonplac1000" & df.clusters[,colnames(df.clusters)==cols[i]]=="no",])),
                                             NoPGC = c(nrow(df.clusters[df.clusters$nonplac1000=="no" & df.clusters[,colnames(df.clusters)==cols[i]]!="no",]),
                                                       nrow(df.clusters[df.clusters$nonplac1000=="no" & df.clusters[,colnames(df.clusters)==cols[i]]=="no",])), 
                                             row.names = c("YesDMR","NoDMR")),
                     placVSnonplac500kb = data.frame(YesPGC = c(nrow(df.clusters[df.clusters$plac500=="plac500" & df.clusters[,colnames(df.clusters)==cols[i]]!="no",]),
                                                                nrow(df.clusters[df.clusters$plac500=="plac500" & df.clusters[,colnames(df.clusters)==cols[i]]=="no",])),
                                                     NoPGC = c(nrow(df.clusters[df.clusters$nonplac500=="nonplac500" & df.clusters[,colnames(df.clusters)==cols[i]]!="no",]),
                                                               nrow(df.clusters[df.clusters$nonplac500=="nonplac500" & df.clusters[,colnames(df.clusters)==cols[i]]=="no",])), 
                                                     row.names = c("YesDMR","NoDMR")),
                     placVSnonplac1Mb = data.frame(YesPGC = c(nrow(df.clusters[df.clusters$plac1000=="plac1000" & df.clusters[,colnames(df.clusters)==cols[i]]!="no",]),
                                                              nrow(df.clusters[df.clusters$plac1000=="plac1000" & df.clusters[,colnames(df.clusters)==cols[i]]=="no",])),
                                                   NoPGC = c(nrow(df.clusters[df.clusters$nonplac1000=="nonplac1000" & df.clusters[,colnames(df.clusters)==cols[i]]!="no",]),
                                                             nrow(df.clusters[df.clusters$nonplac1000=="nonplac1000" & df.clusters[,colnames(df.clusters)==cols[i]]=="no",])), 
                                                   row.names = c("YesDMR","NoDMR")))
}
names(tables) = cols
save(tables, file="/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/DMR/placenta_SNPregions_enrichment_DMRs_tables.rda")
lapply(tables, function(x) lapply(x, function(y) sum(rowSums(y))))


fisher = lapply(tables, function(x) lapply(x,fisher.test))
df = do.call(rbind, Map(cbind, SigGroup = as.list(names(fisher)), lapply(fisher, function(y) 
  do.call(rbind, Map(cbind, SNPgroup = as.list(names(y)), lapply(y, function(x) data.frame(pval = x$p.value, OddsRatio = x$estimate)))))))
df$FDR = p.adjust(df$pval, method = "fdr")
rownames(df) = NULL
df$SigGroup = gsub("Gr1", "G-:N+", df$SigGroup)
df$SigGroup = gsub("Gr2", "G0:N+", df$SigGroup)
df$SigGroup = gsub("Gr3", "G0:N-", df$SigGroup)
df$SigGroup = gsub("Gr4", "G+:N0", df$SigGroup)
df$SigGroup = gsub("Gr5", "G+:N-", df$SigGroup)
df$SigGroup = gsub("Gr6", "G-:N0", df$SigGroup)

write.csv(df,file="/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/DMR/placenta_SNPregions_enrichment_DMRs.csv",quote=F)

df[which(df$FDR<=0.05),]


## How many genes fall in each category?

length(unique(geneMap$gencodeID)) # 60252

ooP5 = findOverlaps(plac500, makeGRangesListFromDataFrame(geneMap))
ooP1 = findOverlaps(plac1000, makeGRangesListFromDataFrame(geneMap))
ooNP5 = findOverlaps(nonplac500, makeGRangesListFromDataFrame(geneMap))
ooNP1 = findOverlaps(nonplac1000, makeGRangesListFromDataFrame(geneMap))

length(unique(subjectHits(ooP5))) # 1612
length(unique(subjectHits(ooP1))) # 2769
length(unique(subjectHits(ooNP5))) # 519
length(unique(subjectHits(ooNP1))) # 1036


## anchor in genes rather than clusters: nearest gene

geneuniverse = na.omit(unique(geneMap$gencodeID))
genes = list(plac500kb = geneMap[subjectHits(ooP5),], plac1Mb = geneMap[subjectHits(ooP1),], nonplac500kb = geneMap[subjectHits(ooNP5),], 
             nonplac1Mb = geneMap[subjectHits(ooNP1),])
genes = lapply(genes, function(x) unique(x$gencodeID))

DMRgr = lapply(DMR, function(x) makeGRangesFromDataFrame(x[which(x$sig=="FWER < 0.05"),], keep.extra.columns = T))
DMRgr = c(DMRgr, list("CellType.N+" = DMRgr$CellType[which(DMRgr$CellType$Dir=="pos" & DMRgr$CellType$sig=="FWER < 0.05"),],
                      "CellType.N-" = DMRgr$CellType[which(DMRgr$CellType$Dir=="neg" & DMRgr$CellType$sig=="FWER < 0.05"),]),
          lapply(dmrs, function(x) makeGRangesFromDataFrame(x[which(x$sig=="FWER < 0.05"),], keep.extra.columns = T)))
sig = lapply(DMRgr, function(x) unique(as.character(x$nearestID)))
notsig = lapply(sig, function(y) geneuniverse[!(geneuniverse %in% y)])

DMRenrich = mapply(function(sig,notsig) lapply(genes, function(x) {
  DE_OVERLAP = c( sum( sig %in% x),sum(!(sig %in% x)))
  NOT_DE_OVERLAP= c(sum(notsig %in% x), sum(!(notsig %in% x)))
  enrich_table = cbind(DE_OVERLAP, NOT_DE_OVERLAP)
  res = fisher.test(enrich_table)
  dat=c(res$p.value, res$estimate)
  names(dat) <- c("P.Value","Odds.Ratio")
  return(dat)
}), sig, notsig,SIMPLIFY =F) 

tog = lapply(sig, function(s) list(placVSnonplac500kb = fisher.test(data.frame(c(sum( s %in% genes$plac500kb),sum(!(s %in% genes$plac500kb))),
                                                                   c(sum( s %in% genes$nonplac500kb),sum(!(s %in% genes$nonplac500kb))))),
                                   placVSnonplac1Mb = fisher.test(data.frame(c(sum( s %in% genes$plac1Mb),sum(!(s %in% genes$plac1Mb))),
                                                                 c(sum( s %in% genes$nonplac1Mb),sum(!(s %in% genes$nonplac1Mb)))))))
tog = lapply(tog, function(x) lapply(x, function(y) data.frame("P.Value"=y$p.value,"Odds.Ratio"=y$estimate)))
DMRenrich = mapply(function(d,t) c(d,t), DMRenrich, tog, SIMPLIFY = F)

DMRenrich = do.call(rbind, Map(cbind, lapply(DMRenrich, function(x) do.call(rbind, Map(cbind, lapply(x, function(y) data.frame(P.Value = y["P.Value"], Odds.Ratio = y["Odds.Ratio"])), 
                                                                                       GeneSet = as.list(names(x))))), Model = as.list(names(DMRenrich))))
DMRenrich$FDR = p.adjust(DMRenrich$P.Value, method = "fdr")
rownames(DMRenrich) = NULL
write.csv(DMRenrich, quote = F, "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/DMR/placenta_SNPregions_enrichment_DMRs_nearestGenes.csv")


## anchor in genes rather than clusters: overlapping genes only

geneuniverse = na.omit(unique(geneMap$gencodeID))
genes = list(plac500kb = geneMap[subjectHits(ooP5),], plac1Mb = geneMap[subjectHits(ooP1),], nonplac500kb = geneMap[subjectHits(ooNP5),], 
             nonplac1Mb = geneMap[subjectHits(ooNP1),])
genes = lapply(genes, function(x) unique(x$gencodeID))

DMRgr = lapply(DMR, function(x) makeGRangesFromDataFrame(x[which(x$sig=="FWER < 0.05" & x$distToGene==0),], keep.extra.columns = T))
DMRgr = c(DMRgr, list("CellType.N+" = DMRgr$CellType[which(DMRgr$CellType$Dir=="pos" & DMRgr$CellType$sig=="FWER < 0.05" & DMRgr$CellType$distToGene==0),],
                      "CellType.N-" = DMRgr$CellType[which(DMRgr$CellType$Dir=="neg" & DMRgr$CellType$sig=="FWER < 0.05" & DMRgr$CellType$distToGene==0),]),
          lapply(dmrs, function(x) makeGRangesFromDataFrame(x[which(x$sig=="FWER < 0.05" & x$distToGene==0),], keep.extra.columns = T)))
sig = lapply(DMRgr, function(x) unique(as.character(x$nearestID)))
notsig = lapply(sig, function(y) geneuniverse[!(geneuniverse %in% y)])

DMRenrich = mapply(function(sig,notsig) lapply(genes, function(x) {
  DE_OVERLAP = c( sum( sig %in% x),sum(!(sig %in% x)))
  NOT_DE_OVERLAP= c(sum(notsig %in% x), sum(!(notsig %in% x)))
  enrich_table = cbind(DE_OVERLAP, NOT_DE_OVERLAP)
  res = fisher.test(enrich_table)
  dat=c(res$p.value, res$estimate)
  names(dat) <- c("P.Value","Odds.Ratio")
  return(dat)
}), sig, notsig,SIMPLIFY =F) 

tog = lapply(sig, function(s) list(placVSnonplac500kb = fisher.test(data.frame(c(sum( s %in% genes$plac500kb),sum(!(s %in% genes$plac500kb))),
                                                                               c(sum( s %in% genes$nonplac500kb),sum(!(s %in% genes$nonplac500kb))))),
                                   placVSnonplac1Mb = fisher.test(data.frame(c(sum( s %in% genes$plac1Mb),sum(!(s %in% genes$plac1Mb))),
                                                                             c(sum( s %in% genes$nonplac1Mb),sum(!(s %in% genes$nonplac1Mb)))))))
tog = lapply(tog, function(x) lapply(x, function(y) data.frame("P.Value"=y$p.value,"Odds.Ratio"=y$estimate)))
DMRenrich = mapply(function(d,t) c(d,t), DMRenrich, tog, SIMPLIFY = F)

DMRenrich = do.call(rbind, Map(cbind, lapply(DMRenrich, function(x) do.call(rbind, Map(cbind, lapply(x, function(y) data.frame(P.Value = y["P.Value"], Odds.Ratio = y["Odds.Ratio"])), 
                                                                                       GeneSet = as.list(names(x))))), Model = as.list(names(DMRenrich))))
DMRenrich$FDR = p.adjust(DMRenrich$P.Value, method = "fdr")
rownames(DMRenrich) = NULL
write.csv(DMRenrich, quote = F, "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/DMR/placenta_SNPregions_enrichment_DMRs_overlappingGenes.csv")


## are there significantly more genes in the placenta loci?

cols = c("CellType","CellType.N+","CellType.N-","Age","Interaction","Gr1","Gr2","Gr3","Gr4","Gr5","Gr6")
f = list(list(plac1000 = glm(plac1000 ~ CellType + Gene, family = "binomial", data=df.clusters),
              nonplac1000 = glm(nonplac1000 ~ CellType + Gene, family = "binomial", data=df.clusters),
              plac500 = glm(plac500 ~ CellType + Gene, family = "binomial", data=df.clusters),
              nonplac500 = glm(nonplac500 ~ CellType + Gene, family = "binomial", data=df.clusters)),
         list(plac1000 = glm(plac1000 ~ CellType.N + Gene, family = "binomial", data=df.clusters),
              nonplac1000 = glm(nonplac1000 ~ CellType.N + Gene, family = "binomial", data=df.clusters),
              plac500 = glm(plac500 ~ CellType.N + Gene, family = "binomial", data=df.clusters),
              nonplac500 = glm(nonplac500 ~ CellType.N + Gene, family = "binomial", data=df.clusters)),
         list(plac1000 = glm(plac1000 ~ CellType.G + Gene, family = "binomial", data=df.clusters),
              nonplac1000 = glm(nonplac1000 ~ CellType.G + Gene, family = "binomial", data=df.clusters),
              plac500 = glm(plac500 ~ CellType.G + Gene, family = "binomial", data=df.clusters),
              nonplac500 = glm(nonplac500 ~ CellType.G + Gene, family = "binomial", data=df.clusters)),
         list(plac1000 = glm(plac1000 ~ Age + Gene, family = "binomial", data=df.clusters),
              nonplac1000 = glm(nonplac1000 ~ Age + Gene, family = "binomial", data=df.clusters),
              plac500 = glm(plac500 ~ Age + Gene, family = "binomial", data=df.clusters),
              nonplac500 = glm(nonplac500 ~ Age + Gene, family = "binomial", data=df.clusters)),
         list(plac1000 = glm(plac1000 ~ Interaction + Gene, family = "binomial", data=df.clusters),
              nonplac1000 = glm(nonplac1000 ~ Interaction + Gene, family = "binomial", data=df.clusters),
              plac500 = glm(plac500 ~ Interaction + Gene, family = "binomial", data=df.clusters),
              nonplac500 = glm(nonplac500 ~ Interaction + Gene, family = "binomial", data=df.clusters)),
         list(plac1000 = glm(plac1000 ~ Gr1 + Gene, family = "binomial", data=df.clusters),
              nonplac1000 = glm(nonplac1000 ~ Gr1 + Gene, family = "binomial", data=df.clusters),
              plac500 = glm(plac500 ~ Gr1 + Gene, family = "binomial", data=df.clusters),
              nonplac500 = glm(nonplac500 ~ Gr1 + Gene, family = "binomial", data=df.clusters)),
         list(plac1000 = glm(plac1000 ~ Gr2 + Gene, family = "binomial", data=df.clusters),
              nonplac1000 = glm(nonplac1000 ~ Gr2 + Gene, family = "binomial", data=df.clusters),
              plac500 = glm(plac500 ~ Gr2 + Gene, family = "binomial", data=df.clusters),
              nonplac500 = glm(nonplac500 ~ Gr2 + Gene, family = "binomial", data=df.clusters)),
         list(plac1000 = glm(plac1000 ~ Gr3 + Gene, family = "binomial", data=df.clusters),
              nonplac1000 = glm(nonplac1000 ~ Gr3 + Gene, family = "binomial", data=df.clusters),
              plac500 = glm(plac500 ~ Gr3 + Gene, family = "binomial", data=df.clusters),
              nonplac500 = glm(nonplac500 ~ Gr3 + Gene, family = "binomial", data=df.clusters)),
         list(plac1000 = glm(plac1000 ~ Gr4 + Gene, family = "binomial", data=df.clusters),
              nonplac1000 = glm(nonplac1000 ~ Gr4 + Gene, family = "binomial", data=df.clusters),
              plac500 = glm(plac500 ~ Gr4 + Gene, family = "binomial", data=df.clusters),
              nonplac500 = glm(nonplac500 ~ Gr4 + Gene, family = "binomial", data=df.clusters)),
         list(plac1000 = glm(plac1000 ~ Gr5 + Gene, family = "binomial", data=df.clusters),
              nonplac1000 = glm(nonplac1000 ~ Gr5 + Gene, family = "binomial", data=df.clusters),
              plac500 = glm(plac500 ~ Gr5 + Gene, family = "binomial", data=df.clusters),
              nonplac500 = glm(nonplac500 ~ Gr5 + Gene, family = "binomial", data=df.clusters)),
         list(plac1000 = glm(plac1000 ~ Gr6 + Gene, family = "binomial", data=df.clusters),
              nonplac1000 = glm(nonplac1000 ~ Gr6 + Gene, family = "binomial", data=df.clusters),
              plac500 = glm(plac500 ~ Gr6 + Gene, family = "binomial", data=df.clusters),
              nonplac500 = glm(nonplac500 ~ Gr6 + Gene, family = "binomial", data=df.clusters)))
df = do.call(rbind, Map(cbind, lapply(f, function(x) data.frame(AdjustedOdds = unlist(lapply(x, function(y) exp(y$coef[2]))), SNPgroup = names(x))), SigGroup = as.list(cols)))
rownames(df) = NULL
df$SigGroup = gsub("Gr1", "G-:N+", df$SigGroup)
df$SigGroup = gsub("Gr2", "G0:N+", df$SigGroup)
df$SigGroup = gsub("Gr3", "G0:N-", df$SigGroup)
df$SigGroup = gsub("Gr4", "G+:N0", df$SigGroup)
df$SigGroup = gsub("Gr5", "G+:N-", df$SigGroup)
df$SigGroup = gsub("Gr6", "G-:N0", df$SigGroup)
write.csv(df,file="/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/DMR/placenta_SNPregions_enrichment_DMRs_geneAdjusted.csv",quote=F)


f.proteincoding = list(list(plac1000 = glm(plac1000 ~ CellType + ProteinCoding, family = "binomial", data=df.clusters),
                            nonplac1000 = glm(nonplac1000 ~ CellType + ProteinCoding, family = "binomial", data=df.clusters),
                            plac500 = glm(plac500 ~ CellType + ProteinCoding, family = "binomial", data=df.clusters),
                            nonplac500 = glm(nonplac500 ~ CellType + ProteinCoding, family = "binomial", data=df.clusters)),
                       list(plac1000 = glm(plac1000 ~ CellType.N + ProteinCoding, family = "binomial", data=df.clusters),
                            nonplac1000 = glm(nonplac1000 ~ CellType.N + ProteinCoding, family = "binomial", data=df.clusters),
                            plac500 = glm(plac500 ~ CellType.N + ProteinCoding, family = "binomial", data=df.clusters),
                            nonplac500 = glm(nonplac500 ~ CellType.N + ProteinCoding, family = "binomial", data=df.clusters)),
                       list(plac1000 = glm(plac1000 ~ CellType.G + ProteinCoding, family = "binomial", data=df.clusters),
                            nonplac1000 = glm(nonplac1000 ~ CellType.G + ProteinCoding, family = "binomial", data=df.clusters),
                            plac500 = glm(plac500 ~ CellType.G + ProteinCoding, family = "binomial", data=df.clusters),
                            nonplac500 = glm(nonplac500 ~ CellType.G + ProteinCoding, family = "binomial", data=df.clusters)),
                       list(plac1000 = glm(plac1000 ~ Age + ProteinCoding, family = "binomial", data=df.clusters),
                            nonplac1000 = glm(nonplac1000 ~ Age + ProteinCoding, family = "binomial", data=df.clusters),
                            plac500 = glm(plac500 ~ Age + ProteinCoding, family = "binomial", data=df.clusters),
                            nonplac500 = glm(nonplac500 ~ Age + ProteinCoding, family = "binomial", data=df.clusters)),
                       list(plac1000 = glm(plac1000 ~ Interaction + ProteinCoding, family = "binomial", data=df.clusters),
                            nonplac1000 = glm(nonplac1000 ~ Interaction + ProteinCoding, family = "binomial", data=df.clusters),
                            plac500 = glm(plac500 ~ Interaction + ProteinCoding, family = "binomial", data=df.clusters),
                            nonplac500 = glm(nonplac500 ~ Interaction + ProteinCoding, family = "binomial", data=df.clusters)),
                       list(plac1000 = glm(plac1000 ~ Gr1 + ProteinCoding, family = "binomial", data=df.clusters),
                            nonplac1000 = glm(nonplac1000 ~ Gr1 + ProteinCoding, family = "binomial", data=df.clusters),
                            plac500 = glm(plac500 ~ Gr1 + ProteinCoding, family = "binomial", data=df.clusters),
                            nonplac500 = glm(nonplac500 ~ Gr1 + ProteinCoding, family = "binomial", data=df.clusters)),
                       list(plac1000 = glm(plac1000 ~ Gr2 + ProteinCoding, family = "binomial", data=df.clusters),
                            nonplac1000 = glm(nonplac1000 ~ Gr2 + ProteinCoding, family = "binomial", data=df.clusters),
                            plac500 = glm(plac500 ~ Gr2 + ProteinCoding, family = "binomial", data=df.clusters),
                            nonplac500 = glm(nonplac500 ~ Gr2 + ProteinCoding, family = "binomial", data=df.clusters)),
                       list(plac1000 = glm(plac1000 ~ Gr3 + ProteinCoding, family = "binomial", data=df.clusters),
                            nonplac1000 = glm(nonplac1000 ~ Gr3 + ProteinCoding, family = "binomial", data=df.clusters),
                            plac500 = glm(plac500 ~ Gr3 + ProteinCoding, family = "binomial", data=df.clusters),
                            nonplac500 = glm(nonplac500 ~ Gr3 + ProteinCoding, family = "binomial", data=df.clusters)),
                       list(plac1000 = glm(plac1000 ~ Gr4 + ProteinCoding, family = "binomial", data=df.clusters),
                            nonplac1000 = glm(nonplac1000 ~ Gr4 + ProteinCoding, family = "binomial", data=df.clusters),
                            plac500 = glm(plac500 ~ Gr4 + ProteinCoding, family = "binomial", data=df.clusters),
                            nonplac500 = glm(nonplac500 ~ Gr4 + ProteinCoding, family = "binomial", data=df.clusters)),
                       list(plac1000 = glm(plac1000 ~ Gr5 + ProteinCoding, family = "binomial", data=df.clusters),
                            nonplac1000 = glm(nonplac1000 ~ Gr5 + ProteinCoding, family = "binomial", data=df.clusters),
                            plac500 = glm(plac500 ~ Gr5 + ProteinCoding, family = "binomial", data=df.clusters),
                            nonplac500 = glm(nonplac500 ~ Gr5 + ProteinCoding, family = "binomial", data=df.clusters)),
                       list(plac1000 = glm(plac1000 ~ Gr6 + ProteinCoding, family = "binomial", data=df.clusters),
                            nonplac1000 = glm(nonplac1000 ~ Gr6 + ProteinCoding, family = "binomial", data=df.clusters),
                            plac500 = glm(plac500 ~ Gr6 + ProteinCoding, family = "binomial", data=df.clusters),
                            nonplac500 = glm(nonplac500 ~ Gr6 + ProteinCoding, family = "binomial", data=df.clusters)))
df = do.call(rbind, Map(cbind, lapply(f, function(x) data.frame(AdjustedOdds = unlist(lapply(x, function(y) exp(y$coef[2]))), SNPgroup = names(x))), SigGroup = as.list(cols)))
rownames(df) = NULL
df$SigGroup = gsub("Gr1", "G-:N+", df$SigGroup)
df$SigGroup = gsub("Gr2", "G0:N+", df$SigGroup)
df$SigGroup = gsub("Gr3", "G0:N-", df$SigGroup)
df$SigGroup = gsub("Gr4", "G+:N0", df$SigGroup)
df$SigGroup = gsub("Gr5", "G+:N-", df$SigGroup)
df$SigGroup = gsub("Gr6", "G-:N0", df$SigGroup)
write.csv(df,file="/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/DMR/placenta_SNPregions_enrichment_DMRs_proteincodinggeneAdjusted.csv",quote=F)


## Associate with LD blocks

pgc = read.delim("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/tableS3_Nature2014_pgc2_loci.txt",as.is=TRUE)
pgc = GRanges(pgc$Position..hg19.)

ooP = findOverlaps(pgc, plac)
plac[c(7,19,26,33,35,40,49,50,55,57)]

