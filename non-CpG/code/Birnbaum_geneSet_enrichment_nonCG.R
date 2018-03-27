library(jaffelab)


load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/non-CpG/CH_object.rda")
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/non-CpG/CHneurons_object.rda")


aej_sets = openxlsx::read.xlsx('/dcl01/lieber/ajaffe/lab/Nicotine/DLPFC/RNAseq/csvs/sets/Supplementary Tables for paper.Birnbaum November 2013.AJP.xlsx')
xx=load("/users/ajaffe/Lieber/Projects/RNAseq/n36/finalCode/gwasResults_lifted.rda")

## create pgc2 loci granges object

pgc = read.delim("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/tableS3_Nature2014_pgc2_loci.txt",as.is=TRUE)
pgc$chr = ss(pgc$Position..hg19.,":")
tmp = ss(pgc$Position..hg19.,":",2)
pgc$start = as.numeric(ss(tmp, "-"))
pgc$end = as.numeric(ss(tmp, "-",2))
pgcGR = GRanges(pgc$chr, IRanges(pgc$start,pgc$end))

geneMapgr = makeGRangesFromDataFrame(geneMap, keep.extra.columns = T)
oo = lapply(c(list(pgcGR), as.list(split(gwasLift, gwasLift$Dx))), function(x) findOverlaps(x, geneMapgr))
PGCgenes = lapply(oo, function(x) geneMap[subjectHits(x),])
names(PGCgenes)[1] = "SCZ PGC2"

## Enrichment in dmCH by cell type, age and interaction

CH = as.data.frame(CHdt)
geneuniverse = na.omit(unique(CH$EntrezID))
sigCH = list("CellType" = na.omit(unique(CH[which(CH$padj.CellType<=0.05),"EntrezID"])),
             "Age" = na.omit(unique(CH[which(CH$padj.Age<=0.01),"EntrezID"])), 
             "Interaction" = na.omit(unique(CH[which(CH$padj.Interaction<=0.05),"EntrezID"])))
notsigCH = lapply(sigCH, function(x) geneuniverse[!(geneuniverse %in% x)])
aej_sets_expressed = aej_sets[which(aej_sets$EntrezGene.ID %in% geneuniverse), ] # drop genes that are not present in the test set
splitSets = c(split(aej_sets_expressed$EntrezGene.ID, aej_sets_expressed$Gene.Set), lapply(PGCgenes, function(x) unique(x$EntrezID)[which(unique(x$EntrezID) %in% geneuniverse)]))
splitSets = lapply(splitSets, unique)

CHenrich = mapply(function(sig,notsig) lapply(splitSets, function(x) {
  DE_OVERLAP = c( sum( sig %in% x),sum(!(sig %in% x)))
  NOT_DE_OVERLAP= c(sum(notsig %in% x), sum(!(notsig %in% x)))
  enrich_table = cbind(DE_OVERLAP, NOT_DE_OVERLAP)
  res = fisher.test(enrich_table)
  dat=c(res$p.value, res$estimate)
  names(dat) <- c("P.Value","Odds Ratio")
  return(dat)
}), sigCH, notsigCH,SIMPLIFY =F) 


## Enrichment in dmCH in neurons by age only

CHneurons = as.data.frame(CHneuronsdt)
geneuniverse = na.omit(unique(CHneurons$EntrezID))
sigCHneurons = list("Age in Neurons: All" = na.omit(unique(CHneurons[which(CHneurons$padj<=0.05),"EntrezID"])), 
                    "Age in Neurons: Increasing" = na.omit(unique(CHneurons[which(CHneurons$padj<=0.05 & CHneurons$Dir=="pos"),"EntrezID"])), 
                    "Age in Neurons: Decreasing" = na.omit(unique(CHneurons[which(CHneurons$padj<=0.05 & CHneurons$Dir=="neg"),"EntrezID"])))
notsigCHneurons = lapply(sigCHneurons, function(x) geneuniverse[!(geneuniverse %in% x)])

CHNeuronsenrich = mapply(function(sig,notsig) lapply(splitSets, function(x) {
  DE_OVERLAP = c( sum( sig %in% x), sum(!(sig %in% x)))
  NOT_DE_OVERLAP= c(sum(notsig %in% x), sum(!(notsig %in% x)))
  enrich_table = cbind(DE_OVERLAP, NOT_DE_OVERLAP)
  res = fisher.test(enrich_table)
  dat=c(res$p.value, res$estimate)
  names(dat) <- c("P.Value","Odds Ratio")
  return(dat)
}), sigCHneurons, notsigCHneurons, SIMPLIFY = F)


CHenrich = do.call(rbind, Map(cbind, lapply(CHenrich, function(x) do.call(rbind, Map(cbind, lapply(x, function(y) data.frame(Odds.Ratio = y[names(y)=="Odds Ratio"], P.value = y[names(y)=="P.Value"])), 
                                                                                     GeneSet = as.list(names(x))))), Model = as.list(names(CHenrich))))
CHNeuronsenrich = do.call(rbind, Map(cbind, lapply(CHNeuronsenrich, function(x) do.call(rbind, Map(cbind, lapply(x, function(y) data.frame(Odds.Ratio = y[names(y)=="Odds Ratio"], P.value = y[names(y)=="P.Value"])), 
                                                                                     GeneSet = as.list(names(x))))), Model = as.list(names(CHNeuronsenrich))))
CHenrich = rbind(CHenrich, CHNeuronsenrich)
CHenrich$FDR = p.adjust(CHenrich$P.value, method = "fdr")
rownames(CHenrich) = NULL

write.csv(CHenrich,file="/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/non-CpG/Birnbaum_PGC2_geneSet_enrichment_CpHs.csv",quote=F)


CHenrich[CHenrich$FDR<=0.05,]
#   Odds.Ratio      P.value      GeneSet                      Model          FDR
#4   0.1571897 4.544958e-09           ID                   CellType 1.908882e-07
#18  0.2357243 7.116369e-10           ID                        Age 5.977750e-08
#25  0.5907584 5.399937e-05     SCZ PGC2                        Age 7.559911e-04
#46  0.2543436 2.912996e-06           ID        Age in Neurons: All 6.117292e-05
#60  0.2674113 3.811205e-06           ID Age in Neurons: Increasing 6.402824e-05
#72  2.1405198 2.224576e-08 ASD DATABASE Age in Neurons: Decreasing 6.228814e-07
#74  0.3418399 9.126282e-05           ID Age in Neurons: Decreasing 1.095154e-03
#75  3.0906668 4.121048e-03          NDD Age in Neurons: Decreasing 4.327100e-02


## Weird about the ID genes being depleted in both increasing and decreasing CpHs

CHNeuronsenrich_table = mapply(function(sig,notsig) lapply(splitSets, function(x) {
  DE_OVERLAP = c( sum( sig %in% x), sum(!(sig %in% x)))
  NOT_DE_OVERLAP= c(sum(notsig %in% x), sum(!(notsig %in% x)))
  enrich_table = cbind(DE_OVERLAP, NOT_DE_OVERLAP)
  return(enrich_table)
}), sigCHneurons, notsigCHneurons, SIMPLIFY = F)

CHNeuronsenrich_table$`Age in Neurons: Decreasing`$ID
CHNeuronsenrich_table$`Age in Neurons: Increasing`$ID

increasingIDgenes = splitSets$ID




