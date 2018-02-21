load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/DMR/DMR_objects.rda")

aej_sets = openxlsx::read.xlsx('/dcl01/lieber/ajaffe/lab/Nicotine/DLPFC/RNAseq/csvs/sets/Supplementary Tables for paper.Birnbaum November 2013.AJP.xlsx')

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
  names(dat) <- c("P.Value","Odds.Ratio")
  return(dat)
}), sig, notsig,splitSets,SIMPLIFY =F) 


DMRenrich = lapply(DMRenrich,t)

DMRenrich = do.call(rbind, Map(cbind, DMRenrich, model = as.list(names(DMRenrich))))
DMRenrich = data.frame(DMRenrich, GeneSet = rownames(DMRenrich))
DMRenrich$P.Value = as.numeric(as.character(DMRenrich$P.Value))
DMRenrich$Odds.Ratio = as.numeric(as.character(DMRenrich$Odds.Ratio))
DMRenrich = rbind(DMRenrich, data.frame(P.Value = c(3.612019e-08,2.071183e-02,4.356683e-04), 
                                        Odds.Ratio = c(1.776787,3.300659,1.898462), 
                                        model = c("CellType","Age","Interaction"), GeneSet = rep.int("PGC2",3)))
DMRenrich$padj = p.adjust(DMRenrich$P.Value, method = "fdr")


write.csv(DMRenrich,file="/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/DMR/Birnbaum_PGC2_geneSet_enrichment_DMRs.csv",quote=F)

DMRenrich[which(DMRenrich$padj<=0.05),]
unique(DMR$Interaction[which(DMR$Interaction$nearestID %in% PGCgenes$gencodeID & DMR$Interaction$distToGene==0),"nearestSymbol"])







