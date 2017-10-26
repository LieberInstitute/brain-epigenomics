load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/non-CpG/CH_object.rda")
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/non-CpG/CHneurons_object.rda")


aej_sets = openxlsx::read.xlsx('/dcl01/lieber/ajaffe/lab/Nicotine/DLPFC/RNAseq/csvs/sets/Supplementary Tables for paper.Birnbaum November 2013.AJP.xlsx')

## Enrichment in dmCH by cell type, age and interaction

geneuniverse = na.omit(unique(CH$EntrezID))
sigCH = list("CellType" = na.omit(unique(CH[which(CH$padj.CellType<=0.01),"EntrezID"])),
             "Age" = na.omit(unique(CH[which(CH$padj.Age<=0.01),"EntrezID"])), 
             "Interaction" = na.omit(unique(CH[which(CH$padj.Interaction<=0.01),"EntrezID"])))
notsigCH = lapply(sigCH, function(x) geneuniverse[!(geneuniverse %in% x)])
aej_sets_expressed = aej_sets[which(aej_sets$EntrezGene.ID %in% geneuniverse), ] # drop genes that are not present in the test set
splitSets = split(aej_sets_expressed, aej_sets_expressed$Gene.Set)

CHenrich = mapply(function(sig,notsig) sapply(splitSets, function(x) {
  DE_OVERLAP = c( sum( sig %in% x$EntrezGene.ID),sum(!(sig %in% x$EntrezGene.ID)))
  NOT_DE_OVERLAP= c(sum(notsig %in% x$EntrezGene.ID), sum(!(notsig %in% x$EntrezGene.ID)))
  enrich_table = cbind(DE_OVERLAP, NOT_DE_OVERLAP)
  res = fisher.test(enrich_table)
  dat=c(res$p.value, res$estimate)
  names(dat) <- c("P.Value","Odds Ratio")
  return(dat)
}), sigCH, notsigCH,SIMPLIFY =F) 

mapply(function(x,y) write.csv(x,file=paste0("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/non-CpG/Birnbaum_geneSet_enrichment_nonCG_by",y,".csv"),
                               quote=F),CHenrich,names(CHenrich))



## Enrichment in dmCH in neurons by age only

geneuniverse = na.omit(unique(CHneurons$EntrezID))
sigCHneurons = na.omit(unique(CHneurons[which(CHneurons$padj<=0.01),"EntrezID"]))
notsigCHneurons = geneuniverse[!(geneuniverse %in% sigCHneurons)]
aej_sets_expressed = aej_sets[which(aej_sets$EntrezGene.ID %in% geneuniverse), ] # drop genes that are not present in the test set
splitSets = split(aej_sets_expressed, aej_sets_expressed$Gene.Set)

CHNeuronsenrich = sapply(splitSets, function(x) {
  DE_OVERLAP = c( sum( sigCHneurons %in% x$EntrezGene.ID),sum(!(sigCHneurons %in% x$EntrezGene.ID)))
  NOT_DE_OVERLAP= c(sum(notsigCHneurons %in% x$EntrezGene.ID), sum(!(notsigCHneurons %in% x$EntrezGene.ID)))
  enrich_table = cbind(DE_OVERLAP, NOT_DE_OVERLAP)
  res = fisher.test(enrich_table)
  dat=c(res$p.value, res$estimate)
  names(dat) <- c("P.Value","Odds Ratio")
  return(dat)
})

write.csv(CHNeuronsenrich,file="/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/non-CpG/Birnbaum_geneSet_enrichment_nonCG_byAge_neuronsOnly.csv",quote=F)