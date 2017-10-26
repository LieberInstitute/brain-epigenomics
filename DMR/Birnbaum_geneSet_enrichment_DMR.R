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
  names(dat) <- c("P.Value","Odds Ratio")
  return(dat)
}), sig, notsig,splitSets,SIMPLIFY =F) 

mapply(function(x,y,z) write.csv(x,file=paste0("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/DMR/",y,"/Birnbaum_geneSet_enrichment_DMR_",z,".csv"),
                                                 quote=F),DMRenrich,c("CellType","Age","CT_Age_Interaction"),names(DMRenrich))