load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/non-CpG/CH_object.rda")
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/non-CpG/CHneurons_object.rda")

aej_sets = openxlsx::read.xlsx('/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/Birnbaum_2013_AJP_Supplementary_table.xlsx')


## Enrichment in dmCH by cell type, age and interaction

CH = as.data.frame(CHdt)
CHneurons = as.data.frame(CHneuronsdt)
geneuniverse = as.character(na.omit(unique(geneMap$EntrezID)))
sig = list("Cell Type" = as.character(unique(CH[CH$distToGene==0 & CH$CT.sig=="FDR < 0.05","EntrezID"])),
           "Age" = as.character(unique(CH[CH$distToGene==0 & CH$Age.sig=="FDR < 0.05","EntrezID"])),
           "Interaction" = as.character(unique(CH[CH$distToGene==0 & CH$Int.sig=="FDR < 0.05","EntrezID"])), 
           "Age in Neurons" = as.character(unique(CHneurons[CHneurons$distToGene==0 & CHneurons$sig=="FDR < 0.05","EntrezID"])))
notsig = lapply(sig, function(y) geneuniverse[!(geneuniverse %in% y)])
aej_sets_expressed = aej_sets[which(aej_sets$EntrezGene.ID %in% geneuniverse), ]
splitSets = split(aej_sets_expressed[aej_sets_expressed$Gene.Set!="SCZ PGC GWAS",], aej_sets_expressed$Gene.Set[aej_sets_expressed$Gene.Set!="SCZ PGC GWAS"])
splitSets = lapply(splitSets, function(x) as.character(unique(x$EntrezGene.ID)))

CHenrich = mapply(function(sig,notsig) lapply(splitSets, function(x) {
  DE_OVERLAP = c( sum( sig %in% x),sum(!(sig %in% x)))
  NOT_DE_OVERLAP= c(sum(notsig %in% x), sum(!(notsig %in% x)))
  enrich_table = cbind(DE_OVERLAP, NOT_DE_OVERLAP)
  res = fisher.test(enrich_table)
  dat=c(res$p.value, res$estimate)
  names(dat) <- c("P.Value","Odds Ratio")
  return(dat)
}), sig, notsig, SIMPLIFY =F) 

CHenrich = do.call(rbind, Map(cbind, Model = as.list(names(CHenrich)), 
                              lapply(CHenrich, function(x) do.call(rbind, Map(cbind, GeneSet = as.list(names(x)), 
                                                                              lapply(x, function(y) data.frame(Odds.Ratio = y[names(y)=="Odds Ratio"], 
                                                                                                               P.value = y[names(y)=="P.Value"])))))))
CHenrich$FDR = p.adjust(CHenrich$P.value, method = "fdr")
rownames(CHenrich) = NULL

write.csv(CHenrich,file="/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/non-CpG/Birnbaum_geneSet_enrichment_CpHs.csv",quote=F)


CHenrich[CHenrich$FDR<=0.05,]
#            Model           GeneSet Odds.Ratio      P.value          FDR
#2       Cell Type      ASD DATABASE  2.6395744 5.063303e-08 3.645578e-07
#4       Cell Type                ID  0.2923644 1.506089e-08 1.355480e-07
#5       Cell Type               NDD  5.5668594 7.095733e-03 2.128720e-02
#9       Cell Type           SCZ SNV  3.3331061 2.266932e-09 4.080477e-08
#11            Age      ASD DATABASE  2.2731621 1.203790e-07 7.222738e-07
#13            Age                ID  0.3397681 9.482268e-07 4.876595e-06
#14            Age               NDD  3.8158873 7.087042e-03 2.128720e-02
#17            Age SCZ Meta-analysis  2.8373908 1.394511e-02 3.861723e-02
#18            Age           SCZ SNV  2.7152725 6.486395e-09 7.783674e-08
#29 Age in Neurons      ASD DATABASE  2.0210391 6.737858e-06 2.695143e-05
#31 Age in Neurons                ID  0.3592397 2.445902e-06 1.100656e-05
#32 Age in Neurons               NDD  4.8281465 3.335997e-03 1.200959e-02
#36 Age in Neurons           SCZ SNV  2.9858497 9.777669e-10 3.519961e-08