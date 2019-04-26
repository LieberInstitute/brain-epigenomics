library(GenomicRanges)

load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/non-CpG/CH_object.rda")
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/non-CpG/CHneurons_object.rda")

aej_sets = openxlsx::read.xlsx('/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/Birnbaum_2013_AJP_Supplementary_table.xlsx')


## Enrichment in dmCH by cell type, age and interaction

CH = as.data.frame(CHdt)
CHneurons = as.data.frame(CHneuronsdt)
geneuniverse = as.character(na.omit(unique(geneMap$EntrezID)))
sig = list("Cell Type (Neurons > Glia)" = as.character(unique(CH[CH$distToGene==0 & CH$CT.sig=="FDR < 0.05" & CH$CT.dir=="pos","EntrezID"])),
           "Cell Type (Glia > Neurons)" = as.character(unique(CH[CH$distToGene==0 & CH$CT.sig=="FDR < 0.05" & CH$CT.dir=="neg","EntrezID"])),
           "Age (Older > Younger)" = as.character(unique(CH[CH$distToGene==0 & CH$Age.sig=="FDR < 0.05" & CH$Age.dir=="pos","EntrezID"])),
           "Age (Younger > Older)" = as.character(unique(CH[CH$distToGene==0 & CH$Age.sig=="FDR < 0.05" & CH$Age.dir=="neg","EntrezID"])),
           "Age in Neurons (Older > Younger)" = as.character(unique(CHneurons[CHneurons$distToGene==0 & CHneurons$sig=="FDR < 0.05" & CHneurons$Dir=="pos","EntrezID"])),
           "Age in Neurons (Younger > Older)" = as.character(unique(CHneurons[CHneurons$distToGene==0 & CHneurons$sig=="FDR < 0.05" & CHneurons$Dir=="neg","EntrezID"])))
elementNROWS(sig)
notsig = lapply(sig, function(y) geneuniverse[!(geneuniverse %in% y)])
aej_sets_expressed = aej_sets[which(aej_sets$EntrezGene.ID %in% geneuniverse), ]
splitSets = split(aej_sets_expressed[!aej_sets_expressed$Gene.Set %in% c("SCZ PGC GWAS","BPAD GWAS","SCZ Meta-analysis"),], 
                  aej_sets_expressed$Gene.Set[!aej_sets_expressed$Gene.Set %in% c("SCZ PGC GWAS","BPAD GWAS","SCZ Meta-analysis")])
splitSets = lapply(splitSets, function(x) as.character(unique(x$EntrezGene.ID)))

CHenrich = mapply(function(sig,notsig) lapply(splitSets, function(x) {
  DE_OVERLAP = c( sum( sig %in% x),sum(!(sig %in% x)))
  NOT_DE_OVERLAP= c(sum(notsig %in% x), sum(!(notsig %in% x)))
  enrich_table = cbind(DE_OVERLAP, NOT_DE_OVERLAP)
  res = fisher.test(enrich_table)
  dat=c(res$p.value, res$estimate, res$conf.int[1], res$conf.int[2])
  names(dat) <- c("pval","OR", "lower","upper")
  return(dat)
}), sig, notsig, SIMPLIFY =F) 

enrich_table = mapply(function(sig,notsig) lapply(splitSets, function(x) {
  DE_OVERLAP = c( sum( sig %in% x),sum(!(sig %in% x)))
  NOT_DE_OVERLAP= c(sum(notsig %in% x), sum(!(notsig %in% x)))
  cbind(DE_OVERLAP, NOT_DE_OVERLAP)
  }), sig, notsig, SIMPLIFY =F) 

CHenrich = do.call(rbind, Map(cbind, Model = as.list(names(CHenrich)), 
                              lapply(CHenrich, function(x) do.call(rbind, Map(cbind, GeneSet = as.list(names(x)), 
                                                                               lapply(x, function(y) data.frame(P.Value = y["pval"], 
                                                                                                                Odds.Ratio = y["OR"],
                                                                                                                Lower = y["lower"],
                                                                                                                Upper = y["upper"])))))))
CHenrich$FDR = p.adjust(CHenrich$P.Value, method = "fdr")
rownames(CHenrich) = NULL

ta = do.call(rbind, lapply(enrich_table, function(x) do.call(rbind, lapply(x, function(y) data.frame(YesGeneSet.YesSigCH = y[1,1], 
                                                                                                     NoGeneSet.YesSigCH = y[2,1],
                                                                                                     YesGeneSet.NoSigCH = y[1,2], 
                                                                                                     NoGeneSet.NoSigCH = y[2,2])))))
rownames(ta) = NULL
df = cbind(CHenrich, ta)
df$C_Context = "CpH"

write.csv(df,file="/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/non-CpG/Birnbaum_geneSet_enrichment_CpHs.csv",quote=F)


df[df$FDR<=0.05,colnames(df) %in% c("Model","GeneSet","Odds.Ratio", "FDR")]
#                              Model           GeneSet Odds.Ratio          FDR
#2        Cell Type (Neurons > Glia)      ASD DATABASE  2.3267216 1.996105e-06
#3        Cell Type (Neurons > Glia)                ID  0.3146172 5.518187e-07
#4        Cell Type (Neurons > Glia)               NDD  4.0352192 1.758225e-02
#7        Cell Type (Neurons > Glia)           SCZ SNV  3.0841059 2.878336e-08
#8        Cell Type (Glia > Neurons)           ASD CNV  1.6002859 8.825119e-03
#9        Cell Type (Glia > Neurons)      ASD DATABASE  3.2851986 5.294648e-17
#10       Cell Type (Glia > Neurons)                ID  0.4745279 7.097966e-03
#11       Cell Type (Glia > Neurons)               NDD  9.1413212 4.955402e-07
#12       Cell Type (Glia > Neurons) Neurodegenerative  2.5958558 3.722698e-03
#14       Cell Type (Glia > Neurons)           SCZ SNV  2.8433179 5.888615e-12
#16            Age (Older > Younger)      ASD DATABASE  2.1159828 2.243289e-06
#17            Age (Older > Younger)                ID  0.3425688 3.270887e-06
#18            Age (Older > Younger)               NDD  4.0435044 8.368437e-03
#21            Age (Older > Younger)           SCZ SNV  2.5170287 1.973036e-07
#23            Age (Younger > Older)      ASD DATABASE  3.2804188 1.432070e-14
#25            Age (Younger > Older)               NDD  5.2496969 3.761386e-05
#28            Age (Younger > Older)           SCZ SNV  2.3866170 5.518187e-07
#30 Age in Neurons (Older > Younger)      ASD DATABASE  2.1078740 3.270887e-06
#31 Age in Neurons (Older > Younger)                ID  0.3665019 1.271120e-05
#32 Age in Neurons (Older > Younger)               NDD  3.7321228 1.332551e-02
#35 Age in Neurons (Older > Younger)           SCZ SNV  2.7509152 4.038519e-08
#37 Age in Neurons (Younger > Older)      ASD DATABASE  2.8667072 1.170154e-12
#39 Age in Neurons (Younger > Older)               NDD  5.8549822 6.626613e-06
#42 Age in Neurons (Younger > Older)           SCZ SNV  2.1788649 2.243289e-06

