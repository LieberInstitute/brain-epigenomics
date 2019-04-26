library(GenomicRanges)

load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/DMR/DMR_objects.rda")

brain_categories = readRDS("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/ldsc/categories.rds")
brain_categories = brain_categories[which(names(brain_categories) %in% c("Celltype_HypoNeuron","Celltype_HypoGlia",
                                                                   "Age_Decreasing", "Age_Increasing",
                                                                   "Gr1_cdDMRs","Gr2_cdDMRs","Gr3_cdDMRs",
                                                                   "Gr4_cdDMRs","Gr5_cdDMRs","Gr6_cdDMRs"))]

brain_categories_df <- data.frame(Category = names(brain_categories),	
                                  Extended = c("Cell Type (Glia > Neuron)", "Cell Type (Neuron > Glia)", 
                                               "Age (Younger > Older)", "Age (Older > Younger)", 
                                               "Group 1 (Decreasing Glial; Increasing Neuronal)", 
                                               "Group 2 (Static Glial; Increasing Neuronal)", 
                                               "Group 3 (Static Glial; Decreasing Neuronal)", 
                                               "Group 4 (Increasing Glial; Static Neuronal)", 
                                               "Group 5 (Increasing Glial; Decreasing Neuronal)", 
                                               "Group 6 (Decreasing Glial; Static Neuronal)"))

DMR = list(DMR$CellType, DMR$CellType, DMR$Age, DMR$Age, DMR$Interaction, DMR$Interaction, DMR$Interaction, 
           DMR$Interaction, DMR$Interaction, DMR$Interaction)

oo = mapply(function(x,y) findOverlaps(x, makeGRangesFromDataFrame(y), type = "equal"), brain_categories, DMR, SIMPLIFY = F)
brCat = mapply(function(x,y) x[subjectHits(y),], DMR, oo, SIMPLIFY = F)
names(brCat) = names(brain_categories)

aej_sets = openxlsx::read.xlsx('/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/Birnbaum_2013_AJP_Supplementary_table.xlsx')


## Enrichment in DMRs by cell type, age and interaction

geneuniverse = na.omit(unique(geneMap$EntrezID))
aej_sets_expressed = aej_sets[which(aej_sets$EntrezGene.ID %in% geneuniverse), ] # drop genes that are not present in the test set
splitSets = split(aej_sets_expressed[!aej_sets_expressed$Gene.Set %in% c("SCZ PGC GWAS","BPAD GWAS","SCZ Meta-analysis"),], 
                  aej_sets_expressed$Gene.Set[!aej_sets_expressed$Gene.Set %in% c("SCZ PGC GWAS","BPAD GWAS","SCZ Meta-analysis")])
sig = lapply(brCat, function(x) as.character(na.omit(unique(x[which(x$distToGene==0),"EntrezID"]))))
notsig = lapply(sig, function(y) as.character(geneuniverse[!(geneuniverse %in% y)]))

DMRenrich = mapply(function(sig,notsig) lapply(splitSets, function(x) {
  DE_OVERLAP = c( sum( sig %in% x$EntrezGene.ID),sum(!(sig %in% x$EntrezGene.ID)))
  NOT_DE_OVERLAP= c(sum(notsig %in% x$EntrezGene.ID), sum(!(notsig %in% x$EntrezGene.ID)))
  enrich_table = cbind(DE_OVERLAP, NOT_DE_OVERLAP)
  res = fisher.test(enrich_table)
  dat=c(res$p.value, res$estimate, res$conf.int[1], res$conf.int[2])
  names(dat) <- c("pval","OR", "lower","upper")
  return(dat)
}), sig, notsig,SIMPLIFY =F) 

enrich_table = mapply(function(sig,notsig) lapply(splitSets, function(x) {
  DE_OVERLAP = c( sum( sig %in% x$EntrezGene.ID),sum(!(sig %in% x$EntrezGene.ID)))
  NOT_DE_OVERLAP= c(sum(notsig %in% x$EntrezGene.ID), sum(!(notsig %in% x$EntrezGene.ID)))
  cbind(DE_OVERLAP, NOT_DE_OVERLAP)
  }), sig, notsig,SIMPLIFY =F) 

DMRenrich = do.call(rbind, Map(cbind, Model = as.list(names(DMRenrich)), 
                               lapply(DMRenrich, function(x) do.call(rbind, Map(cbind, GeneSet = as.list(names(x)), 
                                                                                lapply(x, function(y) data.frame(P.Value = y["pval"], 
                                                                                                                 Odds.Ratio = y["OR"],
                                                                                                                 Lower = y["lower"],
                                                                                                                 Upper = y["upper"])))))))
DMRenrich$FDR = p.adjust(DMRenrich$P.Value, method = "fdr")
rownames(DMRenrich) = NULL
DMRenrich$Extended = brain_categories_df[match(DMRenrich$Model, brain_categories_df$Category),"Extended"]

ta = do.call(rbind, lapply(enrich_table, function(x) do.call(rbind, lapply(x, function(y) data.frame(YesGeneSet.YesDMR = y[1,1], 
                                                                                               NoGeneSet.YesDMR = y[2,1],
                                                                                               YesGeneSet.NoDMR = y[1,2], 
                                                                                               NoGeneSet.NoDMR = y[2,2])))))
rownames(ta) = NULL
df = cbind(DMRenrich, ta)


write.csv(df,file="/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/DMR/Birnbaum_geneSet_enrichment_DMRs.csv",quote=F)

df[which(df$FDR<=0.05),colnames(df) %in% c("Model","GeneSet","Odds.Ratio", "FDR")]
#                 Model      GeneSet Odds.Ratio          FDR
#1  Celltype_HypoNeuron      ASD CNV   2.864093 1.512611e-05
#2  Celltype_HypoNeuron ASD DATABASE   3.898453 3.865049e-15
#3  Celltype_HypoNeuron           ID   2.695148 4.107475e-03
#4  Celltype_HypoNeuron          NDD   6.012221 2.677592e-04
#6  Celltype_HypoNeuron      SCZ CNV   2.225852 1.935033e-02
#7  Celltype_HypoNeuron      SCZ SNV   2.036401 4.107475e-03
#9    Celltype_HypoGlia ASD DATABASE   2.739534 2.009946e-08
#14   Celltype_HypoGlia      SCZ SNV   1.927204 5.011053e-03
#58          Gr5_cdDMRs ASD DATABASE   5.687546 4.107475e-03
#65          Gr6_cdDMRs ASD DATABASE   3.137463 9.332885e-03
#66          Gr6_cdDMRs           ID   4.551188 1.858431e-02
#67          Gr6_cdDMRs          NDD   9.423341 9.332885e-03



## What are some example regions?

sig = lapply(brCat, function(x) as.character(na.omit(unique(x[which(x$distToGene==0),"EntrezID"]))))
notsig = lapply(sig, function(y) as.character(geneuniverse[!(geneuniverse %in% y)]))

DE_OVERLAP = mapply(function(sig,notsig) lapply(splitSets, function(x) sig[sig %in% x$EntrezGene.ID]), sig, notsig,SIMPLIFY =F)
DE_OVERLAP = mapply(function(de,dmr) lapply(de, function(x) dmr[dmr$EntrezID %in% x,]), DE_OVERLAP, brCat, SIMPLIFY = F)
DE_OVERLAP = lapply(DE_OVERLAP, function(x) lapply(x, function(y) y[which(y$distToGene==0),]))
DE_OVERLAP = lapply(DE_OVERLAP, function(x) lapply(x, makeGRangesFromDataFrame))
DE_OVERLAP = lapply(DE_OVERLAP, function(x) lapply(x, reduce))
DE_OVERLAP = lapply(DE_OVERLAP, function(x) reduce(do.call(getMethod(c, "GenomicRanges"), GRangesList(x))))
DE_OVERLAP = reduce(do.call(getMethod(c, "GenomicRanges"), GRangesList(DE_OVERLAP)))
save(DE_OVERLAP, file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/DMR/Birnbaum_PGC2_geneSet_enrichment_interactionClusters_toPLot.rda")


DE_OVERLAP = mapply(function(sig,notsig) lapply(splitSets, function(x) sig[sig %in% x$EntrezGene.ID]), sig, notsig,SIMPLIFY =F)
DE_OVERLAP = mapply(function(de,dmr) lapply(de, function(x) dmr[dmr$EntrezID %in% x,]), DE_OVERLAP, brCat, SIMPLIFY = F)
DE_OVERLAP = lapply(DE_OVERLAP, function(x) lapply(x, function(y) y[which(y$distToGene==0),]))
DE_OVERLAP = lapply(DE_OVERLAP, function(x) x[elementNROWS(x)>0])
DE_OVERLAP = do.call(rbind, Map(cbind, lapply(DE_OVERLAP, function(x) 
              do.call(rbind, Map(cbind, x, GeneSet = as.list(names(x))))), DMRgroup = as.list(names(DE_OVERLAP))))
rownames(DE_OVERLAP) = NULL

write.csv(DE_OVERLAP, "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/DMR/Birnbaum_enrichment_overlappingGenesList.csv")


