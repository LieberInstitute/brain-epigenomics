load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/sorted_nuclear_RNA/DE_limma_results_objects.rda")

aej_sets = openxlsx::read.xlsx('/dcl01/lieber/ajaffe/lab/Nicotine/DLPFC/RNAseq/csvs/sets/Supplementary Tables for paper.Birnbaum November 2013.AJP.xlsx')

# Test sorted nuclear RNA seq for cell type-specific expression pattern enrichment in the gene sets 

nucRNAres$Sym = geneMap[match(nucRNAres$gencodeID, geneMap$gencodeID),"Symbol"]
Tstat_for_cell_type = nucRNAres$Tstat.CellTypeNeuron
splitSets = split(aej_sets, aej_sets$Gene.Set)

genesetEnrich = sapply(splitSets, function(x) {
  ind = rep(FALSE, nrow(nucRNAres))
  ind[x$Gene.Symbol %in% nucRNAres$Sym] = TRUE
  c(up=geneSetTest(ind,Tstat_for_cell_type, alter = "up"), down=geneSetTest(ind,Tstat_for_cell_type, alter = "down"))
})

write.csv(genesetEnrich,file="/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/sorted_nuclear_RNA/Birnbaum_geneSet_enrichment_sortedNuclearRNA_byCellType.csv",quote=F)