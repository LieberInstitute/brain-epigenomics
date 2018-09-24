load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/DMR/DMR_objects.rda")
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/bumphunting/rda/limma_Neuron_CpGs_minCov_3_ageInfo_dmrs.Rdata")

aej_sets = openxlsx::read.xlsx('/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/Birnbaum_2013_AJP_Supplementary_table.xlsx')

dmrs = split(dmrs, dmrs$k6cluster_label)
names(dmrs) = c("Group 1 (G-N+)","Group 2 (G0N+)","Group 3 (G0N-)","Group 4 (G+N0)","Group 5 (G+N-)","Group 6 (G-N0)")
oo = lapply(dmrs, function(x) findOverlaps(x, makeGRangesFromDataFrame(DMR$Interaction)))
dmrs = lapply(oo, function(x) DMR$Interaction[subjectHits(x),])
CT = split(DMR$CellType, DMR$CellType$Dir)
names(CT) = c("Hypomethylated in Neurons", "Hypomethylated in Glia")
DMRgr = lapply(c(CT, DMR[names(DMR)!="CellType"], dmrs), function(x) x[which(x$fwer<=0.05),])


## Enrichment in DMRs by cell type, age and interaction

geneuniverse = na.omit(unique(geneMap$EntrezID))
aej_sets_expressed = aej_sets[which(aej_sets$EntrezGene.ID %in% geneuniverse), ] # drop genes that are not present in the test set
splitSets = split(aej_sets_expressed[aej_sets_expressed$Gene.Set!="SCZ PGC GWAS",], aej_sets_expressed$Gene.Set[aej_sets_expressed$Gene.Set!="SCZ PGC GWAS"])
sig = lapply(DMRgr, function(x) as.character(na.omit(unique(x[which(x$distToGene==0),"EntrezID"]))))
notsig = lapply(sig, function(y) as.character(geneuniverse[!(geneuniverse %in% y)]))

DMRenrich = mapply(function(sig,notsig) lapply(splitSets, function(x) {
  DE_OVERLAP = c( sum( sig %in% x$EntrezGene.ID),sum(!(sig %in% x$EntrezGene.ID)))
  NOT_DE_OVERLAP= c(sum(notsig %in% x$EntrezGene.ID), sum(!(notsig %in% x$EntrezGene.ID)))
  enrich_table = cbind(DE_OVERLAP, NOT_DE_OVERLAP)
  res = fisher.test(enrich_table)
  dat=c(res$p.value, res$estimate)
  names(dat) <- c("P.Value","Odds.Ratio")
  return(dat)
}), sig, notsig,SIMPLIFY =F) 


DMRenrich = do.call(rbind, Map(cbind, Model = as.list(names(DMRenrich)), 
                               lapply(DMRenrich, function(x) do.call(rbind, Map(cbind, GeneSet = as.list(names(x)), 
                                                                                lapply(x, function(y) data.frame(P.Value = y["P.Value"], 
                                                                                                                 Odds.Ratio = y["Odds.Ratio"])))))))
DMRenrich$FDR = p.adjust(DMRenrich$P.Value, method = "fdr")
rownames(DMRenrich) = NULL

write.csv(DMRenrich,file="/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/DMR/Birnbaum_geneSet_enrichment_DMRs.csv",quote=F)

DMRenrich[which(DMRenrich$FDR<=0.05),]
#                        Model      GeneSet      P.Value Odds.Ratio          FDR
#1  Hypomethylated\nin Neurons      ASD CNV 6.482617e-07   2.864093 1.458589e-05
#2  Hypomethylated\nin Neurons ASD DATABASE 5.521499e-17   3.898453 4.969349e-15
#3  Hypomethylated\nin Neurons    BPAD GWAS 4.694118e-03   2.089947 2.640441e-02
#4  Hypomethylated\nin Neurons           ID 4.107475e-04   2.695148 4.107475e-03
#5  Hypomethylated\nin Neurons          NDD 1.530052e-05   6.012221 2.754094e-04
#7  Hypomethylated\nin Neurons      SCZ CNV 3.317200e-03   2.225852 1.990320e-02
#9  Hypomethylated\nin Neurons      SCZ SNV 3.021680e-04   2.036401 3.885017e-03
#11    Hypomethylated\nin Glia ASD DATABASE 5.742703e-10   2.739534 2.584216e-08
#18    Hypomethylated\nin Glia      SCZ SNV 5.726918e-04   1.927204 5.154226e-03
#29                Interaction ASD DATABASE 3.281056e-07   2.970068 9.843169e-06
#31                Interaction           ID 8.755463e-04   3.052202 7.163560e-03
#32                Interaction          NDD 5.095723e-03   4.329763 2.697736e-02
#36                Interaction      SCZ SNV 2.076395e-04   2.434389 3.114592e-03
#74            Group 5 (G+,N-) ASD DATABASE 3.698765e-04   5.687546 4.107475e-03
#83            Group 6 (G-,N0) ASD DATABASE 1.306790e-03   3.137463 9.230326e-03
#85            Group 6 (G-,N0)           ID 2.920391e-03   4.551188 1.877394e-02
#86            Group 6 (G-,N0)          NDD 1.333269e-03   9.423341 9.230326e-03



## What are some example regions?

sig = lapply(dmrs, function(x) na.omit(unique(x[which(x$fwer<=0.05),"EntrezID"])))
notsig = lapply(sig, function(y) geneuniverse[!(geneuniverse %in% y)])

DE_OVERLAP = mapply(function(sig,notsig) lapply(splitSets, function(x) sig[sig %in% x$EntrezGene.ID]), sig, notsig,SIMPLIFY =F)
DE_OVERLAP = mapply(function(de,dmr) lapply(de, function(x) dmr[dmr$EntrezID %in% x,]), DE_OVERLAP, dmrs, SIMPLIFY = F)
DE_OVERLAP = lapply(DE_OVERLAP, function(x) lapply(x, function(y) y[which(y$distToGene==0),]))
DE_OVERLAP = lapply(DE_OVERLAP, function(x) lapply(x, makeGRangesFromDataFrame))
DE_OVERLAP = lapply(DE_OVERLAP, function(x) lapply(x, reduce))
DE_OVERLAP = lapply(DE_OVERLAP, function(x) reduce(do.call(getMethod(c, "GenomicRanges"), GRangesList(x))))
DE_OVERLAP = reduce(do.call(getMethod(c, "GenomicRanges"), GRangesList(DE_OVERLAP)))
save(DE_OVERLAP, file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/DMR/Birnbaum_PGC2_geneSet_enrichment_interactionClusters_toPLot.rda")


DE_OVERLAP = mapply(function(sig,notsig) lapply(splitSets, function(x) sig[sig %in% x$EntrezGene.ID]), sig, notsig,SIMPLIFY =F)
DE_OVERLAP = mapply(function(de,dmr) lapply(de, function(x) dmr[dmr$EntrezID %in% x,]), DE_OVERLAP, dmrs, SIMPLIFY = F)
DE_OVERLAP = lapply(DE_OVERLAP, function(x) lapply(x, function(y) y[which(y$distToGene==0),]))
DE_OVERLAP = lapply(DE_OVERLAP, function(x) x[elementNROWS(x)>0])
DE_OVERLAP = lapply(DE_OVERLAP, function(x) do.call(rbind, Map(cbind, x, GeneSet = as.list(names(x)))))

lapply(DE_OVERLAP, function(x) unique(x[x$nearestSymbol %in% c("CACNA1C", "TCF4","CACNA1B","ARX","NRXN3","HDAC4","AKT3"),"regionID"]))
elementNROWS(lapply(DE_OVERLAP, function(x) unique(x[x$regionID %in% c("chr2:240039881-240040446","chr9:140815338-140816026",
                                                          "chr9:140773700-140777324","chr2:240239500-240240989","chr2:240097796-240098355",
                                                          "chr1:243651653-243652298"),])))

unique(DE_OVERLAP$Gr6[DE_OVERLAP$Gr6$regionID %in% c("chr2:240039881-240040446","chr9:140815338-140816026",
                           "chr9:140773700-140777324","chr2:240239500-240240989","chr2:240097796-240098355",
                           "chr1:243651653-243652298"),])



