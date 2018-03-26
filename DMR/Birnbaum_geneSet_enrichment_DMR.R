load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/DMR/DMR_objects.rda")

aej_sets = openxlsx::read.xlsx('/dcl01/lieber/ajaffe/lab/Nicotine/DLPFC/RNAseq/csvs/sets/Supplementary Tables for paper.Birnbaum November 2013.AJP.xlsx')

## Enrichment in DMRs by cell type, age and interaction

geneuniverse = na.omit(unique(geneMap$EntrezID))
aej_sets_expressed = aej_sets[which(aej_sets$EntrezGene.ID %in% geneuniverse), ] # drop genes that are not present in the test set
splitSets = split(aej_sets_expressed, aej_sets_expressed$Gene.Set)
sig = lapply(DMR, function(x) na.omit(unique(x[which(x$fwer<=0.05),"EntrezID"])))
notsig = lapply(sig, function(y) geneuniverse[!(geneuniverse %in% y)])

DMRenrich = mapply(function(sig,notsig) lapply(splitSets, function(x) {
  DE_OVERLAP = c( sum( sig %in% x$EntrezGene.ID),sum(!(sig %in% x$EntrezGene.ID)))
  NOT_DE_OVERLAP= c(sum(notsig %in% x$EntrezGene.ID), sum(!(notsig %in% x$EntrezGene.ID)))
  enrich_table = cbind(DE_OVERLAP, NOT_DE_OVERLAP)
  res = fisher.test(enrich_table)
  dat=c(res$p.value, res$estimate)
  names(dat) <- c("P.Value","Odds.Ratio")
  return(dat)
}), sig, notsig,SIMPLIFY =F) 


DMRenrich = do.call(rbind, Map(cbind, lapply(DMRenrich, function(x) do.call(rbind, Map(cbind, lapply(x, function(y) data.frame(P.Value = y["P.Value"], Odds.Ratio = y["Odds.Ratio"])), 
                               GeneSet = as.list(names(x))))), Model = as.list(names(DMRenrich))))
DMRenrich = rbind(DMRenrich, data.frame(P.Value = c(1.973096e-06,5.978461e-03,1.541955e-03), 
                                        Odds.Ratio = c(1.770269,3.880647,1.816927), 
                                        GeneSet = rep.int("PGC2",3), Model = c("CellType","Age","Interaction")))


## Add interaction clusters

load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/bumphunting/rda/limma_Neuron_CpGs_minCov_3_ageInfo_dmrs.Rdata")

dmrs = split(dmrs, dmrs$k6cluster_label)
names(dmrs) = c("Gr1","Gr2","Gr3","Gr4","Gr5","Gr6")
oo = lapply(dmrs, function(x) findOverlaps(x, makeGRangesFromDataFrame(DMR$Interaction)))
dmrs = lapply(oo, function(x) DMR$Interaction[subjectHits(x),])

sig = lapply(dmrs, function(x) na.omit(unique(x[which(x$fwer<=0.05),"EntrezID"])))
notsig = lapply(sig, function(y) geneuniverse[!(geneuniverse %in% y)])

DMRenrichClusters = mapply(function(sig,notsig) lapply(splitSets, function(x) {
  DE_OVERLAP = c( sum( sig %in% x$EntrezGene.ID),sum(!(sig %in% x$EntrezGene.ID)))
  NOT_DE_OVERLAP= c(sum(notsig %in% x$EntrezGene.ID), sum(!(notsig %in% x$EntrezGene.ID)))
  enrich_table = cbind(DE_OVERLAP, NOT_DE_OVERLAP)
  res = fisher.test(enrich_table)
  dat=c(res$p.value, res$estimate)
  names(dat) <- c("P.Value","Odds.Ratio")
  return(dat)
}), sig, notsig,SIMPLIFY =F) 


DMRenrichClusters = do.call(rbind, Map(cbind, lapply(DMRenrichClusters, function(x) do.call(rbind, Map(cbind, lapply(x, function(y) data.frame(P.Value = y["P.Value"], Odds.Ratio = y["Odds.Ratio"])), 
                                                                                       GeneSet = as.list(names(x))))), Model = as.list(names(DMRenrichClusters))))
DMRenrich = rbind(DMRenrich, DMRenrichClusters)
DMRenrich$padj = p.adjust(DMRenrich$P.Value, method = "fdr")


write.csv(DMRenrich,file="/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/DMR/Birnbaum_PGC2_geneSet_enrichment_DMRs.csv",quote=F)

DMRenrich[which(DMRenrich$padj<=0.05),]
y = data.frame("Gene Set" = x[which(x$FDR<=0.05),"GeneSet"], Model = x[which(x$FDR<=0.05),"Model"],
               "Odds Ratio" = x[which(x$FDR<=0.05),"Odds.Ratio"], "P Value" = x[which(x$FDR<=0.05),"P.Value"],
               FDR = x[which(x$FDR<=0.05),"FDR"])





unique(DMR$Interaction[which(DMR$Interaction$nearestID %in% PGCgenes$gencodeID & DMR$Interaction$distToGene==0),"nearestSymbol"])


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



