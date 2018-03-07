library(GenomicRanges)
library(clusterProfiler)
library(org.Hs.eg.db)
library(BSgenome.Hsapiens.UCSC.hg19)
library(PWMEnrich.Hsapiens.background)

load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/UMRs_LMRs_annotated.rda")
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/DMR/DMR_objects.rda")
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/PWMEnrich/DMR_PWMEnrich_objects.rda")
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/PWMEnrich/UMR_LMR_PWMEnrich_object.rda")
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/PWMEnrich/filtered_TF_list_expressionInfo.rda")

## Objects for analysis:
  # all_int: PWMEnrich result for all interaction DMRs
  # intergenic_int: PWMEnrich result for interaction DMRs that overlap intergenic regions
  # introns_int: PWMEnrich result for interaction DMRs that overlap introns
  # promoters_int: PWMEnrich result for interaction DMRs that overlap promoters
  # LMR_UMR_pwmenrich: PWMEnrich result for shared and not shared LMR and UMR sequence for All, Prenatal, Postnatal, Neurons, and Glia 
  # all_split: PWMEnrich result for all regions in three models split by direction
  # intergenic_split: PWMEnrich result for intergenic regions in three models split by direction
  # introns_split: PWMEnrich result for introns in three models split by direction
  # promoters_split: PWMEnrich result for promoters in three models split by direction
  # DMR, uDMR, lDMR: lists of annotated methylated regions 
  # geneMap
  # hompd, pd: phenotype tables
  # targettogeneID: map of target ID to gencodeID
  # TFhomRPKM, TFnucres: RNA-seq results for TFs in homogenate (hom) and nuclear RNA (nuc)
  # TFdiff: difference in test statistic between positive and negative beta values of DMRs (all, intergenic, introns, and promoters)
  # 


## collate TF group report results

TF = list(allTF = lapply(all_split, groupReport), promTF = lapply(promoters_split, groupReport),
          intronTF = lapply(introns_split, groupReport), intergenicTF = lapply(intergenic_split, groupReport))
TF = lapply(TF, function(x) lapply(x, as.data.frame))
TF = lapply(TF, function(t) Map(cbind, t, padj = lapply(t, function(x) p.adjust(x$p.value, method = "fdr"))))
TF = lapply(TF, function(t) lapply(t, function(x) x[which(x$target %in% names(targettogeneID)),]))


### Are different genomic features enriched for different TFs?

## Correlate TF enrichment between features

## Are motifs for the same target similarly enriched?




# excitatory neurons (Mo et al): Egr, Ap-1, Neurod2, Rfx1/3/5, and TBR1
# PV neurons: Mafb/g, Mef2a/c/d
# VIP neurons: developmental role TFs such as Dlx, Pou, Sox family, Arx and Vax2, which are also enriched in fetal and glial hypo-DMRs

# Figures 5D (construction of TF-TF regulatory networks) TF A was predicted to regulate TF B when: 
  # (1) TF A was expressed in a cell type (≥30 TPM), (2) TF A had a predicted footprint (FP A) in a cell type-specific ATAC-seq peak, 
  # (3) the ATAC-seq peak was within 10 kb of the TSS for TF B, and (4) TF B was expressed in that cell type (≥30 TPM). 
  # The resulting set of predicted regulatory interactions was visualized as a network (igraph package in R), 
  # omitting TFs with more than 20 connections to ease visualization. To define a pan-neuronal regulatory network, 
  # we identified footprints common to all three cell types that occurred in shared ATAC-seq peaks and did not overlap ubiquitous DNaseI peaks 
  # (peaks occurring in at least 40 out of 53 processed DNaseI-seq samples). The full networks are listed in Table S4.












## annotate all the TFs to their entrez IDS

length(unique(TF$allTF$Age.pos$target)) # 1549 unique targets
noW = unique(TF$allTF$Age.pos$target)[-grep("UW", unique(TF$allTF$Age.pos$target))]
length(noW) # 866 named targets
noW[-which(noW %in% geneMap$Symbol)]

ensIDs = c("HINFP1" = "ENSG00000172273", "MZF1_1-4" = "ENSG00000099326", "Myf" = "ENSG00000111046", "HIF1A::ARNT" = "ENSG00000100644", 
           "MZF1_5-13" = "ENSG00000099326", "MYC::MAX" = "ENSG00000136997", "RXRA::VDR" = "ENSG00000186350", "Trp53" = "ENSG00000141510", 
           "ZNF306" = "ENSG00000189298", "MIZF" = "ENSG00000172273", "EWSR1-FLI1" = "ENSG00000182944", "BHLHB3" = "ENSG00000123095",
           "BHLHB2" = "ENSG00000134107", "TLX1::NFIC" = "ENSG00000107807", "RORA_1" = "ENSG00000069667", "NR1H2::RXRA" = "ENSG00000131408", 
           "Trp73" = "ENSG00000078900", "JUND (var.2)" = "ENSG00000130522", "TAL1::TCF3" = "ENSG00000162367", "RXR::RAR_DR5" = "ENSG00000186350",
           "Pax6" = "ENSG00000007372", "DUX4" = "ENSG00000260596", "ZNF238" = "ENSG00000179456", "SMAD2::SMAD3::SMAD4" = "ENSG00000175387",
           "JUN (var.2)" = "ENSG00000177606", "NFE2::MAF" = "ENSG00000123405", "BATF::JUN" = "ENSG00000156127", "RAXL1" = "ENSG00000173976", 
           "CART1" = "ENSG00000180318", "ZNF435" = "ENSG00000196812", "TAL1::GATA1" = "ENSG00000162367",  "RORA_2" = "ENSG00000069667", 
           "Sox4" = "ENSG00000124766", "RBM8A" = "ENSG00000265241", "POU5F1P1" = "ENSG00000212993", "AP1" = "ENSG00000248394",
           "STAT2::STAT1" = "ENSG00000170581", "Oct-1" = "ENSG00000143190", "SMCR7L" = "ENSG00000100335", "HNF1B" = "ENSG00000275410",
           "GLTPD1" = "ENSG00000224051", "C9orf156" = "ENSG00000136932", "PRKRIR" = "ENSG00000137492", "C19orf40" = "ENSG00000131944", "MEF2BNB-MEF2B" = "ENSG00000254901")
# DUX4 not included in geneMap

info = data.frame(targetID = c(noW[which(noW %in% geneMap$Symbol)], names(ensIDs)), 
                  geneID = c(geneMap[match(noW[which(noW %in% geneMap$Symbol)], geneMap$Symbol), "gencodeID"], 
                             geneMap[match(ensIDs, geneMap$ensemblID), "gencodeID"]),
                  entrezID = c(geneMap[match(noW[which(noW %in% geneMap$Symbol)], geneMap$Symbol), "EntrezID"], 
                               geneMap[match(ensIDs, geneMap$ensemblID), "EntrezID"]))
dim(info) # 870 3
TF = lapply(TF, function(t) Map(cbind, t, entrezID = lapply(t, function(x) info[match(x$target, info$targetID),"entrezID"])))


## Do promoter and intergenic (distal) DMRs correlate in their enrichment for TFs?

corr = mapply(function(p,d) cor.test(p[order(p$id),"p.value"], d[order(d$id),"p.value"]), TF$promTF[names(TF$promTF)!="Age.neg"], TF$intergenicTF, SIMPLIFY =F)
data.frame(rho = unlist(lapply(corr, function(x) x$estimate)), Tstat = unlist(lapply(corr, function(x) x$statistic)), Pval = unlist(lapply(corr, function(x) x$p.value)))
#                          rho    Tstat          Pval
#CellType.pos.cor    0.5061881 28.05657 4.351829e-149
#CellType.neg.cor    0.6670978 42.80493 1.719231e-294
#Age.pos.cor         0.6667358 42.76310 4.646373e-294
#Interaction.pos.cor 0.6813734 44.49950 5.127435e-312
#Interaction.neg.cor 0.5583861 32.17502 1.168557e-187


## Filter significant enrichment

TFsig = unlist(lapply(TF, function(t) lapply(t, function(x) na.omit(unique(x[which(x$padj<=0.05),"entrezID"])))), recursive=F)
TFsig = lapply(TFsig, as.character)

entrez = list(celltype = TFsig[grep("CellType", names(TFsig))], age = TFsig[grep("Age", names(TFsig))],
              interaction = TFsig[grep("Interaction", names(TFsig))])
names(entrez$celltype) = c("Hypermethylated\n(All)", "Hypomethylated\n(All)", "Hypermethylated\n(Promoters)", "Hypomethylated\n(Promoters)",       
                           "Hypermethylated\n(Introns)", "Hypomethylated\n(Introns)", "Hypermethylated\n(Intergenic)", "Hypomethylated\n(Intergenic)")    
names(entrez$age) = c("Hypermethylated\n(All)", "Hypomethylated\n(All)", "Hypermethylated\n(Promoters)", "Hypomethylated\n(Promoters)",       
                      "Hypermethylated\n(Introns)", "Hypermethylated\n(Intergenic)")    
names(entrez$interaction) = c("Hypermethylated\n(All)", "Hypomethylated\n(All)", "Hypermethylated\n(Promoters)", "Hypomethylated\n(Promoters)",       
                              "Hypermethylated\n(Introns)", "Hypomethylated\n(Introns)", "Hypermethylated\n(Intergenic)", "Hypomethylated\n(Intergenic)")    


## Assess enriched terms limiting the gene universe to the terms associated with the master list of TFs

GeneUniverse = na.omit(unique(as.character(info$entrezID)))
length(GeneUniverse)

# Find enriched Pathways via KEGG
keggList = lapply(entrez, function(y) lapply(y, function(x) enrichKEGG(as.character(x), organism="human", universe= GeneUniverse, 
                                                 minGSSize=5, pAdjustMethod="BH", qvalueCutoff=1)))
keggListdf = lapply(keggList, function(y) lapply(y, function(x) as.data.frame(x)))
# Enriched Molecular Function GOs
goList_MF = lapply(entrez, function(y) lapply(y, function(x) enrichGO(as.character(x), ont = "MF", OrgDb = org.Hs.eg.db, 
                                                universe= GeneUniverse, minGSSize=5, pAdjustMethod="BH", qvalueCutoff=1)))
goListdf_MF = lapply(goList_MF, function(y) lapply(y, function(x) as.data.frame(x)))
# Biological Process GO enrichment
goList_BP = lapply(entrez, function(y) lapply(y, function(x) enrichGO(as.character(x), ont = "BP", OrgDb = org.Hs.eg.db, 
                                                universe= GeneUniverse, minGSSize=5, pAdjustMethod="BH", qvalueCutoff=1)))
goListdf_BP = lapply(goList_BP, function(y) lapply(y, function(x) as.data.frame(x)))
# Cellular Compartment GO enrichment
goList_CC = lapply(entrez, function(y) lapply(y, function(x) enrichGO(as.character(x), ont = "CC", OrgDb = org.Hs.eg.db, 
                                                universe= GeneUniverse, minGSSize=5, pAdjustMethod="BH", qvalueCutoff=1)))
goListdf_CC = lapply(goList_CC, function(y) lapply(y, function(x) as.data.frame(x)))
# Disease Ontology
goList_DO = lapply(entrez, function(y) lapply(y, function(x) enrichDO(as.character(x), ont = "DO", universe= GeneUniverse, 
                                                minGSSize=5, pAdjustMethod="BH", qvalueCutoff=1, readable=TRUE)))
goListdf_DO = lapply(goList_DO, function(y) lapply(y, function(x) as.data.frame(x)))

# Compare the enriched terms
# KEGG
compareKegg = lapply(entrez, function(x) compareCluster(x, fun="enrichKEGG", qvalueCutoff = 0.05, pvalueCutoff = 0.05))
# Biological Process
compareBP = lapply(entrez, function(x) compareCluster(x, fun="enrichGO", ont = "BP", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05))
# Molecular Function
compareMF = lapply(entrez, function(x) compareCluster(x, fun="enrichGO",  ont = "MF", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05))
# Cellular Component
compareCC = lapply(entrez, function(x) compareCluster(x, fun="enrichGO",  ont = "CC", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05))
# Disease Ontology
compareDO = lapply(entrez, function(x) compareCluster(x, fun="enrichDO",  ont = "DO", qvalueCutoff = 0.05, pvalueCutoff = 0.05))

# Save
save(keggListdf, goListdf_MF, goListdf_BP, goListdf_CC, goListdf_DO, compareBP, compareMF, compareCC, compareDO, compareKegg,
     file="/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/PWMEnrich/GO.objects.PWMEnrich.rda")

## plot Gene ontology enrichment

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/PWMEnrich/figures/KEGG_MF_CC_DO_cellType_TFs.pdf", width=20,height=20)
plot(compareKegg$celltype,colorBy="p.adjust",  showCategory = 400, title= "KEGG Pathway Enrichment: Cell Type")
plot(compareMF$celltype,colorBy="p.adjust",  showCategory = 400, title= "MF Pathway Enrichment: Cell Type")
plot(compareCC$celltype,colorBy="p.adjust",  showCategory = 400, title= "CC Pathway Enrichment: Cell Type")
plot(compareDO$celltype,colorBy="p.adjust",  showCategory = 400, title= "Disease Ontology Enrichment: Cell Type")
dev.off()
pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/PWMEnrich/figures/BP_cellType_TFs.pdf", width=22,height=90)
plot(compareBP$celltype,colorBy="p.adjust",  showCategory = 400, title= "Biological Process GO Enrichment: Cell Type")
dev.off()

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/PWMEnrich/figures/KEGG_MF_CC_DO_age_TFs.pdf", width=18,height=14)
plot(compareKegg$age,colorBy="p.adjust",  showCategory = 400, title= "KEGG Pathway Enrichment: Age")
plot(compareMF$age,colorBy="p.adjust",  showCategory = 400, title= "MF Pathway Enrichment: Age")
plot(compareCC$age,colorBy="p.adjust",  showCategory = 400, title= "CC Pathway Enrichment: Age")
plot(compareDO$age,colorBy="p.adjust",  showCategory = 400, title= "Disease Ontology Enrichment: Age")
dev.off()
pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/PWMEnrich/figures/BP_age_TFs.pdf", width=16,height=80)
plot(compareBP$age,colorBy="p.adjust",  showCategory = 400, title= "Biological Process GO Enrichment: Age")
dev.off()

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/PWMEnrich/figures/KEGG_MF_CC_DO_interaction_TFs.pdf", width=18,height=20)
plot(compareKegg$interaction,colorBy="p.adjust",  showCategory = 400, title= "KEGG Pathway Enrichment: Interaction")
plot(compareMF$interaction,colorBy="p.adjust",  showCategory = 400, title= "MF Pathway Enrichment: Interaction")
plot(compareCC$interaction,colorBy="p.adjust",  showCategory = 400, title= "CC Pathway Enrichment: Interaction")
plot(compareDO$interaction,colorBy="p.adjust",  showCategory = 400, title= "Disease Ontology Enrichment: Interaction")
dev.off()
pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/PWMEnrich/figures/BP_interaction_TFs.pdf", width=18,height=80)
plot(compareBP$interaction,colorBy="p.adjust",  showCategory = 400, title= "Biological Process GO Enrichment: Interaction")
dev.off()


#### Explore Differential TF binding results

## Explore distribution of statistics 

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/PWMEnrich/figures/TFdiff_distribution.pdf")
for (i in 1:length(TFdiff)) {
h = hist(TFdiff[[i]]$group.bg)
print(h)
}
dev.off()

TFstats = list()
for (i in 1:length(allTF)) {
  TFstats[[i]] = data.frame(max = max(allTF[[i]]$raw.score), min = min(allTF[[i]]$raw.score), mean = mean(allTF[[i]]$raw.score),
                            median = median(allTF[[i]]$raw.score), row.names = names(allTF)[i])
}
TFstats = do.call(rbind, TFstats)
#                        max          min       mean   median
#Age.pos            1160.988 1.453014e-05   2.279963 1.184627
#Age.neg         1222196.263 1.451343e-06 539.750818 1.077923
#Interaction.pos  123159.926 1.790824e-04  56.594015 1.492703
#Interaction.neg  270312.842 1.920342e-05 121.869597 1.232809
#CellType.pos     261998.592 1.350354e-03 117.092000 1.486878
#CellType.neg     282960.419 4.265580e-04 126.183604 1.475019


## Filter results

diffres = lapply(TFdiff, function(x) data.frame(target = TF$allTF$Age.pos$target, id = TF$allTF$Age.pos$id, 
                                                TFdiff = x$group.bg[match(TF$allTF$Age.pos$id, names(x$group.bg))]))
difffilt = lapply(diffres, function(x) x[which(x$target %in% info$targetID & abs(x$TFdiff)>=100),])
elementNROWS(difffilt)
elementNROWS(lapply(difffilt, function(x) unique(x$target)))


## Split by direction (positive mean greater methylation in older samples and neurons)

difffilt = unlist(lapply(difffilt, function(x) split(x, x$TFdiff>0)), recursive = F)
diffentrez = lapply(difffilt, function(x) info[which(info$targetID %in% x$target),"entrezID"])
diffentrez = lapply(diffentrez, function(x) unique(as.character(na.omit(x))))
diffentrez = list(Age = list("Increasing\n(All DMRs)" = diffentrez$allAge.TRUE, "Increasing\n(Promoters)" = diffentrez$promAge.TRUE,
                             "Decreasing\n(All DMRs)" = diffentrez$allAge.FALSE, "Decreasing\n(Promoters)" = diffentrez$promAge.FALSE),
                  CellType = list("Hypermethylated\n(All DMRs)" = diffentrez$allCT.TRUE, "Hypermethylated\n(Promoters)" = diffentrez$promCT.TRUE,
                                  "Hypermethylated\n(Introns)" = diffentrez$intronsCT.TRUE, "Hypermethylated\n(Intergenic)" = diffentrez$interCT.TRUE,
                                  "Hypomethylated\n(All DMRs)" = diffentrez$allCT.FALSE, "Hypomethylated\n(Promoters)" = diffentrez$promCT.FALSE,
                                  "Hypomethylated\n(Introns)" = diffentrez$intronsCT.FALSE, "Hypomethylated\n(Intergenic)" = diffentrez$interCT.FALSE),
                  Interaction = list("Hypermethylated\n(All DMRs)" = diffentrez$allInt.TRUE, "Hypermethylated\n(Promoters)" = diffentrez$promInt.TRUE,
                                     "Hypermethylated\n(Introns)" = diffentrez$intronsInt.TRUE, "Hypermethylated\n(Intergenic)" = diffentrez$interInt.TRUE,
                                     "Hypomethylated\n(All DMRs)" = diffentrez$allInt.FALSE, "Hypomethylated\n(Promoters)" = diffentrez$promInt.FALSE))

# term enrichment
TFdiffKegg = lapply(diffentrez, function(x) compareCluster(x, fun="enrichKEGG", qvalueCutoff = 0.05, pvalueCutoff = 0.05))
TFdiffBP = lapply(diffentrez, function(x) compareCluster(x, fun="enrichGO", ont = "BP", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05))
TFdiffMF = lapply(diffentrez, function(x) compareCluster(x, fun="enrichGO",  ont = "MF", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05))
TFdiffCC = lapply(diffentrez, function(x) compareCluster(x, fun="enrichGO",  ont = "CC", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05))
TFdiffDO = lapply(diffentrez, function(x) compareCluster(x, fun="enrichDO",  ont = "DO", qvalueCutoff = 0.05, pvalueCutoff = 0.05))

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/PWMEnrich/figures/KEGG_BP_MF_CC_DO_statgreaterthan100_diffTFs.pdf", width=20,height=80)
for (i in 1:length(TFdiffKegg)) {
d = plot(TFdiffKegg[[i]],colorBy="p.adjust",  showCategory = 400,title= paste0("KEGG Pathway Enrichment: ",names(TFdiffKegg)[i]))
print(d)
d = plot(TFdiffMF[[i]],colorBy="p.adjust",  showCategory = 400,title= paste0("MF Pathway Enrichment: ",names(TFdiffMF)[i]))
print(d)
d = plot(TFdiffCC[[i]],colorBy="p.adjust",  showCategory = 400,title= paste0("CC Pathway Enrichment: ",names(TFdiffCC)[i]))
print(d)
d = plot(TFdiffBP[[i]],colorBy="p.adjust",  showCategory = 400,title= paste0("Biological Process GO Enrichment: ",names(TFdiffBP)[i]))
print(d)
d = plot(TFdiffDO[[i]],colorBy="p.adjust",  showCategory = 400,title= paste0("Biological Process GO Enrichment: ",names(TFdiffDO)[i]))
print(d)
}
dev.off()




# other tools
# GeneNetworkBuilder: Appliation for discovering direct or indirect targets of transcription factors using ChIP-chip or ChIP-seq, and microarray or RNA-seq gene expression data. Inputting a list of genes of potential targets of one TF from ChIP-chip or ChIP-seq, and the gene expression results, GeneNetworkBuilder generates a regulatory network of the TF.
# http://bioconductor.org/packages/release/bioc/html/GeneNetworkBuilder.html

# RTN: reconstruction of transcriptional networks and analysis of master regulators.
# Circos plots
# HiveR and Gephi plots, http://www.vesnam.com/Rblog/viznets3/
# mfinder: Network motifs detection tool








