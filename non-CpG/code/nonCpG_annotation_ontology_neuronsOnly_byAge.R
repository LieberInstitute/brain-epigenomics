library(GenomicRanges)
library(GenomicFeatures)
library(data.table)
library(ggplot2)
library(clusterProfiler)
require(org.Hs.eg.db)

load("/dcl01/lieber/ajaffe/CellSorting/RNAseq_pipeline/rawCounts_CellSorting_July5_n12.rda")
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/non-CpG/nonCG_highCov_neuronsOnly_pca_pd_methTable.Rdata")
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/non-CpG/limma_exploration_nonCG_highCov_neuronsOnly.Rdata")

### Just Neurons ###

## How many sites are DM? (model: ~ age)

CHneurons = cbind(as.data.frame(gr), lods = ebResults$lods[,2], Tstat = ebResults$t[,2], pval = ebResults$p.value[,2])
CHneurons$padj = p.adjust(CHneurons$pval, method = "fdr")
dim(CHneurons) # 36099584       11
dim(CHneurons[which(CHneurons$padj<=0.05),]) # 4020371      11 


# Annotate editing sites to features in the genome
txdb = loadDb("/dcl01/lieber/ajaffe/Amanda/annotation_objects/gencode.v25lift37.annotation.sqlite")
islands = read.table("/dcl01/lieber/ajaffe/Amanda/annotation_objects/cpgIslandExt.hg19.txt", sep="\t", header = T)
features = list(CDS = cdsBy(txdb, by="tx", use.names=T), Introns = intronsByTranscript(txdb, use.names=T), 
                UTR5 = fiveUTRsByTranscript(txdb, use.names=T), UTR3 = threeUTRsByTranscript(txdb, use.names=T))
features = lapply(features, function(x) unlist(x, recursive = TRUE, use.names = TRUE))
for (i in 1:length(features)){
  tmp = features[[i]]
  tmp$TxID = names(tmp)
  features[[i]] = tmp
}
features = c(features, islands = makeGRangesFromDataFrame(islands, keep.extra.columns = T, start.field = "chromStart", end.field = "chromEnd"),
             promoters = promoters(txdb, upstream=2000, downstream=200))
lapply(features, head)

grCHneurons = makeGRangesFromDataFrame(CHneurons, keep.extra.columns = T)
annotation = lapply(features, function(y) findOverlaps(grCHneurons, y))

grCHneurons$rnum = 1:length(grCHneurons)
grCHneurons$cds = ifelse(grCHneurons$rnum %in% queryHits(annotation[["CDS"]]), "CDS", NA)
grCHneurons$intron = ifelse(grCHneurons$rnum %in% queryHits(annotation[["Introns"]]), "Intron", NA)
grCHneurons$UTR5 = ifelse(grCHneurons$rnum %in% queryHits(annotation[["UTR5"]]), "UTR5", NA)
grCHneurons$UTR3 = ifelse(grCHneurons$rnum %in% queryHits(annotation[["UTR3"]]), "UTR3", NA)
grCHneurons$islands = ifelse(grCHneurons$rnum %in% queryHits(annotation[["islands"]]), "CpG_Island", "non-Island")
grCHneurons$promoter = ifelse(grCHneurons$rnum %in% queryHits(annotation[["promoters"]]), "Promoter", NA)
grCHneurons$anno = paste0(grCHneurons$cds,":",grCHneurons$intron, ":", grCHneurons$UTR5, ":", grCHneurons$UTR3, ":", grCHneurons$promoter)

CHneurons = as.data.frame(grCHneurons)
CHneurons[which(CHneurons$anno == "NA:NA:NA:NA:NA"),"annotation"] = "Other" 
CHneurons[grep("CDS", CHneurons$cds),"annotation"] = "CDS"
CHneurons[which(is.na(CHneurons$annotation) & CHneurons$UTR5 == "UTR5"),"annotation"] = "UTR5"
CHneurons[which(is.na(CHneurons$annotation) & CHneurons$UTR3 == "UTR3"),"annotation"] = "UTR3"
CHneurons[which(is.na(CHneurons$annotation) & CHneurons$intron == "Intron"),"annotation"] = "Intron"
CHneurons[which(is.na(CHneurons$annotation) & CHneurons$promoter == "Promoter"),"annotation"] = "Promoter"

# Mapping cell sites to nearest gene
geneMapGR = makeGRangesFromDataFrame(geneMap, start.field="Start",end.field="End",strand.field="Strand",keep=TRUE)
dA = distanceToNearest(grCHneurons, geneMapGR)
CHneurons$nearestSymbol = geneMapGR$Symbol[subjectHits(dA)]
CHneurons$nearestID = names(geneMapGR)[subjectHits(dA)]
CHneurons$distToGene = mcols(dA)$distance
CHneurons$EntrezID = geneMapGR$EntrezID[subjectHits(dA)]
CHneurons$regionID = paste0(CHneurons$seqnames,":",CHneurons$start,"-", CHneurons$end)
CHneurons$sig = ifelse(CHneurons$padj<=0.05, "FDR < 0.05", "FDR > 0.05")
CHneurons$Dir = ifelse(CHneurons$Tstat<0, "neg", "pos")
dtCHneurons = data.table(CHneurons)
save(CHneurons, file="/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/non-CpG/CHneurons_object.rda")
write.csv(CHneurons[which(CHneurons$padj<=0.001),], quote = F, 
          file="/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/non-CpG/CH_neuronsOnly_FDR_0.001.csv")

### Explore annotation of regions

## how many fall within CpG islands?

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/non-CpG/figures/non-CpG_overap_with_CpG_Islands_neuronsOnly_byAge.pdf")
x = dtCHneurons[,length(unique(regionID)), by = "islands"]
x$perc = round(x$V1/sum(x$V1)*100,2)
ggplot(x, aes(x = islands, y = V1)) + geom_bar(stat = "identity") +
  geom_text(aes( label = paste0(perc,"%")), vjust = -.5) +
  labs(fill="") +
  ylab("Count") + 
  xlab("") +
  ggtitle("non-CpGs Overlapping CpG Islands:\nAll non-CpGs") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))

x = dtCHneurons[sig=="FDR < 0.05",length(unique(regionID)), by = "islands"]
x$perc = round(x$V1/sum(x$V1)*100,2)
ggplot(x, aes(x = islands, y = V1)) + geom_bar(stat = "identity") +
  geom_text(aes( label = paste0(perc,"%")), vjust = -.5) +
  labs(fill="") +
  ylab("Count") + 
  xlab("") +
  ggtitle("non-CpGs Overlapping CpG Islands:\nFDR < 0.05") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))

x = dtCHneurons[,length(unique(regionID)), by = c("islands", "sig")]
x$perc = unlist(c(round(x[1,"V1"]/sum(x[c(1,3),"V1"])*100,2),
                  round(x[2,"V1"]/sum(x[c(2,4),"V1"])*100,2),
                  round(x[3,"V1"]/sum(x[c(1,3),"V1"])*100,2),
                  round(x[4,"V1"]/sum(x[c(2,4),"V1"])*100,2)))
ggplot(x, aes(x = islands, y = V1)) + geom_bar(stat = "identity") +
  geom_text(aes( label = paste0(perc,"%")), vjust = -.5) +
  facet_grid(. ~ sig) +
  labs(fill="") +
  ylab("Count") + 
  xlab("") +
  ggtitle("non-CpGs Overlapping CpG Islands") +
  theme(axis.text.x=element_text(angle=45,hjust=1)) +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
dev.off()

# Is there a relationship between being significantly DM and overlapping a CpG island?

fisher.test(data.frame(c(nrow(CHneurons[which(CHneurons$sig=="FDR < 0.05" & CHneurons$islands=="CpG_Island"),]),
                         nrow(CHneurons[which(CHneurons$sig=="FDR > 0.05" & CHneurons$islands=="CpG_Island"),])),
                       c(nrow(CHneurons[which(CHneurons$sig=="FDR < 0.05" & CHneurons$islands=="non-Island"),]),
                         nrow(CHneurons[which(CHneurons$sig=="FDR > 0.05" & CHneurons$islands=="non-Island"),]))))
# CpG islands are depleted in dmCH
#p-value < 2.2e-16
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.1925635 0.2024278
#sample estimates:
#  odds ratio 
#0.1974567

# assign genomic features

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/non-CpG/figures/non-CpG_annotation_neuronsOnly_byAge.pdf")
x = dtCHneurons[,length(unique(regionID)), by = "annotation"]
x$perc = round(x$V1/sum(x$V1)*100,2)
ggplot(x, aes(x = annotation, y = V1)) + geom_bar(stat = "identity") +
  geom_text(aes( label = paste0(perc,"%")), vjust = -.5) +
  labs(fill="") +
  ylab("Count") + 
  xlab("") +
  ggtitle("non-CpG annotation: All non-CpGs") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))

x = dtCHneurons[sig=="FDR < 0.05",length(unique(regionID)), by = "annotation"]
x$perc = round(x$V1/sum(x$V1)*100,2)
ggplot(x, aes(x = annotation, y = V1)) + geom_bar(stat = "identity") +
  geom_text(aes( label = paste0(perc,"%")), vjust = -.5) +
  labs(fill="") +
  ylab("Count") + 
  xlab("") +
  ggtitle("non-CpG annotation: FDR < 0.05") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))

x = dtCHneurons[,length(unique(regionID)), by = c("annotation", "sig")]
x$perc = NA
x$perc[grep("<", x$sig)] = unlist(c(round(x[grep("<", x$sig),"V1"]/sum(x[grep("<", x$sig),"V1"])*100,2)))
x$perc[grep(">", x$sig)] = unlist(c(round(x[grep(">", x$sig),"V1"]/sum(x[grep(">", x$sig),"V1"])*100,2)))
ggplot(x, aes(x = annotation, y = V1)) + geom_bar(stat = "identity") +
  geom_text(aes( label = paste0(perc,"%")), vjust = -.5) +
  facet_grid(. ~ sig) +
  labs(fill="") +
  ylab("Count") + 
  theme(axis.text.x=element_text(angle=45,hjust=1)) +
  xlab("") +
  ggtitle("non-CpG annotation") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
dev.off()

# Is there a relationship between being significantly DM and overlapping a gene?

fisher.test(data.frame(c(nrow(CHneurons[which(CHneurons$sig=="FDR < 0.05" & CHneurons$distToGene==0),]),
                         nrow(CHneurons[which(CHneurons$sig=="FDR > 0.05" & CHneurons$distToGene==0),])),
                       c(nrow(CHneurons[which(CHneurons$sig=="FDR < 0.05" & CHneurons$distToGene>0),]),
                         nrow(CHneurons[which(CHneurons$sig=="FDR > 0.05" & CHneurons$distToGene>0),]))))
# Genes are underrepresented in dmCH
#p-value < 2.2e-16
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.8328637 0.8367461
#sample estimates:
#  odds ratio 
#0.8347676  

# Is there a relationship between being significantly DM and overlapping a promoter?

fisher.test(data.frame(c(nrow(CHneurons[which(CHneurons$sig=="FDR < 0.05" & CHneurons$annotation=="Promoter"),]),
                         nrow(CHneurons[which(CHneurons$sig=="FDR > 0.05" & CHneurons$annotation=="Promoter"),])),
                       c(nrow(CHneurons[which(CHneurons$sig=="FDR < 0.05" & CHneurons$annotation!="Promoter"),]),
                         nrow(CHneurons[which(CHneurons$sig=="FDR > 0.05" & CHneurons$annotation!="Promoter"),]))))
# promoters are underrepresented in dmCH
#p-value < 2.2e-16
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.9166502 0.9308900
#sample estimates:
#  odds ratio 
#0.9237535

# Is there a relationship between being significantly DM and overlapping a gene and/or promoter?

fisher.test(data.frame(c(nrow(CHneurons[which(CHneurons$sig=="FDR < 0.05" & CHneurons$annotation!="Other"),]),
                         nrow(CHneurons[which(CHneurons$sig=="FDR > 0.05" & CHneurons$annotation!="Other"),])),
                       c(nrow(CHneurons[which(CHneurons$sig=="FDR < 0.05" & CHneurons$annotation=="Other"),]),
                         nrow(CHneurons[which(CHneurons$sig=="FDR > 0.05" & CHneurons$annotation=="Other"),]))))
# genes and promoters together are underrepresented in dmCH
#p-value < 2.2e-16
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.8311978 0.8350066
#sample estimates:
#  odds ratio 
#0.8331474


### Gene Ontology
entrez = list(All = dtCHneurons[sig=="FDR < 0.05",list(na.omit(EntrezID)),], 
              GenesPlusPromoters = dtCHneurons[annotation != "Other" & sig=="FDR < 0.05",list(na.omit(EntrezID)),],
              Genes = dtCHneurons[distToGene==0 & sig=="FDR < 0.05",list(na.omit(EntrezID)),],
              Promoters = dtCHneurons[annotation == "Promoter" & sig=="FDR < 0.05",list(na.omit(EntrezID)),])
entrez.dir = list(All.pos = dtCHneurons[Tstat>0 & sig=="FDR < 0.05",list(na.omit(EntrezID)),], 
                  GenesPlusPromoters.pos = dtCHneurons[Tstat>0 & annotation != "Other" & sig=="FDR < 0.05",list(na.omit(EntrezID)),],
                  Genes.pos = dtCHneurons[Tstat>0 & distToGene==0 & sig=="FDR < 0.05",list(na.omit(EntrezID)),],
                  Promoters.pos = dtCHneurons[Tstat>0 & annotation == "Promoter" & sig=="FDR < 0.05",list(na.omit(EntrezID)),],
                  All.neg = dtCHneurons[Tstat<0 & sig=="FDR < 0.05",list(na.omit(EntrezID)),], 
                  GenesPlusPromoters.neg = dtCHneurons[Tstat<0 & annotation != "Other" & sig=="FDR < 0.05",list(na.omit(EntrezID)),],
                  Genes.neg = dtCHneurons[Tstat<0 & distToGene==0 & sig=="FDR < 0.05",list(na.omit(EntrezID)),],
                  Promoters.neg = dtCHneurons[Tstat<0 & annotation == "Promoter" & sig=="FDR < 0.05",list(na.omit(EntrezID)),])
entrez = lapply(entrez, function(x) as.character(unique(x$V1)))              
entrez.dir = lapply(entrez.dir, function(x) as.character(unique(x$V1)))       

GeneUniverse = dtCHneurons[,list(na.omit(EntrezID)),]
GeneUniverse = as.character(unique(GeneUniverse$V1))

# Find enriched Pathways via KEGG
elementNROWS(entrez)
keggList = lapply(entrez, function(x) enrichKEGG(as.character(x), organism="human", universe= GeneUniverse, 
                                                   minGSSize=5, pAdjustMethod="BH", qvalueCutoff=1))
keggListdf = lapply(keggList, function(x) as.data.frame(x))
keggList.dir = lapply(entrez.dir, function(x) enrichKEGG(as.character(x), organism="human", universe= GeneUniverse, 
                                                   minGSSize=5, pAdjustMethod="BH", qvalueCutoff=1))
keggList.dir.df = lapply(keggList.dir, function(x) as.data.frame(x))

# Enriched Molecular Function GOs
goList_MF = lapply(entrez, function(x) enrichGO(as.character(x), ont = "MF", OrgDb = org.Hs.eg.db, 
                                                  universe= GeneUniverse, minGSSize=5, pAdjustMethod="BH",
                                                  qvalueCutoff=1))
goListdf_MF = lapply(goList_MF, function(x) as.data.frame(x))
goList_MF.dir = lapply(entrez.dir, function(x) enrichGO(as.character(x), ont = "MF", OrgDb = org.Hs.eg.db, 
                                                universe= GeneUniverse, minGSSize=5, pAdjustMethod="BH",
                                                qvalueCutoff=1))
goListdf_MF.dir = lapply(goList_MF.dir, function(x) as.data.frame(x))

# Biological Process GO enrichment
goList_BP = lapply(entrez, function(x) enrichGO(as.character(x), ont = "BP", OrgDb = org.Hs.eg.db, 
                                                universe= GeneUniverse, minGSSize=5, pAdjustMethod="BH",
                                                qvalueCutoff=1))
goListdf_BP = lapply(goList_BP, function(x) as.data.frame(x))
goList_BP.dir = lapply(entrez.dir, function(x) enrichGO(as.character(x), ont = "BP", OrgDb = org.Hs.eg.db, 
                                                        universe= GeneUniverse, minGSSize=5, pAdjustMethod="BH",
                                                        qvalueCutoff=1))
goListdf_BP.dir = lapply(goList_BP.dir, function(x) as.data.frame(x))

# Cellular Compartment GO enrichment
goList_CC = lapply(entrez, function(x) enrichGO(as.character(x), ont = "CC", OrgDb = org.Hs.eg.db, 
                                                universe= GeneUniverse, minGSSize=5, pAdjustMethod="BH",
                                                qvalueCutoff=1))
goListdf_CC = lapply(goList_CC, function(x) as.data.frame(x))
goList_CC.dir = lapply(entrez.dir, function(x) enrichGO(as.character(x), ont = "CC", OrgDb = org.Hs.eg.db, 
                                                        universe= GeneUniverse, minGSSize=5, pAdjustMethod="BH",
                                                        qvalueCutoff=1))
goListdf_CC.dir = lapply(goList_CC.dir, function(x) as.data.frame(x))

# Disease Ontology
goList_DO = lapply(entrez, function(x) enrichDO(as.character(x), ont = "DO", universe= GeneUniverse, 
                                                  minGSSize=5, pAdjustMethod="BH", qvalueCutoff=1, readable=TRUE))
goListdf_DO = lapply(goList_DO, function(x) as.data.frame(x))
goList_DO.dir = lapply(entrez.dir, function(x) enrichDO(as.character(x), ont = "DO", universe= GeneUniverse, 
                                                minGSSize=5, pAdjustMethod="BH", qvalueCutoff=1, readable=TRUE))
goListdf_DO.dir = lapply(goList_DO.dir, function(x) as.data.frame(x))

# Compare the enriched terms between 7 groups
# KEGG
compareKegg = compareCluster(entrez, fun="enrichKEGG", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
compareKegg.dir = compareCluster(entrez.dir, fun="enrichKEGG", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
# Biological Process
compareBP = compareCluster(entrez, fun="enrichGO", ont = "BP", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
compareBP.dir = compareCluster(entrez.dir, fun="enrichGO",  ont = "BP", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
# Molecular Function
compareMF = compareCluster(entrez, fun="enrichGO",  ont = "MF", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
compareMF.dir = compareCluster(entrez.dir, fun="enrichGO",  ont = "MF", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
# Cellular Component
compareCC = compareCluster(entrez, fun="enrichGO",  ont = "CC", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
compareCC.dir = compareCluster(entrez.dir, fun="enrichGO",  ont = "CC", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
# Disease Ontology
compareDO = compareCluster(entrez, fun="enrichDO",  ont = "DO", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
compareDO.dir = compareCluster(entrez.dir, fun="enrichDO",  ont = "DO", qvalueCutoff = 0.05, pvalueCutoff = 0.05)

# save, write to csv
save(compareKegg, compareKegg.dir, compareBP, compareBP.dir, compareMF, compareMF.dir, compareCC, compareCC.dir, compareDO, compareDO.dir,
     keggList, keggList.dir, goList_BP, goList_BP.dir, goList_MF, goList_MF.dir, goList_CC, goList_CC.dir, goList_DO, goList_DO.dir,
     keggListdf, keggList.dir.df, goListdf_BP, goListdf_BP.dir, goListdf_MF, goListdf_MF.dir, goListdf_CC, goListdf_CC.dir, goListdf_DO, goListdf_DO.dir,
     file="/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/non-CpG/non-CpG_KEGG_GO_DO_objects_neuronsOnly_byAge.rda")

# plot compared results
pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/non-CpG/figures/non-CpG_KEGG_GO_DO_plots_neuronsOnly_byAge.pdf", height = 20, width = 20)
plot(compareKegg, colorBy="p.adjust", showCategory = 45, title= "KEGG Pathway Enrichment")
plot(compareBP, colorBy="p.adjust", showCategory = 45, title= "Biological Process GO Enrichment")
plot(compareMF, colorBy="p.adjust", showCategory = 45, title= "Molecular Function GO Enrichment")
plot(compareCC, colorBy="p.adjust", showCategory = 45, title= "Cellular Compartment GO Enrichment")
plot(compareDO, colorBy="p.adjust", showCategory = 30, title= "Disease Ontology Enrichment")
plot(compareKegg.dir, colorBy="p.adjust", showCategory = 45, title= "KEGG Pathway Enrichment")
plot(compareBP.dir, colorBy="p.adjust", showCategory = 45, title= "Biological Process GO Enrichment")
plot(compareMF.dir, colorBy="p.adjust", showCategory = 45, title= "Molecular Function GO Enrichment")
plot(compareCC.dir, colorBy="p.adjust", showCategory = 45, title= "Cellular Compartment GO Enrichment")
plot(compareDO.dir, colorBy="p.adjust", showCategory = 30, title= "Disease Ontology Enrichment")
dev.off()