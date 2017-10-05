library(GenomicRanges)
library(GenomicFeatures)
library(data.table)
library(ggplot2)
library(clusterProfiler)
require(org.Hs.eg.db)

load("/dcl01/lieber/ajaffe/CellSorting/RNAseq_pipeline/rawCounts_CellSorting_July5_n12.rda")
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/bumphunting/bumps_bsseqSmooth_Neuron_cell_250_perm.Rdata")


## How many regions are DMRs? (model: ~ age + cell)

cell = bumps[[1]]
dim(cell) # 241924     14
dim(cell[which(cell$fwer<=0.05),]) # 11179    14


# Annotate editing sites to features in the genome
txdb = loadDb("/dcl01/lieber/ajaffe/lab/brain-epigenomics/annotation_objects/gencode.v25lift37.annotation.sqlite")
islands = read.table("/dcl01/lieber/ajaffe/lab/brain-epigenomics/annotation_objects/cpgIslandExt.hg19.txt", sep="\t", header = T)
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

grcell = makeGRangesFromDataFrame(cell, keep.extra.columns = T)
annotation = lapply(features, function(y) findOverlaps(grcell, y))

grcell$rnum = 1:length(grcell)
grcell$cds = ifelse(grcell$rnum %in% queryHits(annotation[["CDS"]]), "CDS", NA)
grcell$intron = ifelse(grcell$rnum %in% queryHits(annotation[["Introns"]]), "Intron", NA)
grcell$UTR5 = ifelse(grcell$rnum %in% queryHits(annotation[["UTR5"]]), "UTR5", NA)
grcell$UTR3 = ifelse(grcell$rnum %in% queryHits(annotation[["UTR3"]]), "UTR3", NA)
grcell$islands = ifelse(grcell$rnum %in% queryHits(annotation[["islands"]]), "CpG_Island", "non-Island")
grcell$promoter = ifelse(grcell$rnum %in% queryHits(annotation[["promoters"]]), "Promoter", NA)
grcell$anno = paste0(grcell$cds,":",grcell$intron, ":", grcell$UTR5, ":", grcell$UTR3, ":", grcell$promoter)

cell = as.data.frame(grcell)
cell[which(cell$anno == "NA:NA:NA:NA:NA"),"annotation"] = "Other" 
cell[grep("CDS", cell$cds),"annotation"] = "CDS"
cell[which(is.na(cell$annotation) & cell$UTR5 == "UTR5"),"annotation"] = "UTR5"
cell[which(is.na(cell$annotation) & cell$UTR3 == "UTR3"),"annotation"] = "UTR3"
cell[which(is.na(cell$annotation) & cell$intron == "Intron"),"annotation"] = "Intron"
cell[which(is.na(cell$annotation) & cell$promoter == "Promoter"),"annotation"] = "Promoter"

# Mapping cell sites to nearest gene
geneMapGR = makeGRangesFromDataFrame(geneMap, start.field="Start",end.field="End",strand.field="Strand",keep=TRUE)
dA = distanceToNearest(grcell, geneMapGR)
cell$nearestSymbol = geneMapGR$Symbol[subjectHits(dA)]
cell$nearestID = names(geneMapGR)[subjectHits(dA)]
cell$distToGene = mcols(dA)$distance
cell$EntrezID = geneMapGR$EntrezID[subjectHits(dA)]
cell$regionID = paste0(cell$seqnames,":",cell$start,"-", cell$end)
cell$sig = ifelse(cell$fwer<=0.05, "FWER < 0.05", "FWER > 0.05")
cell$Dir = ifelse(cell$value<0, "neg", "pos")
dtcell = data.table(cell)


### Explore annotation of regions

## how many fall within CpG islands?

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/DMR/CellType/DMR_overap_with_CpG_Islands_byCellType.pdf")
x = dtcell[,length(unique(regionID)), by = "islands"]
x$perc = round(x$V1/sum(x$V1)*100,2)
ggplot(x, aes(x = islands, y = V1)) + geom_bar(stat = "identity") +
  geom_text(aes( label = paste0(perc,"%")), vjust = -.5) +
  labs(fill="") +
  ylab("Count") + 
  xlab("") +
  ggtitle("DMRs Overlapping CpG Islands: All DMRs") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))

x = dtcell[sig=="FWER < 0.05",length(unique(regionID)), by = "islands"]
x$perc = round(x$V1/sum(x$V1)*100,2)
ggplot(x, aes(x = islands, y = V1)) + geom_bar(stat = "identity") +
  geom_text(aes( label = paste0(perc,"%")), vjust = -.5) +
  labs(fill="") +
  ylab("Count") + 
  xlab("") +
  ggtitle("DMRs Overlapping CpG Islands: FWER < 0.05") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))

x = dtcell[,length(unique(regionID)), by = c("islands", "sig")]
x$perc = unlist(c(round(x[1:2,"V1"]/sum(x[1:2,"V1"])*100,2),
           round(x[3:4,"V1"]/sum(x[3:4,"V1"])*100,2)))
ggplot(x, aes(x = islands, y = V1)) + geom_bar(stat = "identity") +
  geom_text(aes( label = paste0(perc,"%")), vjust = -.5) +
  facet_grid(. ~ sig) +
  labs(fill="") +
  ylab("Count") + 
  xlab("") +
  ggtitle("DMRs Overlapping CpG Islands") +
  theme(axis.text.x=element_text(angle=45,hjust=1)) +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
dev.off()

# Is there a relationship between being significantly DM and overlapping a CpG island?

fisher.test(data.frame(c(nrow(cell[which(cell$sig=="FWER < 0.05" & cell$islands=="CpG_Island"),]),
                         nrow(cell[which(cell$sig=="FWER > 0.05" & cell$islands=="CpG_Island"),])),
                       c(nrow(cell[which(cell$sig=="FWER < 0.05" & cell$islands=="non-Island"),]),
                         nrow(cell[which(cell$sig=="FWER > 0.05" & cell$islands=="non-Island"),]))))
# CpG islands are overrepresented in DMRs
#p-value < 2.2e-16
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  5.070749 5.818562
#sample estimates:
#  odds ratio 
#5.432825


# assign genomic features

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/DMR/CellType/DMR_annotation_byCellType.pdf")
x = dtcell[,length(unique(regionID)), by = "annotation"]
x$perc = round(x$V1/sum(x$V1)*100,2)
ggplot(x, aes(x = annotation, y = V1)) + geom_bar(stat = "identity") +
  geom_text(aes( label = paste0(perc,"%")), vjust = -.5) +
  labs(fill="") +
  ylab("Count") + 
  xlab("") +
  ggtitle("DMR annotation: All DMRs") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))

x = dtcell[sig=="FWER < 0.05",length(unique(regionID)), by = "annotation"]
x$perc = round(x$V1/sum(x$V1)*100,2)
ggplot(x, aes(x = annotation, y = V1)) + geom_bar(stat = "identity") +
  geom_text(aes( label = paste0(perc,"%")), vjust = -.5) +
  labs(fill="") +
  ylab("Count") + 
  xlab("") +
  ggtitle("DMR annotation: FWER < 0.05") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))

x = dtcell[,length(unique(regionID)), by = c("annotation", "sig")]
x$perc = unlist(c(round(x[1:6,"V1"]/sum(x[1:6,"V1"])*100,2),
                  round(x[7:12,"V1"]/sum(x[7:12,"V1"])*100,2)))
ggplot(x, aes(x = annotation, y = V1)) + geom_bar(stat = "identity") +
  geom_text(aes( label = paste0(perc,"%")), vjust = -.5) +
  facet_grid(. ~ sig) +
  labs(fill="") +
  ylab("Count") + 
  theme(axis.text.x=element_text(angle=45,hjust=1)) +
  xlab("") +
  ggtitle("DMR annotation") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
dev.off()

# Is there a relationship between being significantly DM and overlapping a gene?

fisher.test(data.frame(c(nrow(cell[which(cell$sig=="FWER < 0.05" & cell$distToGene==0),]),
                         nrow(cell[which(cell$sig=="FWER > 0.05" & cell$distToGene==0),])),
                       c(nrow(cell[which(cell$sig=="FWER < 0.05" & cell$distToGene>0),]),
                         nrow(cell[which(cell$sig=="FWER > 0.05" & cell$distToGene>0),]))))
# Genes are overrepresented in DMRs
#p-value < 2.2e-16
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  2.268516 2.480824
#sample estimates:
#  odds ratio 
#2.371966 

# Is there a relationship between being significantly DM and overlapping a promoter?

fisher.test(data.frame(c(nrow(cell[which(cell$sig=="FWER < 0.05" & cell$annotation=="Promoter"),]),
                         nrow(cell[which(cell$sig=="FWER > 0.05" & cell$annotation=="Promoter"),])),
                       c(nrow(cell[which(cell$sig=="FWER < 0.05" & cell$annotation!="Promoter"),]),
                         nrow(cell[which(cell$sig=="FWER > 0.05" & cell$annotation!="Promoter"),]))))
# promoters by themselves aren't overrepresented in DMRs
#p-value = 0.5655
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.931030 1.132771
#sample estimates:
#  odds ratio 
#1.028018 

# Is there a relationship between being significantly DM and overlapping a gene and/or promoter?

fisher.test(data.frame(c(nrow(cell[which(cell$sig=="FWER < 0.05" & cell$annotation!="Other"),]),
                         nrow(cell[which(cell$sig=="FWER > 0.05" & cell$annotation!="Other"),])),
                       c(nrow(cell[which(cell$sig=="FWER < 0.05" & cell$annotation=="Other"),]),
                         nrow(cell[which(cell$sig=="FWER > 0.05" & cell$annotation=="Other"),]))))
# genes and promoters together are overrepresented in DMRs
#p-value < 2.2e-16
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.3801622 0.4175264
#sample estimates:
#  odds ratio 
#0.3984746


### Gene Ontology
entrez = list(All = dtcell[sig=="FWER < 0.05",list(na.omit(EntrezID)),], 
              GenesPlusPromoters = dtcell[annotation != "Other" & sig=="FWER < 0.05",list(na.omit(EntrezID)),],
              Genes = dtcell[distToGene==0 & sig=="FWER < 0.05",list(na.omit(EntrezID)),],
              Promoters = dtcell[annotation == "Promoter" & sig=="FWER < 0.05",list(na.omit(EntrezID)),])
entrez.dir = list(All.pos = dtcell[value>0 & sig=="FWER < 0.05",list(na.omit(EntrezID)),], 
                  GenesPlusPromoters.pos = dtcell[value>0 & annotation != "Other" & sig=="FWER < 0.05",list(na.omit(EntrezID)),],
                  Genes.pos = dtcell[value>0 & distToGene==0 & sig=="FWER < 0.05",list(na.omit(EntrezID)),],
                  Promoters.pos = dtcell[value>0 & annotation == "Promoter" & sig=="FWER < 0.05",list(na.omit(EntrezID)),],
                  All.neg = dtcell[value<0 & sig=="FWER < 0.05",list(na.omit(EntrezID)),], 
                  GenesPlusPromoters.neg = dtcell[value<0 & annotation != "Other" & sig=="FWER < 0.05",list(na.omit(EntrezID)),],
                  Genes.neg = dtcell[value<0 & distToGene==0 & sig=="FWER < 0.05",list(na.omit(EntrezID)),],
                  Promoters.neg = dtcell[value<0 & annotation == "Promoter" & sig=="FWER < 0.05",list(na.omit(EntrezID)),])
entrez = lapply(entrez, function(x) as.character(unique(x$V1)))              
entrez.dir = lapply(entrez.dir, function(x) as.character(unique(x$V1)))       

GeneUniverse = dtcell[,list(na.omit(EntrezID)),]
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
     file="/dcl01/lieber/ajaffe/lab/brain-epigenomics/DMR/CellType/KEGG_GO_DO_objects.rda")

# plot compared results
pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/DMR/CellType/DMR_KEGG_GO_DO_plots_byCellType.pdf", height = 20, width = 20)
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