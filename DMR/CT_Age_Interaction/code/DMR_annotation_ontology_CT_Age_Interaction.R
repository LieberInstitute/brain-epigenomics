library(GenomicRanges)
library(GenomicFeatures)
library(data.table)
library(ggplot2)
library(clusterProfiler)
require(org.Hs.eg.db)

load("/dcl01/lieber/ajaffe/CellSorting/RNAseq_pipeline/rawCounts_CellSorting_July5_n12.rda")
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/bumphunting/bumps_bsseqSmooth_Neuron_interaction_250_perm.Rdata")


## How many regions are DMRs? (model: ~ age + cell + age:cell)

interaction = bumps[[1]]
dim(interaction) # 26282     14
dim(interaction[which(interaction$fwer<=0.05),]) # 2178    14


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

grinteraction = makeGRangesFromDataFrame(interaction, keep.extra.columns = T)
annotation = lapply(features, function(y) findOverlaps(grinteraction, y))

grinteraction$rnum = 1:length(grinteraction)
grinteraction$cds = ifelse(grinteraction$rnum %in% queryHits(annotation[["CDS"]]), "CDS", NA)
grinteraction$intron = ifelse(grinteraction$rnum %in% queryHits(annotation[["Introns"]]), "Intron", NA)
grinteraction$UTR5 = ifelse(grinteraction$rnum %in% queryHits(annotation[["UTR5"]]), "5'UTR", NA)
grinteraction$UTR3 = ifelse(grinteraction$rnum %in% queryHits(annotation[["UTR3"]]), "3'UTR", NA)
grinteraction$islands = ifelse(grinteraction$rnum %in% queryHits(annotation[["islands"]]), "CpG Island", "Non-Island")
grinteraction$promoter = ifelse(grinteraction$rnum %in% queryHits(annotation[["promoters"]]), "Promoter", NA)
grinteraction$anno = paste0(grinteraction$cds,":",grinteraction$intron, ":", grinteraction$UTR5, ":", grinteraction$UTR3, ":", grinteraction$promoter)

interaction = as.data.frame(grinteraction)
interaction[which(interaction$anno == "NA:NA:NA:NA:NA"),"annotation"] = "Intergenic" 
interaction[grep("CDS", interaction$cds),"annotation"] = "CDS"
interaction[which(is.na(interaction$annotation) & interaction$UTR5 == "5'UTR"),"annotation"] = "5'UTR"
interaction[which(is.na(interaction$annotation) & interaction$UTR3 == "3'UTR"),"annotation"] = "3'UTR"
interaction[which(is.na(interaction$annotation) & interaction$intron == "Intron"),"annotation"] = "Intron"
interaction[which(is.na(interaction$annotation) & interaction$promoter == "Promoter"),"annotation"] = "Promoter"

# Mapping cell sites to nearest gene
geneMapGR = makeGRangesFromDataFrame(geneMap, start.field="Start",end.field="End",strand.field="Strand",keep=TRUE)
dA = distanceToNearest(grinteraction, geneMapGR)
interaction$nearestSymbol = geneMapGR$Symbol[subjectHits(dA)]
interaction$nearestID = names(geneMapGR)[subjectHits(dA)]
interaction$distToGene = mcols(dA)$distance
interaction$EntrezID = geneMapGR$EntrezID[subjectHits(dA)]
interaction$regionID = paste0(interaction$seqnames,":",interaction$start,"-", interaction$end)
interaction$sig = ifelse(interaction$fwer<=0.05, "FWER < 0.05", "FWER > 0.05")
interaction$Dir = ifelse(interaction$value<0, "neg", "pos")
dtinteraction = data.table(interaction)


### Explore annotation of regions

## how many fall within CpG islands?

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/DMR/CT_Age_Interaction/figures/DMR_overap_with_CpG_Islands_CT_Age_Interaction.pdf",width=8.5)
x = dtinteraction[,length(unique(regionID)), by = "islands"]
x$perc = round(x$V1/sum(x$V1)*100,2)
ggplot(x, aes(x = islands, y = V1)) + geom_bar(stat = "identity") +
  geom_text(aes( label = paste0(perc,"%")), vjust = -.5) +
  labs(fill="") +
  ylab("Count") + 
  xlab("") +
  ggtitle("DMRs Overlapping CpG Islands: All DMRs") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))

x = dtinteraction[sig=="FWER < 0.05",length(unique(regionID)), by = "islands"]
x$perc = round(x$V1/sum(x$V1)*100,2)
ggplot(x, aes(x = islands, y = V1)) + geom_bar(stat = "identity") +
  geom_text(aes( label = paste0(perc,"%")), vjust = -.5) +
  labs(fill="") +
  ylab("Count") + 
  xlab("") +
  ggtitle("DMRs Overlapping CpG Islands: FWER < 0.05") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))

x = dtinteraction[,length(unique(regionID)), by = c("islands", "sig")]
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

fisher.test(data.frame(c(nrow(interaction[which(interaction$sig=="FWER < 0.05" & interaction$islands=="CpG Island"),]),
                         nrow(interaction[which(interaction$sig=="FWER > 0.05" & interaction$islands=="CpG Island"),])),
                       c(nrow(interaction[which(interaction$sig=="FWER < 0.05" & interaction$islands=="Non-Island"),]),
                         nrow(interaction[which(interaction$sig=="FWER > 0.05" & interaction$islands=="Non-Island"),]))))
# CpG islands are overrepresented in DMRs
#p-value < 2.2e-16
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  2.864436 3.825320
#sample estimates:
#  odds ratio 
#3.314157


# assign genomic features

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/DMR/CT_Age_Interaction/figures/DMR_annotation_CT_Age_Interaction.pdf")
x = dtinteraction[,length(unique(regionID)), by = "annotation"]
x$perc = round(x$V1/sum(x$V1)*100,2)
ggplot(x, aes(x = annotation, y = V1)) + geom_bar(stat = "identity") +
  geom_text(aes( label = paste0(perc,"%")), vjust = -.5) +
  labs(fill="") +
  ylab("Count") + 
  xlab("") +
  ggtitle("DMR annotation: All DMRs") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))

x = dtinteraction[sig=="FWER < 0.05",length(unique(regionID)), by = "annotation"]
x$perc = round(x$V1/sum(x$V1)*100,2)
ggplot(x, aes(x = annotation, y = V1)) + geom_bar(stat = "identity") +
  geom_text(aes( label = paste0(perc,"%")), vjust = -.5) +
  labs(fill="") +
  ylab("Count") + 
  xlab("") +
  ggtitle("DMR annotation: FWER < 0.05") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))

x = dtinteraction[,length(unique(regionID)), by = c("annotation", "sig")]
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

fisher.test(data.frame(c(nrow(interaction[which(interaction$sig=="FWER < 0.05" & interaction$distToGene==0),]),
                         nrow(interaction[which(interaction$sig=="FWER > 0.05" & interaction$distToGene==0),])),
                       c(nrow(interaction[which(interaction$sig=="FWER < 0.05" & interaction$distToGene>0),]),
                         nrow(interaction[which(interaction$sig=="FWER > 0.05" & interaction$distToGene>0),]))))
# Genes are overrepresented in DMRs
#p-value = 2.511e-11
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  1.274575 1.574807
#sample estimates:
#  odds ratio 
#1.415813 

# Is there a relationship between being significantly DM and overlapping a promoter?

fisher.test(data.frame(c(nrow(interaction[which(interaction$sig=="FWER < 0.05" & interaction$annotation=="Promoter"),]),
                         nrow(interaction[which(interaction$sig=="FWER > 0.05" & interaction$annotation=="Promoter"),])),
                       c(nrow(interaction[which(interaction$sig=="FWER < 0.05" & interaction$annotation!="Promoter"),]),
                         nrow(interaction[which(interaction$sig=="FWER > 0.05" & interaction$annotation!="Promoter"),]))))
# promoters by themselves aren't overrepresented in DMRs
#p-value = 0.26
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.9219585 1.3248631
#sample estimates:
#  odds ratio 
#1.108489 

# Is there a relationship between being significantly DM and overlapping a gene and/or promoter?

fisher.test(data.frame(c(nrow(interaction[which(interaction$sig=="FWER < 0.05" & interaction$annotation!="Intergenic"),]),
                         nrow(interaction[which(interaction$sig=="FWER > 0.05" & interaction$annotation!="Intergenic"),])),
                       c(nrow(interaction[which(interaction$sig=="FWER < 0.05" & interaction$annotation=="Intergenic"),]),
                         nrow(interaction[which(interaction$sig=="FWER > 0.05" & interaction$annotation=="Intergenic"),]))))
# genes and promoters together are overrepresented in DMRs
#p-value = 1.759e-14
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  1.373659 1.734670
#sample estimates:
#  odds ratio 
#1.542175 


### Gene Ontology
entrez = list(All = dtinteraction[sig=="FWER < 0.05",list(na.omit(EntrezID)),], 
              GenesPlusPromoters = dtinteraction[annotation != "Intergenic" & sig=="FWER < 0.05",list(na.omit(EntrezID)),],
              Genes = dtinteraction[distToGene==0 & sig=="FWER < 0.05",list(na.omit(EntrezID)),],
              Promoters = dtinteraction[annotation == "Promoter" & sig=="FWER < 0.05",list(na.omit(EntrezID)),])
entrez.dir = list(All.pos = dtinteraction[value>0 & sig=="FWER < 0.05",list(na.omit(EntrezID)),], 
                  GenesPlusPromoters.pos = dtinteraction[value>0 & annotation != "Intergenic" & sig=="FWER < 0.05",list(na.omit(EntrezID)),],
                  Genes.pos = dtinteraction[value>0 & distToGene==0 & sig=="FWER < 0.05",list(na.omit(EntrezID)),],
                  Promoters.pos = dtinteraction[value>0 & annotation == "Promoter" & sig=="FWER < 0.05",list(na.omit(EntrezID)),],
                  All.neg = dtinteraction[value<0 & sig=="FWER < 0.05",list(na.omit(EntrezID)),], 
                  GenesPlusPromoters.neg = dtinteraction[value<0 & annotation != "Intergenic" & sig=="FWER < 0.05",list(na.omit(EntrezID)),],
                  Genes.neg = dtinteraction[value<0 & distToGene==0 & sig=="FWER < 0.05",list(na.omit(EntrezID)),],
                  Promoters.neg = dtinteraction[value<0 & annotation == "Promoter" & sig=="FWER < 0.05",list(na.omit(EntrezID)),])
entrez = lapply(entrez, function(x) as.character(unique(x$V1)))              
entrez.dir = lapply(entrez.dir, function(x) as.character(unique(x$V1)))       

GeneUniverse = dtinteraction[,list(na.omit(EntrezID)),]
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
save(compareKegg, compareKegg.dir, compareBP, compareBP.dir, compareMF, compareMF.dir, compareCC, compareCC.dir,
     keggList, keggList.dir, goList_BP, goList_BP.dir, goList_MF, goList_MF.dir, goList_CC, goList_CC.dir, goList_DO, goList_DO.dir,
     keggListdf, keggList.dir.df, goListdf_BP, goListdf_BP.dir, goListdf_MF, goListdf_MF.dir, goListdf_CC, goListdf_CC.dir, goListdf_DO, goListdf_DO.dir,
     file="/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/DMR/CT_Age_Interaction/DMR_KEGG_GO_DO_objects_CT_Age_Interaction.rda")

# plot compared results
pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/DMR/CT_Age_Interaction/figures/DMR_KEGG_GO_DO_plots_CT_Age_Interaction.pdf", height = 20, width = 20)
plot(compareKegg, colorBy="p.adjust", showCategory = 45, title= "KEGG Pathway Enrichment")
plot(compareBP, colorBy="p.adjust", showCategory = 45, title= "Biological Process GO Enrichment")
plot(compareMF, colorBy="p.adjust", showCategory = 45, title= "Molecular Function GO Enrichment")
plot(compareCC, colorBy="p.adjust", showCategory = 45, title= "Cellular Compartment GO Enrichment")
plot(compareKegg.dir, colorBy="p.adjust", showCategory = 45, title= "KEGG Pathway Enrichment")
plot(compareBP.dir, colorBy="p.adjust", showCategory = 45, title= "Biological Process GO Enrichment")
plot(compareMF.dir, colorBy="p.adjust", showCategory = 45, title= "Molecular Function GO Enrichment")
plot(compareCC.dir, colorBy="p.adjust", showCategory = 45, title= "Cellular Compartment GO Enrichment")
dev.off()


### GO analysis of interaction DMRs that don't overlap cell type DMRs
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/DMR/DMR_objects.rda")

DMRgr = lapply(DMR, function(x) makeGRangesFromDataFrame(x[which(x$fwer<=0.05),], keep.extra.columns=T))
CellTypexInt = findOverlaps(DMRgr$CellType, DMRgr$Interaction)

int = DMR$Interaction[which(DMR$Interaction$fwer<=0.05),]
int = int[-subjectHits(CellTypexInt),]

### Gene Ontology
entrez = list(GenesPlusPromoters = na.omit(int[int$annotation != "Intergenic","EntrezID"]),
              Genes = na.omit(int[int$distToGene==0,"EntrezID"]),
              Promoters = na.omit(int[int$annotation == "Promoter","EntrezID"]))
entrez.dir = list(GenesPlusPromoters.pos = na.omit(int[int$value>0 & int$annotation != "Intergenic","EntrezID"]),
                  Genes.pos = na.omit(int[int$value>0 & int$distToGene==0,"EntrezID"]),
                  Promoters.pos = na.omit(int[int$value>0 & int$annotation == "Promoter","EntrezID"]),
                  GenesPlusPromoters.neg = na.omit(int[int$value<0 & int$annotation != "Intergenic","EntrezID"]),
                  Genes.neg = na.omit(int[int$value<0 & int$distToGene==0,"EntrezID"]),
                  Promoters.neg = na.omit(int[int$value<0 & int$annotation == "Promoter","EntrezID"]))
entrez = lapply(entrez, function(x) as.character(unique(x)))              
entrez.dir = lapply(entrez.dir, function(x) as.character(unique(x)))       

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
save(compareKegg, compareKegg.dir, compareBP, compareBP.dir, compareMF, compareMF.dir, compareCC, compareCC.dir,compareDO, compareDO.dir,
     file="/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/DMR/CT_Age_Interaction/DMR_KEGG_GO_DO_objects_Interaction_noCToverlap.rda")

# plot compared results
pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/DMR/CT_Age_Interaction/figures/DMR_KEGG_GO_DO_plots_Interaction_noCToverlap.pdf", height = 20, width = 20)
plot(compareKegg, colorBy="p.adjust", showCategory = 45, title= "KEGG Pathway Enrichment")
plot(compareBP, colorBy="p.adjust", showCategory = 45, title= "Biological Process GO Enrichment")
plot(compareMF, colorBy="p.adjust", showCategory = 45, title= "Molecular Function GO Enrichment")
plot(compareCC, colorBy="p.adjust", showCategory = 45, title= "Cellular Compartment GO Enrichment")
plot(compareDO, colorBy="p.adjust", showCategory = 45, title= "Disease Ontology Enrichment")
plot(compareKegg.dir, colorBy="p.adjust", showCategory = 45, title= "KEGG Pathway Enrichment")
plot(compareBP.dir, colorBy="p.adjust", showCategory = 45, title= "Biological Process GO Enrichment")
plot(compareMF.dir, colorBy="p.adjust", showCategory = 45, title= "Molecular Function GO Enrichment")
plot(compareCC.dir, colorBy="p.adjust", showCategory = 45, title= "Cellular Compartment GO Enrichment")
plot(compareDO.dir, colorBy="p.adjust", showCategory = 45, title= "Disease Ontology Enrichment")
dev.off()

