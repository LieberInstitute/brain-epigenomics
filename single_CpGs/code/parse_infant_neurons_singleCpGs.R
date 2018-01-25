library(bsseq)
library('limma')
library(pheatmap)
library(RColorBrewer)
library(genefilter)
library(data.table)
library(ggplot2)
library(clusterProfiler)
require(org.Hs.eg.db)


load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/single_CpGs/singleCpG_object.rda")


## load data
load('/dcl01/lieber/ajaffe/lab/brain-epigenomics/bumphunting/BSobj_bsseqSmooth_Neuron_minCov_3.Rdata')
postpd = pData(BSobj)
postpd$Race[postpd$Race== "CAUC "] <- 'CAUC'
postpd$Sex[postpd$Sex == " M"] <- 'M'
postpd$RIN <- as.numeric(gsub(" ", "", postpd$RIN))
postpd$pg.DNA.nuclei.input <- as.numeric(postpd$pg.DNA.nuclei.input)
postpd$Reads <- as.numeric(postpd$Reads)
postpd$Percent.GreaterThan.Q30 <- as.numeric(postpd$Percent.GreaterThan.Q30)

## get mean meth per DMR
postmeth <- getMeth(BSobj, type = 'raw')

load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/bumphunting/BSobj_bsseqSmooth_Neuron_minCov_3_prenatal.Rdata")
prenmeth = getMeth(BSobj, type = 'raw')
prenpd = pData(BSobj)

meth = data.frame(postmeth, prenmeth)
meth = as.matrix(meth)


### All CpGs
sampleDists <- dist(t(meth))
sampleDistMatrix <- as.matrix(sampleDists)
colnames(sampleDistMatrix) = rownames(sampleDistMatrix) = c(paste(
		postpd $Cell.Type, postpd $Age.Bin, postpd $Working.Num, sep = ":"), paste0("Prenatal:",prenpd$Brain.Num))

### interaction
MethInt = meth[which(CpG$padj.Interaction<=0.05),]
sampleDistsInt <- dist(t(MethInt))
sampleDistMatrixInt <- as.matrix(sampleDistsInt)
colnames(sampleDistMatrixInt) = rownames(sampleDistMatrixInt) = c(paste(
		postpd $Cell.Type, postpd $Age.Bin, postpd $Working.Num, sep = ":"), paste0("Prenatal:",prenpd$Brain.Num))

### age
MethAge = meth[which(CpG$padj.Age<=0.05),]
sampleDistsAge <- dist(t(MethAge))
sampleDistMatrixAge <- as.matrix(sampleDistsAge)
colnames(sampleDistMatrixAge) = rownames(sampleDistMatrixAge) = c(paste(
		postpd $Cell.Type, postpd $Age.Bin, postpd $Working.Num, sep = ":"), paste0("Prenatal:",prenpd$Brain.Num))

### cell type
MethCT = meth[which(CpG$padj.CellType<=0.05),]
sampleDistsCT <- dist(t(MethCT))
sampleDistMatrixCT <- as.matrix(sampleDistsCT)
colnames(sampleDistMatrixCT) = rownames(sampleDistMatrixCT) = c(paste(
		postpd $Cell.Type, postpd $Age.Bin, postpd $Working.Num, sep = ":"), paste0("Prenatal:",prenpd$Brain.Num))
		
		
pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/single_CpGs/figures/singleCpGs_heatmap_euclidean_dist_3models_plusFetals.pdf",width=10,height=10)
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix, clustering_distance_rows = sampleDists, clustering_distance_cols = sampleDists,
         col=colors, main = "All CpGs - Euclidean Distance")
pheatmap(sampleDistMatrixInt,clustering_distance_rows=sampleDistsInt, clustering_distance_cols=sampleDistsInt,
         col=colors, main = "90,227 Interaction CpGs - Euclidean Distance")
pheatmap(sampleDistMatrixAge,clustering_distance_rows=sampleDistsAge, clustering_distance_cols=sampleDistsAge,
         col=colors, main = "536,164 Age CpGs - Euclidean Distance")
pheatmap(sampleDistMatrixCT, clustering_distance_rows=sampleDistsCT, clustering_distance_cols=sampleDistsCT,
         col=colors, main = "4,824,804 Cell Type CpGs - Euclidean Distance")
dev.off()


## which CpGs show the association?
postpd$Neonate = ifelse(postpd$Age < 1, "Neonate", "Older") 
gIndexes = list(YoungNeuronVsGlia = which(postpd$Age < 1 | postpd$Cell.Type == "Glia"),
	YoungVsOldNeuron = which(postpd$Cell.Type == "Neuron"))

tt_YoungNeuronVsGlia_Int = rowttests(MethInt[,gIndexes$YoungNeuronVsGlia],
	factor(postpd$Cell.Type[gIndexes$YoungNeuronVsGlia])) # neg is glia effect
tt_YoungVsOldNeuron_Int = rowttests(MethInt[,gIndexes$YoungVsOldNeuron],
	factor(postpd$Neonate[gIndexes$YoungVsOldNeuron])) # neg is neonate effect
plot(tt_YoungNeuronVsGlia_Int$statistic, tt_YoungVsOldNeuron_Int$statistic)
plot(-log10(tt_YoungNeuronVsGlia_Int$p.value), -log10(tt_YoungVsOldNeuron_Int$p.value))

table(tt_YoungNeuronVsGlia_Int$p.value < 0.01, tt_YoungVsOldNeuron_Int$p.value < 0.01,
	dnn = c("YoungNeuronVsGlia", "YoungVsOldNeuron"))
#                 YoungVsOldNeuron
#YoungNeuronVsGlia FALSE  TRUE
#            FALSE 42788 17656
#            TRUE  21318  8173


## Identify the CpG sites in each quadrant of the table
intCpG = CpG[which(CpG$padj.Interaction<=0.05),]
tableCells = split(intCpG, paste0(tt_YoungNeuronVsGlia_Int$p.value < 0.01, 
	"-",tt_YoungVsOldNeuron_Int$p.value < 0.01))
tableCells = tableCells[-grep("NA", names(tableCells))] # drop the missing ones

## Plot distribution of annotated features
tableCells = Map(cbind, tableCells, quadrant = list("FALSE-FALSE", "FALSE-TRUE", "TRUE-FALSE", "TRUE-TRUE"))
tableCellsdt = data.table(do.call(rbind, tableCells))


# how many fall within CpG islands?
pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/single_CpGs/figures/singleCpGs_overlap_with_CpG_Islands_infant_neurons_parsing_byQuadrant.pdf",width=16)
x = tableCellsdt[,length(unique(regionID)), by = c("islands", "quadrant")] 
x$perc = c(round(x$V1[1:2]/sum(x$V1[1:2])*100,2),round(x$V1[3:4]/sum(x$V1[3:4])*100,2),
           round(x$V1[5:6]/sum(x$V1[5:6])*100,2),round(x$V1[7:8]/sum(x$V1[7:8])*100,2))
ggplot(x, aes(x = quadrant, y = V1, fill=islands)) + geom_bar(stat = "identity") +
  geom_text(aes( label = paste0(perc,"%")), vjust = -.5) +
  labs(fill="") +
  ylab("Count") + 
  xlab("") +
  ggtitle("CpGs Overlapping CpG Islands By Quadrant") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
dev.off()


# anntotate CpGs w/in each list
pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/single_CpGs/figures/singleCpGs_annotation_infant_neurons_parsing_byQuadrant.pdf",height=18,width=10)
x = tableCellsdt[,length(unique(regionID)), by = c("annotation", "quadrant")]
x$perc = unlist(c(round(x[1:6,"V1"]/sum(x[1:6,"V1"])*100,2),
                  round(x[7:12,"V1"]/sum(x[7:12,"V1"])*100,2),
                  round(x[13:18,"V1"]/sum(x[13:18,"V1"])*100,2),
                  round(x[19:24,"V1"]/sum(x[19:24,"V1"])*100,2)))
ggplot(x, aes(x = quadrant, y = V1, fill=annotation)) + geom_bar(stat = "identity") + # fix so that percentages aren't all bunched up
  geom_text(aes(x = quadrant, y = V1, label = paste0(perc,"%"), group = annotation),
            position = position_stack(vjust = .5)) +
  labs(fill="") +
  ylab("Count") + 
  xlab("") +
  ggtitle("CpG annotation by Quadrant") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
dev.off()


## Identify the Entrez IDs in each group
entrez = lapply(tableCells, function(x) as.character(unique(na.omit(x[which(x$annotation != "Intergenic"),"EntrezID"]))))
entrez.dir = list("GenesPlusPromoters.pos\n" = lapply(tableCells, function(x) as.character(unique(na.omit(x[which(x$Tstat.Interaction > 0 & x$annotation != "Intergenic"),"EntrezID"])))),
                  "GenesPlusPromoters.neg\n" = lapply(tableCells, function(x) as.character(unique(na.omit(x[which(x$Tstat.Interaction < 0 & x$annotation != "Intergenic"),"EntrezID"])))))
entrez.dir = unlist(entrez.dir, recursive=F)

GeneUniverse = as.character(unique(na.omit(tableCellsdt$EntrezID)))

# Find enriched Pathways via KEGG
elementNROWS(entrez)
keggList = lapply(entrez, function(x) enrichKEGG(as.character(x), organism="human", universe= GeneUniverse, 
                                                 minGSSize=5, pAdjustMethod="BH", qvalueCutoff=1))
keggList.dir = lapply(entrez.dir, function(x) enrichKEGG(as.character(x), organism="human", universe= GeneUniverse, 
                                                         minGSSize=5, pAdjustMethod="BH", qvalueCutoff=1))
# Enriched Molecular Function GOs
goList_MF = lapply(entrez, function(x) enrichGO(as.character(x), ont = "MF", OrgDb = org.Hs.eg.db, 
                                                universe= GeneUniverse, minGSSize=5, pAdjustMethod="BH",
                                                qvalueCutoff=1))
goList_MF.dir = lapply(entrez.dir, function(x) enrichGO(as.character(x), ont = "MF", OrgDb = org.Hs.eg.db, 
                                                        universe= GeneUniverse, minGSSize=5, pAdjustMethod="BH",
                                                        qvalueCutoff=1))
# Biological Process GO enrichment
goList_BP = lapply(entrez, function(x) enrichGO(as.character(x), ont = "BP", OrgDb = org.Hs.eg.db, 
                                                universe= GeneUniverse, minGSSize=5, pAdjustMethod="BH",
                                                qvalueCutoff=1))
goList_BP.dir = lapply(entrez.dir, function(x) enrichGO(as.character(x), ont = "BP", OrgDb = org.Hs.eg.db, 
                                                        universe= GeneUniverse, minGSSize=5, pAdjustMethod="BH",
                                                        qvalueCutoff=1))
# Cellular Compartment GO enrichment
goList_CC = lapply(entrez, function(x) enrichGO(as.character(x), ont = "CC", OrgDb = org.Hs.eg.db, 
                                                universe= GeneUniverse, minGSSize=5, pAdjustMethod="BH",
                                                qvalueCutoff=1))
goList_CC.dir = lapply(entrez.dir, function(x) enrichGO(as.character(x), ont = "CC", OrgDb = org.Hs.eg.db, 
                                                        universe= GeneUniverse, minGSSize=5, pAdjustMethod="BH",
                                                        qvalueCutoff=1))
# Disease Ontology
goList_DO = lapply(entrez, function(x) enrichDO(as.character(x), ont = "DO", universe= GeneUniverse, 
                                                minGSSize=5, pAdjustMethod="BH", qvalueCutoff=1, readable=TRUE))
goList_DO.dir = lapply(entrez.dir, function(x) enrichDO(as.character(x), ont = "DO", universe= GeneUniverse, 
                                                        minGSSize=5, pAdjustMethod="BH", qvalueCutoff=1, readable=TRUE))

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

# save
save(compareKegg, compareKegg.dir, compareBP, compareBP.dir, compareMF, compareMF.dir, compareCC, compareCC.dir, compareDO.dir,
     keggList, keggList.dir, goList_BP, goList_BP.dir, goList_MF, goList_MF.dir, goList_CC, goList_CC.dir, goList_DO, goList_DO.dir,
     file="/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/single_CpGs/singleCpGs_KEGG_GO_DO_objects_infant_neuron_parsing_byQuadrant.rda")

# plot the GO results
pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/single_CpGs/figures/singleCpGs_KEGG_GO_DO_plots_infant_neuron_parsing_byQuadrant.pdf", height = 24, width = 24)
plot(compareKegg, colorBy="p.adjust", showCategory = 45, title= "KEGG Pathway Enrichment")
plot(compareBP, colorBy="p.adjust", showCategory = 45, title= "Biological Process GO Enrichment")
plot(compareMF, colorBy="p.adjust", showCategory = 45, title= "Molecular Function GO Enrichment")
plot(compareCC, colorBy="p.adjust", showCategory = 45, title= "Cellular Compartment GO Enrichment")
plot(compareKegg.dir, colorBy="p.adjust", showCategory = 45, title= "KEGG Pathway Enrichment")
plot(compareBP.dir, colorBy="p.adjust", showCategory = 45, title= "Biological Process GO Enrichment")
plot(compareMF.dir, colorBy="p.adjust", showCategory = 45, title= "Molecular Function GO Enrichment")
plot(compareCC.dir, colorBy="p.adjust", showCategory = 45, title= "Cellular Compartment GO Enrichment")
plot(compareDO.dir, colorBy="p.adjust", showCategory = 45, title= "Disease Ontology Enrichment")
dev.off()
