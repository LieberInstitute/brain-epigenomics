##
library(bsseq)
library(RColorBrewer)
library(pheatmap)
library(genefilter)
library(data.table)
library(ggplot2)
library(clusterProfiler)
require(org.Hs.eg.db)


## load dmrs
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/bumphunting/bumps_bsseqSmooth_Neuron_interaction_250_perm.Rdata")
sigInt = bumps$table[bumps$table$fwer < 0.05,]
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/bumphunting/bumps_bsseqSmooth_Neuron_age_250_perm.Rdata")
sigAge = bumps$table[bumps$table$fwer < 0.05,]
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/bumphunting/bumps_bsseqSmooth_Neuron_cell_250_perm.Rdata")
sigCT = bumps$table[bumps$table$fwer < 0.05,]
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/DMR/DMR_objects.rda")

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

meth = dataframe(postmeth, prenmeth)
meth = as.matrix(meth)

### interaction
topIndsInt = mapply(function(s,e) s:e, sigInt$indexStart, sigInt$indexEnd)
meanMethInt = sapply(topIndsInt, function(ii) colMeans(t(t(meth[ii,]))))
meanMethInt = do.call("rbind", meanMethInt)

sampleDistsInt <- dist(t(meanMethInt))
sampleDistMatrixInt <- as.matrix(sampleDistsInt)
colnames(sampleDistMatrixInt) = rownames(sampleDistMatrixInt) = c(paste(
		pd$Cell.Type, pd$Age.Bin, pd$Working.Num, sep = ":"), paste0("Prenatal:",prenpd$Brain.Num))

### age
topIndsAge = mapply(function(s,e) s:e, sigAge$indexStart, sigAge$indexEnd)
meanMethAge = sapply(topIndsAge, function(ii) colMeans(t(t(meth[ii,]))))
meanMethAge = do.call("rbind", meanMethAge)

sampleDistsAge <- dist(t(meanMethAge))
sampleDistMatrixAge <- as.matrix(sampleDistsAge)
colnames(sampleDistMatrixAge) = rownames(sampleDistMatrixAge) = c(paste(
		pd$Cell.Type, pd$Age.Bin, pd$Working.Num, sep = ":"), paste0("Prenatal:",prenpd$Brain.Num))
		
### cell type
topIndsCT = mapply(function(s,e) s:e, sigCT$indexStart, sigCT$indexEnd)
meanMethCT = sapply(topIndsCT, function(ii) colMeans(t(t(meth[ii,]))))
meanMethCT = do.call("rbind", meanMethCT)

sampleDistsCT <- dist(t(meanMethCT))
sampleDistMatrixCT <- as.matrix(sampleDistsCT)
colnames(sampleDistMatrixCT) = rownames(sampleDistMatrixCT) = c(paste(
		pd$Cell.Type, pd$Age.Bin, pd$Working.Num, sep = ":"), paste0("Prenatal:",prenpd$Brain.Num))

		
pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/DMR/figures/heatmap_euclidean_dist_3models_includingFetals.pdf",width=10,height=10)
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrixInt,clustering_distance_rows=sampleDistsInt, clustering_distance_cols=sampleDistsInt,
         col=colors, main = "2178 Interaction DMRs - Euclidean Distance")
pheatmap(sampleDistMatrixAge,clustering_distance_rows=sampleDistsAge, clustering_distance_cols=sampleDistsAge,
         col=colors, main = "129 Age DMRs - Euclidean Distance")
pheatmap(sampleDistMatrixCT, clustering_distance_rows=sampleDistsCT, clustering_distance_cols=sampleDistsCT,
         col=colors, main = "11179 Cell Type DMRs - Euclidean Distance")
dev.off()

## which DMRs show the association?
postpd$Neonate = ifelse(postpd$Age < 1, "Neonate", "Older") 
gIndexes = list(YoungNeuronVsGlia = which(postpd$Age < 1 | postpd$Cell.Type == "Glia"),
	YoungVsOldNeuron = which(postpd$Cell.Type == "Neuron"))

tt_YoungNeuronVsGlia_Int = rowttests(meanMethInt[,gIndexes$YoungNeuronVsGlia],
	factor(postpd$Cell.Type[gIndexes$YoungNeuronVsGlia])) # neg is glia effect
tt_YoungVsOldNeuron_Int = rowttests(meanMethInt[,gIndexes$YoungVsOldNeuron],
	factor(postpd$Neonate[gIndexes$YoungVsOldNeuron])) # neg is neonate effect
plot(tt_YoungNeuronVsGlia_Int$statistic, tt_YoungVsOldNeuron_Int$statistic)
plot(-log10(tt_YoungNeuronVsGlia_Int$p.value), -log10(tt_YoungVsOldNeuron_Int$p.value))

table(tt_YoungNeuronVsGlia_Int$p.value < 0.01, tt_YoungVsOldNeuron_Int$p.value < 0.01,
	dnn = c("YoungNeuronVsGlia", "YoungVsOldNeuron"))
	

## Identify the regions in each quadrant of the table
tableCells = split(sigInt, paste0(tt_YoungNeuronVsGlia_Int$p.value < 0.01, 
	"-",tt_YoungVsOldNeuron_Int$p.value < 0.01))
tableCells = tableCells[names(tableCells)!="NA-NA"] # drop the missing one
tableCells = Map(cbind, tableCells,regionID = lapply(tableCells, function(x) paste0(x$chr, ":", x$start,"-", x$end)))

## Plot distribution of annotated features
tableCells = Map(cbind, lapply(tableCells, function(x) DMR$CellType[match(x$regionID,DMR$Interaction$regionID),]), quadrant = list("FALSE-FALSE", "FALSE-TRUE", "TRUE-FALSE", "TRUE-TRUE"))
tableCells = do.call(rbind, tableCells)
tableCellsdt = data.table(tableCells)

# how many fall within CpG islands?
pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/DMR/figures/DMR_overap_with_CpG_Islands_infant_neurons_parsing_byQuadrant.pdf",width=16)
x = tableCellsdt[,length(unique(regionID)), by = c("islands", "quadrant")] 
x$perc = c(round(x$V1[1:2]/sum(x$V1[1:2])*100,2),round(x$V1[3:4]/sum(x$V1[3:4])*100,2),
           round(x$V1[5:6]/sum(x$V1[5:6])*100,2),round(x$V1[7:8]/sum(x$V1[7:8])*100,2))
ggplot(x, aes(x = quadrant, y = V1, fill=islands)) + geom_bar(stat = "identity") +
  geom_text(aes( label = paste0(perc,"%")), vjust = -.5) +
  labs(fill="") +
  ylab("Count") + 
  xlab("") +
  ggtitle("DMRs Overlapping CpG Islands By Quadrant") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
dev.off()


# anntotate DMRs w/in each list
pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/DMR/figures/DMR_annotation_infant_neurons_parsing_byQuadrant.pdf",height=18,width=10)
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
  ggtitle("DMR annotation by Quadrant") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
dev.off()


## Identify the Entrez IDs in each group
entrez = split(tableCells[which(tableCells$annotation != "Intergenic" & tableCells$sig=="FWER < 0.05"),"EntrezID"], tableCells$quadrant) # DMRs overlapping genes and promoters
entrez = lapply(entrez, function(x) as.character(unique(na.omit(x))))
entrez.dir = list(GenesPlusPromoters.pos = tableCells[which(tableCells$value>0 & tableCells$annotation != "Intergenic" & 
                                                                    tableCells$sig=="FWER < 0.05"),],
                  GenesPlusPromoters.neg = tableCells[which(tableCells$value<0 & tableCells$annotation != "Intergenic" & 
                                                                    tableCells$sig=="FWER < 0.05"),])
entrez.dir = lapply(entrez.dir, function(x) split(x,x$quadrant))
entrez.dir = lapply(entrez.dir, function(x) lapply(x, function(y) as.character(unique(na.omit(y$EntrezID)))))
entrez.dir = unlist(entrez.dir, recursive=F)

GeneUniverse = as.character(unique(na.omit(tableCells$EntrezID)))

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

# save
save(compareKegg, compareKegg.dir, compareBP, compareBP.dir, compareMF, compareMF.dir, compareCC, compareCC.dir, compareDO.dir,
     keggList, keggList.dir, goList_BP, goList_BP.dir, goList_MF, goList_MF.dir, goList_CC, goList_CC.dir, goList_DO, goList_DO.dir,
     keggListdf, keggList.dir.df, goListdf_BP, goListdf_BP.dir, goListdf_MF, goListdf_MF.dir, goListdf_CC, goListdf_CC.dir, goListdf_DO, goListdf_DO.dir,
     file="/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/DMR/DMR_KEGG_GO_DO_objects_infant_neuron_parsing_byQuadrant.rda")

# plot the GO results
pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/DMR/figures/DMR_KEGG_GO_DO_plots_infant_neuron_parsing_byQuadrant.pdf", height = 20, width = 24)
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
