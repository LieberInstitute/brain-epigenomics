library(PWMEnrich)
library(GenomicRanges)
library(PWMEnrich.Hsapiens.background)
library(BSgenome.Hsapiens.UCSC.hg19)
library(pheatmap)
library(RColorBrewer)
library(rafalib)
library(pvclust)
library(clusterProfiler)
library(org.Hs.eg.db)


load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/DMR/DMR_objects.rda")
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/PWMEnrich/DMR_PWMEnrich_objects.rda")
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/PWMEnrich/filtered_TF_list_expressionInfo.rda")
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/bumphunting/rda/limma_Neuron_CpGs_minCov_3_ageInfo_dmrs.Rdata")


## Cluster Interaction DMRs by TF enrichment to see how the groups materialize

## Extract matrices of pvalue enrichment for each DMR

groupint = groupReport(all_int)
groupint = as.data.frame(groupint)
groupint = cbind(groupint, padj = p.adjust(groupint$p.value, method = "fdr"))
groupint = groupint[which(groupint$target %in% names(PostnataltargettogeneID)),]
groupint = groupint[order(groupint$id),]

int = list()
for (i in 1:length(all_int$sequences)) {
  int[[i]] = sequenceReport(all_int, seq.id=i)
}
int = lapply(int, function(y) as.data.frame(y)[order(y$id),])
int = lapply(int, function(x) x[which(x$target %in% names(PostnataltargettogeneID)),])
int = do.call(cbind, lapply(int, function(z) p.adjust(z$p.value, method = "fdr")))
rownames(int) = groupint$target
colnames(int) = names(all_int$sequences)

dim(int)
thresh = int >= 0.9
table(rowSums(thresh==TRUE)==ncol(thresh))
#FALSE 
#1157 
table(colSums(thresh==TRUE)==nrow(thresh))
#FALSE  TRUE 
#1764   366 

lMatInt = -log10(int)


##  cluster by DMR (bootstrapped)

hc.int <- pvclust(t(lMatInt), method.hclust="ward", method.dist="euclidean")
hc_cut.int = lapply(hc.int, function(x) cutree(x, k= 10))

save(hc.int, file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/PWMEnrich/hc.int_bootstrapped.rda")

## cluster by TF (bootstrapped)

hct.int <- pvclust(lMatInt, method.hclust="ward", method.dist="euclidean")
hct_cut.int = lapply(hct.int, function(x) cutree(x, k= 10))

##  cluster by TF (not bootstrapped)

int.ward.TF = hclust(dist(lMatInt, method = "euclidean"), method="ward") 
int.ward.TF.cut = cutree(int.ward.TF, k= 6)

## cluster by DMR (not bootstrapped)

dmr.lMatInt = t(lMatInt)
dmrs$regionID = paste0(as.character(seqnames(dmrs)), ":",as.character(start(dmrs)), "-",as.character(end(dmrs)))
dmr.colors = dmrs$k6cluster[match(rownames(dmr.lMatInt), dmrs$regionID)]
names(dmr.colors) = dmrs$regionID[match(rownames(dmr.lMatInt), dmrs$regionID)]

int.ward.DMR = hclust(dist(dmr.lMatInt, method = "euclidean"), method="ward") 
int.ward.DMR.cut = cutree(int.ward.DMR, k= 6)


pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/PWMEnrich/figures/Interaction_DMR_TFenrichment_cluster.pdf", h = 12, w = 60)
palette(brewer.pal(12,"Paired"))
myplclust(int.ward.TF, lab.col=int.ward.TF.cut, xlab="", hang=0.05, cex=1.1, main = "Cluster by TF")
myplclust(int.ward.DMR, lab.col=int.ward.DMR.cut, xlab="", hang=0.05, cex=1.1, main = "Cluster by DMR")
palette(brewer.pal(8, 'Dark2'))
myplclust(int.ward.DMR, lab.col=dmr.colors, xlab="", hang=0.05, cex=1.1, main = "Cluster by DMR, Split by kmeans cluster")
dev.off()


# split TFs by cluster
sIndexes = splitit(hc2_cut)


## PCA

pca = prcomp(lMatInt)
pcaVars <- jaffelab::getPcaVars(pca)
names(pcaVars) <- paste0('PC', seq_len(length(pcaVars)))

pcat = prcomp(t(lMatInt))
pcatVars <- jaffelab::getPcaVars(pcat)
names(pcatVars) <- paste0('PC', seq_len(length(pcatVars)))


pdf('/dcl01/lieber/ajaffe/lab/brain-epigenomics/PWMEnrich/figures/pca_TF_Interaction_DMRs.pdf')
barplot(pcaVars[1:10], col = '#377EB8', ylab = 'Percent of Variance Explained',
        main = "PCA lMatInt")
  
plot(pca$x[, 1] ~ pca$x[, 2],
     ylab = paste0('PC1: ', pcaVars[1], '% of Var Explained'),
     xlab = paste0('PC2: ', pcaVars[2], '% of Var Explained'), pch = 19,
     main = "PCA lMatInt")
  
barplot(pcatVars[1:10], col = '#377EB8', ylab = 'Percent of Variance Explained',
        main = "PCA t(lMatInt)")
  
plot(pcat$x[, 1] ~ pcat$x[, 2],
     ylab = paste0('PC1: ', pcatVars[1], '% of Var Explained'),
     xlab = paste0('PC2: ', pcatVars[2], '% of Var Explained'), pch = 19,
     main = "PCA t(lMatInt)")
dev.off()


pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/PWMEnrich/figures/Interaction_DMR_TFenrichment_heatmap.pdf",h=16,w=8)
sampleDists <- dist(t(lMatInt))
sampleDistMatrix <- as.matrix(sampleDists)
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists,
         col=colors, main="Euclidean Distance Between\nInteraction DMRs by TF Enrichment")
sampleDists <- dist(lMatInt)
sampleDistMatrix <- as.matrix(sampleDists)
pheatmap(sampleDistMatrix,clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists,
         col=colors, main="Euclidean Distance Between\nTFs by Enrichment in Interaction DMRs")
dev.off()
# too large to be informative


## Pull out significant TFs per kmeans cluster and compare

intgroups = lapply(int.kmeans, groupReport)
intgroups = lapply(intgroups, as.data.frame)
intgroups = lapply(intgroups, function(x) x[which(x$target %in% names(PostnataltargettogeneID)),])
intgroups = Map(cbind, intgroups, padj = lapply(intgroups, function(x) p.adjust(x$p.value, method = "fdr")))
intgroups = lapply(intgroups, function(x) x[order(x$id),])

inttargets = lapply(intgroups, function(x) unique(as.character(x[which(x$padj<=0.01),"target"])))
inttargets = inttargets[elementNROWS(inttargets)>0]
venn.diagram(inttargets, "/dcl01/lieber/ajaffe/lab/brain-epigenomics/PWMEnrich/figures/Interaction_DMR_TF_venn_6kmeans.jpeg", 
             main="Interaction Kmeans Cluster TF Enrichment Overlap (FDR<0.01)",
             col = "transparent",
             fill = c("lightpink2","cornflowerblue", "olivedrab2", "khaki1", "plum2"),
             cat.col = c("palevioletred4", "darkblue", "olivedrab4", "lightgoldenrod4", "orchid4"),
             alpha = 0.50, fontfamily = "Arial", fontface = "bold", cat.fontfamily = "Arial", margin=0.2)

VennInt = calculate.overlap(inttargets)
entrez = list(all = VennInt$a31, "1:G-N+" = VennInt$a1, "2:G0N+" = VennInt$a2, "3:G0N-" = VennInt$a3, 
              "5:G+N-" = VennInt$a4, "6:G-N0" = VennInt$a5)
entrez = lapply(entrez, function(x) unlist(strsplit(x, "::", fixed=T)))
entrez = lapply(entrez, function(x) unique(na.omit(as.character(geneMap[which(geneMap$Symbol %in% x),"EntrezID"]))))


## Assess enriched terms limiting the gene universe to the terms associated with the master list of TFs

GeneUniverse = unique(na.omit(as.character(geneMap[which(geneMap$gencodeID %in% PostnataltargettogeneID), "EntrezID"])))
length(GeneUniverse) # 616

# Find enriched pathways and processes

keggList = lapply(entrez, function(x) enrichKEGG(x, organism="human", universe= GeneUniverse, minGSSize=5, 
                                                                      pAdjustMethod="BH", qvalueCutoff=1))
goList_MF = lapply(entrez, function(x) enrichGO(x, ont = "MF", OrgDb = org.Hs.eg.db, universe= GeneUniverse, 
                                                                     minGSSize=5, pAdjustMethod="BH", qvalueCutoff=1))
goList_BP = lapply(entrez, function(x) enrichGO(x, ont = "BP", OrgDb = org.Hs.eg.db,universe= GeneUniverse, 
                                                                     minGSSize=5, pAdjustMethod="BH", qvalueCutoff=1))
goList_CC = lapply(entrez, function(x) enrichGO(x, ont = "CC", OrgDb = org.Hs.eg.db, universe= GeneUniverse, 
                                                                     minGSSize=5, pAdjustMethod="BH", qvalueCutoff=1))
goList_DO = lapply(entrez, function(x) enrichDO(x, ont = "DO", universe= GeneUniverse, minGSSize=5, 
                                                                     pAdjustMethod="BH", qvalueCutoff=1, readable=TRUE))

# Compare the enriched terms

compareKegg = compareCluster(entrez, fun="enrichKEGG", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
compareBP = compareCluster(entrez, fun="enrichGO", ont = "BP", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
compareMF = compareCluster(entrez, fun="enrichGO",  ont = "MF", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
compareCC = compareCluster(entrez, fun="enrichGO",  ont = "CC", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
compareDO = compareCluster(entrez, fun="enrichDO",  ont = "DO", qvalueCutoff = 0.05, pvalueCutoff = 0.05)

save(keggList, goList_MF, goList_BP, goList_CC, goList_DO, compareBP, compareMF, compareCC, compareKegg, compareDO,
     file="/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/PWMEnrich/interaction_DMR_kmeans_GO.objects.PWMEnrich.rda")


pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/PWMEnrich/figures/KEGG_MF_CC_DO_interaction_DMR_kmeans_TFs.pdf", width=14,height=70)
plot(compareKegg, colorBy= "p.adjust",  showCategory = 700, title= "KEGG Pathway Enrichment")
plot(compareMF, colorBy= "p.adjust",  showCategory = 700, title= "MF Pathway Enrichment")
plot(compareCC, colorBy= "p.adjust",  showCategory = 700, title= "CC Pathway Enrichment")
plot(compareDO, colorBy= "p.adjust",  showCategory = 700, title= "Disease Ontology Enrichment")
plot(compareBP, colorBy= "p.adjust",  showCategory = 700, title= "Biological Process GO")
dev.off()


### Cluster Cell Type and overall age DMRs by TF enrichment to see how the groups materialize

## Extract matrices of pvalue enrichment for each DMR

ct.age = list()
for (i in 1:length(all_split[-grep("Interaction", names(all_split))])) {
  ct.age[[i]] = list(vector("list", length(all_split[-grep("Interaction", names(all_split))][[i]]$sequences)))
  for (j in 1:length(all_split[-grep("Interaction", names(all_split))][[i]]$sequences)) {
    ct.age[[i]][[j]] = sequenceReport(all_split[-grep("Interaction", names(all_split))][[i]], seq.id=j)
  }
}
names(ct.age) = names(all_split[-grep("Interaction", names(all_split))])
ct.age = lapply(ct.age, function(x) lapply(x, function(y) as.data.frame(y)[order(y$id),]))
ct.age = lapply(ct.age, function(t) lapply(t, function(x) x[which(x$target %in% names(targettogeneID)),]))
ct.age.Mat = lapply(ct.age, function(y) do.call(cbind, lapply(y, function(z) p.adjust(z$p.value, method = "fdr"))))
for (i in 1:length(ct.age.Mat)) { rownames(ct.age.Mat[[i]]) = ct.age[[1]][[1]]$target }
for (i in 1:length(ct.age.Mat)) { 
  colnames(ct.age.Mat[[i]]) = names(all_split[-grep("Interaction", names(all_split))][[i]]$sequences) }
ct.age.gr = lapply(all_split[-grep("Interaction", names(all_split))], function(x) GRanges(names(x$sequences)))
elementNROWS(ct.age.gr)
lMat = lapply(ct.age.Mat, function(x) -log10(x))
names(lMat) = c("Age (+Beta)", "Age (-Beta)", "Cell Type (+Beta)", "Cell Type (-Beta)")


## Filter TFs to those that are enriched

TF = list(allTF = lapply(all_split, groupReport), promTF = lapply(promoters_split, groupReport),
          intronTF = lapply(introns_split, groupReport), intergenicTF = lapply(intergenic_split, groupReport))
TF = lapply(TF, function(x) lapply(x, as.data.frame))
TF = lapply(TF, function(t) Map(cbind, t, padj = lapply(t, function(x) p.adjust(x$p.value, method = "fdr"))))
TF = lapply(TF, function(t) lapply(t, function(x) x[which(x$target %in% names(targettogeneID)),]))
TF = lapply(TF, function(t) lapply(t, function(x) x[order(x$id),]))
TF = lapply(TF, function(t) t[order(names(t))])
TF = unlist(TF, recursive = F)
TF = TF[-grep("Interaction", names(TF))]
sigTF = lapply(TF, function(x) x[which(x$padj<=0.01),])
sigTarg = list(Age.pos = unique(unlist(lapply(sigTF[grep("Age.pos", names(sigTF))], function(x) x$target))), 
               Age.neg = unique(unlist(lapply(sigTF[grep("Age.neg", names(sigTF))], function(x) x$target))), 
               CellType.pos = unique(unlist(lapply(sigTF[grep("CellType.pos", names(sigTF))], function(x) x$target))), 
               CellType.neg = unique(unlist(lapply(sigTF[grep("CellType.neg", names(sigTF))], function(x) x$target))))

lMat = list("Age (+Beta)" = lMat[["Age (+Beta)"]][which(rownames(lMat[["Age (+Beta)"]]) %in% sigTarg$Age.pos),],
            "Age (-Beta)" = lMat[["Age (-Beta)"]][which(rownames(lMat[["Age (-Beta)"]]) %in% sigTarg$Age.neg),],
            "Cell Type (+Beta)" = lMat[["Cell Type (+Beta)"]][which(rownames(lMat[["Cell Type (+Beta)"]]) %in% sigTarg$CellType.pos),],
            "Cell Type (-Beta)" = lMat[["Cell Type (-Beta)"]][which(rownames(lMat[["Cell Type (-Beta)"]]) %in% sigTarg$CellType.neg),])


##  cluster by DMR

hcT.ctage = lapply(lMat, function(x) hclust(dist(t(x)), method="ward"))
cutT.ctage = lapply(hcT.ctage, function(x) cutree(x, k= 6))

## Cluster by TF

hc.ctage = lapply(lMat, function(x) hclust(dist(x), method="ward"))
cut.ctage = lapply(hc.ctage, function(x) cutree(x, k= 8))


pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/PWMEnrich/figures/CellType_Age_DMR_TFenrichment_cluster.pdf", h = 10, w = 60)
palette(brewer.pal(12,"Paired"))
mapply(function(hc, hc_cut) myplclust(hc, lab.col=hc_cut,xlab="",hang=0.05,cex=1.1, main = "Cluster by TF"), 
       hc.ctage, cut.ctage, SIMPLIFY = F)
mapply(function(hc, hc_cut) myplclust(hc, lab.col=hc_cut,xlab="",hang=0.05,cex=1.1, main = "Cluster by DMR"), 
       hcT.ctage, cutT.ctage, SIMPLIFY = F)
dev.off()


## PCA

pca = lapply(lMat, prcomp)
pcaVars <- lapply(pca, function(x) jaffelab::getPcaVars(x))
for (i in 1:length(pcaVars)) { names(pcaVars[[i]]) <- paste0('PC', seq_len(length(pcaVars[[i]]))) }

pcat = lapply(lMat, function(x) prcomp(t(x)))
pcatVars <- lapply(pcat, function(x) jaffelab::getPcaVars(x))
for (i in 1:length(pcatVars)) { names(pcatVars[[i]]) <- paste0('PC', seq_len(length(pcatVars[[i]]))) }


pdf('/dcl01/lieber/ajaffe/lab/brain-epigenomics/PWMEnrich/figures/pca_TF_CellType_Age.pdf')
for (i in 1:length(pcaVars)) {
  barplot(pcaVars[[i]][1:10], col = '#377EB8', ylab = 'Percent of Variance Explained',
          main = names(pcaVars)[i])
  
  plot(pca[[i]]$x[, 1] ~ pca[[i]]$x[, 2],
       ylab = paste0('PC1: ', pcaVars[[i]][1], '% of Var Explained'),
       xlab = paste0('PC2: ', pcaVars[[i]][2], '% of Var Explained'), pch = 19,
       main = names(pca)[i])
  
  barplot(pcatVars[[i]][1:10], col = '#377EB8', ylab = 'Percent of Variance Explained',
          main = names(pcatVars)[i])
  
  plot(pcat[[i]]$x[, 1] ~ pcat[[i]]$x[, 2],
       ylab = paste0('PC1: ', pcatVars[[i]][1], '% of Var Explained'),
       xlab = paste0('PC2: ', pcatVars[[i]][2], '% of Var Explained'), pch = 19,
       main = names(pcat)[i])
}
dev.off()



## co-occuring motifs in Cell Type and Age DMRs


