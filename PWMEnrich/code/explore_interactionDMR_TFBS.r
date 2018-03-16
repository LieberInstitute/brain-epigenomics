library(PWMEnrich)
library(GenomicRanges)
library(PWMEnrich.Hsapiens.background)
library(BSgenome.Hsapiens.UCSC.hg19)
library(pheatmap)
library(RColorBrewer)
library(rafalib)
library(pvclust)

load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/DMR/DMR_objects.rda")
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/PWMEnrich/DMR_PWMEnrich_objects.rda")
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/PWMEnrich/filtered_TF_list_expressionInfo.rda")


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





## Cluster Interaction DMRs by TF enrichment to see how the groups materialize

load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/bumphunting/rda/limma_Neuron_CpGs_minCov_3_ageInfo_dmrs.Rdata")


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


##  cluster by DMR

hc.int <- pvclust(t(lMatInt), method.hclust="ward", method.dist="euclidean")
hc_cut.int = lapply(hc.int, function(x) cutree(x, k= 10))

## cluster by TF

hct.int <- pvclust(lMatInt, method.hclust="ward", method.dist="euclidean")
hct_cut.int = lapply(hct.int, function(x) cutree(x, k= 10))



pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/PWMEnrich/figures/Interaction_DMR_TFenrichment_cluster.pdf", h = 8, w = 16)
palette(brewer.pal(12,"Paired"))
mapply(function(hc, hc_cut) myplclust(hc, lab.col=hc_cut,xlab="",hang=0.05,cex=1.1), hc.ctage, hc_cut.ctage, SIMPLIFY = F)
dev.off()

# split TFs by cluster
sIndexes = splitit(hc2_cut)
nullgenes =  read.delim("/users/ajaffe/Lieber/Projects/450k/grant/ref_gene_hg19.txt", header=T,as.is=T)
genes = sig$annotation
goListByTF = lapply(sIndexes, function(x) dogo(genes[x],nullgenes[,2])[,-8])


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




pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/PWMEnrich/figures/Interaction_DMR_TFenrichment_heatmap.pdf",h=16,w=8)
sampleDists <- dist(t(lMat))
sampleDistMatrix <- as.matrix(sampleDists)
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists,
         col=colors, main="Euclidean Distance Between\nInteraction DMRs by TF Enrichment")
sampleDists <- dist(lMat)
sampleDistMatrix <- as.matrix(sampleDists)
pheatmap(sampleDistMatrix,clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists,
         col=colors, main="Euclidean Distance Between\nTFs by Enrichment in Interaction DMRs")
dev.off()
# too large to be informative




pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/PWMEnrich/figures/Interaction_DMR_TFenrichment_cluster.pdf", h = 5, w = 14)
##  cluster by TF
hc = hclust(dist(t(lMat)))
hc_cut = cutree(hc, h= 40)
palette(brewer.pal(12,"Paired"))
myplclust(hc, lab.col=hc_cut,xlab="",hang=0.05,cex=1.1)
# cluster by DMR
hc2 = hclust(dist(lMat))
hc2_cut = cutree(hc2, k= 20)
myplclust(hc2, lab.col=hc2_cut,xlab="",hang=0)
dev.off()

### go by cluster
sIndexes = splitit(hc2_cut)
nullgenes =  read.delim("/users/ajaffe/Lieber/Projects/450k/grant/ref_gene_hg19.txt", header=T,as.is=T)
genes = sig$annotation
goListByTF = lapply(sIndexes, function(x) dogo(genes[x],nullgenes[,2])[,-8])

