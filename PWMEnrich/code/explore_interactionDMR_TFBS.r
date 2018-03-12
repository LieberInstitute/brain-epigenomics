library(PWMEnrich)
library(GenomicRanges)
library(PWMEnrich.Hsapiens.background)
library(BSgenome.Hsapiens.UCSC.hg19)
library(pheatmap)
library(RColorBrewer)

load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/DMR/DMR_objects.rda")
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/PWMEnrich/DMR_PWMEnrich_objects.rda")
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/PWMEnrich/filtered_TF_list_expressionInfo.rda")


## Cluster Cell Type and overall age DMRs by TF enrichment to see how the groups materialize

ct.age = list()
for (i in 1:length(all_split[-grep("Interaction", names(all_split))])) {
  for (j in 1:length(all_split[-grep("Interaction", names(all_split))][[i]]$sequences)) {
    ct.age[[i]] = list(vector("list", length(all_split[-grep("Interaction", names(all_split))][[i]]$sequences)))
    ct.age[[i]][[j]] = sequenceReport(all_split[-grep("Interaction", names(all_split))][[i]], seq.id=j)
  }
}

ct.age = lapply(ct.age, function(x) lapply(x, function(y) as.data.frame(y)[order(y$id),]))
ct.age = lapply(ct.age, function(t) lapply(t, function(x) x[which(x$target %in% names(targettogeneID)),]))
ct.age.Mat = lapply(ct.age, function(y) lapply(y, do.call(cbind, lapply(y, p.adjust(y$p.value, method = "fdr")))))
rownames(ct.age.Mat) = ct.age[[1]][[1]]$target

for (i in 1:length(all_split[-grep("Interaction", names(all_split))])) { 
  colnames(ct.age.Mat[[i]]) = names(all_split[-grep("Interaction", names(all_split))][[i]]$sequences) }
ct.age.gr = lapply(all_split[-grep("Interaction", names(all_split))], function(x) GRanges(names(x$sequences)))
elementNROWS(ct.age.gr)

lMat = -log10(pvalMat)

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



pca = prcomp(lMat)
pcat = prcomp(t(lMat))
ggplot(data = d, aes_string(x = "PC2", y = "PC3", color = "group")) + 
  geom_point(size = 3) + 
  ylab(paste0("PC3: ", round(percentVar[3] * 100), "% variance")) + 
  xlab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))

pca= prcomp(t(zMat))
km = kmeans(t(zMat), 5)




## Cluster Interaction DMRs by TF enrichment to see how the groups materialize


int = list()
for (i in 1:length(all_int$sequences)) {
  all_int_seqRep[[i]] = sequenceReport(all_int, seq.id=i)
}
targ = all_int_seqRep[[1]]$target    
pvalMat = do.call(cbind, lapply(all_int_seqRep, function(x) x$p.value[match(targ, x$target)]))
rownames(pvalMat) = targ
colnames(pvalMat) = names(all_int$sequences)

lMat = -log10(pvalMat)

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



pca = prcomp(lMat)
pcat = prcomp(t(lMat))
ggplot(data = d, aes_string(x = "PC2", y = "PC3", color = "group")) + 
  geom_point(size = 3) + 
  ylab(paste0("PC3: ", round(percentVar[3] * 100), "% variance")) + 
  xlab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))

pca= prcomp(t(zMat))
km = kmeans(t(zMat), 5)