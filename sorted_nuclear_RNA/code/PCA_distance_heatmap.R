library(DESeq2)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)


load("/dcl01/lieber/ajaffe/CellSorting/RNAseq_pipeline/rawCounts_CellSorting_July5_n12.rda")

# process counts for PCA
metrics[grep("12", metrics$SampleID),"Age"] = "Neonate"
metrics[grep("13", metrics$SampleID),"Age"] = "Teenager"
metrics[grep("14", metrics$SampleID),"Age"] = "Toddler"
metrics$CellType = ifelse(metrics$NeuN=="NeuN_Minus", "Glia", "Neuron")

dds = DESeqDataSetFromMatrix(countData = geneCounts, colData = metrics, design = ~ Prep + Age + CellType)
rlog.dds  = rlog(dds)
rlog = rlog(geneCounts)
save(rlog,rlog.dds,file="/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/sorted_nuclear_RNA/rlog_transformed_dds_sorted_nuclear_RNA.rda")
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/sorted_nuclear_RNA/rlog_transformed_dds_sorted_nuclear_RNA.rda")

### PCA Functions ###

plotPCA <- function (x, intgroup = "condition", ntop = 500, returnData = FALSE) 
{
  rv <- rowVars(assay(x))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, 
                                                     length(rv)))]
  pca <- prcomp(t(assay(x)[select, ]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  if (!all(intgroup %in% names(colData(x)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  intgroup.df <- as.data.frame(colData(x)[, intgroup, drop = FALSE])
  group <- factor(apply(intgroup.df, 1, paste, collapse = " : "))
  d <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], group = group, 
                  intgroup.df, names = colnames(x))
  if (returnData) {
    attr(d, "percentVar") <- percentVar[1:2]
    return(d)
  }
  ggplot(data = d, aes_string(x = "PC1", y = "PC2", color = "group")) + 
    geom_point(size = 3) + 
    xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) +
    ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) +
    theme(title = element_text(size = 20)) +
    theme(text = element_text(size = 20))
}

plotPCA2 <- function (x, intgroup = "condition", ntop = 500, returnData = FALSE) 
{
  rv <- rowVars(assay(x))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, 
                                                     length(rv)))]
  pca <- prcomp(t(assay(x)[select, ]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  if (!all(intgroup %in% names(colData(x)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  intgroup.df <- as.data.frame(colData(x)[, intgroup, drop = FALSE])
  group <- factor(apply(intgroup.df, 1, paste, collapse = " : "))
  d <- data.frame(PC2 = pca$x[, 2], PC3 = pca$x[, 3], group = group, 
                  intgroup.df, names = colnames(x))
  if (returnData) {
    attr(d, "percentVar") <- percentVar[1:2]
    return(d)
  }
  ggplot(data = d, aes_string(x = "PC2", y = "PC3", color = "group")) + 
    geom_point(size = 3) + 
    ylab(paste0("PC3: ", round(percentVar[3] * 100), "% variance")) + 
    xlab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) +
    theme(title = element_text(size = 20)) +
    theme(text = element_text(size = 20))
}

plotPCA3 <- function (x, intgroup = "condition", ntop = 500, returnData = FALSE) 
{
  rv <- rowVars(assay(x))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, 
                                                     length(rv)))]
  pca <- prcomp(t(assay(x)[select, ]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  if (!all(intgroup %in% names(colData(x)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  intgroup.df <- as.data.frame(colData(x)[, intgroup, drop = FALSE])
  group <- factor(apply(intgroup.df, 1, paste, collapse = " : "))
  d <- data.frame(PC3 = pca$x[, 3], PC4 = pca$x[, 4], group = group, 
                  intgroup.df, names = colnames(x))
  if (returnData) {
    attr(d, "percentVar") <- percentVar[1:2]
    return(d)
  }
  ggplot(data = d, aes_string(x = "PC3", y = "PC4", color = "group")) + 
    geom_point(size = 3) + 
    xlab(paste0("PC3: ", round(percentVar[3] * 100), "% variance")) + 
    ylab(paste0("PC4: ", round(percentVar[4] * 100), "% variance")) +
    theme(title = element_text(size = 20)) +
    theme(text = element_text(size = 20))
}

plotPCA4 <- function (x, intgroup = "condition", ntop = 500, returnData = FALSE) 
{
  rv <- rowVars(assay(x))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, 
                                                     length(rv)))]
  pca <- prcomp(t(assay(x)[select, ]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  if (!all(intgroup %in% names(colData(x)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  intgroup.df <- as.data.frame(colData(x)[, intgroup, drop = FALSE])
  group <- factor(apply(intgroup.df, 1, paste, collapse = " : "))
  d <- data.frame(PC4 = pca$x[, 4], PC5 = pca$x[, 5], group = group, 
                  intgroup.df, names = colnames(x))
  if (returnData) {
    attr(d, "percentVar") <- percentVar[1:2]
    return(d)
  }
  ggplot(data = d, aes_string(x = "PC4", y = "PC5", color = "group")) + 
    geom_point(size = 3) + 
    xlab(paste0("PC4: ", round(percentVar[4] * 100), "% variance")) +
    ylab(paste0("PC5: ", round(percentVar[5] * 100), "% variance")) +
    theme(title = element_text(size = 20)) +
    theme(text = element_text(size = 20))
}

### Run PCA
require(genefilter)

plotPCA(rlog.dds, intgroup = "CellType")
plotPCA(rlog.dds, intgroup = "Prep")
plotPCA2(rlog.dds, intgroup = "Age")
plotPCA3(rlog.dds, intgroup = "Age")


# Euclidean distance between samples heatmap
sampleDists <- dist(t(rlog))
sampleDistMatrix <- as.matrix(sampleDists)
head(sampleDistMatrix)
colnames(metrics)
match(colnames(geneCounts),metrics$SAMPLE_ID)
colnames(sampleDistMatrix) = rownames(sampleDistMatrix) = paste0(metrics$CellType,":",metrics$Age,":",metrics$Prep)
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists,
         col=colors, main="Euclidean Distance Between Samples")

sampleDists <- dist(t(rlog.down))
sampleDistMatrix <- as.matrix(sampleDists)
colnames(sampleDistMatrix) = rownames(sampleDistMatrix) = pd$WorkingID
pheatmap(sampleDistMatrix,clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists,
         col=colors, main="Euclidean Distance Between Samples (Downsampled)")



