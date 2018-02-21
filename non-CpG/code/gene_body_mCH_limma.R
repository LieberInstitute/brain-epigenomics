library(data.table)
library(ggplot2)
library(plyr)
library('bumphunter')
library(GenomicRanges)
library('devtools')
library('limma')
library('jaffelab')
library('shinycsv') # must download
library('RColorBrewer')
library(pheatmap)


load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/non-CpG/CH_object.rda")
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/non-CpG/CHneurons_object.rda")
load('/dcl01/lieber/ajaffe/lab/brain-epigenomics/bsseq/bsobj_by_chr/allChrs_postNatal_cleaned_nonCG_noHomogenate_highCov.Rdata')


## Get methylation
meth = getMeth(BSobj, type = "raw")
dim(meth)
gr = granges(BSobj)


#### Find mean mCH per gene and the regions flanking the gene

## split CH into gene bodies
geneMapgr = makeGRangesFromDataFrame(geneMap, keep.extra.columns = T)
mean(width(geneMapgr)) # 29389.18
ov = findOverlaps(gr, geneMapgr)
rIndexesIn = split(queryHits(ov), subjectHits(ov))

## get flanking regions 
leftShift = makeGRangesFromDataFrame(data.frame(seqnames = geneMap$Chr, start = (geneMap$Start-30000), end = geneMap$Start-1))
rightShift = makeGRangesFromDataFrame(data.frame(seqnames = geneMap$Chr, start = (geneMap$End+1), end = geneMap$End+30000))

## match up
ooLeft = findOverlaps(leftShift, gr)
ooRight = findOverlaps(rightShift, gr)
ooOut = rbind(as.matrix(ooLeft), as.matrix(ooRight))
rIndexesOut = split(ooOut[,"subjectHits"], ooOut[,"queryHits"])

## calculate mean mCH within each gene
inGene = sapply(rIndexesIn, function(ii) colMeans(t(t(meth[ii,]))))
length(inGene) # 40563
names(inGene) = geneMap[as.numeric(names(inGene)),"gencodeID"]
inGene = do.call("rbind", inGene)

## calculate flanking mCH
outGene= sapply(rIndexesOut, function(ii) colMeans(t(t(meth[ii,]))))
length(outGene) # 58907
names(outGene) = geneMap[as.numeric(names(outGene)),"gencodeID"]
outGene = do.call("rbind", outGene)
outGene = outGene[which(rownames(outGene) %in% rownames(inGene)),]
InGene = inGene[which(rownames(inGene) %in% rownames(outGene)),]


## differences in all genes

t.test(rowMeans(InGene), rowMeans(outGene))
#t = -23.239, df = 63707, p-value < 2.2e-16
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -0.004306001 -0.003636154
#sample estimates:
#  mean of x  mean of y 
#0.04424734 0.04821842 

df = data.frame(methDiff = -1*rowMeans(InGene - outGene))
pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/non-CpG/figures/nonCpG_meth_in.VS.out_GeneBodies.pdf",w=12)
ggplot(df, aes(methDiff)) + theme_classic() + geom_histogram() +
  labs(fill="") +
  ylab("") + 
  xlab("(Outside Gene) - (Inside Gene)") + xlim(-0.4,0.4) +
  ggtitle("mCH Within VS. Outside Gene Bodies") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20), legend.title=element_blank()) + theme(legend.position="bottom")
dev.off()


## Make PCA of mean gene body mCH 

pcs <- prcomp(t(inGene))
pcaVars <- getPcaVars(pcs)
names(pcaVars) <- paste0('PC', seq_len(length(pcaVars)))

pdf('/dcl01/lieber/ajaffe/lab/brain-epigenomics/non-CpG/figures/pca_mean_mCH_perGeneBody.pdf')
barplot(pcaVars[1:10], col = '#377EB8', ylab = 'Percent of Variance Explained')

plot(pcs$x[, 1] ~ pcs$x[, 2],
     ylab = paste0('PC1: ', pcaVars[1], '% of Var Explained'),
     xlab = paste0('PC2: ', pcaVars[2], '% of Var Explained'), pch = 19)

to_plot <- c('Cell.Type', 'Age', 'Age.Bin', 'PMI', 'Sex', 'Race', 'RIN', 'pH', 'Proportion.of.Neurons', 'yield.nuclei.mg.tissue', 'pg.DNA.nuclei.input', 'X260.280.DNA', 'Library.Technician', 'Flowcell', 'Reads', 'Percent.GreaterThan.Q30', 'avg.Cov', 'cov.sDev', 'Percent.Duplication', 'total_num_trimmed_reads', 'total_num_untrimmed_reads', 'alignment.efficiency')
to_plot <- which(colnames(pd) %in% to_plot)
names(to_plot) <- colnames(pd)[to_plot]

for(pc in 1:4) {
  mapply(function(id, name) {
    plot_twoway(y = pcs$x[, pc], x = pd[, id],
                yvar = paste0('PC', pc, ': ', pcaVars[pc], '% of Var Explained'),
                xvar = name, color = '#377EB8', pal = 'Set1')
    return(invisible(NULL))
  }, to_plot, names(to_plot))
}
dev.off()

pdf('/dcl01/lieber/ajaffe/lab/brain-epigenomics/non-CpG/figures/pca_mean_mCH_perGeneBody_PC1vsAge.pdf')
palette(brewer.pal(8,"Dark2"))
par(mar=c(5,6,2,2), cex.axis=2, cex.lab=2)
plot(pcs$x[, 1] ~ pd$Age, cex=2, pch = 21, bg=factor(pd$Cell.Type),
     ylab = paste0('PC1: ', pcaVars[1], '% of Var Explained'),  xlab = "Age")
legend("right", levels(factor(pd$Cell.Type)), pch = 15, col=1:2, cex=1.4)
dev.off()


## Heatmap of Euclidean Distance

sampleDists <- dist(t(inGene))
sampleDistMatrix <- as.matrix(sampleDists)
head(sampleDistMatrix)
colnames(sampleDistMatrix) = rownames(sampleDistMatrix) = paste(pd[match(colnames(sampleDistMatrix),pd$Data.ID),"Cell.Type"],
                                                                pd[match(colnames(sampleDistMatrix),pd$Data.ID),"Age.Bin"], sep = ":")

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/non-CpG/figures/heatmap_mean_mCH_perGeneBody.pdf")
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists,
         col=colors, main = "Mean Gene Body mCH - Euclidean Distance")
dev.off()

## fit linear models

models <- list(
  'cell' = with(pd, model.matrix(~ Cell.Type + Age)),
  'age' = with(pd, model.matrix(~ Age + Cell.Type)),
  'interaction' = with(pd, model.matrix(~ Age * Cell.Type)))

fits <- lapply(models, function(mod) { lmFit(inGene, design = mod) })
coefs <- c(2, 2, 4)
names(coefs) <- names(fits)

coef_interest <- mapply(function(f, coef) { f$coefficients[, coef] }, fits, coefs)
summary(coef_interest)
summary(abs(coef_interest))

ebList <- lapply(fits, ebayes)

## Explore neurons only with limma

inGene.NeuronsOnly = sapply(rIndexesIn, function(ii) colMeans(t(t(meth[ii,which(colnames(meth) %in% pd[which(pd$Cell.Type=="Neuron"),"Data.ID"])]))))
length(inGene.NeuronsOnly) # 40563
names(inGene.NeuronsOnly) = geneMap[as.numeric(names(inGene.NeuronsOnly)),"gencodeID"]
inGene.NeuronsOnly = do.call("rbind", inGene.NeuronsOnly)

model <- with(pd[which(pd$Cell.Type=="Neuron"),], model.matrix(~ Age))

fitNO <- lmFit(inGene.NeuronsOnly, design = model)

coef_interest <- fitNO$coefficients[,2]
summary(coef_interest)
summary(abs(coef_interest))

ebResultsNO <- ebayes(fitNO)

CHgene = data.frame(geneID = rownames(inGene), lods.CellType = ebList$cell$lods[,"Cell.TypeNeuron"], 
                    Tstat.CellType = ebList$cell$t[,"Cell.TypeNeuron"], 
                    pval.CellType = ebList$cell$p.value[,"Cell.TypeNeuron"],
                    padj.CellType = p.adjust(ebList$cell$p.value[,"Cell.TypeNeuron"], method = "fdr"),
                    lods.Age = ebList$age$lods[,"Age"], 
                    Tstat.Age = ebList$age$t[,"Age"], 
                    pval.Age = ebList$age$p.value[,"Age"], 
                    padj.Age = p.adjust(ebList$age$p.value[,"Age"], method = "fdr"),
                    lods.Interaction = ebList$interaction$lods[,"Age:Cell.TypeNeuron"], 
                    Tstat.Interaction = ebList$interaction$t[,"Age:Cell.TypeNeuron"], 
                    pval.Interaction = ebList$interaction$p.value[,"Age:Cell.TypeNeuron"], 
                    padj.Interaction = p.adjust(ebList$interaction$p.value[,"Age:Cell.TypeNeuron"], method = "fdr"),
                    lods.AgeNeuron = ebResultsNO$lods[,"Age"], 
                    Tstat.AgeNeuron = ebResultsNO$t[,"Age"], 
                    pval.AgeNeuron = ebResultsNO$p.value[,"Age"], 
                    padj.AgeNeuron = p.adjust(ebResultsNO$p.value[,"Age"], method = "fdr"))

save(fitNO, ebResultsNO, fits, ebList, inGene, inGene.NeuronsOnly, pd, CHgene,
     file = '/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/non-CpG/limma_exploration_mean_mCH_inGeneBody.Rdata')


## What is the distribution of mCH within each gene?

InCH = sapply(rIndexesIn, function(ii) meth[ii,])
length(InCH) # 40563
names(InCH) = geneMap[as.numeric(names(InCH)),"gencodeID"]
InCH = do.call("rbind", InCH)
df = as.data.frame(gr)
perc = list()
names = as.numeric(names(rIndexesIn))
for (i in 1:length(rIndexesIn)) {
  perc[[i]] = (df[rIndexesIn[[i]],"start"]-geneMap[names[i],"Start"]) / (geneMap[names[i],"End"]-geneMap[names[i],"Start"])
}
df = Map(data.frame, perc = perc, InCH)
df[which(lapply(df, ncol)!=33)]
df = do.call(rbind, df[which(lapply(df, ncol)==33)])


## Plot the position of the CH within the gene by the methylation by sample

corr = list()
for (i in 2:ncol(df)) {
corr[[i]] = cor.test(x = df[,1], y = df[,i])
}

x = do.call(rbind, Map(cbind, lapply(corr[2:33], function(y) 
  data.frame(pval = y$p.value, estimate = y$estimate)), id = as.list(colnames(df)[2:length(colnames(df))])))
x$Age.Bin = pd[match(x$id, pd$Data.ID),"Age.Bin"]
x$Cell.Type = pd[match(x$id, pd$Data.ID),"Cell.Type"]
write.csv(x, quote = F, file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/non-CpG/nonCG_geneBody_position_correlation.csv")