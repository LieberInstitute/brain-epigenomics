library(GenomicRanges)
library(pheatmap)
library(RColorBrewer)


load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/processed_beta_values_plusMap.rda")
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/DMR/DMR_objects.rda")


## Find mean of all regions identified in the "Interaction" DMR analsysis

Int_gr = makeGRangesFromDataFrame(DMR$Interaction, keep.extra.columns = T)
oo = findOverlaps(Int_gr, gr)
rIndexes = split(subjectHits(oo), queryHits(oo))

meanCov_inDMR = sapply(rIndexes, function(ii) colMeans(t(t(meth[ii,]))))
names(meanCov_inDMR) = Int_gr$regionID
meanCov_inDMR = do.call("rbind", meanCov_inDMR)
match(colnames(meanCov_inDMR),pd$Data.ID)
meansig = meanCov_inDMR[which(rownames(meanCov_inDMR) %in% DMR$Interaction[which(DMR$Interaction$sig == "FWER < 0.05"),"regionID"]),]


## Calculate Euclidean distance between samples

sampleDists <- dist(t(meanCov_inDMR))
sampleDistMatrix <- as.matrix(sampleDists)
head(sampleDistMatrix)
sigDists <- dist(t(meansig))
sigDistMatrix <- as.matrix(sigDists)
head(sigDistMatrix)


## Plot heatmap of distances between samples

colors <- colorRampPalette(rev(brewer.pal(9, "Blues")) )(255)
pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/DMR/CT_Age_Interaction/figures/Euclidean_Distance_btwnSamples_DMR_byInt.pdf")
pheatmap(sampleDistMatrix,clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists,
         col=colors, main="Euclidean Distance Between Samples")
colnames(sampleDistMatrix) = rownames(sampleDistMatrix) = paste(pd$Cell.Type, pd$Age.Bin, pd$Working.Num, sep = ":")
pheatmap(sampleDistMatrix,clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists,
         col=colors, main="Euclidean Distance Between Samples")
pheatmap(sigDistMatrix,clustering_distance_rows= sigDists, clustering_distance_cols= sigDists,
         col=colors, main="Euclidean Distance Between Samples: FWER < 0.05")
colnames(sigDistMatrix) = rownames(sigDistMatrix) = paste(pd$Cell.Type, pd$Age.Bin, pd$Working.Num, sep = ":")
pheatmap(sigDistMatrix,clustering_distance_rows= sigDists, clustering_distance_cols= sigDists,
         col=colors, main="Euclidean Distance Between Samples: FWER < 0.05")
dev.off()


## Split mean DMR table by group

meanbygroup = list(Glia = meansig[,which(colnames(meansig) %in% pd[which(pd$Cell.Type=="Glia"),"Data.ID"])],
				   Neonatal.Neurons = meansig[,which(colnames(meansig) %in% pd[which(pd$Cell.Type=="Neuron" & pd$Age.Bin=="Neonate"),"Data.ID"])],
				   Other.Neurons = meansig[,which(colnames(meansig) %in% pd[which(pd$Cell.Type=="Neuron" & pd$Age.Bin!="Neonate"),"Data.ID"])])
meanbygroup = lapply(meanbygroup, function(x) split(x, seq(nrow(x))))



## Calculate T statistics for each DMR between the groups

tstats = list(GliavsNeonate = mapply(function(x, y) t.test(x, y), meanbygroup$Glia, meanbygroup$Neonatal.Neurons),
			  NeuronsvsNeonate = mapply(function(x, y) t.test(x, y), meanbygroup$Other.Neurons, meanbygroup$Neonatal.Neurons))









