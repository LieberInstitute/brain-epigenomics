##
library(bsseq)
library(pheatmap)
library(genefilter)

## load dmrs
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/bumphunting/bumps_bsseqSmooth_Neuron_interaction_250_perm.Rdata")
sigInt = bumps$table[bumps$table$fwer < 0.05,]
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/bumphunting/bumps_bsseqSmooth_Neuron_age_250_perm.Rdata")
sigAge = bumps$table[bumps$table$fwer < 0.05,]

## load data
load('/dcl01/lieber/ajaffe/lab/brain-epigenomics/bumphunting/BSobj_bsseqSmooth_Neuron_minCov_3.Rdata')
pd = pData(BSobj)
## get mean meth per DMR
meth <- getMeth(BSobj, type = 'raw')

### interaction
topIndsInt = mapply(function(s,e) s:e, sigInt$indexStart, sigInt$indexEnd)
meanMethInt = sapply(topIndsInt, function(ii) colMeans(t(t(meth[ii,]))))
meanMethInt = do.call("rbind", meanMethInt)

sampleDistsInt <- dist(t(meanMethInt))
sampleDistMatrixInt <- as.matrix(sampleDistsInt)
colnames(sampleDistMatrixInt) = rownames(sampleDistMatrixInt) = paste(
		pd$Cell.Type, pd$Age.Bin, pd$Working.Num, sep = ":")

### age
topIndsAge = mapply(function(s,e) s:e, sigAge$indexStart, sigAge$indexEnd)
meanMethAge = sapply(topIndsAge, function(ii) colMeans(t(t(meth[ii,]))))
meanMethAge = do.call("rbind", meanMethAge)

sampleDistsAge <- dist(t(meanMethAge))
sampleDistMatrixAge <- as.matrix(sampleDistsAge)
colnames(sampleDistMatrixAge) = rownames(sampleDistMatrixAge) = paste(
		pd$Cell.Type, pd$Age.Bin, pd$Working.Num, sep = ":")
		
pdf("heatmap_euclidean_dist.pdf")
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrixInt,clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists,
         col=colors, main = "2178 Interaction DMRs - Euclidean Dist")
pheatmap(sampleDistMatrixAge,clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists,
         col=colors, main = "129 Age DMRs - Euclidean Dist")
dev.off()

## which DMRs show the association?
pd$Neonate = ifelse(pd$Age < 1, "Neonate", "Older") 
gIndexes = list(YoungNeuronVsGlia = which(pd$Age < 1 | pd$Cell.Type == "Glia"),
	YoungVsOldNeuron = which(pd$Cell.Type == "Neuron"))

tt_YoungNeuronVsGlia_Int = rowttests(meanMethInt[,gIndexes$YoungNeuronVsGlia],
	factor(pd$Cell.Type[gIndexes$YoungNeuronVsGlia])) # neg is glia effect
tt_YoungVsOldNeuron_Int = rowttests(meanMethInt[,gIndexes$YoungVsOldNeuron],
	factor(pd$Neonate[gIndexes$YoungVsOldNeuron])) # neg is neonate effect
plot(tt_YoungNeuronVsGlia_Int$statistic, tt_YoungVsOldNeuron_Int$statistic)
plot(-log10(tt_YoungNeuronVsGlia_Int$p.value), -log10(tt_YoungVsOldNeuron_Int$p.value))

table(tt_YoungNeuronVsGlia_Int$p.value < 0.01, tt_YoungVsOldNeuron_Int$p.value < 0.01,
	dnn = c("YoungNeuronVsGlia", "YoungVsOldNeuron"))
	

## AP : do GO on each quadrant of that table
tableCells = split(sigInt, paste0(tt_YoungNeuronVsGlia_Int$p.value < 0.01, 
	"-",tt_YoungVsOldNeuron_Int$p.value < 0.01))
tableCells = tableCells[names(tableCells)!="NA-NA"] # drop the missing one

## anntotate DMRs w/in each list
## run clusterProfiler within each element