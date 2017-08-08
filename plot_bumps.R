###
library('bsseq')
library('bumphunter')
library('doParallel')

# load data
# load("bumphunting/BSobj_Neuron.Rdata")

# load bumps
load("bumphunting/bumps_Neuron_interaction_0.Rdata")

# # prep data
# gr <- granges(BSobj)
# meth <- getMeth(BSobj, type = 'raw')
# cl = clusterMaker( chr = as.character(seqnames(gr)), 
	# pos = start(gr),   maxGap = 1000)
# gr$cluster = cl
# save(meth,gr, pd, 	compress=TRUE,
	# file="processed_beta_values_plusMap.rda")
load("processed_beta_values_plusMap.rda")

## tale of top loci
tab = bumps$table[1:100,]
## top CpGs
topInds = mapply(function(s,e) s:e, tab$indexStart, tab$indexEnd)

meanMeth = sapply(topInds, function(ii) colMeans(t(t(meth[ii,]))))
meanMeth = do.call("rbind", meanMeth)

for(i in 1:nrow(meanMeth)) {
	plot(meanMeth[i,] ~ pd$Age, pch = 21,
		bg = factor(pd$Cell.Type),ylim = c(0,1),
		ylab = "DNAm Level")
	abline(lm(meanMeth[i,] ~ pd$Age, subset = pd$Cell.Type =="Glia"),
		col = 1)
	abline(lm(meanMeth[i,] ~ pd$Age, subset = pd$Cell.Type =="Neuron"),
		col = 2)
}