library('bsseq')
library('bumphunter')
library('devtools')
library('getopt')
library('doParallel')
library('limma')
library('jaffelab')
library('shinycsv')

## Load data
load(paste0('bumphunting/BSobj_Neuron_topMillionVar.Rdata'))

## extract pheno
pd <- pData(BSobj_top)

## Fix pheno data
pd$Race[pd$Race == "CAUC "] <- 'CAUC'
pd$Sex[pd$Sex == " M"] <- 'M'
pd$RIN <- as.numeric(gsub(' ', '', pd$RIN))
pd$pg.DNA.nuclei.input <- as.numeric(pd$pg.DNA.nuclei.input)
pd$Data.Yield. <- as.numeric(pd$Data.Yield.)
pd$Reads <- as.numeric(pd$Reads)
pd$Percent.GreaterThan.Q30 <- as.numeric(pd$Percent.GreaterThan.Q30)

## Get methylation info
meth <- getMeth(BSobj_top, type = 'raw')

pcs <- prcomp(t(meth))
pcaVars <- getPcaVars(pcs)
names(pcaVars) <- paste0('PC', seq_len(length(pcaVars)))

pdf('plots/pca_topMillionVar_PC1vsAge.pdf')
palette(brewer.pal(8,"Dark2"))
par(mar=c(5,6,2,2), cex.axis=2,cex.lab=2)
plot(pcs$x[, 1] ~ pd$Age, cex=2, pch = 21, bg=factor(pd$Cell.Type),
    ylab = paste0('PC1: ', pcaVars[1], '% of Var Explained'),  xlab = "Age")
legend("right", levels(factor(pd$Cell.Type)), pch = 15, col=1:2,cex=1.4)
dev.off()

## load limma output
load('bumphunting/limma_exploration_topMillionVar.Rdata')

ebList = lapply(fits, ebayes)
plot(ebList[[1]]$t[,2:3])
plot(ebList[[3]]$t[,c(2,4)])

### examples to plot ###
pdf("plots/cellType_DMPs.pdf")
par(mar=c(5,6,2,2),cex.axis=2,cex.lab=2)
sigIndex_cell = order(ebList[[1]]$p[,2])[1:50]
for(i in seq(along=sigIndex_cell)) {
	ii = sigIndex_cell[i]
	plot(meth[ii,] ~ pd$Age, pch = 21, cex=1.5,
		bg=factor(pd$Cell.Type), ylab="DNAm Level",xlab="Age")
	legend("right", paste0("p=",signif(ebList[[1]]$p[ii,2], 3)),cex=1.5)
	legend("left", levels(factor(pd$Cell.Type)), 
		pch = 15, col = 1:2,cex=1.5)
}
dev.off()
		
pdf("plots/overallAge_DMPs.pdf")
par(mar=c(5,6,2,2),cex.axis=2,cex.lab=2)
sigIndex_age = order(ebList[[2]]$p[,2])[1:50]
for(i in seq(along=sigIndex_age)) {
	ii = sigIndex_age[i]
	plot(meth[ii,] ~ pd$Age, pch = 21, cex=1.5,
		bg=factor(pd$Cell.Type), ylab="DNAm Level",xlab="Age")
	ll = ifelse(ebList[[2]]$t[ii,2] > 0, "bottomright", "topright")
	legend(ll, paste0("p=",signif(ebList[[2]]$p[ii,2], 3)),cex=1.5)
	ll2 = ifelse(ebList[[2]]$t[ii,2] > 0, "topleft", "bottomleft")
	legend(ll2, levels(factor(pd$Cell.Type)), 
		pch = 15, col = 1:2,cex=1.5)
}
dev.off()		

pdf("plots/interactionCelltypeAge_DMPs.pdf")
par(mar=c(5,6,2,2),cex.axis=2,cex.lab=2)
sigIndex_int = order(ebList[[3]]$p[,4])[1:50]
for(i in seq(along=sigIndex_int)) {
	ii = sigIndex_int[i]
	plot(meth[ii,] ~ pd$Age, pch = 21, cex=1.5,
		bg=factor(pd$Cell.Type), ylab="DNAm Level",xlab="Age")
	legend("right", paste0("p=",signif(ebList[[3]]$p[ii,4], 3)),cex=1.5)
	legend("left", levels(factor(pd$Cell.Type)), pch = 15, col = 1:2,cex=1.5)
}
dev.off()
		
		## Reproducibility info
proc.time()
message(Sys.time())
options(width = 120)
session_info()
