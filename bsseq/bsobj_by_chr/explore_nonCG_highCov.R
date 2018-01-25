library('bsseq')
library('devtools')
library('limma')
library('jaffelab')
library('shinycsv')
library('RColorBrewer')

## Load data
load('allChrs_postNatal_cleaned_nonCG_noHomogenate_highCov.Rdata')

## Size of the data
dim(BSobj)

## extract pheno
pd <- pData(BSobj)

## Fix pheno data
pd$Race[pd$Race == "CAUC "] <- 'CAUC'
pd$Sex[pd$Sex == " M"] <- 'M'
pd$RIN <- as.numeric(gsub(' ', '', pd$RIN))
pd$pg.DNA.nuclei.input <- as.numeric(pd$pg.DNA.nuclei.input)
pd$Reads <- as.numeric(pd$Reads)
pd$Percent.GreaterThan.Q30 <- as.numeric(pd$Percent.GreaterThan.Q30)

## Get methylation info
meth <- getMeth(BSobj, type = 'raw')
rm(BSobj)
meth.g0 <- meth > 0
meth.filt <- rowSums(meth.g0) >= 5
meth.tab <- table(meth.filt)
meth.tab
round(meth.tab / sum(meth.tab) * 100, 2)
meth <- meth[meth.filt, ]
rm(meth.g0, meth.filt)

pcs <- prcomp(t(meth))
pcaVars <- getPcaVars(pcs)
names(pcaVars) <- paste0('PC', seq_len(length(pcaVars)))

pdf('pdf/pca_nonCG_highCov.pdf')
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

pdf('pdf/pca_nonCG_highCov_PC1vsAge.pdf')
palette(brewer.pal(8,"Dark2"))
par(mar=c(5,6,2,2), cex.axis=2, cex.lab=2)
plot(pcs$x[, 1] ~ pd$Age, cex=2, pch = 21, bg=factor(pd$Cell.Type),
    ylab = paste0('PC1: ', pcaVars[1], '% of Var Explained'),  xlab = "Age")
legend("right", levels(factor(pd$Cell.Type)), pch = 15, col=1:2, cex=1.4)
dev.off()

## Explore with limma
models <- list(
    'cell' = with(pd, model.matrix(~ Cell.Type + Age)),
    'age' = with(pd, model.matrix(~ Age + Cell.Type)),
    'interaction' = with(pd, model.matrix(~ Age * Cell.Type)))

fits <- lapply(models, function(mod) {
    lmFit(meth, design = mod)
})
coefs <- c(2, 2, 4)
names(coefs) <- names(fits)

coef_interest <- mapply(function(f, coef) {
    f$coefficients[, coef]
}, fits, coefs)
summary(coef_interest)
summary(abs(coef_interest))

ebList <- lapply(fits, ebayes)

save(fits, models, coef_interest, ebList, file = 'limma_exploration_nonCG_highCov.Rdata')
rm(fits, models, coef_interest)



### examples to plot ###
pdf("pdf/cellType_DMPs_nonCG_highCov.pdf")
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
		
pdf("pdf/overallAge_DMPs_nonCG_highCov.pdf")
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

pdf("pdf/interactionCelltypeAge_DMPs_nonCG_highCov.pdf")
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

xx <- mapply(function(coef, fit, model) {
    print(paste('Info for model:', model))
    med <- table(fit$p.value[, coef] < 1e-04)
    print("Cutoff: 1e-04")
    print(med)
    print(round(med / sum(med) * 100, 2))
    high <- table(fit$p.value[, coef] < 1e-06)
    print("Cutoff: 1e-06")
    print(high)
    print(round(high / sum(high) * 100, 2))
    return(invisible(NULL))
}, coefs, ebList, names(ebList))

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
