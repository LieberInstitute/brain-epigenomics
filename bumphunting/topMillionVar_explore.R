library('bsseq')
library('bumphunter')
library('devtools')
library('getopt')
library('doParallel')
library('limma')
library('jaffelab')
library('shinycsv')

## Specify parameters
spec <- matrix(c(
    'subset', 's', 1, 'character', 'Either: Homogenate or Neuron.',
#    'chromosome', 'c', 1, 'character', 'One of chr1 to chr22, chrX, chrY or chrM'
	'help' , 'h', 0, 'logical', 'Display help'
), byrow=TRUE, ncol=5)
opt <- getopt(spec)

## if help was asked for print a friendly message
## and exit with a non-zero error code
if (!is.null(opt$help)) {
	cat(getopt(spec, usage=TRUE))
	q(status=1)
}

## For testing
if(FALSE) {
    opt <- list('subset' = 'Neuron')
}

## Check inputs
stopifnot(opt$model %in% c('cell', 'age', 'interaction'))
stopifnot(opt$subset %in% c('all', 'Homogenate', 'Neuron'))
if(opt$subset == 'all') stop("Subset = 'all' is not supported.")

## Load data
load(paste0('BSobj_', opt$subset, '_topMillionVar.Rdata'))

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

pdf('pca_topMillionVar.pdf')
barplot(pcaVars[1:10], col = '#377EB8', ylab = 'Percent of Variance Explained')

plot(pcs$x[, 1] ~ pcs$x[, 2],
    ylab = paste0('PC1: ', pcaVars[1], '% of Var Explained'),
    xlab = paste0('PC2: ', pcaVars[2], '% of Var Explained'), pch = 19)

to_plot <- c('Cell.Type', 'Age', 'Age.Bin', 'PMI', 'Sex', 'Race', 'RIN', 'pH', 'Proportion.of.Neurons', 'yield.nuclei.mg.tissue', 'pg.DNA.nuclei.input', 'X260.280.DNA', 'Library.Technician', 'Flowcell', 'Data.Yield.', 'Reads', 'Percent.GreaterThan.Q30', 'avg.Cov', 'cov.sDev', 'Percent.Duplication', 'total_num_trimmed_reads', 'total_num_untrimmed_reads', 'alignment.efficiency')
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

save(fits, models, coef_interest, file = 'limma_exploration_topMillionVar.Rdata')

## Reproducibility info
proc.time()
message(Sys.time())
options(width = 120)
session_info()
