library('bsseq')
library('bumphunter')
library('devtools')
library('getopt')
library('RColorBrewer')

## Specify parameters
spec <- matrix(c(
	'model', 'm', 1, 'character', 'Either: cell, age or interaction',
    'subset', 's', 1, 'character', 'Either: Homogenate or Neuron.',
    'cores', 't', 1, 'integer', 'Number of cores to use',
#    'chromosome', 'c', 1, 'character', 'One of chr1 to chr22, chrX, chrY or chrM'
    'permutations', 'p', 1, 'integer', 'Number of permutations to run',
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
    opt <- list('model' = 'cell', 'subset' = 'Neuron',
        permutations = 0)
}

## Check inputs
inputFile <- paste0('bumps_bsseqSmooth_', opt$subset, '_', opt$model,
    '_', opt$permutations, '.Rdata')
stopifnot(file.exists(inputFile))

load(inputFile)

load("../processed_beta_values_plusMap.rda")

## tale of top loci
tab = bumps$table[1:100,]
## top CpGs
topInds = mapply(function(s,e) s:e, tab$indexStart, tab$indexEnd)

meanMeth = lapply(topInds, function(ii) colMeans(t(t(meth[ii,]))))
meanMeth = do.call("rbind", meanMeth)

pdfFile <- paste0('bumps_bsseqSmooth_', opt$subset, '_', opt$model,
    '_', opt$permutations, '.pdf')

pdf(pdfFile)
palette(brewer.pal(8,"Dark2"))
for(i in 1:nrow(meanMeth)) {
	plot(meanMeth[i,] ~ pd$Age, pch = 21,
		bg = factor(pd$Cell.Type),ylim = c(0,1),
		ylab = "DNAm Level", xlab = 'Age')
	abline(lm(meanMeth[i,] ~ pd$Age, subset = pd$Cell.Type =="Glia"),
		col = 1)
	abline(lm(meanMeth[i,] ~ pd$Age, subset = pd$Cell.Type =="Neuron"),
		col = 2)
    legend("right", levels(factor(pd$Cell.Type)), pch = 15, col=1:2, cex=1.4)
}
dev.off()

## Reproducibility info
proc.time()
message(Sys.time())
options(width = 120)
session_info()
