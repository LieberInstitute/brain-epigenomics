library('bsseq')
library('bumphunter')
library('devtools')
library('getopt')
library('RColorBrewer')
library('ggplot2')

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


pdf(paste0('bumps_bsseqSmooth_', opt$subset, '_', opt$model,
    '_', opt$permutations, '_CpG_in_DMR.pdf'))
for(i in seq_len(100)) {
    matplot(meth[topInds[[i]], ], pch = 20, bg = factor(pd$Cell.Type), ylim = c(0, 1), ylab = 'DNAm Level', xlab = 'CpG in DMR', col = factor(pd$Cell.Type))
}
dev.off()



pdf(paste0('bumps_bsseqSmooth_', opt$subset, '_', opt$model,
    '_', opt$permutations, '_boxplot_by_age.pdf'))
for(i in seq_len(100)) {
    df <-  data.frame(
        Meth = as.vector(t(meth[topInds[[i]], ])),
        Age = rep(pd$Age, each = length(topInds[[i]])),
        cell = rep(pd$Cell.Type, each = length(topInds[[i]])),
        dmr = rep(seq_len(length(topInds[[i]])), each = nrow(pd))
    )
    g <- ggplot(data = df, aes(x = factor(Age), y = Meth, fill = cell)) + geom_boxplot() + theme_bw() + scale_fill_brewer(palette = 'Dark2') + xlab('Age') + ylab('DNAm Level')
    print(g)
}
dev.off()

## Reproducibility info
proc.time()
message(Sys.time())
options(width = 120)
session_info()
