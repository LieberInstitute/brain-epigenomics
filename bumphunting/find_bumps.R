library('bsseq')
library('bumphunter')
library('devtools')
library('getopt')
library('doparallel')

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
    opt <- list('model' = 'cell', 'subset' = 'Homogenate', cores = 1, permutations = 5)
}

## Check inputs
stopifnot(opt$model %in% c('cell', 'age', 'interaction'))
stopifnot(opt$subset %in% c('all', 'Homogenate', 'Neuron'))

## Load data
load("/dcl01/lieber/WGBS/LIBD_Data/bsseqObj/bsseqObj_postNatal_cleaned_CpGonly.rda")

## extract pheno
pd <- pData(BSobj)

## Keep only a given cell type if specified
if(opt$subset != 'all') {
    BSobj <- BSobj[, which(pd$Cell.Type %in% c(opt$subset, 'Glia'))]
    pd <- pData(BSobj)
}

print(paste('Number of samples used:', nrow(pd)))

coef <- 2
if(opt$model == 'cell') {
    design <- with(pd, model.matrix(~ Cell.Type + Age))
    cut <- 0.1
} else if (opt$model == 'age') {
    design <- with(pd, model.matrix(~ Age + Cell.Type))
    cut <- 0.01
} else if (opt$model == 'interaction') {
    design <- with(pd, model.matrix(~ Age * Cell.Type))
    coef <- grep(':', colnames(design))
    cut <- 0.1 * 0.01
}

## Get chr coordinates and methylation values
gr <- granges(BSobj)
#gr <- gr[1:1000]
meth <- getMeth(BSobj, type = 'raw')
#meth <- getMeth(BSobj[1:1000,], type = 'raw')

## Free some memory
rm(BSobj)

## Define parallel environment
registerDoParallel(cores = opt$cores)

## Run bumphunter
bumps <- bumphunter(object = meth,
    design = design, chr = as.character(seqnames(gr)), pos = start(gr),
    maxGap = 300, B = opt$permutations, coef = coef, cutoff = cut)

## Save results
save(bumps, file = paste0('bumps_', opt$subset, '_', opt$model, '_',
    opt$permutations, '.Rdata'))

## Reproducibility info
proc.time()
message(Sys.time())
options(width = 120)
session_info()
