library('bsseq')
library('bumphunter')
library('devtools')
library('getopt')
library('doParallel')

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
    opt <- list('model' = 'cell', 'subset' = 'Neuron', cores = 1,
        permutations = 0)
}

## Check inputs
stopifnot(opt$model %in% c('cell', 'age', 'interaction'))
stopifnot(opt$subset %in% c('all', 'Homogenate', 'Neuron'))
if(opt$subset == 'all') stop("Subset = 'all' is not supported.")


if(!file.exists(paste0('BSobj_', opt$subset '.Rdata'))) {
    ## Load data
    load("/dcl01/lieber/WGBS/LIBD_Data/bsseqObj/bsseqObj_postNatal_cleaned_CpGonly.rda")

    ## extract pheno
    pd <- pData(BSobj)

    ## Keep only a given cell type if specified
    if(opt$subset != 'all') {
        BSobj <- BSobj[, which(pd$Cell.Type %in% c(opt$subset, 'Glia'))]
        pd <- pData(BSobj)
    }

    ## Filter low coverage bases
    cov <- getCoverage(BSobj, type = 'Cov')
    cov.ge1 <- cov >= 1
    cov.filt <- rowSums(cov.ge1) == ncol(cov)
    print("Number of bases filtered")
    table(cov.filt)
#    FALSE     TRUE
#  3655383 24562065
    BSobj <- BSobj[cov.filt, ]
    rm(cov, cov.ge1, cov.filt)
    
    save(BSobj, file = paste0('BSobj_', opt$subset, '.Rdata'))
    
    ## Get the top million most variable CpGs
    rvar <- genefilter::rowVars(getMeth(BSobj, type = 'raw'))
    rvar.filt <- which(order(rvar, decreasing = TRUE) %in% seq_len(1e6))
    BSobj_top <- BSobj[rvar.filt, ]
    save(BSobj_top, file = paste0('BSobj_', opt$subset, '_topMillionVar.Rdata'))
    rm(rvar, rvar.filt)
    
} else {
    load(paste0('BSobj_', opt$subset, '.Rdata'))
    pd <- pData(BSobj)
}


print(paste('Number of samples used:', nrow(pd)))

coef <- 2
if(opt$model == 'cell') {
    design <- with(pd, model.matrix(~ Cell.Type + Age))
    cut <- 0.1
} else if (opt$model == 'age') {
    design <- with(pd, model.matrix(~ Age + Cell.Type))
    cut <- 0.005
} else if (opt$model == 'interaction') {
    design <- with(pd, model.matrix(~ Age * Cell.Type))
    coef <- grep(':', colnames(design))
    cut <- 0.009
}

## Get chr coordinates and methylation values
gr <- granges(BSobj)
#gr <- gr[1:1e6]
## This codes leads to NaN values when not filtering
meth <- getMeth(BSobj, type = 'raw')
#meth <- getMeth(BSobj[1:1e6,], type = 'raw')
## Use this code when not filtering
#meth <- getCoverage(BSobj, type = 'M')  / (getCoverage(BSobj, type = 'Cov') + 1e-05)
stopifnot(all(is.finite(range(meth))))


## Free some memory
rm(BSobj)

## Define parallel environment
registerDoParallel(cores = opt$cores)

## Run bumphunter
bumps <- bumphunter(object = meth,
    design = design, chr = as.character(seqnames(gr)), pos = start(gr),
    maxGap = 1000, B = opt$permutations, coef = coef, cutoff = cut,
    nullMethod = 'bootstrap', smooth = TRUE, smoothFunction = locfitByCluster,
    useWeights = TRUE, bpSpan = 500)

## Save results
save(bumps, pd, file = paste0('bumps_', opt$subset, '_', opt$model, '_',
    opt$permutations, '.Rdata'))

## Reproducibility info
proc.time()
message(Sys.time())
options(width = 120)
session_info()
