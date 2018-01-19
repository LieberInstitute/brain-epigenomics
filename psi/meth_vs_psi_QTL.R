## Partially based on /users/ajaffe/Lieber/Projects/450k/scalpVsDura/eQTL.R
library('SGSeq')
library('bsseq')
library('MatrixEQTL')
library('getopt')

## Specify parameters
spec <- matrix(c(
    'cpg', 'c', 1, 'logical', 'Either TRUE for CpG or FALSE for nonCpG',
	'model', 'm', 1, 'character', 'Either: cell, age or interaction',
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
    opt <- list('model' = 'interaction', 'cpg' = TRUE)
}

## Load data
load('rda/sgvc10_postnatal_polyA.Rdata', verbose = TRUE)
cpg <- ifelse(opt$cpg, 'CpG', 'nonCpG')

load_dmp <- function(is_cpg) {
    load(paste0('rda/DMP_', cpg, '_', opt$model, '.Rdata'), verbose = TRUE)
    if(is_cpg) {
        DMP <- DMP_CpG
    } else {
        DMP <- DMP_nonCpG
    }
    ## Keep Neurons only
    DMP <- DMP[, colData(DMP)$Cell.Type == 'Neuron']
    return(DMP)
}
DMP <- load_dmp(opt$cpg)

## For matching the brain ids
getid <- function(x) {
    as.integer(gsub('Br', '', x))
}

## Match and subset appropriately
message(paste(Sys.time(), 'subsetting the data to use'))
m  <- match(getid(colData(sgvc10)$BrNum), getid(colData(DMP)$Brain.ID))
sgvc10 <- sgvc10[, which(!is.na(m))]
colnames(sgvc10) <- paste0('Br', getid(colData(sgvc10)$BrNum))
DMP <- DMP[, m[!is.na(m)]]
colnames(DMP) <- paste0('Br', getid(colData(DMP)$Brain.ID))
# rse <- rse[, which(!is.na(m))]
# colnames(rse) <- paste0('Br', getid(colData(rse)$BrNum))

print('DMP info dimensions')
dim(DMP)
print('PSI info dimensions')
dim(sgvc10)

## Get methylation
message(paste(Sys.time(), 'preparing methylation info'))
meth <- SlicedData$new( getMeth(DMP, type = 'raw') )
meth$fileSliceSize <- 2000

methpos <- data.frame(
    cpg = paste0(cpg, '_', as.character(seqnames(DMP)), '_', start(DMP)),
    chr = as.character(seqnames(DMP)),
    pos = start(DMP),
    stringsAsFactors = FALSE
)

## Get PSI
message(paste(Sys.time(), 'preparing psi info'))
psi <- SlicedData$new(variantFreq(sgvc10))
psi$fileSliceSize <- 2000

psipos <- data.frame(
    spliceid = paste0('splice_', seq_len(nrow(sgvc10))),
    chr = as.character(seqnames(unlist(range(rowRanges(sgvc10))))),
    start = min(start(rowRanges(sgvc10))),
    end = max(end(rowRanges(sgvc10))),
    stringsAsFactors = FALSE
)

message(paste(Sys.time(), 'running MatrixEQTL'))
me <- Matrix_eQTL_main(snps = meth, gene = psi, 
    output_file_name.cis = paste0('.', cpg, '_', opt$model, '.txt'), # invis file, temporary
    pvOutputThreshold = 0, pvOutputThreshold.cis = 1e-5, 
	useModel = modelLINEAR,
	snpspos = methpos, genepos = genepos, cisDist = 1e4)

message(paste(Sys.time(), 'saving MatrixEQTL results'))
dir.create('rda', showWarnings = FALSE)
save(me, file = paste0('rda/me_', cpg, '_', opt$model, '.Rdata'))

## Reproducibility info
proc.time()
message(Sys.time())
options(width = 120)
session_info()
