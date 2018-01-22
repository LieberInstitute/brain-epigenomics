## Partially based on /users/ajaffe/Lieber/Projects/450k/scalpVsDura/eQTL.R
library('SGSeq')
library('bsseq')
library('MatrixEQTL')
library('getopt')
library('devtools')

## Specify parameters
spec <- matrix(c(
    'cpg', 'c', 1, 'logical', 'Either TRUE for CpG or FALSE for nonCpG',
    'feature', 'f', 1, 'character', 'Either: gene, exon, jx or psi',
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
    opt <- list('cpg' = FALSE, 'feature' = 'gene')
    opt <- list('cpg' = FALSE, 'feature' = 'exon')
    opt <- list('cpg' = FALSE, 'feature' = 'jx')
    opt <- list('cpg' = FALSE, 'feature' = 'psi')
}

stopifnot(opt$feature %in% c('psi', 'gene', 'exon', 'jx'))
cpg <- ifelse(opt$cpg, 'CpG', 'nonCpG')



## Load methylation data
load_dmp <- function(is_cpg) {
    if(is_cpg) {
        load('/dcl01/lieber/ajaffe/lab/brain-epigenomics/bumphunting/BSobj_bsseqSmooth_Neuron_minCov_3.Rdata', verbose = TRUE)
    } else {
        load('/dcl01/lieber/ajaffe/lab/brain-epigenomics/bsseq/bsobj_by_chr/allChrs_postNatal_cleaned_nonCG_noHomogenate_highCov.Rdata', verbose = TRUE)
    }
    ## Keep Neurons only
    BSobj <- BSobj[, colData(BSobj)$Cell.Type == 'Neuron']
    return(BSobj)
}
BSobj <- load_dmp(opt$cpg)

## Load expr data
load_expr <- function(type) {
    if(type == 'psi') {
        load('/dcl01/lieber/ajaffe/lab/brain-epigenomics/psi/rda/sgvc10_postnatal_polyA.Rdata', verbose = TRUE)
        expr <- sgvc10
    } else {
        load(paste0('rda/expr_', opt$feature, '.Rdata'), verbose = TRUE)
    }
    
    ## Subset to postnatal and PolyA only
    expr <- expr[, colData(expr)$Age > 0 & colData(expr)$Experiment == 'PolyA']
    return(expr)
}
expr <- load_expr(opt$feature)

## For matching the brain ids
getid <- function(x) {
    as.integer(gsub('Br', '', x))
}

## Match and subset appropriately
message(paste(Sys.time(), 'subsetting the data to use'))
m  <- match(getid(colData(expr)$BrNum), getid(colData(BSobj)$Brain.ID))
expr <- expr[, which(!is.na(m))]
colnames(expr) <- paste0('Br', getid(colData(expr)$BrNum))
BSobj <- BSobj[, m[!is.na(m)]]
colnames(BSobj) <- paste0('Br', getid(colData(BSobj)$Brain.ID))

print('BSobj info dimensions')
dim(BSobj)
print('Expression info dimensions')
dim(expr)

## Load meth vs expr results
message(paste(Sys.time(), 'loading and annotating me vs expr results'))
load(paste0('rda/me_', cpg, '_', opt$feature, '.Rdata'), verbose = TRUE)

## Annotate meth vs expr results
snp_id <- function(x) { as.integer(gsub('row', '', x)) }
me_annotated <- list(
    'eqtls' = DataFrame(me$cis$eqtls),
    'cinfo' = rowRanges(BSobj)[snp_id(me$cis$eqtls$snps), ],
    'exprinfo' = rowRanges(expr)[match(me$cis$eqtls$gene, names(rowRanges(expr))), ]
)

## Save the results
dir.create('rda', showWarnings = FALSE)
save(me_annotated, file = paste0('rda/me_annotated_', cpg, '_', opt$feature,
    '.Rdata'))

## Reproducibility info
proc.time()
message(Sys.time())
options(width = 120)
session_info()
