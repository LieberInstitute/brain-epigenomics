## Partially based on /users/ajaffe/Lieber/Projects/450k/scalpVsDura/eQTL.R
library('SGSeq')
library('bsseq')
library('MatrixEQTL')
library('getopt')
library('jaffelab')
library('recount')
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

## Load data
load_expr <- function(type) {
    if(type == 'psi') {
        load('rda/sgvc10_postnatal_polyA.Rdata', verbose = TRUE)
        expr <- sgvc10
    } else if (type == 'gene') {
        load('/dcl01/lieber/ajaffe/lab/brain-epigenomics/brainseq_pipeline/polyA_unstranded/rse_gene_polyA_dlpfc_n41.Rdata', verbose = TRUE)
        assays(rse_gene)$norm <- recount::getRPKM(rse_gene, 'Length')
        ## RPKM
        expr <- rse_gene
    } else if (type == 'exon') {
        load('/dcl01/lieber/ajaffe/lab/brain-epigenomics/brainseq_pipeline/polyA_unstranded/rse_exon_polyA_dlpfc_n41.Rdata', verbose = TRUE)
        assays(rse_exon)$norm <- recount::getRPKM(rse_exon, 'Length')
        ## RPKM
        expr <- rse_exon
    } else if (type == 'jx') {
        if(file.exists(paste0('rda/expr_', opt$feature, '.Rdata'))) {
            ## Ran interactively and saved the results on 2018-01-19
            load(paste0('rda/expr_', opt$feature, '.Rdata'), verbose = TRUE)
            return(expr)
        }
        load('/dcl01/lieber/ajaffe/lab/brain-epigenomics/brainseq_pipeline/polyA_unstranded/rse_jx_polyA_dlpfc_n41.Rdata', verbose = TRUE)
        rowRanges(rse_jx)$Length <- 100 / 8
        ## RP80m
        assays(rse_jx)$norm <- recount::getRPKM(rse_jx, 'Length')
        expr <- rse_jx
    }
    
    ## Subset to postnatal and PolyA only
    expr <- expr[, colData(expr)$Age > 0 & colData(expr)$Experiment == 'PolyA']
    
    ## Drop low expressed features
    if(type != 'psi') {
        dir.create('pdf', showWarnings = FALSE)
        pdf(paste0('pdf/suggested_expr_cutoffs_', tolower(opt$feature),
            '.pdf'), width = 12)
        cuts <- expression_cutoff(assays(expr)$norm, seed = 20180119)
        message(paste(cuts, collapse = ' '))
        cut <- max(cuts)
        dev.off()
        
        meanExpr <- rowMeans(assays(expr)$norm)
        rowRanges(expr)$meanExprs <- meanExpr
        rowRanges(expr)$passExprsCut <- meanExpr > cut
        dir.create('rda', showWarnings = FALSE)
        save(expr, file = paste0('rda/expr_', opt$feature, '_unfiltered.Rdata'))
        expr <- expr[rowRanges(expr)$passExprsCut]
        save(expr, file = paste0('rda/expr_', opt$feature, '.Rdata'))
    }
        
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

## Get methylation
message(paste(Sys.time(), 'preparing methylation info'))
meth <- SlicedData$new( getMeth(BSobj, type = 'raw') )
meth$fileSliceSize <- 2000

methpos <- data.frame(
    cname = paste0('row', seq_len(nrow(BSobj))),
    chr = as.character(seqnames(BSobj)),
    pos = start(BSobj),
    stringsAsFactors = FALSE
)

## Get PSI
message(paste(Sys.time(), 'preparing expression info'))
if(opt$feature == 'psi') {
    exprinfo <- SlicedData$new(variantFreq(expr))
} else {
    exprinfo <- SlicedData$new(log2(assays(expr)$norm + 1))
}
exprinfo$fileSliceSize <- 2000

get_exprpos <- function(type) {
    if(type == 'psi') {
        exprpos <- data.frame(
            spliceid = paste0('splice_', seq_len(nrow(expr))),
            chr = as.character(seqnames(unlist(range(rowRanges(expr))))),
            start = min(start(rowRanges(expr))),
            end = max(end(rowRanges(expr))),
            stringsAsFactors = FALSE
        )
    } else if (type == 'gene') {
        exprpos <- data.frame(
            spliceid = rowRanges(expr)$gencodeID,
            chr = as.character(seqnames(rowRanges(expr))),
            start = start(rowRanges(expr)),
            end = end(rowRanges(expr)),
            stringsAsFactors = FALSE
        )
    } else if (type %in% c('exon', 'jx')) {
        exprpos <- data.frame(
            spliceid = names(rowRanges(expr)),
            chr = as.character(seqnames(rowRanges(expr))),
            start = start(rowRanges(expr)),
            end = end(rowRanges(expr)),
            stringsAsFactors = FALSE
        )
    }
    return(exprpos)
}
exprpos <- get_exprpos(opt$feature)

message(paste(Sys.time(), 'running MatrixEQTL'))
me <- Matrix_eQTL_main(snps = meth, gene = exprinfo, 
    output_file_name.cis = paste0('.', cpg, '_', opt$feature,
        '.txt'), # invis file, temporary
    pvOutputThreshold = 0, pvOutputThreshold.cis = 1e-5, 
	useModel = modelLINEAR,
	snpspos = methpos, genepos = exprpos, cisDist = 1000)

message(paste(Sys.time(), 'saving MatrixEQTL results'))
dir.create('rda', showWarnings = FALSE)
save(me, file = paste0('rda/me_', cpg, '_', opt$feature, '.Rdata'))

## Reproducibility info
proc.time()
message(Sys.time())
options(width = 120)
session_info()
