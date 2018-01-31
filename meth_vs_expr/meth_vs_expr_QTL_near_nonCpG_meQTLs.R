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
    opt <- list('cpg' = FALSE, 'feature' = 'psi')
    opt <- list('cpg' = TRUE, 'feature' = 'gene')
    opt <- list('cpg' = TRUE, 'feature' = 'exon')
    opt <- list('cpg' = TRUE, 'feature' = 'jx')
    opt <- list('cpg' = TRUE, 'feature' = 'psi')
}

stopifnot(opt$feature %in% c('psi', 'gene', 'exon', 'jx'))
cpg <- ifelse(opt$cpg, 'CpG', 'nonCpG')



## Load methylation data
load_dmp <- function(is_cpg) {
    if(is_cpg) {
        load('/dcl01/lieber/ajaffe/lab/brain-epigenomics/bumphunting/BSobj_bsseqSmooth_Neuron_minCov_3.Rdata', verbose = TRUE)
        rowRanges(BSobj)$c_context <- Rle('CG')
        rowRanges(BSobj)$trinucleotide_context <- Rle('CpG')
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
        load('/dcl01/lieber/ajaffe/lab/brain-epigenomics/psi/rda/sgvc10_postnatal_polyA.Rdata', verbose = TRUE)
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
        if(!file.exists(paste0('rda/expr_', opt$feature, '_unfiltered.Rdata'))) {
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
        }
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

## Subset to around a 1kb window from the nonCpG meQTl results
load_meqtl <- function() {
    load(paste0('rda/me_annotated_FDR5_nonCpG_', opt$feature,
        '.Rdata'), verbose = TRUE)
    return(me_annotated)
}
meqtl <- load_meqtl()
me_ov <- findOverlaps(resize(rowRanges(meqtl$expr), width(rowRanges(meqtl$expr)) + 2000, fix = 'center'), expr)

message(paste(Sys.time(), 'keeping the following percent of genes'))
round(length(unique(subjectHits(me_ov))) / nrow(expr) * 100, 2)
expr <- expr[sort(unique(subjectHits(me_ov))), ]

cp_ov <- findOverlaps(resize(rowRanges(expr), width(rowRanges(expr)) + 2000, fix = 'center'), BSobj)

message(paste(Sys.time(), 'keeping the following percent of CpGs'))
round(length(unique(subjectHits(cp_ov))) / nrow(BSobj) * 100, 2)
BSobj <- BSobj[sort(unique(subjectHits(cp_ov))), ]

if(FALSE) {
    ## For debugging
    BSobj_raw <- BSobj
    expr_raw <- expr
    
    BSobj <- BSobj_raw
    expr <- expr_raw
        
    BSobj <- BSobj[seqnames(rowRanges(BSobj)) %in% c('chr1', 'chr6', 'chr10', 'chr21'), ]
    expr <- expr[seqnames(rowRanges(expr)) %in% c('chr1', 'chr6', 'chr10', 'chr21'), ]
}

print('BSobj info dimensions')
dim(BSobj)
print('Expression info dimensions')
dim(expr)

## Get methylation
message(paste(Sys.time(), 'preparing methylation info'))
meth <- SlicedData$new( getMeth(BSobj, type = 'raw') )
meth$fileSliceSize <- ifelse(opt$cpg, 300, 2000)

methpos <- data.frame(
    cname = paste0('row', seq_len(nrow(BSobj))),
    chr = as.character(seqnames(BSobj)),
    pos = start(BSobj),
    stringsAsFactors = FALSE
)

## Get expression or PSI info
message(paste(Sys.time(), 'preparing expression info'))
if(opt$feature == 'psi') {
    exprinfo <- SlicedData$new(variantFreq(expr))
} else {
    exprinfo <- SlicedData$new(log2(assays(expr)$norm + 1))
}
exprinfo$fileSliceSize <- ifelse(opt$cpg, 300, 2000)

get_exprpos <- function(type) {
    if(type == 'psi') {
        exprpos <- data.frame(
            spliceid = rownames(exprinfo),
            chr = as.character(seqnames(unlist(range(rowRanges(expr))))),
            start = min(start(rowRanges(expr))),
            end = max(end(rowRanges(expr))),
            stringsAsFactors = FALSE
        )
    } else {
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
        '_near_nonCpG_meQTLs.txt'), # invis file, temporary
    pvOutputThreshold = 0, pvOutputThreshold.cis = 0.01, 
	useModel = modelLINEAR,
	snpspos = methpos, genepos = exprpos, cisDist = 1000)
## Check warnings if any are present
warnings()
    
message(paste(Sys.time(), 'saving MatrixEQTL results'))
dir.create('rda', showWarnings = FALSE)
save(me, file = paste0('rda/me_', cpg, '_', opt$feature, '_near_nonCpG_meQTLs.Rdata'))

print('neqtls and ntests')
me$cis$neqtls
me$cis$ntests

print('Number of meQTLs without row ids')
table(!is.na(me$cis$eqtls$snps))
me$cis$eqtls <- me$cis$eqtls[!is.na(me$cis$eqtls$snps), ]

print('FDR distribution')
summary(me$cis$eqtls$FDR)

## Keep only FDR <5%
print('meQTLs less than 0.05 FDR (will annotate those)')
table(me$cis$eqtls$FDR < 0.05)
me$cis$eqtls <- me$cis$eqtls[me$cis$eqtls$FDR < 0.05, ]

print('Number of finite statistics (for meQTLs with FDR < 5%)')
table(is.finite(me$cis$eqtls$statistic))

me$cis$eqtls <- me$cis$eqtls[is.finite(me$cis$eqtls$statistic), ]

## Annotate meth vs expr results
snp_id <- function(x) { as.integer(gsub('row', '', x)) }
if(opt$feature == 'psi') {
    me_annotated <- list(
        'eqtls' = DataFrame(me$cis$eqtls),
        'meth' = BSobj[snp_id(me$cis$eqtls$snps), ],
        'expr' = expr[snp_id(me$cis$eqtls$gene), ]
    )
} else {
    me_annotated <- list(
        'eqtls' = DataFrame(me$cis$eqtls),
        'meth' = BSobj[snp_id(me$cis$eqtls$snps), ],
        'expr' = expr[match(me$cis$eqtls$gene, names(rowRanges(expr))), ]
    )
}
save(me_annotated,
    file = paste0('rda/me_annotated_FDR5_', cpg, '_', opt$feature,
    '_near_nonCpG_meQTLs.Rdata'))
    
## Explore annotated mEQTLs with FDR < 5%
print('Trinucleotide vs C context')
addmargins(table('trinucleotide_context' = as.vector(rowRanges(me_annotated$meth)$trinucleotide_context), 'c_context' = as.vector(rowRanges(me_annotated$meth)$c_context)))

print('Trinucleotide context vs C strand')
addmargins(table('trinucleotide_context' = as.vector(rowRanges(me_annotated$meth)$trinucleotide_context), 'C strand' = as.vector(strand(me_annotated$meth))))


print('Number of meQTLs per feature')
table(table(me_annotated$eqtls$gene))

print('Number of features per C pos')
table(table(me_annotated$eqtls$snps))

print('C strand vs feature strand')
if(opt$feature == 'psi') {
    print(addmargins(table('C strand' = as.vector(strand(me_annotated$meth)), 'feature strand' = as.vector(strand(unlist(range(rowRanges(me_annotated$expr))))))))
    print(chisq.test(table('C strand' = as.vector(strand(me_annotated$meth)), 'feature strand' = as.vector(strand(unlist(range(rowRanges(me_annotated$expr))))))))
} else {
    print(addmargins(table('C strand' = as.vector(strand(me_annotated$meth)), 'feature strand' = as.vector(strand(me_annotated$expr)))))
    print(chisq.test(table('C strand' = as.vector(strand(me_annotated$meth)), 'feature strand' = as.vector(strand(me_annotated$expr)))))
}   


if(opt$feature %in% c('gene', 'exon')) {
    print('meQTLs by gene type (all, then by unique gene)')
    print(table(rowRanges(me_annotated$expr)$gene_type))
    print(table(rowRanges(me_annotated$expr)$gene_type[!duplicated(names(rowRanges(me_annotated$expr)))]))
} else if (opt$feature == 'jx') {
    print('meQTls by jx class (all, then by unique jx)')
    print(table(rowRanges(me_annotated$expr)$Class))
    print(table(rowRanges(me_annotated$expr)$Class[!duplicated(names(rowRanges(me_annotated$expr)))]))
} else if (opt$feature == 'psi') {
    print('meQTLs by variant type (all, then by unique splicing event)')
    print(table(paste(mcols(rowRanges(me_annotated$expr))$variantType, collapse = ' ')))
    print(table(paste(mcols(rowRanges(me_annotated$expr))$variantType[!duplicated(mcols(rowRanges(me_annotated$expr))$variantName)], collapse = ' ')))
}

## Make scatter plots of the top 100
ylab <- ifelse(opt$feature == 'psi', 'PSI', ifelse(opt$feature == 'jx', 'log2 (RP80M + 1)', 'log2 (RPKM + 1)'))
y <- if(opt$feature == 'psi') variantFreq(me_annotated$expr) else log2(assays(me_annotated$expr)$norm + 1)

plotting_code <- function(i) {
    if(opt$feature == 'psi') {
        main <- paste(opt$feature, mcols(me_annotated$expr)$variantName[i], 'FDR', signif(me_annotated$eqtls$FDR[i], 3))
    } else {
        main <- paste(opt$feature, me_annotated$eqtls$gene[i], 'FDR', signif(me_annotated$eqtls$FDR[i], 3))
    }
    
    plot(x = jitter(getMeth(me_annotated$meth[i, ], type = 'raw'), 0.05), y = jitter(y[i, ], 0.05), xlab = 'Methylation', ylab = ylab, main = main, sub = paste(as.vector(seqnames(rowRanges(me_annotated$meth)[i])), start(rowRanges(me_annotated$meth)[i]), as.vector(strand(rowRanges(me_annotated$meth)[i])), as.vector(rowRanges(me_annotated$meth)$c_context[i])))
}

dir.create('pdf', showWarnings = FALSE)
pdf(paste0('pdf/top100_FDR5_', cpg, '_', opt$feature, '_near_nonCpG_meQTLs.pdf'))
for(i in seq_len(100)) {
    plotting_code(i)
}
dev.off()


meth_n <- rowSums(getMeth(me_annotated$meth, type = 'raw') > 0)
message(paste(Sys.time(), 'meth_n summary data (raw, percent, cumulative percent)'))
table(meth_n)
round(table(meth_n) / length(meth_n) * 100, 2)
cumsum(round(table(meth_n) / length(meth_n) * 100, 2))

if(opt$feature == 'psi') {
    expr_delta <- apply(variantFreq(me_annotated$expr), 1, function(x) diff(range(x)) )
} else {
    expr_delta <- apply(log2(assays(me_annotated$expr)$norm + 1), 1, function(x) diff(range(x)) )
}

message(paste(Sys.time(), 'expression delta summary'))
summary(expr_delta)

addmargins(table('Expr Delta >= 0.1' = expr_delta >= 0.1, 'Meth N >= 4' = meth_n >= 4))
round(addmargins(table('Expr Delta >= 0.1' = expr_delta >= 0.1, 'Meth N >= 4' = meth_n >= 4)) / length(meth_n) * 100, 2)

## Make scatter plots of the top 100 with at least 4 samples with non-zero meth
if(length(which(meth_n >= 4)) > 0) {
    pdf(paste0('pdf/top100_FDR5_min_Meth4_', cpg, '_', opt$feature, '_near_nonCpG_meQTLs.pdf'))
    for(i in head(which(meth_n >= 4), 100)) {
        plotting_code(i)
    }
    dev.off()
} else {
    message(paste(Sys.time(), 'found no meQTLs with meth_n >= 4'))
}

## Make scatter plots of the top 100 with at least 4 samples with non-zero meth and a expr change of at least 0.1
if(length(which(meth_n >= 4 & expr_delta >= 0.1)) > 0) {
    pdf(paste0('pdf/top100_FDR5_min_Meth4_exprDelta0.1_', cpg, '_', opt$feature, '_near_nonCpG_meQTLs.pdf'))
    for(i in head(which(meth_n >= 4 & expr_delta >= 0.1), 100)) {
        plotting_code(i)
    }
    dev.off()
} else {
    message(paste(Sys.time(), 'found no meQTLs with meth_n >= 4 & expr_delta > 0.1'))
}



## Reproducibility info
proc.time()
message(Sys.time())
options(width = 120)
session_info()
