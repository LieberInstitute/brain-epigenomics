library('SGSeq')
library('bsseq')
library('getopt')
library('devtools')
library('limma')
library('gplots')
library('VennDiagram')
library('RColorBrewer')
library('clusterProfiler')
library('org.Hs.eg.db')
library('ggplot2')

## For testing
opt <- list('feature' = 'psi')

stopifnot(opt$feature %in% c('psi', 'gene', 'exon', 'jxleft', 'jxright'))
cpgs <- c('CpG', 'nonCpG', 'CpGmarg')

## For the CpGmarg
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

## Load data
load_expr <- function(type) {
    ## For jxright and jxleft
    type <- gsub('left|right', '', type)
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
        if(!file.exists(paste0('rda/expr_', type, '_unfiltered.Rdata'))) {
            dir.create('pdf', showWarnings = FALSE)
            pdf(paste0('pdf/suggested_expr_cutoffs_', tolower(type),
                       '.pdf'), width = 12)
            cuts <- expression_cutoff(assays(expr)$norm, seed = 20180119)
            message(paste(cuts, collapse = ' '))
            cut <- max(cuts)
            dev.off()

            meanExpr <- rowMeans(assays(expr)$norm)
            rowRanges(expr)$meanExprs <- meanExpr
            rowRanges(expr)$passExprsCut <- meanExpr > cut
            dir.create('rda', showWarnings = FALSE)
            save(expr, file = paste0('rda/expr_', type, '_unfiltered.Rdata'))
        }
    }

    return(expr)
}



## For matching the brain ids
getid <- function(x) {
    as.integer(gsub('Br', '', x))
}

## Subset to around a 1kb window from the nonCpG meQTl results
load_meqtl <- function() {
    load(paste0('rda/me_annotated_FDR5_nonCpG_', gsub('left|right', '',
                                                      opt$feature), '.Rdata'), verbose = TRUE)
    return(me_annotated)
}

message(paste(Sys.time(), 'loading the mres object'))
load(paste0('rda/meqtl_mres_', opt$feature, '_using_near.Rdata'), verbose = TRUE)

## Filter further to keep only those with meth_n >= 11
## and protein coding


gene <- load_expr('gene')

## Use the gene-level info to simplify things
for(cp in names(mres)) {
    ## Simplify rowRanges() for the psi info just to make the rest of the
    ## code work
    newgr <- unlist(range(rowRanges(mres[[cp]]$expr)))
    ## Match using gene id (just first one in cases with more than 1)
    mgr <- match(sapply(mcols(rowRanges(mres[[cp]]$expr))$geneName, '[[', 1), rowRanges(gene)$gencodeID)
    mcols(newgr) <- mcols(gene)[mgr, ]
    names(newgr) <- names(gene)[mgr]
    mres[[cp]]$eqtls$gene <- names(newgr)

    ## Instead of using the gene ranges, simply filter by them

    pc <- newgr$gene_type == 'protein_coding'
    meth11 <- mres[[cp]]$eqtls$meth_n >= 11
    meth11_all <- mres[[cp]]$eqtls$meth_all_n >= 11
    message(paste(Sys.time(), 'Filtering by meth_n >= 11 and meth_all_n >= 11 and protein_coding genes for set', cp))
    pc_meth_tab <- addmargins(table('protein coding' = pc, 'meth_n >= 11 & meth_all_n >= 11' = meth11 & meth11_all, useNA = 'ifany'))
    print(pc_meth_tab)
    print(round(pc_meth_tab / max(pc_meth_tab) * 100, 2))
    mres[[cp]]$eqtls <- mres[[cp]]$eqtls[which(pc & meth11 & meth11_all), ]
    mres[[cp]]$expr <- mres[[cp]]$expr[which(pc & meth11 & meth11_all), ]
    mres[[cp]]$meth <- mres[[cp]]$meth[which(pc & meth11 & meth11_all), ]
    rm(pc, meth11, meth11_all, pc_meth_tab, newgr)
}
rm(cp)

save(mres, file = paste0('rda/meqtl_mres_', opt$feature,
                         '_using_near_meth11_proteincoding_intactSGSeq.Rdata'))

## Reproducibility info
proc.time()
message(Sys.time())
options(width = 120)
session_info()

# Session info ----------------------------------------------------------------------------------------------------------
#     setting  value
# version  R version 3.4.3 Patched (2018-01-20 r74142)
# system   x86_64, linux-gnu
# ui       X11
# language (EN)
# collate  en_US.UTF-8
# tz       US/Eastern
# date     2018-03-16
#
# Packages --------------------------------------------------------------------------------------------------------------
#     package              * version   date       source
# acepack                1.4.1     2016-10-29 CRAN (R 3.4.1)
# AnnotationDbi        * 1.40.0    2017-11-29 Bioconductor
# assertthat             0.2.0     2017-04-11 CRAN (R 3.4.1)
# backports              1.1.2     2017-12-13 CRAN (R 3.4.2)
# base                 * 3.4.3     2018-01-20 local
# base64enc              0.1-3     2015-07-28 CRAN (R 3.4.1)
# bindr                  0.1       2016-11-13 CRAN (R 3.4.1)
# bindrcpp               0.2       2017-06-17 CRAN (R 3.4.1)
# Biobase              * 2.38.0    2017-11-07 Bioconductor
# BiocGenerics         * 0.24.0    2017-11-29 Bioconductor
# BiocParallel           1.12.0    2017-11-29 Bioconductor
# biomaRt                2.34.2    2018-02-17 Bioconductor
# Biostrings           * 2.46.0    2017-11-29 Bioconductor
# bit                    1.1-12    2014-04-09 CRAN (R 3.4.1)
# bit64                  0.9-7     2017-05-08 CRAN (R 3.4.1)
# bitops                 1.0-6     2013-08-17 CRAN (R 3.4.1)
# blob                   1.1.0     2017-06-17 CRAN (R 3.4.1)
# BSgenome               1.46.0    2017-11-29 Bioconductor
# bsseq                * 1.14.0    2017-11-07 Bioconductor
# bumphunter             1.20.0    2017-11-29 Bioconductor
# caTools                1.17.1    2014-09-10 CRAN (R 3.4.1)
# checkmate              1.8.5     2017-10-24 CRAN (R 3.4.2)
# cluster                2.0.6     2017-03-10 CRAN (R 3.4.3)
# clusterProfiler      * 3.6.0     2017-11-29 Bioconductor
# codetools              0.2-15    2016-10-05 CRAN (R 3.4.3)
# colorout             * 1.2-0     2018-02-19 Github (jalvesaq/colorout@2f01173)
# colorspace             1.3-2     2016-12-14 CRAN (R 3.4.1)
# compiler               3.4.3     2018-01-20 local
# data.table             1.10.4-3  2017-10-27 CRAN (R 3.4.2)
# datasets             * 3.4.3     2018-01-20 local
# DBI                    0.8       2018-03-02 CRAN (R 3.4.3)
# DelayedArray         * 0.4.1     2017-11-07 Bioconductor
# derfinder              1.12.6    2018-02-17 Bioconductor
# derfinderHelper        1.12.0    2017-11-29 Bioconductor
# devtools             * 1.13.4    2017-11-09 CRAN (R 3.4.2)
# digest                 0.6.15    2018-01-28 cran (@0.6.15)
# DO.db                  2.9       2017-08-10 Bioconductor
# doRNG                  1.6.6     2017-04-10 CRAN (R 3.4.1)
# DOSE                 * 3.4.0     2017-11-29 Bioconductor
# downloader             0.4       2015-07-09 CRAN (R 3.4.1)
# dplyr                  0.7.4     2017-09-28 CRAN (R 3.4.1)
# fastmatch              1.1-0     2017-01-28 CRAN (R 3.4.1)
# fgsea                  1.4.1     2018-02-17 Bioconductor
# foreach                1.4.4     2017-12-12 CRAN (R 3.4.2)
# foreign                0.8-69    2017-06-22 CRAN (R 3.4.3)
# Formula                1.2-2     2017-07-10 CRAN (R 3.4.1)
# futile.logger        * 1.4.3     2016-07-10 CRAN (R 3.4.1)
# futile.options         1.0.0     2010-04-06 CRAN (R 3.4.1)
# gdata                  2.18.0    2017-06-06 CRAN (R 3.4.1)
# GenomeInfoDb         * 1.14.0    2017-11-29 Bioconductor
# GenomeInfoDbData       1.0.0     2018-01-09 Bioconductor
# GenomicAlignments      1.14.1    2017-11-29 Bioconductor
# GenomicFeatures        1.30.3    2018-02-17 Bioconductor
# GenomicFiles           1.14.0    2017-11-29 Bioconductor
# GenomicRanges        * 1.30.2    2018-02-17 Bioconductor
# GEOquery               2.46.15   2018-03-06 Bioconductor
# getopt               * 1.20.2    2018-02-16 CRAN (R 3.4.3)
# ggplot2              * 2.2.1     2016-12-30 CRAN (R 3.4.1)
# glue                   1.2.0     2017-10-29 CRAN (R 3.4.2)
# GO.db                  3.5.0     2017-11-29 Bioconductor
# GOSemSim               2.4.1     2018-02-17 Bioconductor
# gplots               * 3.0.1     2016-03-30 CRAN (R 3.4.1)
# graphics             * 3.4.3     2018-01-20 local
# grDevices            * 3.4.3     2018-01-20 local
# grid                 * 3.4.3     2018-01-20 local
# gridExtra              2.3       2017-09-09 CRAN (R 3.4.1)
# gtable                 0.2.0     2016-02-26 CRAN (R 3.4.1)
# gtools                 3.5.0     2015-05-29 CRAN (R 3.4.1)
# Hmisc                  4.1-1     2018-01-03 CRAN (R 3.4.2)
# hms                    0.4.1     2018-01-24 CRAN (R 3.4.3)
# htmlTable              1.11.2    2018-01-20 CRAN (R 3.4.3)
# htmltools              0.3.6     2017-04-28 CRAN (R 3.4.1)
# htmlwidgets            0.9       2017-07-10 CRAN (R 3.4.1)
# httpuv                 1.3.6.2   2018-03-02 CRAN (R 3.4.3)
# httr                   1.3.1     2017-08-20 CRAN (R 3.4.1)
# igraph                 1.1.2     2017-07-21 CRAN (R 3.4.1)
# IRanges              * 2.12.0    2017-11-29 Bioconductor
# iterators              1.0.9     2017-12-12 CRAN (R 3.4.2)
# jsonlite               1.5       2017-06-01 CRAN (R 3.4.1)
# KernSmooth             2.23-15   2015-06-29 CRAN (R 3.4.3)
# knitr                  1.20      2018-02-20 CRAN (R 3.4.3)
# lambda.r               1.2       2017-09-16 CRAN (R 3.4.1)
# lattice                0.20-35   2017-03-25 CRAN (R 3.4.3)
# latticeExtra           0.6-28    2016-02-09 CRAN (R 3.4.1)
# lazyeval               0.2.1     2017-10-29 CRAN (R 3.4.2)
# limma                * 3.34.8    2018-02-17 Bioconductor
# locfit                 1.5-9.1   2013-04-20 CRAN (R 3.4.1)
# magrittr               1.5       2014-11-22 CRAN (R 3.4.1)
# Matrix                 1.2-12    2017-11-30 CRAN (R 3.4.3)
# matrixStats          * 0.53.1    2018-02-11 CRAN (R 3.4.3)
# memoise                1.1.0     2017-04-21 CRAN (R 3.4.1)
# methods              * 3.4.3     2018-01-20 local
# munsell                0.4.3     2016-02-13 CRAN (R 3.4.1)
# nnet                   7.3-12    2016-02-02 CRAN (R 3.4.3)
# org.Hs.eg.db         * 3.5.0     2017-11-29 Bioconductor
# parallel             * 3.4.3     2018-01-20 local
# permute                0.9-4     2016-09-09 CRAN (R 3.4.1)
# pillar                 1.1.0     2018-01-14 CRAN (R 3.4.2)
# pkgconfig              2.0.1     2017-03-21 CRAN (R 3.4.1)
# pkgmaker               0.22      2014-05-14 CRAN (R 3.4.1)
# plyr                   1.8.4     2016-06-08 CRAN (R 3.4.1)
# png                    0.1-7     2013-12-03 CRAN (R 3.4.1)
# prettyunits            1.0.2     2015-07-13 CRAN (R 3.4.1)
# progress               1.1.2     2016-12-14 CRAN (R 3.4.1)
# purrr                  0.2.4     2017-10-18 CRAN (R 3.4.2)
# qvalue                 2.10.0    2017-11-29 Bioconductor
# R.methodsS3            1.7.1     2016-02-16 CRAN (R 3.4.1)
# R.oo                   1.21.0    2016-11-01 CRAN (R 3.4.1)
# R.utils                2.6.0     2017-11-05 CRAN (R 3.4.2)
# R6                     2.2.2     2017-06-17 CRAN (R 3.4.1)
# RColorBrewer         * 1.1-2     2014-12-07 CRAN (R 3.4.1)
# Rcpp                   0.12.14   2017-11-23 CRAN (R 3.4.2)
# RCurl                  1.95-4.10 2018-01-04 CRAN (R 3.4.2)
# readr                  1.1.1     2017-05-16 CRAN (R 3.4.1)
# recount                1.4.4     2018-02-17 Bioconductor
# registry               0.5       2017-12-03 CRAN (R 3.4.2)
# rentrez                1.2.0     2018-02-12 CRAN (R 3.4.3)
# reshape2               1.4.3     2017-12-11 CRAN (R 3.4.2)
# rlang                  0.1.6     2017-12-21 CRAN (R 3.4.2)
# rmote                * 0.3.4     2018-02-16 deltarho (R 3.4.3)
# RMySQL                 0.10.14   2018-02-26 CRAN (R 3.4.3)
# rngtools               1.2.4     2014-03-06 CRAN (R 3.4.1)
# rpart                  4.1-12    2018-01-12 CRAN (R 3.4.3)
# Rsamtools            * 1.30.0    2017-11-29 Bioconductor
# RSQLite                2.0       2017-06-19 CRAN (R 3.4.1)
# rstudioapi             0.7       2017-09-07 CRAN (R 3.4.1)
# rtracklayer            1.38.3    2018-02-17 Bioconductor
# RUnit                  0.4.31    2015-11-06 CRAN (R 3.4.1)
# rvcheck                0.0.9     2017-07-10 CRAN (R 3.4.2)
# S4Vectors            * 0.16.0    2017-11-29 Bioconductor
# scales                 0.5.0     2017-08-24 CRAN (R 3.4.1)
# servr                  0.8       2017-11-06 CRAN (R 3.4.3)
# SGSeq                * 1.12.0    2018-01-16 Bioconductor
# splines                3.4.3     2018-01-20 local
# stats                * 3.4.3     2018-01-20 local
# stats4               * 3.4.3     2018-01-20 local
# stringi                1.1.6     2017-11-17 CRAN (R 3.4.2)
# stringr                1.2.0     2017-02-18 CRAN (R 3.4.1)
# SummarizedExperiment * 1.8.1     2018-01-09 Bioconductor
# survival               2.41-3    2017-04-04 CRAN (R 3.4.3)
# tibble                 1.4.1     2017-12-25 CRAN (R 3.4.2)
# tidyr                  0.8.0     2018-01-29 CRAN (R 3.4.3)
# tools                  3.4.3     2018-01-20 local
# utils                * 3.4.3     2018-01-20 local
# VariantAnnotation      1.24.5    2018-01-16 Bioconductor
# VennDiagram          * 1.6.18    2017-11-21 CRAN (R 3.4.2)
# withr                  2.1.1     2017-12-19 CRAN (R 3.4.2)
# XML                    3.98-1.10 2018-02-19 CRAN (R 3.4.3)
# xml2                   1.2.0     2018-01-24 CRAN (R 3.4.3)
# xtable                 1.8-2     2016-02-05 CRAN (R 3.4.1)
# XVector              * 0.18.0    2017-11-29 Bioconductor
# zlibbioc               1.24.0    2017-11-07 Bioconductor
