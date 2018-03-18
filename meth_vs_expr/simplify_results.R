library('SGSeq')
library('bsseq')
library('devtools')

features <- c('gene', 'exon', 'psi')

meth_data <- lapply(features, function(feature) {
    if(feature == 'psi') {
        f <- paste0('rda/meqtl_mres_', feature, '_using_near_meth11_proteincoding_intactSGSeq.Rdata')
    } else {
        f <- paste0('rda/meqtl_mres_', feature, '_using_near_meth11_proteincoding.Rdata')
    }
    message(paste(Sys.time(), 'loading', f))
    load(f, verbose = TRUE)

    message(paste(Sys.time(), 'processing mres'))
    res <- lapply(mres[c('CpG', 'nonCpG')], function(x) {
        y <- x[c('expr', 'meth')]
        y$ids <- data.frame('c_id' = as.character(x$eqtls$snps),
            'feature_id' = as.character(x$eqtls$gene),
            'FDR' = x$eqtls$FDR,
            stringsAsFactors = FALSE)
        return(y)
    })
    return(res)
})
names(meth_data) <- features

print(object.size(meth_data), units = 'Mb')
# 31.4 Mb

save(meth_data, file = 'rda/meth_data.Rdata')
system('ls -lh rda/meth_data.Rdata')
# 25 Mb

meth_summary <- do.call(rbind, lapply(features, function(feature) {
    f <-  paste0('rda/meqtl_data_by_venn_', feature, '_using_near.Rdata')
    message(paste(Sys.time(), 'loading', f))
    load(f, verbose = TRUE)

    message(paste(Sys.time(), 'processing data_by_venn'))
    ## Keep only FDR < 0.05
    res <- data_by_venn[data_by_venn$FDR < 0.05, ]
    res$typeref[res$typeref == 'CpGmarg'] <- 'CpG'

    colnamesmap <- c(
        'snps' = 'c_id',
        'gene' = 'feature_id',
        'statistic' = 'meth_statistic',
        'pvalue' = 'meth_pvalue',
        'FDR' = 'meth_FDR',
        'beta' = 'meth_coefficient',
        'meth_n' = 'n_samples_with_meth_not0',
        'meth_all_n' = 'n_samples_with_meth_not1',
        'expr_delta' = 'expr_delta',
        'Estimate' = 'age_coefficient',
        'Std. Error' = 'age_std_error',
        't value' = 'age_statistic',
        'Pr(>|t|)' = 'age_pvalue',
        'ageFDR' = 'age_FDR',
        'agemethFDR' = 'meth_adjust_age_FDR',
        'noage' = 'meth_adjust_age_FDR_less5percent',
        'promoter_geneid' = 'promoter_geneid',
        'promoter_present' = 'promoter_present',
        'promoter_gencodeID' = 'promoter_gencodeID',
        'body_geneid' = 'body_geneid',
        'body_present' = 'body_present',
        'body_gencodeID' = 'body_gencodeID',
        'flanking_geneid' = 'flanking_geneid',
        'flanking_present' = 'flanking_present',
        'flanking_gencodeID' = 'flanking_gencodeID',
        'vset' = 'venn_set',
        'typeref' = 'meth_type'
    )
    colnames(res) <- colnamesmap[colnames(res)]

    res$feature <- feature
    return(res)
}))


## Add Neuron or glia info
venn_summ <- lapply(features, function(feature) {
    f <-  paste0('rda/meqtl_data_venn_summ_', feature, '_using_near.Rdata')
    message(paste(Sys.time(), 'loading', f))
    load(f, verbose = TRUE)
    return(data_venn_summ)
})
names(venn_summ) <- features

meth_summary$gene_type <- NA
for(i in features) {
    j <- which(meth_summary$feature == i)
    m <- match(meth_summary$feature_id[j], venn_summ[[i]]$gene)
    print(table(is.na(m)))
    meth_summary$gene_type[j[!is.na(m)]] <- venn_summ[[i]]$gtype[m[!is.na(m)]]
}

## Add id mapping
meth_summary$i <- NA
for(i in unique(meth_summary$feature)) {
    message(paste(Sys.time(), 'processing', i))
    for(j in unique(meth_summary$meth_type)) {
        message(paste(Sys.time(), 'processing', j))
        k <- which(meth_summary$feature == i & meth_summary$meth_type == j)
        meth_summary$i[k] <- mapply(function(g, s) {
            which(meth_data[[i]][[j]]$ids$feature_id == g & meth_data[[i]][[j]]$ids$c_id == s)
        }, meth_summary$feature_id[k], meth_summary$c_id[k])
    }
}
meth_summary$i <- NumericList(meth_summary$i)

print(object.size(meth_summary), units = 'Mb')
# 64.4 Mb
save(meth_summary, file = 'rda/meth_summary.Rdata')
system('ls -lh rda/meth_summary.Rdata')
# 9.8 Mb


meth_df <- as.data.frame(meth_summary)
for(i in which(sapply(meth_df, class) == 'AsIs')) {
    meth_df[, i] <- paste(meth_summary[, i], collapse = ',')
}

print(object.size(meth_df), units = 'Mb')
# 64.2 Mb
save(meth_df, file = 'rda/meth_df.Rdata')
system('ls -lh rda/meth_df.Rdata')
# 9.2 M

## summarize TF results
load('/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/PWMEnrich/spliceRegions_PWMEnrich_objects_exonAdjustedBackgrounds.rda', verbose = TRUE)

# info <- score; pval_info <- pval; idvar <- 'exonID'
# info <- scoreC; pval_info <- pvalC; idvar <- 'C_ID'

summarize_tf <- function(info, pval_info, idvar, motifs = c('MECP2', 'CTCF')) {
    s <- info[which(info$motifID %in% motifs), ]
    p <- pval_info[which(pval_info$motifID %in% motifs), ]

    stopifnot(identical(rownames(s), rownames(p)))
    stopifnot(identical(s$motifID, p$motifID))
    stopifnot(identical(s$Target, p$Target))
    stopifnot(identical(as.character(s$motifID), as.character(s$Target)))

    if(idvar == 'C_ID') {
        s$C_ID <- as.character(s$C_ID)
        res <- GRanges(
            sapply(strsplit(s$C_ID, ':'), '[', 1),
            IRanges(
                start = as.numeric(sapply(strsplit(gsub('.*:', '', s$C_ID), '-'), '[', 1)),
                end = as.numeric(sapply(strsplit(gsub('.*:', '', s$C_ID), '-'), '[', 2))
            ),
            motif = s$motifID,
            direction = ifelse(s$Direction == 'Positive Association', 'up', 'down'),
            c_in_exon =  s$C_in_exon == 'C in Exon',
            raw_score = s$value,
            FDR = p$FDR,
            i = NA
        )
        ov <- findOverlaps(res, rowRanges(meth_data$exon$CpG$meth))
        res$i[s$Context == 'CpG'] <- split(subjectHits(ov), queryHits(ov))
        ov <- findOverlaps(res, rowRanges(meth_data$exon$nonCpG$meth))
        res$i[s$Context == 'nonCpG'] <- split(subjectHits(ov), queryHits(ov))
        res$i <- IntegerList(res$i)
 #       names(res) <- rownames(s)
    } else {
        res <- data.frame(
            motif = s$motifID,
            direction = ifelse(s$Direction == 'Positive Association', 'up', 'down'),
            c_in_exon = s$C_in_exon == 'C in Exon',
            raw_score = s$value,
            FDR = p$value,
            feature_id = s$exonID,
            i = NA
        )
        res$i[s$Context == 'CpG'] <- match(res$feature_id[s$Context == 'CpG'], rowRanges(meth_data$exon$CpG$expr)$exon_libdID )
        res$i[s$Context == 'nonCpG'] <- match(res$feature_id[s$Context == 'nonCpG'], rowRanges(meth_data$exon$nonCpG$expr)$exon_libdID )
 #       rownames(res) <- rownames(s)
    }

    final <- split(res, s$Context)
    return(final)
}

tf_data <- mapply(summarize_tf, list(score, scoreC), list(pval, pvalC), c('exonID', 'C_ID'))
names(tf_data) <- c('exon', 'cbase')

lapply(tf_data, function(x) lapply(x, head))

print(object.size(tf_data), units = 'Mb')
# 4.9Mb
save(tf_data, file = 'rda/tf_data.Rdata')
system('ls -lh rda/tf_data.Rdata')
# 1.2M

tools::md5sum(c('rda/meth_data.Rdata', 'rda/meth_summary.Rdata', 'rda/meth_df.Rdata', 'rda/tf_data.Rdata'))
# rda/meth_data.Rdata             rda/meth_summary.Rdata
# "c27b0f98c853be4a1e250665bae6c293" "060efade56edae016daef466f691e5fe"
# rda/meth_df.Rdata                  rda/tf_data.Rdata
# "16f274ba4843a423e854d9a5846c2705" "b061076ae8c39e8ab040f0dc37f52ccf"


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
# AnnotationDbi          1.40.0    2017-11-29 Bioconductor
# assertthat             0.2.0     2017-04-11 CRAN (R 3.4.1)
# base                 * 3.4.3     2018-01-20 local
# Biobase              * 2.38.0    2017-11-07 Bioconductor
# BiocGenerics         * 0.24.0    2017-11-29 Bioconductor
# BiocParallel           1.12.0    2017-11-29 Bioconductor
# biomaRt                2.34.2    2018-02-17 Bioconductor
# Biostrings           * 2.46.0    2017-11-29 Bioconductor
# bit                    1.1-12    2014-04-09 CRAN (R 3.4.1)
# bit64                  0.9-7     2017-05-08 CRAN (R 3.4.1)
# bitops                 1.0-6     2013-08-17 CRAN (R 3.4.1)
# blob                   1.1.0     2017-06-17 CRAN (R 3.4.1)
# bsseq                * 1.14.0    2017-11-07 Bioconductor
# colorout             * 1.2-0     2018-02-19 Github (jalvesaq/colorout@2f01173)
# colorspace             1.3-2     2016-12-14 CRAN (R 3.4.1)
# compiler               3.4.3     2018-01-20 local
# data.table             1.10.4-3  2017-10-27 CRAN (R 3.4.2)
# datasets             * 3.4.3     2018-01-20 local
# DBI                    0.8       2018-03-02 CRAN (R 3.4.3)
# DelayedArray         * 0.4.1     2017-11-07 Bioconductor
# devtools             * 1.13.4    2017-11-09 CRAN (R 3.4.2)
# digest                 0.6.15    2018-01-28 cran (@0.6.15)
# GenomeInfoDb         * 1.14.0    2017-11-29 Bioconductor
# GenomeInfoDbData       1.0.0     2018-01-09 Bioconductor
# GenomicAlignments      1.14.1    2017-11-29 Bioconductor
# GenomicFeatures        1.30.3    2018-02-17 Bioconductor
# GenomicRanges        * 1.30.2    2018-02-17 Bioconductor
# ggplot2                2.2.1     2016-12-30 CRAN (R 3.4.1)
# graphics             * 3.4.3     2018-01-20 local
# grDevices            * 3.4.3     2018-01-20 local
# grid                   3.4.3     2018-01-20 local
# gtable                 0.2.0     2016-02-26 CRAN (R 3.4.1)
# gtools                 3.5.0     2015-05-29 CRAN (R 3.4.1)
# htmltools              0.3.6     2017-04-28 CRAN (R 3.4.1)
# htmlwidgets            0.9       2017-07-10 CRAN (R 3.4.1)
# httpuv                 1.3.6.2   2018-03-02 CRAN (R 3.4.3)
# httr                   1.3.1     2017-08-20 CRAN (R 3.4.1)
# igraph                 1.1.2     2017-07-21 CRAN (R 3.4.1)
# IRanges              * 2.12.0    2017-11-29 Bioconductor
# lattice                0.20-35   2017-03-25 CRAN (R 3.4.3)
# lazyeval               0.2.1     2017-10-29 CRAN (R 3.4.2)
# limma                  3.34.8    2018-02-17 Bioconductor
# locfit                 1.5-9.1   2013-04-20 CRAN (R 3.4.1)
# magrittr               1.5       2014-11-22 CRAN (R 3.4.1)
# Matrix                 1.2-12    2017-11-30 CRAN (R 3.4.3)
# matrixStats          * 0.53.1    2018-02-11 CRAN (R 3.4.3)
# memoise                1.1.0     2017-04-21 CRAN (R 3.4.1)
# methods              * 3.4.3     2018-01-20 local
# munsell                0.4.3     2016-02-13 CRAN (R 3.4.1)
# parallel             * 3.4.3     2018-01-20 local
# permute                0.9-4     2016-09-09 CRAN (R 3.4.1)
# pillar                 1.1.0     2018-01-14 CRAN (R 3.4.2)
# pkgconfig              2.0.1     2017-03-21 CRAN (R 3.4.1)
# plyr                   1.8.4     2016-06-08 CRAN (R 3.4.1)
# png                    0.1-7     2013-12-03 CRAN (R 3.4.1)
# prettyunits            1.0.2     2015-07-13 CRAN (R 3.4.1)
# progress               1.1.2     2016-12-14 CRAN (R 3.4.1)
# R.methodsS3            1.7.1     2016-02-16 CRAN (R 3.4.1)
# R.oo                   1.21.0    2016-11-01 CRAN (R 3.4.1)
# R.utils                2.6.0     2017-11-05 CRAN (R 3.4.2)
# R6                     2.2.2     2017-06-17 CRAN (R 3.4.1)
# Rcpp                   0.12.14   2017-11-23 CRAN (R 3.4.2)
# RCurl                  1.95-4.10 2018-01-04 CRAN (R 3.4.2)
# rlang                  0.1.6     2017-12-21 CRAN (R 3.4.2)
# rmote                * 0.3.4     2018-02-16 deltarho (R 3.4.3)
# RMySQL                 0.10.14   2018-02-26 CRAN (R 3.4.3)
# Rsamtools            * 1.30.0    2017-11-29 Bioconductor
# RSQLite                2.0       2017-06-19 CRAN (R 3.4.1)
# rtracklayer            1.38.3    2018-02-17 Bioconductor
# RUnit                  0.4.31    2015-11-06 CRAN (R 3.4.1)
# S4Vectors            * 0.16.0    2017-11-29 Bioconductor
# scales                 0.5.0     2017-08-24 CRAN (R 3.4.1)
# servr                  0.8       2017-11-06 CRAN (R 3.4.3)
# SGSeq                * 1.12.0    2018-01-16 Bioconductor
# stats                * 3.4.3     2018-01-20 local
# stats4               * 3.4.3     2018-01-20 local
# stringi                1.1.6     2017-11-17 CRAN (R 3.4.2)
# stringr                1.2.0     2017-02-18 CRAN (R 3.4.1)
# SummarizedExperiment * 1.8.1     2018-01-09 Bioconductor
# tibble                 1.4.1     2017-12-25 CRAN (R 3.4.2)
# tools                  3.4.3     2018-01-20 local
# utils                * 3.4.3     2018-01-20 local
# withr                  2.1.1     2017-12-19 CRAN (R 3.4.2)
# XML                    3.98-1.10 2018-02-19 CRAN (R 3.4.3)
# XVector              * 0.18.0    2017-11-29 Bioconductor
# zlibbioc               1.24.0    2017-11-07 Bioconductor
