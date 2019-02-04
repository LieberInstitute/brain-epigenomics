library('bsseq')
library('GenomicRanges')
library('bumphunter')
library('sessioninfo')
library('SummarizedExperiment')

load('/dcl01/lieber/ajaffe/lab/brain-epigenomics/bumphunting/BSobj_bsseqSmooth_Neuron_minCov_3.Rdata', verbose = TRUE)

pd <- pData(BSobj)

## Add lambda info
## adapted code from 
## https://github.com/LieberInstitute/brain-epigenomics/blob/master/lambda/check_lambda.R
lambda <- read.csv('/dcl01/lieber/WGBS/LIBD_Data/lambda_alignment/MasterList_lambda.csv')
m_lambda <- match(pd$Data.ID, lambda$WGC.ID)
stopifnot(!any(is.na(m_lambda)))
pd$avg_conv_efficiency <- colData(BSobj)$avg_conv_efficiency <- lambda$avg_conv_efficiency[m_lambda]
pd$ratio_trimmed <- colData(BSobj)$ratio_trimmed <- pd$total_num_trimmed_reads/pd$total_num_untrimmed_reads

models <- c('interaction', 'cell', 'age')
# model <- models[2]
rse_mean_meth <- lapply(models, function(model) {
    f <- paste0('/dcl01/lieber/ajaffe/lab/brain-epigenomics/bumphunting/bumps_bsseqSmooth_Neuron_', model, '_250_perm.Rdata')
    message(paste(Sys.time(), 'loading', f))
    load(f, verbose = TRUE)

    ## Build a GRanges object from the table output of bumphunter
    bumps_gr <- GRanges(bumps$table[bumps$table$fwer < 0.05, ])

    ## Next find which bases are part of the DMRs
    idx <- countOverlaps(rowRanges(BSobj), bumps_gr) > 0

    ## Get the methylation only for the bases that are in the DMRs
    ## Use the smoothed values, since that's what we used for
    ## finding the DMRs as in
    ## https://github.com/LieberInstitute/brain-epigenomics/blob/master/bumphunting/find_bumps_bsseqSmooth.R#L124
    message(paste(Sys.time(), 'extracting the methylation values'))
    meth <- getMeth(BSobj[idx, ], type = 'smooth')


    message(paste(Sys.time(), 'computing means across the CpGs in the DMRs'))
    ov <- findOverlaps(bumps_gr, rowRanges(BSobj[idx, ]))
    ov_dmr <- split(subjectHits(ov), queryHits(ov))

    ## Finally compute the mean methylation for the bases in the DMR for the glia samples
    dmr_mean_meth <- do.call(rbind, lapply(ov_dmr, function(dmr_bases) {
        colMeans(meth[dmr_bases, , drop = FALSE])
    }))

    stopifnot(length(bumps_gr) == nrow(dmr_mean_meth))
    rownames(dmr_mean_meth) <- names(bumps_gr)
    final <- SummarizedExperiment(assays = list('mean_meth' = dmr_mean_meth), rowRanges = bumps_gr, colData = pData(BSobj))
    return(final)
})
names(rse_mean_meth) <- models
save(rse_mean_meth, file = 'rdas/rse_mean_meth.Rdata')


pd <- colData(rse_mean_meth[[1]])
covs <- c('Race', 'Sex', 'PMI', 'pH', 'avg.Cov', 'Percent.Duplication', 'alignment.efficiency', 'total_num_trimmed_reads', 'total_num_untrimmed_reads', 'ratio_trimmed', 'avg_conv_efficiency')

dmr_cov <- lapply(models, function(model) {
    
    result <- do.call(cbind, lapply(covs, function(cov) {
    # cov <- covs[1]
        coef <- 2
        if(model == 'cell') {
            design <- model.matrix(~ ., data = pd[, c('Cell.Type', 'Age', cov)])
        } else if (model == 'age') {
            design <- model.matrix(~ ., data = pd[, c('Age', 'Cell.Type', cov)])
        } else if (model == 'interaction') {
            design <- model.matrix(~ ~ Age * Cell.Type + ., data = pd[, c('Age', 'Cell.Type', cov)])
            coef <- grep(':', colnames(design))
        }
        
        message(paste(
            Sys.time(), 'extracting the coefficient for',
            colnames(design)[coef], 'after adjusting for',
            paste(colnames(design)[-c(1, coef)], collapse = ', ')
        ))
        
        rawBeta <- bumphunter:::.getEstimate(mat = assays(rse_mean_meth[[model]])$mean_meth, design = design, coef = coef, full = FALSE)
        
        ## If cov = NULL, with model = 'cell'
        ## then
        # > table(rowRanges(rse_mean_meth[[model]])$value - rawBeta[, 1] < 1e-06)
        #
        #  TRUE
        # 11179
        ## this is based on https://github.com/rafalab/bumphunter/blob/master/R/bumphunter.R#L104-L105
        ## which is the piece of code that gets run when we compute the DMRs with
        ## https://github.com/LieberInstitute/brain-epigenomics/blob/master/bumphunting/find_bumps_bsseqSmooth.R#L137-L141
        
        colnames(rawBeta) <- cov
        return(rawBeta)

    }))
    
    ## Add the original beta
    result <- cbind(result, 'original' = rowRanges(rse_mean_meth[[model]])$value)
    return(result)
})
names(dmr_cov) <- models
save(dmr_cov, file = 'rdas/dmr_cov.Rdata')

sapply(dmr_cov, dim)
# [1,]        2178 11179 129
# [2,]          12    12  12
head(dmr_cov$interaction)
#              Race         Sex         PMI          pH     avg.Cov
# 8307   0.01590095  0.01587469  0.01626441  0.01632416  0.01655542
# 8306   0.01596907  0.01589940  0.01641649  0.01643977  0.01652040
# 5886   0.01353393  0.01314937  0.01382809  0.01397798  0.01353401
# 26011 -0.01586735 -0.01582416 -0.01606460 -0.01490312 -0.01511349
# 19972 -0.01798947 -0.01804704 -0.01811076 -0.01699670 -0.01710066
# 22563 -0.01682682 -0.01662047 -0.01663662 -0.01592537 -0.01561785
#       Percent.Duplication alignment.efficiency total_num_trimmed_reads
# 8307           0.01617616           0.01698547              0.01659117
# 8306           0.01622861           0.01718333              0.01656379
# 5886           0.01359882           0.01462163              0.01314812
# 26011         -0.01516003          -0.01473250             -0.01503604
# 19972         -0.01724533          -0.01681139             -0.01721672
# 22563         -0.01620978          -0.01593187             -0.01595356
#       total_num_untrimmed_reads ratio_trimmed avg_conv_efficiency    original
# 8307                 0.01640517    0.01608386          0.01572145  0.01634005
# 8306                 0.01641839    0.01624832          0.01557559  0.01637278
# 5886                 0.01372937    0.01200527          0.01260341  0.01379249
# 26011               -0.01509268   -0.01572078         -0.01432257 -0.01515396
# 19972               -0.01726911   -0.01722871         -0.01641404 -0.01730097
# 22563               -0.01603825   -0.01625742         -0.01593754 -0.01608374

## Make all the plots
dir.create('pdf', showWarnings = FALSE)
dmr_cov_ratio <- lapply(models, function(model) {
    pdf(paste0('pdf/dmr_covariates_model_', model, '.pdf'), useDingbats = FALSE)
    par(cex.lab=1.2, cex.axis=1.2)
    result <- do.call(cbind, lapply(covs, function(cov) {
        dat_or <- dmr_cov[[model]][, 'original']
        dat_cov <- dmr_cov[[model]][, cov]
        dat_ratio <- dat_cov / dat_or
        dat_tab <- addmargins(table(
            sign(dat_or),
            sign(dat_cov),
            dnn = c('Original', cov)
        ))
        perc_agree <- round(sum(diag(dat_tab)[-3]) / diag(dat_tab)[3] * 100, 2)
        # print(dat_tab)
        
        plot(x = dat_or, y = dat_cov,
            xlab = paste('Original beta:', model, 'model'),
            ylab = paste('Adjusted beta for', cov),
            sub = paste0('Agreement by sign: ', perc_agree, '%'),
            pch = 20, col = scales::alpha('black', 1/5)
        )
        abline(a = 0, b = 1, col = 'red')
        plot(x = dat_or, y = dat_ratio,
            xlab = paste('Original beta:', model, 'model'),
            ylab = paste('Ratio of betas:', cov, 'adj. / original'),
            pch = 20, col = scales::alpha('black', 1/5)
        )
        abline(h = 1, col = 'red')
        
        res <- data.frame(dat_ratio)
        colnames(res) <- cov
        return(res)
    }))
    
    par(mar=c(16, 5, 5, 2), cex.main = 1.3, cex.lab=1.2, cex.axis=1.2)
    boxplot(result, las = 2,
        ylab = 'Ratio of betas: adjusted / original',
        main = paste('Model:', model))
    abline(h = 1, col = 'red')
    
    dev.off()
    return(result)
})
names(dmr_cov_ratio) <- models
save(dmr_cov_ratio, file = 'rdas/dmr_cov_ratio.Rdata')

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()

# ─ Session info ───────────────────────────────────────────────────────────────────────────────────────────────────────
#  setting  value
#  version  R version 3.5.1 Patched (2018-10-29 r75535)
#  os       Red Hat Enterprise Linux Server release 6.9 (Santiago)
#  system   x86_64, linux-gnu
#  ui       X11
#  language (EN)
#  collate  en_US.UTF-8
#  ctype    en_US.UTF-8
#  tz       US/Eastern
#  date     2019-02-04
#
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
#  package              * version   date       lib source
#  AnnotationDbi          1.44.0    2018-10-30 [1] Bioconductor
#  assertthat             0.2.0     2017-04-11 [2] CRAN (R 3.5.0)
#  bibtex                 0.4.2     2017-06-30 [1] CRAN (R 3.5.0)
#  bindr                  0.1.1     2018-03-13 [1] CRAN (R 3.5.0)
#  bindrcpp               0.2.2     2018-03-29 [1] CRAN (R 3.5.0)
#  Biobase              * 2.42.0    2018-10-30 [2] Bioconductor
#  BiocGenerics         * 0.28.0    2018-10-30 [1] Bioconductor
#  BiocParallel         * 1.16.5    2019-01-04 [1] Bioconductor
#  biomaRt                2.38.0    2018-10-30 [1] Bioconductor
#  Biostrings             2.50.2    2019-01-03 [1] Bioconductor
#  bit                    1.1-14    2018-05-29 [2] CRAN (R 3.5.0)
#  bit64                  0.9-7     2017-05-08 [2] CRAN (R 3.5.0)
#  bitops                 1.0-6     2013-08-17 [2] CRAN (R 3.5.0)
#  blob                   1.1.1     2018-03-25 [2] CRAN (R 3.5.0)
#  BSgenome               1.50.0    2018-10-30 [1] Bioconductor
#  bsseq                * 1.18.0    2018-10-30 [1] Bioconductor
#  bumphunter           * 1.24.5    2018-12-01 [1] Bioconductor
#  cli                    1.0.1     2018-09-25 [1] CRAN (R 3.5.1)
#  codetools              0.2-15    2016-10-05 [3] CRAN (R 3.5.1)
#  colorout             * 1.2-0     2018-05-02 [1] Github (jalvesaq/colorout@c42088d)
#  colorspace             1.4-0     2019-01-13 [2] CRAN (R 3.5.1)
#  crayon                 1.3.4     2017-09-16 [1] CRAN (R 3.5.0)
#  data.table             1.12.0    2019-01-13 [1] CRAN (R 3.5.1)
#  DBI                    1.0.0     2018-05-02 [2] CRAN (R 3.5.0)
#  DelayedArray         * 0.8.0     2018-10-30 [2] Bioconductor
#  DelayedMatrixStats     1.4.0     2018-10-30 [2] Bioconductor
#  digest                 0.6.18    2018-10-10 [1] CRAN (R 3.5.1)
#  doRNG                  1.7.1     2018-06-22 [2] CRAN (R 3.5.1)
#  dplyr                  0.7.8     2018-11-10 [1] CRAN (R 3.5.1)
#  foreach              * 1.4.4     2017-12-12 [2] CRAN (R 3.5.0)
#  GenomeInfoDb         * 1.18.1    2018-11-12 [1] Bioconductor
#  GenomeInfoDbData       1.2.0     2018-11-02 [2] Bioconductor
#  GenomicAlignments      1.18.1    2019-01-04 [1] Bioconductor
#  GenomicFeatures        1.34.1    2018-11-03 [1] Bioconductor
#  GenomicRanges        * 1.34.0    2018-10-30 [1] Bioconductor
#  ggplot2                3.1.0     2018-10-25 [1] CRAN (R 3.5.1)
#  glue                   1.3.0     2018-07-17 [1] CRAN (R 3.5.1)
#  gtable                 0.2.0     2016-02-26 [2] CRAN (R 3.5.0)
#  gtools                 3.8.1     2018-06-26 [2] CRAN (R 3.5.1)
#  HDF5Array              1.10.1    2018-12-05 [1] Bioconductor
#  hms                    0.4.2     2018-03-10 [2] CRAN (R 3.5.0)
#  htmltools              0.3.6     2017-04-28 [2] CRAN (R 3.5.0)
#  htmlwidgets            1.3       2018-09-30 [1] CRAN (R 3.5.1)
#  httpuv                 1.4.5.1   2018-12-18 [2] CRAN (R 3.5.1)
#  httr                   1.4.0     2018-12-11 [1] CRAN (R 3.5.1)
#  IRanges              * 2.16.0    2018-10-30 [1] Bioconductor
#  iterators            * 1.0.10    2018-07-13 [2] CRAN (R 3.5.1)
#  later                  0.7.5     2018-09-18 [2] CRAN (R 3.5.1)
#  lattice                0.20-38   2018-11-04 [3] CRAN (R 3.5.1)
#  lazyeval               0.2.1     2017-10-29 [2] CRAN (R 3.5.0)
#  limma                  3.38.3    2018-12-02 [1] Bioconductor
#  locfit               * 1.5-9.1   2013-04-20 [2] CRAN (R 3.5.0)
#  magrittr               1.5       2014-11-22 [1] CRAN (R 3.5.0)
#  Matrix                 1.2-15    2018-11-01 [3] CRAN (R 3.5.1)
#  matrixStats          * 0.54.0    2018-07-23 [1] CRAN (R 3.5.1)
#  memoise                1.1.0     2017-04-21 [2] CRAN (R 3.5.0)
#  mime                   0.6       2018-10-05 [1] CRAN (R 3.5.1)
#  munsell                0.5.0     2018-06-12 [2] CRAN (R 3.5.0)
#  permute                0.9-4     2016-09-09 [2] CRAN (R 3.5.0)
#  pillar                 1.3.1     2018-12-15 [1] CRAN (R 3.5.1)
#  pkgconfig              2.0.2     2018-08-16 [1] CRAN (R 3.5.1)
#  pkgmaker               0.27      2018-05-25 [2] CRAN (R 3.5.0)
#  plyr                   1.8.4     2016-06-08 [2] CRAN (R 3.5.0)
#  png                    0.1-7     2013-12-03 [2] CRAN (R 3.5.0)
#  prettyunits            1.0.2     2015-07-13 [1] CRAN (R 3.5.0)
#  progress               1.2.0     2018-06-14 [1] CRAN (R 3.5.1)
#  promises               1.0.1     2018-04-13 [2] CRAN (R 3.5.0)
#  purrr                  0.2.5     2018-05-29 [2] CRAN (R 3.5.0)
#  R.methodsS3            1.7.1     2016-02-16 [2] CRAN (R 3.5.0)
#  R.oo                   1.22.0    2018-04-22 [2] CRAN (R 3.5.0)
#  R.utils                2.7.0     2018-08-27 [1] CRAN (R 3.5.1)
#  R6                     2.3.0     2018-10-04 [2] CRAN (R 3.5.1)
#  Rcpp                   1.0.0     2018-11-07 [1] CRAN (R 3.5.1)
#  RCurl                  1.95-4.11 2018-07-15 [2] CRAN (R 3.5.1)
#  registry               0.5       2017-12-03 [2] CRAN (R 3.5.0)
#  rhdf5                  2.26.2    2019-01-02 [2] Bioconductor
#  Rhdf5lib               1.4.2     2018-12-03 [1] Bioconductor
#  rlang                  0.3.1     2019-01-08 [1] CRAN (R 3.5.1)
#  rmote                * 0.3.4     2018-05-02 [1] deltarho (R 3.5.0)
#  rngtools               1.3.1     2018-05-15 [2] CRAN (R 3.5.0)
#  Rsamtools              1.34.0    2018-10-30 [1] Bioconductor
#  RSQLite                2.1.1     2018-05-06 [2] CRAN (R 3.5.0)
#  rtracklayer            1.42.1    2018-11-22 [1] Bioconductor
#  S4Vectors            * 0.20.1    2018-11-09 [1] Bioconductor
#  scales                 1.0.0     2018-08-09 [2] CRAN (R 3.5.1)
#  servr                  0.11      2018-10-23 [1] CRAN (R 3.5.1)
#  sessioninfo          * 1.1.1     2018-11-05 [1] CRAN (R 3.5.1)
#  stringi                1.2.4     2018-07-20 [2] CRAN (R 3.5.1)
#  stringr                1.3.1     2018-05-10 [1] CRAN (R 3.5.0)
#  SummarizedExperiment * 1.12.0    2018-10-30 [1] Bioconductor
#  tibble                 2.0.1     2019-01-12 [1] CRAN (R 3.5.1)
#  tidyselect             0.2.5     2018-10-11 [2] CRAN (R 3.5.1)
#  withr                  2.1.2     2018-03-15 [2] CRAN (R 3.5.0)
#  xfun                   0.4       2018-10-23 [1] CRAN (R 3.5.1)
#  XML                    3.98-1.16 2018-08-19 [2] CRAN (R 3.5.1)
#  xtable                 1.8-3     2018-08-29 [2] CRAN (R 3.5.1)
#  XVector                0.22.0    2018-10-30 [1] Bioconductor
#  zlibbioc               1.28.0    2018-10-30 [2] Bioconductor
#
# [1] /users/lcollado/R/x86_64-pc-linux-gnu-library/3.5.x
# [2] /jhpce/shared/jhpce/core/conda/miniconda-3/envs/svnR-3.5.x/R/3.5.x/lib64/R/site-library
# [3] /jhpce/shared/jhpce/core/conda/miniconda-3/envs/svnR-3.5.x/R/3.5.x/lib64/R/library
