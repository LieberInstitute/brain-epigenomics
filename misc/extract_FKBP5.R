library('bsseq')
library('GenomicRanges')
library('devtools')

###### CpG

## Load CpG data
load(file.path('/dcl01/lieber/ajaffe/lab/brain-epigenomics', 'bumphunting',
    'BSobj_bsseqSmooth_Neuron_minCov_3.Rdata'))
## These are the CpGs with a minimum coverage of 3 in all samples.
dim(BSobj)
# [1] 18664892       32
   
regions <- GRanges(c('chr6'),
   IRanges(
       c(35541362),
       c(35696397)
   )
)
names(regions) <- c('FKBP5')

regions_long <- resize(regions, width(regions) + 2 * 20000, fix = 'center')

# > regions_long
# GRanges object with 1 range and 0 metadata columns:
#         seqnames               ranges strand
#            <Rle>            <IRanges>  <Rle>
#   FKBP5     chr6 [35521362, 35716397]      *
#   -------
#   seqinfo: 1 sequence from an unspecified genome; no seqlengths

## subset CpG data
idx_cpg <- countOverlaps(rowRanges(BSobj), regions_long) > 0
BSobj_CpG <- BSobj[idx_cpg, ]

# > dim(BSobj_CpG)
# [1] 1884   32

##
models <- c('interaction', 'cell', 'age')
dmrs <- lapply(models, function(x) {
    load(file.path('/dcl01/lieber/ajaffe/lab/brain-epigenomics', 'bumphunting',
        paste0('bumps_bsseqSmooth_Neuron_', x, '_250_perm.Rdata')))
    return(bumps)
})
names(dmrs) <- models

extract_info <- function(dmr, name, idx) {
    res <- data.frame(
        coef = dmr$coef[idx, ],
        pvaluesMarginal = dmr$pvaluesMarginal[idx]
    )
    colnames(res) <- paste0(name, '_', colnames(res))
    return(res)    
}

mcols(rowRanges(BSobj_CpG)) <- do.call(cbind, mapply(extract_info, dmrs, models, MoreArgs = list(idx = idx_cpg), SIMPLIFY = FALSE, USE.NAMES = FALSE))

save(BSobj_CpG, file = 'rdas/BSobj_CpG.Rdata')

extract_dmrs <- function(dmr, name) {
    dmr <- GRanges(dmr$table)
    res <- dmr[countOverlaps(dmr, regions_long) > 0]
    res$model <- name
    return(res)    
}

DMRs_CpG <- GRangesList(mapply(extract_dmrs, dmrs, models, SIMPLIFY = FALSE))

save(DMRs_CpG, file = 'rdas/DMRs_CpG.Rdata')


## Load matching limma output and BSobj data
load(file.path('/dcl01/lieber/ajaffe/lab/brain-epigenomics', 'bumphunting',
    'BSobj_Neuron.Rdata'))
## These are CpGs with minimum coverage of 1 in all samples
dim(BSobj)
# [1] 24562065       32
load(file.path('/dcl01/lieber/ajaffe/lab/brain-epigenomics', 'bumphunting',
    'limma_exploration_Neuron.Rdata'))
    
idx_limma_cpg <- countOverlaps(rowRanges(BSobj), regions_long) > 0
sum(idx_limma_cpg)
# [1] 2371
extract_limma <- function(limma, idx) {
    
    res <- limma
    res$coefficients <- res$coefficients[idx, ]
    res$df.residual <- res$df.residual[idx]
    res$sigma <- res$sigma[idx]
    res$stdev.unscaled <- res$stdev.unscaled[idx, ]
    res$Amean <- res$Amean[idx]
    ## Add the genome information too
    res$genome_position <- rowRanges(BSobj)[idx]
    return(res)
}

limma_CpG <- lapply(fits, extract_limma, idx = idx_limma_cpg)[c('age', 'interaction')]
save(limma_CpG, file = 'rdas/limma_CpG.Rdata')


###### Non-CpG


## Load non-CpG data
load(file.path('/dcl01/lieber/ajaffe/lab/brain-epigenomics',
    'bsseq/bsobj_by_chr',
    'limma_exploration_nonCG_highCov.Rdata'))
load(file.path('/dcl01/lieber/ajaffe/lab/brain-epigenomics',
    'bsseq/bsobj_by_chr',
    'allChrs_postNatal_cleaned_nonCG_noHomogenate_highCov.Rdata'))
## These are the non-CpGs with a minimum coverage of 5 in all samples.
dim(BSobj)
# [1] 58109566       32

## Reproduce methylation filter used for the limma exploration:
## at least 5 samples have to have methylation greater than 0
meth.g0 <- getMeth(BSobj, type = 'raw') > 0
meth.filt <- rowSums(meth.g0) >= 5
BSobj <- BSobj[meth.filt, ]
dim(BSobj)
# [1] 40818742       32
idx_limma_non_cpg <- countOverlaps(rowRanges(BSobj), regions_long) > 0
sum(idx_limma_non_cpg)
# [1] 3055
limma_nonCpG <- lapply(fits, extract_limma, idx = idx_limma_non_cpg)[c('age', 'interaction')]
save(limma_nonCpG, file = 'rdas/limma_nonCpG.Rdata')

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()

# Session info ----------------------------------------------------------------------------------------------------------
#  setting  value
#  version  R version 3.3.3 Patched (2017-03-15 r72696)
#  system   x86_64, linux-gnu
#  ui       X11
#  language (EN)
#  collate  en_US.UTF-8
#  tz       <NA>
#  date     2017-10-05
#
# Packages --------------------------------------------------------------------------------------------------------------
#  package              * version  date       source
#  base                 * 3.3.3    2017-05-18 local
#  Biobase              * 2.34.0   2016-10-20 Bioconductor
#  BiocGenerics         * 0.20.0   2017-06-05 Bioconductor
#  bitops                 1.0-6    2013-08-17 CRAN (R 3.3.0)
#  bsseq                * 1.10.0   2016-10-20 Bioconductor
#  colorout             * 1.1-2    2016-07-22 Github (jalvesaq/colorout@6d84420)
#  colorspace             1.3-2    2016-12-14 CRAN (R 3.3.1)
#  data.table             1.10.4   2017-02-01 CRAN (R 3.3.1)
#  datasets             * 3.3.3    2017-05-18 local
#  devtools             * 1.13.3   2017-08-02 CRAN (R 3.3.3)
#  digest                 0.6.12   2017-01-27 CRAN (R 3.3.1)
#  GenomeInfoDb         * 1.10.3   2017-02-14 Bioconductor
#  GenomicRanges        * 1.26.4   2017-06-05 Bioconductor
#  graphics             * 3.3.3    2017-05-18 local
#  grDevices            * 3.3.3    2017-05-18 local
#  grid                   3.3.3    2017-05-18 local
#  gtools                 3.5.0    2015-05-29 CRAN (R 3.3.0)
#  IRanges              * 2.8.2    2017-06-05 Bioconductor
#  lattice                0.20-34  2016-09-06 CRAN (R 3.3.3)
#  limma                * 3.30.13  2017-03-21 Bioconductor
#  locfit                 1.5-9.1  2013-04-20 CRAN (R 3.3.0)
#  Matrix                 1.2-8    2017-01-20 CRAN (R 3.3.3)
#  matrixStats            0.52.2   2017-04-14 CRAN (R 3.3.1)
#  memoise                1.1.0    2017-04-21 CRAN (R 3.3.1)
#  methods              * 3.3.3    2017-05-18 local
#  munsell                0.4.3    2016-02-13 CRAN (R 3.3.0)
#  parallel             * 3.3.3    2017-05-18 local
#  permute                0.9-4    2016-09-09 CRAN (R 3.3.1)
#  plyr                   1.8.4    2016-06-08 CRAN (R 3.3.1)
#  R.methodsS3            1.7.1    2016-02-16 CRAN (R 3.3.0)
#  R.oo                   1.21.0   2016-11-01 CRAN (R 3.3.1)
#  R.utils                2.5.0    2016-11-07 CRAN (R 3.3.1)
#  Rcpp                   0.12.12  2017-07-15 CRAN (R 3.3.3)
#  RCurl                  1.95-4.8 2016-03-01 CRAN (R 3.3.0)
#  S4Vectors            * 0.12.2   2017-06-05 Bioconductor
#  scales                 0.5.0    2017-08-24 CRAN (R 3.3.3)
#  stats                * 3.3.3    2017-05-18 local
#  stats4               * 3.3.3    2017-05-18 local
#  SummarizedExperiment * 1.4.0    2017-06-05 Bioconductor
#  tools                  3.3.3    2017-05-18 local
#  utils                * 3.3.3    2017-05-18 local
#  withr                  2.0.0    2017-07-28 CRAN (R 3.3.3)
#  XVector                0.14.1   2017-03-21 Bioconductor
#  zlibbioc               1.20.0   2016-10-20 Bioconductor
