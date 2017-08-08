library('bsseq')
library('GenomicRanges')
library('devtools')

## CpGs first
load('gr_cps_minCov_3.Rdata')
load('BSobj_lister_minCov_3.Rdata')

gr <- rowRanges(BSobj)

Cov <- M <- matrix(0, nrow = length(gr_cpgs), ncol = ncol(BSobj))
colnames(Cov) <- colnames(M) <- colnames(BSobj)

ov <- findOverlaps(gr_cpgs, rowRanges(BSobj))
M[queryHits(ov), ] <- assays(BSobj)$M[subjectHits(ov), ]
Cov[queryHits(ov), ] <- assays(BSobj)$Cov[subjectHits(ov), ]

BSobj <- BSseq(M = M, Cov = Cov, gr = gr_cpgs)
save(BSobj, file = 'BSobj_matched_lister_minCov_3.Rdata')

rm(BSobj, gr, Cov, M, ov, gr_cpgs)


## nonCG next
load('gr_all_highCov.Rdata')
load('allChrs_lister_nonCG_highCov.Rdata')

gr <- rowRanges(BSobj)

Cov <- M <- matrix(0, nrow = length(gr_all_highCov), ncol = ncol(BSobj))
colnames(Cov) <- colnames(M) <- colnames(BSobj)

ov <- findOverlaps(gr_all_highCov, rowRanges(BSobj))
M[queryHits(ov), ] <- assays(BSobj)$M[subjectHits(ov), ]
Cov[queryHits(ov), ] <- assays(BSobj)$Cov[subjectHits(ov), ]

BSobj <- BSseq(M = M, Cov = Cov, gr = gr_all_highCov)
save(BSobj, file = 'allChrs_matched_lister_nonCG_highCov.Rdata')

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
#  date     2017-08-08
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
#  devtools             * 1.13.2   2017-06-02 CRAN (R 3.3.3)
#  digest                 0.6.12   2017-01-27 CRAN (R 3.3.1)
#  fcuk                 * 0.1.21   2017-07-08 CRAN (R 3.3.3)
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
#  magrittr               1.5      2014-11-22 CRAN (R 3.3.1)
#  Matrix                 1.2-8    2017-01-20 CRAN (R 3.3.3)
#  matrixStats            0.52.2   2017-04-14 CRAN (R 3.3.1)
#  memoise                1.1.0    2017-04-21 CRAN (R 3.3.1)
#  methods              * 3.3.3    2017-05-18 local
#  munsell                0.4.3    2016-02-13 CRAN (R 3.3.0)
#  parallel             * 3.3.3    2017-05-18 local
#  permute                0.9-4    2016-09-09 CRAN (R 3.3.1)
#  plyr                   1.8.4    2016-06-08 CRAN (R 3.3.1)
#  purrr                  0.2.2.2  2017-05-11 CRAN (R 3.3.1)
#  R.methodsS3            1.7.1    2016-02-16 CRAN (R 3.3.0)
#  R.oo                   1.21.0   2016-11-01 CRAN (R 3.3.1)
#  R.utils                2.5.0    2016-11-07 CRAN (R 3.3.1)
#  Rcpp                   0.12.11  2017-05-22 CRAN (R 3.3.3)
#  RCurl                  1.95-4.8 2016-03-01 CRAN (R 3.3.0)
#  rlang                  0.1.1    2017-05-18 CRAN (R 3.3.1)
#  S4Vectors            * 0.12.2   2017-06-05 Bioconductor
#  scales                 0.4.1    2016-11-09 CRAN (R 3.3.1)
#  stats                * 3.3.3    2017-05-18 local
#  stats4               * 3.3.3    2017-05-18 local
#  stringdist             0.9.4.4  2016-12-16 cran (@0.9.4.4)
#  SummarizedExperiment * 1.4.0    2017-06-05 Bioconductor
#  tibble                 1.3.3    2017-05-28 CRAN (R 3.3.3)
#  tools                  3.3.3    2017-05-18 local
#  utils                * 3.3.3    2017-05-18 local
#  withr                  1.0.2    2016-06-20 CRAN (R 3.3.1)
#  XVector                0.14.1   2017-03-21 Bioconductor
#  zlibbioc               1.20.0   2016-10-20 Bioconductor