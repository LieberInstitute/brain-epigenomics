library('GenomicRanges')
library('devtools')

models <- c('age', 'cell', 'interaction')
dmrs <- lapply(models, function(x) {
    load(paste0('/dcl01/lieber/ajaffe/lab/brain-epigenomics/bumphunting/', 
        'bumps_bsseqSmooth_Neuron_', x, '_250_perm.Rdata'), verbose = TRUE)    
    bumps_gr <- GRanges(bumps$table[bumps$table$fwer < 0.05, ])
})
names(dmrs) <- models
sapply(dmrs, length)
sapply(dmrs, function(x) sum(width(x))) / 1e6

# > sapply(dmrs, length)
#         age        cell interaction
#         129       11179        2178
# > sapply(dmrs, function(x) sum(width(x))) / 1e6
#         age        cell interaction
#    0.407420   31.114074    2.995936
   
## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()

# Session info ----------------------------------------------------------------------------------------------------------
#  setting  value
#  version  R version 3.5.0 Patched (2018-04-30 r74679)
#  system   x86_64, linux-gnu
#  ui       X11
#  language (EN)
#  collate  en_US.UTF-8
#  tz       US/Eastern
#  date     2018-09-12
#
# Packages --------------------------------------------------------------------------------------------------------------
#  package          * version   date       source
#  assertthat         0.2.0     2017-04-11 CRAN (R 3.5.0)
#  base             * 3.5.0     2018-05-02 local
#  bindr              0.1.1     2018-03-13 CRAN (R 3.5.0)
#  bindrcpp           0.2.2     2018-03-29 CRAN (R 3.5.0)
#  BiocGenerics     * 0.26.0    2018-05-03 Bioconductor
#  bitops             1.0-6     2013-08-17 CRAN (R 3.5.0)
#  colorout         * 1.2-0     2018-05-02 Github (jalvesaq/colorout@c42088d)
#  colorspace         1.3-2     2016-12-14 CRAN (R 3.5.0)
#  compiler           3.5.0     2018-05-02 local
#  crayon             1.3.4     2017-09-16 CRAN (R 3.5.0)
#  datasets         * 3.5.0     2018-05-02 local
#  devtools         * 1.13.6    2018-06-27 CRAN (R 3.5.0)
#  digest             0.6.15    2018-01-28 CRAN (R 3.5.0)
#  dplyr              0.7.6     2018-06-29 CRAN (R 3.5.0)
#  GenomeInfoDb     * 1.16.0    2018-05-03 Bioconductor
#  GenomeInfoDbData   1.1.0     2018-04-17 Bioconductor
#  GenomicRanges    * 1.32.6    2018-07-20 Bioconductor
#  ggplot2            3.0.0     2018-07-03 CRAN (R 3.5.0)
#  glue               1.3.0     2018-07-17 CRAN (R 3.5.0)
#  graphics         * 3.5.0     2018-05-02 local
#  grDevices        * 3.5.0     2018-05-02 local
#  grid               3.5.0     2018-05-02 local
#  gtable             0.2.0     2016-02-26 CRAN (R 3.5.0)
#  htmltools          0.3.6     2017-04-28 CRAN (R 3.5.0)
#  htmlwidgets        1.2       2018-04-19 CRAN (R 3.5.0)
#  httpuv             1.4.5     2018-07-19 CRAN (R 3.5.0)
#  IRanges          * 2.14.10   2018-05-17 Bioconductor
#  later              0.7.4     2018-08-31 CRAN (R 3.5.0)
#  lattice            0.20-35   2017-03-25 CRAN (R 3.5.0)
#  lazyeval           0.2.1     2017-10-29 CRAN (R 3.5.0)
#  magrittr           1.5       2014-11-22 CRAN (R 3.5.0)
#  memoise            1.1.0     2017-04-21 CRAN (R 3.5.0)
#  methods          * 3.5.0     2018-05-02 local
#  munsell            0.5.0     2018-06-12 CRAN (R 3.5.0)
#  parallel         * 3.5.0     2018-05-02 local
#  pillar             1.3.0     2018-07-14 CRAN (R 3.5.0)
#  pkgconfig          2.0.1     2017-03-21 CRAN (R 3.5.0)
#  plyr               1.8.4     2016-06-08 CRAN (R 3.5.0)
#  png                0.1-7     2013-12-03 CRAN (R 3.5.0)
#  promises           1.0.1     2018-04-13 CRAN (R 3.5.0)
#  purrr              0.2.5     2018-05-29 CRAN (R 3.5.0)
#  R6                 2.2.2     2017-06-17 CRAN (R 3.5.0)
#  Rcpp               0.12.18   2018-07-23 CRAN (R 3.5.0)
#  RCurl              1.95-4.11 2018-07-15 CRAN (R 3.5.0)
#  rlang              0.2.1     2018-05-30 cran (@0.2.1)
#  rmote            * 0.3.4     2018-05-02 deltarho (R 3.5.0)
#  S4Vectors        * 0.18.3    2018-06-13 Bioconductor
#  scales             1.0.0     2018-08-09 CRAN (R 3.5.0)
#  servr              0.10      2018-05-30 CRAN (R 3.5.0)
#  stats            * 3.5.0     2018-05-02 local
#  stats4           * 3.5.0     2018-05-02 local
#  tibble             1.4.2     2018-01-22 CRAN (R 3.5.0)
#  tidyselect         0.2.4     2018-02-26 CRAN (R 3.5.0)
#  tools              3.5.0     2018-05-02 local
#  utils            * 3.5.0     2018-05-02 local
#  withr              2.1.2     2018-03-15 CRAN (R 3.5.0)
#  xfun               0.3       2018-07-06 CRAN (R 3.5.0)
#  XVector            0.20.0    2018-05-03 Bioconductor
#  zlibbioc           1.26.0    2018-05-02 Bioconductor
