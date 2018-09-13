#* aic/bic for Amanda /dcl01/lieber/ajaffe/lab/brain-epigenomics/bumphunting/pdf/glia_vs_neuron_mean_coef.pdf

library('GenomicRanges')
library('devtools')

load('rda/limma_Neuron_CpGs_minCov_3_ageInfo_dmrs.Rdata', verbose = TRUE)
kdf <- data.frame('glia' = dmrs$age_glia_coef_mean,
    'neuron' = dmrs$age_neuron_coef_mean)
set.seed(20180201)
k2 <- kmeans(kdf, 2, nstart = 100)
k4 <- kmeans(kdf, 4, nstart = 100)
k6 <- kmeans(kdf, 6, nstart = 100)
k8 <- kmeans(kdf, 8, nstart = 100)

## From https://stackoverflow.com/questions/15839774/how-to-calculate-bic-for-k-means-clustering-in-r
## which refers to http://sherrytowers.com/2013/10/24/k-means-clustering/
## and is also implemented in https://rdrr.io/rforge/kmeansstep/man/kmeansAIC.html
kmeansAIC = function(fit){

m = ncol(fit$centers)
n = length(fit$cluster)
k = nrow(fit$centers)
D = fit$tot.withinss
return(data.frame(AIC = D + 2*m*k,
                  BIC = D + log(n)*m*k))
}

ks <- list(k2 = k2, k4 = k4, k6 = k6, k8 = k8)
k_aic <- do.call(rbind, lapply(ks, kmeansAIC))
#          AIC       BIC
# k2  8.137345  30.88199
# k4 16.048332  61.53763
# k6 24.029804  92.26375
# k8 32.021970 123.00057

set.seed(20180201)
k_s <- c(1, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 30, 40, 50)
k_many <- lapply(k_s, function(x) {
    message(paste(Sys.time(), 'processing k =', x))
    kmeans(kdf, x, nstart = 100)
})
names(k_many) <- paste0('k', k_s)
k_many_aic <- do.call(rbind, lapply(k_many, kmeansAIC))
k_many_aic
# > k_many_aic
#            AIC       BIC
# k1    4.342512  15.71484
# k2    8.137345  30.88199
# k4   16.048332  61.53763
# k6   24.029804  92.26375
# k8   32.021970 123.00057
# k10  40.018555 153.74180
# k12  48.015503 184.48340
# k14  56.013346 215.22589
# k16  64.011768 245.96896
# k18  72.010672 276.71251
# k20  80.009725 307.45622
# k30 120.006859 461.17660
# k40 160.005190 614.89817
# k50 200.004220 768.62045

kdf_cent <- kdf
kdf_cent$glia <- (kdf_cent$glia - mean(kdf_cent$glia)) / sd(kdf_cent$glia)
kdf_cent$neuron <- (kdf_cent$neuron - mean(kdf_cent$neuron)) / sd(kdf_cent$neuron)
set.seed(20180201)
k_many_cent <- lapply(k_s, function(x) {
    message(paste(Sys.time(), 'processing k =', x))
    kmeans(kdf_cent, x, nstart = 100)
})
names(k_many_cent) <- paste0('k', k_s)
k_many_aic_cent <- do.call(rbind, lapply(k_many_cent, kmeansAIC))
k_many_aic_cent
# > k_many_aic_cent
#           AIC       BIC
# k1  4358.0000 4369.3723
# k2  1864.8178 1887.5625
# k4   661.1956  706.6849
# k6   421.0191  489.2530
# k8   317.7236  408.7022
# k10  283.8193  397.5425
# k12  253.7665  390.2344
# k14  231.9630  391.1756
# k16  219.0774  401.0346
# k18  211.7903  416.4922
# k20  207.5166  434.9631
# k30  208.7303  549.9001
# k40  228.6892  683.5822
# k50  255.7656  824.3818


k_sum <- k_many_aic_cent
colnames(k_sum) <- paste0(colnames(k_sum), '_cent')
k_sum <- cbind(k_many_aic, k_sum)
k_sum$k <- as.integer(gsub('k', '', rownames(k_sum)))

library(ggplot2)

pdf('pdf/AIC_BIC_kmeans.pdf', useDingbats = FALSE)
ggplot(k_sum, aes(x = k, y = AIC)) + geom_point() + geom_line() + theme_bw(base_size = 30)
ggplot(k_sum, aes(x = k, y = AIC_cent)) + geom_point() + geom_line() + theme_bw(base_size = 30)
ggplot(k_sum, aes(x = k, y = BIC)) + geom_point() + geom_line() + theme_bw(base_size = 30)
ggplot(k_sum, aes(x = k, y = BIC_cent)) + geom_point() + geom_line() + theme_bw(base_size = 30)
dev.off()

## k = 6 cluster labels don't match but that's ok
identical(k_many$k6$cluster, k6$cluster)
# [1] FALSE
## since the numbers per group do match
all(tapply(k_many$k6$cluster, k6$cluster, table) - table(k6$cluster) == 0)
# [1] TRUE

## Most k = 6 (centered) vs k = 6 (raw) groups overlap by quite a bit
tapply(k_many_cent$k6$cluster, k_many$k6$cluster, table)

# $`1`
#
#  5  6
# 94  2
#
# $`2`
#
#   2   4   5   6
#   1   4   2 163
#
# $`3`
#
#   1   2
# 416   7
#
# $`4`
#
#   4   6
# 361   1
#
# $`5`
#
#   1   2   3
#  27 534   6
#
# $`6`
#
#   2   3
#  18 542

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
#  ggplot2          * 3.0.0     2018-07-03 CRAN (R 3.5.0)
#  glue               1.3.0     2018-07-17 CRAN (R 3.5.0)
#  graphics         * 3.5.0     2018-05-02 local
#  grDevices        * 3.5.0     2018-05-02 local
#  grid               3.5.0     2018-05-02 local
#  gtable             0.2.0     2016-02-26 CRAN (R 3.5.0)
#  htmltools          0.3.6     2017-04-28 CRAN (R 3.5.0)
#  htmlwidgets        1.2       2018-04-19 CRAN (R 3.5.0)
#  httpuv             1.4.5     2018-07-19 CRAN (R 3.5.0)
#  IRanges          * 2.14.10   2018-05-17 Bioconductor
#  labeling           0.3       2014-08-23 CRAN (R 3.5.0)
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
#
