library('devtools')
library('ggplot2')

files <- dir(pattern = '^auto_long_chr')

load_auto <- function(f) {
    chr <- gsub('_.*', '', gsub('auto_long_', '', f))
    message(paste(Sys.time(), 'loading', f))
    load(f)
    auto_long$chr <- chr
    return(auto_long)
}

if(!file.exists('auto_long_combined.Rdata')) {
    auto_long <- do.call(rbind, lapply(files, load_auto))
    save(auto_long, file = 'auto_long_combined.Rdata')
} else {
    load('auto_long_combined.Rdata')
}
dim(auto_long)
<<<<<<< HEAD
=======
# [1] 7090476       8
>>>>>>> c342a3e8636f5c1dab8031343885425bc57c00f7

## This never finished running:
# context_summary <- lapply(unique(auto_long)$context, function(context) {
#     sub <- auto_long[auto_long$context == context, ]
#     tapply(sub$acf, sub$lag, summary)
# })
# context_summary
# save(context_summary, file = 'autocorrelation_summary_by_context.Rdata')

<<<<<<< HEAD
auto_long$lag <- as.factor(auto_long$lag)

png('autocorrelation_by_context.png', width = 480 * 2)
ggplot(auto_long, aes(x = lag, y = acf)) + geom_boxplot() + facet_grid(. ~ context) + theme_bw(base_size = 14)
dev.off()

png('autocorrelation_by_context_abs.png', width = 480 * 2)
ggplot(auto_long, aes(x = lag, y = abs(acf))) + geom_boxplot() + facet_grid(. ~ context) + theme_bw(base_size = 14)
dev.off()

=======
dir.create('pdf', showWarnings = FALSE)
pdf('pdf/autocorrelation_by_context.pdf', width = 7 * 2, height = 7 * 2)
ggplot(auto_long, aes(x = lag, y = acf_neuron)) + geom_boxplot() + facet_grid(. ~ context) + theme_bw(base_size = 30) + ggtitle('Neuron') + ylab('Auto correlation') + ylim(c(-1, 1))

ggplot(auto_long, aes(x = lag, y = acf_glia)) + geom_boxplot() + facet_grid(. ~ context) + theme_bw(base_size = 30) + ggtitle('Glia') + ylab('Auto correlation')  + ylim(c(-1, 1))
dev.off()

pdf('pdf/autocorrelation_by_context_abs.pdf', width = 7 * 2, height = 7 * 2)
ggplot(auto_long, aes(x = lag, y = abs(acf_neuron))) + geom_boxplot() + facet_grid(. ~ context) + theme_bw(base_size = 30) + ggtitle('Neuron') + ylab('Absolute auto correlation') + ylim(c(0, 1))

ggplot(auto_long, aes(x = lag, y = abs(acf_glia))) + geom_boxplot() + facet_grid(. ~ context) + theme_bw(base_size = 30) + ggtitle('Glia') + ylab('Absolute auto correlation') + ylim(c(0, 1))
dev.off()


## chr21 files
files <- dir('rda', pattern = '^custom_acf_chr')
load_custom <- function(f) {
    chr <- gsub('_.*', '', gsub('custom_acf_', '', f))
    message(paste(Sys.time(), 'loading', f))
    load(file.path('rda', f))
    auto_long$chr <- chr
    return(auto_long)
}
auto_long <- do.call(rbind, lapply(files, load_custom))
auto_long$lag <- as.factor(auto_long$lag)

## Add cell type info
library('bsseq')
load('../../BSobj_subsets/rda/DMP_nonCpG_interaction.Rdata')
auto_long$cell <- colData(DMP_nonCpG)$Cell.Type[match(auto_long$sample, colData(DMP_nonCpG)$Data.ID)]

dim(auto_long)
# [1] 2850944       8

pdf('pdf/autocorrelation_by_context_by_sample_chr21.pdf', width = 7 * 2, height = 7 * 3)
ggplot(auto_long, aes(x = lag, y = acf)) + geom_boxplot() + facet_grid(cell ~ context) + theme_bw(base_size = 30) + ggtitle('By sample (chr21 only)') + ylab('Auto correlation') + ylim(c(-1, 1))
dev.off()

pdf('pdf/autocorrelation_by_context_abs_by_sample_chr21.pdf', width = 7 * 2, height = 7 * 3)
ggplot(auto_long, aes(x = lag, y = abs(acf))) + geom_boxplot() + facet_grid(cell ~ context) + theme_bw(base_size = 30) + ggtitle('By sample (chr21 only)') + ylab('Absolute auto correlation') + ylim(c(0, 1))
dev.off()

>>>>>>> c342a3e8636f5c1dab8031343885425bc57c00f7
## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
<<<<<<< HEAD
=======

# Session info ----------------------------------------------------------------------------------------------------------
#  setting  value
#  version  R version 3.4.3 Patched (2018-01-20 r74142)
#  system   x86_64, linux-gnu
#  ui       X11
#  language (EN)
#  collate  en_US.UTF-8
#  tz       America/New_York
#  date     2018-01-26
#
# Packages --------------------------------------------------------------------------------------------------------------
#  package              * version   date       source
#  base                 * 3.4.3     2018-01-20 local
#  Biobase              * 2.38.0    2017-11-07 Bioconductor
#  BiocGenerics         * 0.24.0    2017-11-29 Bioconductor
#  bitops                 1.0-6     2013-08-17 CRAN (R 3.4.1)
#  bsseq                * 1.14.0    2017-11-07 Bioconductor
#  colorout             * 1.1-2     2017-08-10 Github (jalvesaq/colorout@020a14d)
#  colorspace             1.3-2     2016-12-14 CRAN (R 3.4.1)
#  compiler               3.4.3     2018-01-20 local
#  data.table             1.10.4-3  2017-10-27 CRAN (R 3.4.2)
#  datasets             * 3.4.3     2018-01-20 local
#  DelayedArray         * 0.4.1     2017-11-07 Bioconductor
#  devtools             * 1.13.4    2017-11-09 CRAN (R 3.4.2)
#  digest                 0.6.14    2018-01-14 CRAN (R 3.4.2)
#  GenomeInfoDb         * 1.14.0    2017-11-29 Bioconductor
#  GenomeInfoDbData       1.0.0     2018-01-09 Bioconductor
#  GenomicRanges        * 1.30.1    2018-01-09 Bioconductor
#  ggplot2              * 2.2.1     2016-12-30 CRAN (R 3.4.1)
#  graphics             * 3.4.3     2018-01-20 local
#  grDevices            * 3.4.3     2018-01-20 local
#  grid                   3.4.3     2018-01-20 local
#  gtable                 0.2.0     2016-02-26 CRAN (R 3.4.1)
#  gtools                 3.5.0     2015-05-29 CRAN (R 3.4.1)
#  IRanges              * 2.12.0    2017-11-29 Bioconductor
#  labeling               0.3       2014-08-23 CRAN (R 3.4.1)
#  lattice                0.20-35   2017-03-25 CRAN (R 3.4.3)
#  lazyeval               0.2.1     2017-10-29 CRAN (R 3.4.2)
#  limma                  3.34.5    2018-01-16 Bioconductor
#  locfit                 1.5-9.1   2013-04-20 CRAN (R 3.4.1)
#  magrittr               1.5       2014-11-22 CRAN (R 3.4.1)
#  Matrix                 1.2-12    2017-11-30 CRAN (R 3.4.3)
#  matrixStats          * 0.52.2    2017-04-14 CRAN (R 3.4.1)
#  memoise                1.1.0     2017-04-21 CRAN (R 3.4.1)
#  methods              * 3.4.3     2018-01-20 local
#  munsell                0.4.3     2016-02-13 CRAN (R 3.4.1)
#  parallel             * 3.4.3     2018-01-20 local
#  permute                0.9-4     2016-09-09 CRAN (R 3.4.1)
#  pillar                 1.1.0     2018-01-14 CRAN (R 3.4.2)
#  plyr                   1.8.4     2016-06-08 CRAN (R 3.4.1)
#  R.methodsS3            1.7.1     2016-02-16 CRAN (R 3.4.1)
#  R.oo                   1.21.0    2016-11-01 CRAN (R 3.4.1)
#  R.utils                2.6.0     2017-11-05 CRAN (R 3.4.2)
#  Rcpp                   0.12.14   2017-11-23 CRAN (R 3.4.2)
#  RCurl                  1.95-4.10 2018-01-04 CRAN (R 3.4.2)
#  reshape2               1.4.3     2017-12-11 CRAN (R 3.4.2)
#  rlang                  0.1.6     2017-12-21 CRAN (R 3.4.2)
#  rstudioapi             0.7       2017-09-07 CRAN (R 3.4.1)
#  S4Vectors            * 0.16.0    2017-11-29 Bioconductor
#  scales                 0.5.0     2017-08-24 CRAN (R 3.4.1)
#  stats                * 3.4.3     2018-01-20 local
#  stats4               * 3.4.3     2018-01-20 local
#  stringi                1.1.6     2017-11-17 CRAN (R 3.4.2)
#  stringr                1.2.0     2017-02-18 CRAN (R 3.4.1)
#  SummarizedExperiment * 1.8.1     2018-01-09 Bioconductor
#  tibble                 1.4.1     2017-12-25 CRAN (R 3.4.2)
#  tools                  3.4.3     2018-01-20 local
#  utils                * 3.4.3     2018-01-20 local
#  withr                  2.1.1     2017-12-19 CRAN (R 3.4.2)
#  XVector                0.18.0    2017-11-29 Bioconductor
#  zlibbioc               1.24.0    2017-11-07 Bioconductor
>>>>>>> c342a3e8636f5c1dab8031343885425bc57c00f7
