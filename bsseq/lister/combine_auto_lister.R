library('devtools')
library('ggplot2')

files <- dir('rda', pattern = '^auto_long_chr')

load_auto <- function(f) {
    chr <- gsub('_.*', '', gsub('auto_long_', '', f))
    message(paste(Sys.time(), 'loading', f))
    load(file.path('rda', f))
    auto_long$chr <- chr
    return(auto_long)
}

dir.create('rda', showWarnings = FALSE)
if(!file.exists('rda/auto_long_combined.Rdata')) {
    auto_long <- do.call(rbind, lapply(files, load_auto))
    auto_long$lag <- as.factor(auto_long$lag)
    save(auto_long, file = 'rda/auto_long_combined.Rdata')
} else {
    load('rda/auto_long_combined.Rdata', verbose = TRUE)
}
dim(auto_long)
# [1] 7043476       8

dir.create('pdf', showWarnings = FALSE)
pdf('pdf/autocorrelation_by_context.pdf', width = 7 * 2, height = 7 * 2)
ggplot(auto_long, aes(x = lag, y = acf_neuron)) + geom_boxplot() + facet_grid(. ~ context) + theme_bw(base_size = 30) + ggtitle('Neuron') + ylab('Auto correlation') + ylim(c(-1, 1))

ggplot(auto_long, aes(x = lag, y = acf_glia)) + geom_boxplot() + facet_grid(. ~ context) + theme_bw(base_size = 30) + ggtitle('Glia') + ylab('Auto correlation')  + ylim(c(-1, 1))
dev.off()

pdf('pdf/autocorrelation_by_context_abs.pdf', width = 7 * 2, height = 7 * 2)
ggplot(auto_long, aes(x = lag, y = abs(acf_neuron))) + geom_boxplot() + facet_grid(. ~ context) + theme_bw(base_size = 30) + ggtitle('Neuron') + ylab('Absolute auto correlation') + ylim(c(0, 1))

ggplot(auto_long, aes(x = lag, y = abs(acf_glia))) + geom_boxplot() + facet_grid(. ~ context) + theme_bw(base_size = 30) + ggtitle('Glia') + ylab('Absolute auto correlation') + ylim(c(0, 1))
dev.off()

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()

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
#  package    * version date       source
#  base       * 3.4.3   2018-01-20 local
#  colorout   * 1.1-2   2017-08-10 Github (jalvesaq/colorout@020a14d)
#  colorspace   1.3-2   2016-12-14 CRAN (R 3.4.1)
#  compiler     3.4.3   2018-01-20 local
#  datasets   * 3.4.3   2018-01-20 local
#  devtools   * 1.13.4  2017-11-09 CRAN (R 3.4.2)
#  digest       0.6.14  2018-01-14 CRAN (R 3.4.2)
#  ggplot2    * 2.2.1   2016-12-30 CRAN (R 3.4.1)
#  graphics   * 3.4.3   2018-01-20 local
#  grDevices  * 3.4.3   2018-01-20 local
#  grid         3.4.3   2018-01-20 local
#  gtable       0.2.0   2016-02-26 CRAN (R 3.4.1)
#  labeling     0.3     2014-08-23 CRAN (R 3.4.1)
#  lazyeval     0.2.1   2017-10-29 CRAN (R 3.4.2)
#  magrittr     1.5     2014-11-22 CRAN (R 3.4.1)
#  memoise      1.1.0   2017-04-21 CRAN (R 3.4.1)
#  methods    * 3.4.3   2018-01-20 local
#  munsell      0.4.3   2016-02-13 CRAN (R 3.4.1)
#  pillar       1.1.0   2018-01-14 CRAN (R 3.4.2)
#  plyr         1.8.4   2016-06-08 CRAN (R 3.4.1)
#  Rcpp         0.12.14 2017-11-23 CRAN (R 3.4.2)
#  reshape2     1.4.3   2017-12-11 CRAN (R 3.4.2)
#  rlang        0.1.6   2017-12-21 CRAN (R 3.4.2)
#  scales       0.5.0   2017-08-24 CRAN (R 3.4.1)
#  stats      * 3.4.3   2018-01-20 local
#  stringi      1.1.6   2017-11-17 CRAN (R 3.4.2)
#  stringr      1.2.0   2017-02-18 CRAN (R 3.4.1)
#  tibble       1.4.1   2017-12-25 CRAN (R 3.4.2)
#  tools        3.4.3   2018-01-20 local
#  utils      * 3.4.3   2018-01-20 local
#  withr        2.1.1   2017-12-19 CRAN (R 3.4.2)
