## Based on https://github.com/LieberInstitute/brain-epigenomics/blob/master/bsseq/bsobj_by_chr/combine_auto.R

library('devtools')
library('ggplot2')

files <- dir('rda', pattern = '^auto_long', full.names = TRUE)

load_auto <- function(f) {
    message(paste(Sys.time(), 'loading', f))
    load(f)
    return(auto_long)
}

if(!file.exists('rda/auto_long_combined.Rdata')) {
    auto_long <- do.call(rbind, lapply(files, load_auto))
    # auto_long$context[auto_long$context == 'nonCG'] <- 'CpH'
    # auto_long$context[auto_long$context == 'CG'] <- 'CpG'
    auto_long$lag <- as.factor(auto_long$lag)
    save(auto_long, file = 'rda/auto_long_combined.Rdata')
} else {
    load('rda/auto_long_combined.Rdata')
}
dim(auto_long)
# [1] 234584      9

head(auto_long)
#   acf_neuron  acf_glia lag clus clus_n clus_range context model avg_distance
# 1  0.4220721 0.4926867   1    1   1513      19172     all   age     12.67923
# 2  0.3547252 0.4194819   2    1   1513      19172     all   age     25.36664
# 3  0.3060108 0.3613006   3    1   1513      19172     all   age     38.06821
# 4  0.2283096 0.2979926   4    1   1513      19172     all   age     50.78396
# 5  0.5084426 0.4796444   1    2    638       9007     all   age     14.13815
# 6  0.3721603 0.3722272   2    2    638       9007     all   age     28.11321

auto_long[head(which(auto_long$avg_distance < 0)), ]
#          acf_neuron    acf_glia lag clus clus_n clus_range context model
# 101732 -0.003217755 -0.06663796   4   10     15       1216     CHG   age
# 101761 -0.012045398 -0.05505406   1   18     36        945     CHG   age
# 101762  0.119849050  0.07904448   2   18     36        945     CHG   age
# 101763 -0.047964536 -0.04692480   3   18     36        945     CHG   age
# 101764 -0.021825900 -0.02983349   4   18     36        945     CHG   age
# 101837 -0.068108113 -0.02281280   1   37     13       1145     CHG   age
#        avg_distance
# 101732    -4.090909
# 101761    -2.114286
# 101762    -5.911765
# 101763   -19.090909
# 101764   -34.437500
# 101837   -12.000000

## This never finished running:
# mod_context_summary <- lapply(unique(auto_long)$model, function(model) {
#     lapply(unique(auto_long)$context, function(context) {
#          sub <- auto_long[auto_long$context == context & auto_long$model == model, ]
#          tapply(sub$acf, sub$lag, summary)
#      })
# })
# mod_context_summary
# save(mod_context_summary, file = 'autocorrelation_summary_by_context_and_model.Rdata')

dir.create('pdf', showWarnings = FALSE)

pdf('pdf/autocorrelation_by_context_and_model.pdf', width = 7 * 2, height = 7 * 3)
ggplot(auto_long, aes(x = lag, y = acf_neuron)) + geom_boxplot() + facet_grid(model ~ context) + theme_bw(base_size = 30) + ggtitle('Neuron') + ylab('Auto correlation') + ylim(c(-1, 1))

ggplot(auto_long, aes(x = lag, y = acf_glia)) + geom_boxplot() + facet_grid(model ~ context) + theme_bw(base_size = 30) + ggtitle('Glia') + ylab('Auto correlation')  + ylim(c(-1, 1))
dev.off()

pdf('pdf/autocorrelation_by_context_and_model_abs.pdf', width = 7 * 2, height = 7 * 3)
ggplot(auto_long, aes(x = lag, y = abs(acf_neuron))) + geom_boxplot() + facet_grid(model ~ context) + theme_bw(base_size = 30) + ggtitle('Neuron') + ylab('Absolute auto correlation') + ylim(c(0, 1))

ggplot(auto_long, aes(x = lag, y = abs(acf_glia))) + geom_boxplot() + facet_grid(model ~ context) + theme_bw(base_size = 30) + ggtitle('Glia') + ylab('Absolute auto correlation') + ylim(c(0, 1))
dev.off()

pdf('pdf/autocorrelation_by_context_and_model_int_subset.pdf', width = 7, height = 7, useDingbats = FALSE)
sub <- subset(auto_long, context %in% c('all', 'CG', 'nonCG') & model == 'interaction')
sub$context <- c('all' = 'C', 'CG' = 'CpG', 'nonCG' = 'CpH')[sub$context]

ggplot(sub, aes(x = lag, y = acf_neuron)) + geom_boxplot() + facet_grid(~ context) + theme_bw(base_size = 30) + ggtitle('Neuron') + ylab('Auto correlation') + ylim(c(-1, 1))
ggplot(sub, aes(x = lag, y = acf_glia)) + geom_boxplot() + facet_grid(~ context) + theme_bw(base_size = 30) + ggtitle('Glia') + ylab('Auto correlation') + ylim(c(-1, 1))
dev.off()



ggplot(auto_long, aes(x = lag, y = avg_distance)) + geom_boxplot() + facet_grid(model ~ context) + theme_bw(base_size = 30) + ylab('Avg base-pair distance')


ggplot(auto_long, aes(x = lag, y = log10(avg_distance + 1))) + geom_boxplot() + facet_grid(model ~ context) + theme_bw(base_size = 30) + ylab('Avg base-pair distance')
ggplot(auto_long, aes(x = lag, y = avg_distance)) + geom_boxplot() + facet_grid(model ~ context) + theme_bw(base_size = 30) + ylab('Avg base-pair distance') + scale_y_log10()

summary(auto_long$avg_distance)
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
# -2236.00    22.42    50.77    94.67   116.08  4538.00

ggplot(sub, aes(x = lag, y = avg_distance)) + geom_boxplot() + facet_grid(~ context) + theme_bw(base_size = 30) + ylab('Avg base-pair distance')

ggplot(sub, aes(x = lag, y = avg_distance)) + geom_boxplot() + facet_grid(~ context) + theme_bw(base_size = 30) + ylab('Avg base-pair distance')  + scale_y_log10()

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
#
# Session info ----------------------------------------------------------------------------------------------------------
#  setting  value
#  version  R version 3.4.3 Patched (2018-01-20 r74142)
#  system   x86_64, linux-gnu
#  ui       X11
#  language (EN)
#  collate  en_US.UTF-8
#  tz       US/Eastern
#  date     2018-03-13
#
# Packages --------------------------------------------------------------------------------------------------------------
#  package     * version date       source
#  base        * 3.4.3   2018-01-20 local
#  colorout    * 1.2-0   2018-02-19 Github (jalvesaq/colorout@2f01173)
#  colorspace    1.3-2   2016-12-14 CRAN (R 3.4.1)
#  compiler      3.4.3   2018-01-20 local
#  datasets    * 3.4.3   2018-01-20 local
#  devtools    * 1.13.4  2017-11-09 CRAN (R 3.4.2)
#  digest        0.6.15  2018-01-28 cran (@0.6.15)
#  ggplot2     * 2.2.1   2016-12-30 CRAN (R 3.4.1)
#  graphics    * 3.4.3   2018-01-20 local
#  grDevices   * 3.4.3   2018-01-20 local
#  grid          3.4.3   2018-01-20 local
#  gtable        0.2.0   2016-02-26 CRAN (R 3.4.1)
#  htmltools     0.3.6   2017-04-28 CRAN (R 3.4.1)
#  htmlwidgets   0.9     2017-07-10 CRAN (R 3.4.1)
#  httpuv        1.3.6.2 2018-03-02 CRAN (R 3.4.3)
#  labeling      0.3     2014-08-23 CRAN (R 3.4.1)
#  lattice       0.20-35 2017-03-25 CRAN (R 3.4.3)
#  lazyeval      0.2.1   2017-10-29 CRAN (R 3.4.2)
#  magrittr      1.5     2014-11-22 CRAN (R 3.4.1)
#  memoise       1.1.0   2017-04-21 CRAN (R 3.4.1)
#  methods     * 3.4.3   2018-01-20 local
#  munsell       0.4.3   2016-02-13 CRAN (R 3.4.1)
#  pillar        1.1.0   2018-01-14 CRAN (R 3.4.2)
#  plyr          1.8.4   2016-06-08 CRAN (R 3.4.1)
#  png           0.1-7   2013-12-03 CRAN (R 3.4.1)
#  Rcpp          0.12.14 2017-11-23 CRAN (R 3.4.2)
#  reshape2      1.4.3   2017-12-11 CRAN (R 3.4.2)
#  rlang         0.1.6   2017-12-21 CRAN (R 3.4.2)
#  rmote       * 0.3.4   2018-02-16 deltarho (R 3.4.3)
#  scales        0.5.0   2017-08-24 CRAN (R 3.4.1)
#  servr         0.8     2017-11-06 CRAN (R 3.4.3)
#  stats       * 3.4.3   2018-01-20 local
#  stringi       1.1.6   2017-11-17 CRAN (R 3.4.2)
#  stringr       1.2.0   2017-02-18 CRAN (R 3.4.1)
#  tibble        1.4.1   2017-12-25 CRAN (R 3.4.2)
#  tools         3.4.3   2018-01-20 local
#  utils       * 3.4.3   2018-01-20 local
#  withr         2.1.1   2017-12-19 CRAN (R 3.4.2)