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
# [1] acf_neuron   acf_glia     lag          clus         clus_n
# [6] clus_range   context      model        avg_distance
# <0 rows> (or 0-length row.names)

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


pdf('pdf/bp_distance_by_lag.pdf', width = 10, height = 10, useDingbats = FALSE)
ggplot(auto_long, aes(x = lag, y = avg_distance)) + geom_boxplot() + facet_grid(model ~ context) + theme_bw(base_size = 30) + ylab('Avg base-pair distance')
ggplot(auto_long, aes(x = lag, y = avg_distance)) + geom_boxplot() + facet_grid(model ~ context) + theme_bw(base_size = 30) + ylab('Avg base-pair distance') + scale_y_log10()

summary(auto_long$avg_distance)
#  Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
# 1.333   31.637   64.437  114.359  136.147 4538.000

ggplot(sub, aes(x = lag, y = avg_distance)) + geom_boxplot() + facet_grid(~ context) + theme_bw(base_size = 30) + ylab('Avg base-pair distance')
ggplot(sub, aes(x = lag, y = avg_distance)) + geom_boxplot() + facet_grid(~ context) + theme_bw(base_size = 30) + ylab('Avg base-pair distance')  + scale_y_log10()
dev.off()

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
#  date     2019-02-11
#
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
#  package     * version date       lib source
#  assertthat    0.2.0   2017-04-11 [2] CRAN (R 3.5.0)
#  backports     1.1.3   2018-12-14 [2] CRAN (R 3.5.1)
#  bindr         0.1.1   2018-03-13 [1] CRAN (R 3.5.0)
#  bindrcpp      0.2.2   2018-03-29 [1] CRAN (R 3.5.0)
#  callr         3.1.1   2018-12-21 [2] CRAN (R 3.5.1)
#  cli           1.0.1   2018-09-25 [1] CRAN (R 3.5.1)
#  colorout    * 1.2-0   2018-05-02 [1] Github (jalvesaq/colorout@c42088d)
#  colorspace    1.4-0   2019-01-13 [2] CRAN (R 3.5.1)
#  crayon        1.3.4   2017-09-16 [1] CRAN (R 3.5.0)
#  desc          1.2.0   2018-05-01 [2] CRAN (R 3.5.1)
#  devtools    * 2.0.1   2018-10-26 [1] CRAN (R 3.5.1)
#  digest        0.6.18  2018-10-10 [1] CRAN (R 3.5.1)
#  dplyr         0.7.8   2018-11-10 [1] CRAN (R 3.5.1)
#  fs            1.2.6   2018-08-23 [2] CRAN (R 3.5.1)
#  ggplot2     * 3.1.0   2018-10-25 [1] CRAN (R 3.5.1)
#  glue          1.3.0   2018-07-17 [1] CRAN (R 3.5.1)
#  gtable        0.2.0   2016-02-26 [2] CRAN (R 3.5.0)
#  htmltools     0.3.6   2017-04-28 [2] CRAN (R 3.5.0)
#  htmlwidgets   1.3     2018-09-30 [1] CRAN (R 3.5.1)
#  httpuv        1.4.5.1 2018-12-18 [2] CRAN (R 3.5.1)
#  labeling      0.3     2014-08-23 [2] CRAN (R 3.5.0)
#  later         0.7.5   2018-09-18 [2] CRAN (R 3.5.1)
#  lattice       0.20-38 2018-11-04 [3] CRAN (R 3.5.1)
#  lazyeval      0.2.1   2017-10-29 [2] CRAN (R 3.5.0)
#  magrittr      1.5     2014-11-22 [1] CRAN (R 3.5.0)
#  memoise       1.1.0   2017-04-21 [2] CRAN (R 3.5.0)
#  mime          0.6     2018-10-05 [1] CRAN (R 3.5.1)
#  munsell       0.5.0   2018-06-12 [2] CRAN (R 3.5.0)
#  pillar        1.3.1   2018-12-15 [1] CRAN (R 3.5.1)
#  pkgbuild      1.0.2   2018-10-16 [2] CRAN (R 3.5.1)
#  pkgconfig     2.0.2   2018-08-16 [1] CRAN (R 3.5.1)
#  pkgload       1.0.2   2018-10-29 [2] CRAN (R 3.5.1)
#  plyr          1.8.4   2016-06-08 [2] CRAN (R 3.5.0)
#  png           0.1-7   2013-12-03 [2] CRAN (R 3.5.0)
#  prettyunits   1.0.2   2015-07-13 [1] CRAN (R 3.5.0)
#  processx      3.2.1   2018-12-05 [1] CRAN (R 3.5.1)
#  promises      1.0.1   2018-04-13 [2] CRAN (R 3.5.0)
#  ps            1.3.0   2018-12-21 [2] CRAN (R 3.5.1)
#  purrr         0.2.5   2018-05-29 [2] CRAN (R 3.5.0)
#  R6            2.3.0   2018-10-04 [2] CRAN (R 3.5.1)
#  Rcpp          1.0.0   2018-11-07 [1] CRAN (R 3.5.1)
#  remotes       2.0.2   2018-10-30 [1] CRAN (R 3.5.1)
#  reshape2      1.4.3   2017-12-11 [2] CRAN (R 3.5.0)
#  rlang         0.3.1   2019-01-08 [1] CRAN (R 3.5.1)
#  rmote       * 0.3.4   2018-05-02 [1] deltarho (R 3.5.0)
#  rprojroot     1.3-2   2018-01-03 [2] CRAN (R 3.5.0)
#  scales        1.0.0   2018-08-09 [2] CRAN (R 3.5.1)
#  servr         0.11    2018-10-23 [1] CRAN (R 3.5.1)
#  sessioninfo   1.1.1   2018-11-05 [1] CRAN (R 3.5.1)
#  stringi       1.2.4   2018-07-20 [2] CRAN (R 3.5.1)
#  stringr       1.3.1   2018-05-10 [1] CRAN (R 3.5.0)
#  testthat      2.0.1   2018-10-13 [1] CRAN (R 3.5.1)
#  tibble        2.0.1   2019-01-12 [1] CRAN (R 3.5.1)
#  tidyselect    0.2.5   2018-10-11 [2] CRAN (R 3.5.1)
#  usethis     * 1.4.0   2018-08-14 [2] CRAN (R 3.5.1)
#  withr         2.1.2   2018-03-15 [2] CRAN (R 3.5.0)
#  xfun          0.4     2018-10-23 [1] CRAN (R 3.5.1)
#
# [1] /users/lcollado/R/x86_64-pc-linux-gnu-library/3.5.x
# [2] /jhpce/shared/jhpce/core/conda/miniconda-3/envs/svnR-3.5.x/R/3.5.x/lib64/R/site-library
# [3] /jhpce/shared/jhpce/core/conda/miniconda-3/envs/svnR-3.5.x/R/3.5.x/lib64/R/library
