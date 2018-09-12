library('devtools')
library('bsseq')

load("/dcl01/lieber/WGBS/LIBD_Data/bsseqObj/bsseqObj_postNatal_cleaned_CpGonly.rda", verbose = TRUE)
pd <- pData(BSobj)

pd_neuron_glia <- subset(pd, Cell.Type %in% c('Neuron', 'Glia'))
pd_homogenate <- subset(pd, Cell.Type %in% 'Homogenate')

lambda <- read.csv('/dcl01/lieber/WGBS/LIBD_Data/lambda_alignment/MasterList_lambda.csv')
stopifnot(all(pd$Data.ID %in% lambda$WGC.ID))

## Use the same names for the cell types as in the DMR code
cell_info <- ifelse(lambda$Cell.Type == 'H', 'Homogenate', ifelse(lambda$Cell.Type == 'NeuN-', 'Glia', 'Neuron'))
table(lambda$Cell.Type)
#  H NeuN- NeuN+
# 42     8    28
table(cell_info)
# Glia Homogenate     Neuron
#    8         42         28

## Just match to the info on the pd we used (that csv file looks outdated?)
m <- match(lambda$WGC.ID, pd$Data.ID)
lambda$Cell.Type <- pd$Cell.Type[m]
table(lambda$Cell.Type)
# Glia Homogenate     Neuron
#          8         23         24
lambda$Age[is.na(lambda$Cell.Type)]
#  [1] -0.421917 -0.383561 -0.498630 -0.421917 -0.441095 -0.402739 -0.402739
#  [8] -0.402739 -0.402739 -0.402739 15.160000  2.490000 -0.479452 -0.421917
# [15] -0.460273 -0.421917 -0.498630 16.790000 -0.421917 -0.402739 -0.498630
# [22] -0.498630 -0.498630
cell_info[is.na(lambda$Cell.Type)]
#  [1] "Homogenate" "Homogenate" "Homogenate" "Homogenate" "Homogenate"
#  [6] "Homogenate" "Homogenate" "Homogenate" "Homogenate" "Homogenate"
# [11] "Neuron"     "Neuron"     "Homogenate" "Homogenate" "Homogenate"
# [16] "Homogenate" "Homogenate" "Neuron"     "Homogenate" "Homogenate"
# [21] "Homogenate" "Homogenate" "Homogenate"

pd$Age[m[!is.na(m)]]
#  [1] 18.15  2.49  0.33  5.34 19.89  8.25 16.79  0.20  8.25  0.20 14.01 15.48
# [13] 13.02 17.22 22.71 13.02 15.15 21.01 21.01 20.77 20.77 15.48 13.69 13.69
# [25] 14.62 14.62 18.11 18.11  0.20  0.25  0.25 18.05 14.01  2.71  2.71  0.36
# [37]  0.36 19.89 13.02 17.22 17.22 15.15 22.58 22.58 18.05 22.71 22.71 19.89
# [49] 18.15  0.33  5.34  2.49  8.25  2.71  5.34
lambda$Age[!is.na(m)]
#  [1] 18.15  2.49  0.33  5.34 19.89  8.25 16.79  0.20  8.25  0.20 14.01 15.48
# [13] 13.02 17.22 22.71 13.02 15.15 21.01 21.01 20.77 20.77 15.48 13.61 13.61
# [25] 14.62 14.62 18.11 18.11  0.20  0.25  0.25 14.01  2.71  2.71  0.36  0.36
# [37]  0.36 19.89 13.02 17.22 17.22 15.15 22.58 22.58 18.05 22.71 22.71 19.89
# [49] 18.15  0.33  5.34  2.49  8.25  2.71  5.34
lambda$Age[!is.na(m)] - pd$Age[m[!is.na(m)]]
#  [1]   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00
# [11]   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00
# [21]   0.00   0.00  -0.08  -0.08   0.00   0.00   0.00   0.00   0.00   0.00
# [31]   0.00  -4.04 -11.30   0.00  -2.35   0.00   0.00   0.00   0.00   0.00
# [41]   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00
# [51]   0.00   0.00   0.00   0.00   0.00

## Age is not matching, replace by info on the pd we used
age_info <- lambda$Age
lambda$Age <- pd$Age[m]

lambda_neuron_glia <- subset(lambda, WGC.ID %in% pd_neuron_glia$Data.ID & Cell.Type != 'Homogenate')
lambda_homogenate <- subset(lambda, WGC.ID %in% pd_homogenate$Data.ID & Cell.Type == 'Homogenate')

stopifnot(nrow(lambda_neuron_glia) == nrow(pd_neuron_glia))
stopifnot(nrow(lambda_homogenate) == nrow(pd_homogenate))

## DMR models
fits <- list(
    'cell' = lm(avg_conv_efficiency ~ Cell.Type + Age, data = lambda_neuron_glia),
    'cell_only' = lm(avg_conv_efficiency ~ Cell.Type, data = lambda_neuron_glia),
    'age' = lm(avg_conv_efficiency ~ Age + Cell.Type, data = lambda_neuron_glia),
    'age_only' = lm(avg_conv_efficiency ~ Age, data = lambda_neuron_glia),
    'interaction' = lm(avg_conv_efficiency ~ Age * Cell.Type, data = lambda_neuron_glia),
    'homogenateage' = lm(avg_conv_efficiency ~ Age, data = lambda_homogenate)
)

fits
# $cell
#
# Call:
# lm(formula = avg_conv_efficiency ~ Cell.Type + Age, data = lambda_neuron_glia)
#
# Coefficients:
#     (Intercept)  Cell.TypeNeuron              Age
#       98.757026        -0.068254        -0.004883
#
#
# $cell_only
#
# Call:
# lm(formula = avg_conv_efficiency ~ Cell.Type, data = lambda_neuron_glia)
#
# Coefficients:
#     (Intercept)  Cell.TypeNeuron
#          98.703           -0.075
#
#
# $age
#
# Call:
# lm(formula = avg_conv_efficiency ~ Age + Cell.Type, data = lambda_neuron_glia)
#
# Coefficients:
#     (Intercept)              Age  Cell.TypeNeuron
#       98.757026        -0.004883        -0.068254
#
#
# $age_only
#
# Call:
# lm(formula = avg_conv_efficiency ~ Age, data = lambda_neuron_glia)
#
# Coefficients:
# (Intercept)          Age
#   98.709436    -0.005178
#
#
# $interaction
#
# Call:
# lm(formula = avg_conv_efficiency ~ Age * Cell.Type,
#     data = lambda_neuron_glia)
#
# Coefficients:
#         (Intercept)                  Age      Cell.TypeNeuron  Age:Cell.TypeNeuron
#            98.84588             -0.01284             -0.19105              0.01066
#
#
# $homogenateage
#
# Call:
# lm(formula = avg_conv_efficiency ~ Age, data = lambda_homogenate)
#
# Coefficients:
# (Intercept)          Age
#   98.677095    -0.002015
#


lapply(fits, summary)
# $cell
#
# Call:
# lm(formula = avg_conv_efficiency ~ Cell.Type + Age, data = lambda_neuron_glia)
#
# Residuals:
#      Min       1Q   Median       3Q      Max
# -0.52780 -0.11850 -0.02207  0.12033  0.49321
#
# Coefficients:
#                  Estimate Std. Error t value Pr(>|t|)
# (Intercept)     98.757026   0.106875 924.041   <2e-16 ***
# Cell.TypeNeuron -0.068254   0.100287  -0.681    0.502
# Age             -0.004883   0.005609  -0.870    0.391
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# Residual standard error: 0.2449 on 29 degrees of freedom
# Multiple R-squared:  0.04354,    Adjusted R-squared:  -0.02242
# F-statistic: 0.6601 on 2 and 29 DF,  p-value: 0.5244
#
#
# $cell_only
#
# Call:
# lm(formula = avg_conv_efficiency ~ Cell.Type, data = lambda_neuron_glia)
#
# Residuals:
#     Min      1Q  Median      3Q     Max
# -0.4675 -0.1275 -0.0375  0.1325  0.4725
#
# Coefficients:
#                 Estimate Std. Error  t value Pr(>|t|)
# (Intercept)     98.70250    0.08624 1144.497   <2e-16 ***
# Cell.TypeNeuron -0.07500    0.09958   -0.753    0.457
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# Residual standard error: 0.2439 on 30 degrees of freedom
# Multiple R-squared:  0.01856,    Adjusted R-squared:  -0.01416
# F-statistic: 0.5672 on 1 and 30 DF,  p-value: 0.4572
#
#
# $age
#
# Call:
# lm(formula = avg_conv_efficiency ~ Age + Cell.Type, data = lambda_neuron_glia)
#
# Residuals:
#      Min       1Q   Median       3Q      Max
# -0.52780 -0.11850 -0.02207  0.12033  0.49321
#
# Coefficients:
#                  Estimate Std. Error t value Pr(>|t|)
# (Intercept)     98.757026   0.106875 924.041   <2e-16 ***
# Age             -0.004883   0.005609  -0.870    0.391
# Cell.TypeNeuron -0.068254   0.100287  -0.681    0.502
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# Residual standard error: 0.2449 on 29 degrees of freedom
# Multiple R-squared:  0.04354,    Adjusted R-squared:  -0.02242
# F-statistic: 0.6601 on 2 and 29 DF,  p-value: 0.5244
#
#
# $age_only
#
# Call:
# lm(formula = avg_conv_efficiency ~ Age, data = lambda_neuron_glia)
#
# Residuals:
#      Min       1Q   Median       3Q      Max
# -0.54840 -0.12071 -0.01397  0.10264  0.47750
#
# Coefficients:
#              Estimate Std. Error  t value Pr(>|t|)
# (Intercept) 98.709436   0.080099 1232.338   <2e-16 ***
# Age         -0.005178   0.005542   -0.934    0.358
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# Residual standard error: 0.2427 on 30 degrees of freedom
# Multiple R-squared:  0.02827,    Adjusted R-squared:  -0.004124
# F-statistic: 0.8727 on 1 and 30 DF,  p-value: 0.3577
#
#
# $interaction
#
# Call:
# lm(formula = avg_conv_efficiency ~ Age * Cell.Type,
#     data = lambda_neuron_glia)
#
# Residuals:
#      Min       1Q   Median       3Q      Max
# -0.49439 -0.10793 -0.02607  0.12899  0.48173
#
# Coefficients:
#                     Estimate Std. Error t value Pr(>|t|)
# (Intercept)         98.84588    0.15240 648.614   <2e-16 ***
# Age                 -0.01284    0.01120  -1.146    0.261
# Cell.TypeNeuron     -0.19105    0.18019  -1.060    0.298
# Age:Cell.TypeNeuron  0.01066    0.01296   0.822    0.418
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# Residual standard error: 0.2463 on 28 degrees of freedom
# Multiple R-squared:  0.0661,    Adjusted R-squared:  -0.03396
# F-statistic: 0.6606 on 3 and 28 DF,  p-value: 0.5832
#
#
# $homogenateage
#
# Call:
# lm(formula = avg_conv_efficiency ~ Age, data = lambda_homogenate)
#
# Residuals:
#      Min       1Q   Median       3Q      Max
# -0.38590 -0.12211 -0.02051  0.09607  0.38525
#
# Coefficients:
#              Estimate Std. Error  t value Pr(>|t|)
# (Intercept) 98.677095   0.073932 1334.693   <2e-16 ***
# Age         -0.002015   0.005055   -0.399    0.694
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# Residual standard error: 0.1894 on 21 degrees of freedom
# Multiple R-squared:  0.007514,    Adjusted R-squared:  -0.03975
# F-statistic: 0.159 on 1 and 21 DF,  p-value: 0.6941
#

coefs <- rep(2, length(fits))
names(coefs) <- names(fits)
coefs['interaction'] <- 4

coef_table <- mapply(function(fit, coef) { fit$coefficients[coef, ] }, lapply(fits, summary), coefs)
coef_table
#                   cell   cell_only          age     age_only interaction homogenateage
# Estimate   -0.06825393 -0.07500000 -0.004882557 -0.005177597  0.01066138  -0.002015474
# Std. Error  0.10028712  0.09958246  0.005609448  0.005542410  0.01296382   0.005054774
# t value    -0.68058525 -0.75314467 -0.870416698 -0.934177998  0.82239536  -0.398726784
# Pr(>|t|)    0.50153210  0.45723305  0.391220966  0.357670832  0.41780069   0.694119870
save(age_info, cell_info, pd, lambda, lambda_neuron_glia, lambda_homogenate, pd_neuron_glia, pd_homogenate, m, fits, coefs, coef_table, file = 'lambda_checks.Rdata')

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
#  date     2018-09-11
#
# Packages --------------------------------------------------------------------------------------------------------------
#  package              * version   date       source
#  assertthat             0.2.0     2017-04-11 CRAN (R 3.5.0)
#  base                 * 3.5.0     2018-05-02 local
#  bindr                  0.1.1     2018-03-13 CRAN (R 3.5.0)
#  bindrcpp               0.2.2     2018-03-29 CRAN (R 3.5.0)
#  Biobase              * 2.40.0    2018-05-02 Bioconductor
#  BiocGenerics         * 0.26.0    2018-05-03 Bioconductor
#  BiocParallel         * 1.14.2    2018-07-08 Bioconductor
#  bitops                 1.0-6     2013-08-17 CRAN (R 3.5.0)
#  bsseq                * 1.16.1    2018-06-14 Bioconductor
#  colorout             * 1.2-0     2018-05-02 Github (jalvesaq/colorout@c42088d)
#  colorspace             1.3-2     2016-12-14 CRAN (R 3.5.0)
#  compiler               3.5.0     2018-05-02 local
#  crayon                 1.3.4     2017-09-16 CRAN (R 3.5.0)
#  data.table             1.11.4    2018-05-27 cran (@1.11.4)
#  datasets             * 3.5.0     2018-05-02 local
#  DelayedArray         * 0.6.5     2018-08-15 Bioconductor
#  DelayedMatrixStats     1.2.0     2018-05-02 Bioconductor
#  devtools             * 1.13.6    2018-06-27 CRAN (R 3.5.0)
#  digest                 0.6.15    2018-01-28 CRAN (R 3.5.0)
#  dplyr                  0.7.6     2018-06-29 CRAN (R 3.5.0)
#  GenomeInfoDb         * 1.16.0    2018-05-03 Bioconductor
#  GenomeInfoDbData       1.1.0     2018-04-17 Bioconductor
#  GenomicRanges        * 1.32.6    2018-07-20 Bioconductor
#  ggplot2                3.0.0     2018-07-03 CRAN (R 3.5.0)
#  glue                   1.3.0     2018-07-17 CRAN (R 3.5.0)
#  graphics             * 3.5.0     2018-05-02 local
#  grDevices            * 3.5.0     2018-05-02 local
#  grid                   3.5.0     2018-05-02 local
#  gtable                 0.2.0     2016-02-26 CRAN (R 3.5.0)
#  gtools                 3.8.1     2018-06-26 CRAN (R 3.5.0)
#  HDF5Array              1.8.1     2018-06-21 Bioconductor
#  htmltools              0.3.6     2017-04-28 CRAN (R 3.5.0)
#  htmlwidgets            1.2       2018-04-19 CRAN (R 3.5.0)
#  httpuv                 1.4.5     2018-07-19 CRAN (R 3.5.0)
#  IRanges              * 2.14.10   2018-05-17 Bioconductor
#  later                  0.7.4     2018-08-31 CRAN (R 3.5.0)
#  lattice                0.20-35   2017-03-25 CRAN (R 3.5.0)
#  lazyeval               0.2.1     2017-10-29 CRAN (R 3.5.0)
#  limma                  3.36.2    2018-06-21 Bioconductor
#  locfit                 1.5-9.1   2013-04-20 CRAN (R 3.5.0)
#  magrittr               1.5       2014-11-22 CRAN (R 3.5.0)
#  Matrix                 1.2-14    2018-04-13 CRAN (R 3.5.0)
#  matrixStats          * 0.54.0    2018-07-23 CRAN (R 3.5.0)
#  memoise                1.1.0     2017-04-21 CRAN (R 3.5.0)
#  methods              * 3.5.0     2018-05-02 local
#  munsell                0.5.0     2018-06-12 CRAN (R 3.5.0)
#  parallel             * 3.5.0     2018-05-02 local
#  permute                0.9-4     2016-09-09 CRAN (R 3.5.0)
#  pillar                 1.3.0     2018-07-14 CRAN (R 3.5.0)
#  pkgconfig              2.0.1     2017-03-21 CRAN (R 3.5.0)
#  plyr                   1.8.4     2016-06-08 CRAN (R 3.5.0)
#  png                    0.1-7     2013-12-03 CRAN (R 3.5.0)
#  promises               1.0.1     2018-04-13 CRAN (R 3.5.0)
#  purrr                  0.2.5     2018-05-29 CRAN (R 3.5.0)
#  R.methodsS3            1.7.1     2016-02-16 CRAN (R 3.5.0)
#  R.oo                   1.22.0    2018-04-22 CRAN (R 3.5.0)
#  R.utils                2.6.0     2017-11-05 CRAN (R 3.5.0)
#  R6                     2.2.2     2017-06-17 CRAN (R 3.5.0)
#  Rcpp                   0.12.18   2018-07-23 CRAN (R 3.5.0)
#  RCurl                  1.95-4.11 2018-07-15 CRAN (R 3.5.0)
#  rhdf5                  2.24.0    2018-05-02 Bioconductor
#  Rhdf5lib               1.2.1     2018-09-11 Bioconductor
#  rlang                  0.2.1     2018-05-30 cran (@0.2.1)
#  rmote                * 0.3.4     2018-05-02 deltarho (R 3.5.0)
#  rstudioapi             0.7       2017-09-07 CRAN (R 3.5.0)
#  S4Vectors            * 0.18.3    2018-06-13 Bioconductor
#  scales                 1.0.0     2018-08-09 CRAN (R 3.5.0)
#  servr                  0.10      2018-05-30 CRAN (R 3.5.0)
#  stats                * 3.5.0     2018-05-02 local
#  stats4               * 3.5.0     2018-05-02 local
#  SummarizedExperiment * 1.10.1    2018-05-17 Bioconductor
#  tibble                 1.4.2     2018-01-22 CRAN (R 3.5.0)
#  tidyselect             0.2.4     2018-02-26 CRAN (R 3.5.0)
#  tools                  3.5.0     2018-05-02 local
#  utils                * 3.5.0     2018-05-02 local
#  withr                  2.1.2     2018-03-15 CRAN (R 3.5.0)
#  xfun                   0.3       2018-07-06 CRAN (R 3.5.0)
#  XVector                0.20.0    2018-05-03 Bioconductor
#  zlibbioc               1.26.0    2018-05-02 Bioconductor
