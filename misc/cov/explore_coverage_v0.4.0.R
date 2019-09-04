# module load conda_R/3.5.x
library('jaffelab')
library('readr')
library('sessioninfo')
library('ggplot2')
library('SummarizedExperiment')
library('tidyr')
library('broom')
library('RColorBrewer')

## For output files
dir.create('rda', showWarnings = FALSE)

## bamcount reader
read_bamcount <- function(f) {
    readr::read_tsv(f, col_names = c('cat', 'n'))$n
}


#### Exploratory code

## Get ids and genome size
ids <- readLines('WGC_IDs_subset.txt')
hg19 <- readr::read_tsv(
    '/dcl01/lieber/ajaffe/lab/brain-epigenomics/misc/cov/test/human.hg19.genome',
    col_names = c('chr', 'length')
)
hg19.size <- sum(hg19$length[hg19$chr %in% paste0('chr', c(1:22, 'X', 'Y', 'M'))])

## Read in fastqc data
load('rda/cov_fastqc.Rdata', verbose = TRUE)

## Read in bamcount output data
get_bamcount_data <- function(id) {
    message(paste(Sys.time(), 'processing sample', id))
    
    ## Find the files
    types <- c('duplicatesRemoved', 'Marked_duplicates')
    bamcount_files <- file.path(
        '/dcl01/lieber/ajaffe/lab/brain-epigenomics/misc/cov/auc_files_v0.4.0',
        paste0(id, '_', types, '.auc.tsv')
    )
    stopifnot(all(file.exists(bamcount_files)))
    
    ## Read the data in
    auc_bam <- sapply(bamcount_files, read_bamcount)
    
    res <- data.frame(
        sample = id,
        type = types,
        auc = auc_bam,
        stringsAsFactors = FALSE
    )
    rownames(res) <- NULL
    return(res)
}

cov_bamcount <- do.call(rbind, lapply(ids, get_bamcount_data))
save(cov_bamcount, file = 'rda/cov_bamcount_v0.4.0.Rdata')

## Merge both types
cov_merged <- rbind(cov_fastqc, cov_bamcount)
cov_merged$stage <- 'raw'

## Check that it's all looking good
stopifnot(length(unique(cov_merged$type)) == 9)
stopifnot(!any(cov_merged$auc == 0))

cov_summarized <- do.call(rbind, lapply(split(cov_merged, cov_merged$sample), function(dat) {
    
    types <- c('initial', 'post-trim', 'mapped', 'no-dups')
    res <- data.frame(
        sample = unique(dat$sample),
        type = factor(types, levels = types),
        auc = c(
            sum(dat$auc[dat$type %in% c('combined_R1', 'combined_R2')]),
            sum(dat$auc[grepl('flash|unpaired', dat$type)]),
            dat$auc[dat$type == 'Marked_duplicates'],
            dat$auc[dat$type == 'duplicatesRemoved']
            ),
        stage = 'summarized',
        stringsAsFactors = FALSE
    )
    return(res)
}))
rownames(cov_summarized) <- NULL

## Compute the genome (hg19) coverage
cov_merged$coverage <- cov_merged$auc / hg19.size
cov_summarized$coverage <- cov_summarized$auc / hg19.size

save(cov_merged, file = 'rda/cov_merged_v0.4.0.Rdata')

## Load pheno data
load('/dcl01/lieber/ajaffe/lab/brain-epigenomics/misc/rdas/rse_mean_meth.Rdata', verbose = FALSE)
pd <- colData(rse_mean_meth[[1]])


## Wide format
combine_topd <- function(covdat) {
    pd_new <- pd
    for(var in c('auc', 'coverage')) {
        info_cov <- tidyr::spread(covdat[, c('type', 'sample', var)], type, var)
        m <- match(rownames(pd), info_cov$sample)
        to_add <- info_cov[m, -1]
        colnames(to_add) <- paste0(var, '_', gsub('-', '_', colnames(to_add)))
        pd_new <- cbind(pd_new, to_add)
    }
    return(pd_new)
}
pd_exp <- combine_topd(cov_summarized)


save(pd_exp, file = 'rda/pd_exp_v0.4.0.Rdata')

## Long format
combine_pd <- function(covdat) {
    m <- match(covdat$sample, rownames(pd))
    res <- cbind(covdat, pd[m, ])
    return(as.data.frame(res))
}
cov_summarized <- combine_pd(cov_summarized)

## Add age group based on 
## https://github.com/LieberInstitute/brain-epigenomics/blob/master/bumphunting/plot_bumps_bsseqSmooth.R#L170-L174
cov_summarized$age_group <- factor(
    ifelse(cov_summarized$Age < 0, 'Prenatal',
    ifelse(cov_summarized$Age < 1, 'Infant',
    ifelse(cov_summarized$Age <= 12, 'Child',
    ifelse(cov_summarized$Age <= 17, 'Teen', 'Adult')))),
    levels = c('Infant', 'Child', 'Teen', 'Adult', 'Prenatal'))


cov_summarized$age_group_cell <- factor(paste0(cov_summarized$age_group, '_',
    cov_summarized$Cell.Type),
    levels = c(paste0(rep(levels(cov_summarized$age_group)[1:4], each = 2),
    '_', c('Glia', 'Neuron')), 'Prenatal_H'))
cov_summarized$col <- c(brewer.pal(8, "Paired"), 'grey50')[c(5:6, 7:8, 3:4, 1:2, 9)][cov_summarized$age_group_cell]
    
save(cov_summarized, file = 'rda/cov_summarized_v0.4.0.Rdata')

# summary(lm(coverage_initial ~ Age, data = pd_exp))
model_fits <- list(
    'age' = list(
        'initial' = lm(coverage_initial ~ Age + Cell.Type, data = pd_exp),
        'post_trim' = lm(coverage_post_trim ~ Age + Cell.Type, data = pd_exp),
        'mapped' = lm(coverage_mapped ~ Age + Cell.Type, data = pd_exp),
        'no_dups' = lm(coverage_no_dups ~ Age + Cell.Type, data = pd_exp)
    ),
    'cell' = list(
        'initial' = lm(coverage_initial ~ Cell.Type + Age, data = pd_exp),
        'post_trim' = lm(coverage_post_trim ~ Cell.Type + Age, data = pd_exp),
        'mapped' = lm(coverage_mapped ~ Cell.Type + Age, data = pd_exp),
        'no_dups' = lm(coverage_no_dups ~ Cell.Type + Age, data = pd_exp)
    ),
    'interaction' = list(
        'initial' = lm(coverage_initial ~ Age * Cell.Type, data = pd_exp),
        'post_trim' = lm(coverage_post_trim ~ Age * Cell.Type, data = pd_exp),
        'mapped' = lm(coverage_mapped ~ Age * Cell.Type, data = pd_exp),
        'no_dups' = lm(coverage_no_dups ~ Age * Cell.Type, data = pd_exp)
    )
)

## Tidy the model outputs
# lapply(model_fits, function(fits) {
#     lapply(fits, tidy)
# })

model_summ <- mapply(function(mods, coef) {
    sapply(mods, function(fit) {
        tidy(fit)$p.value[coef]
    })
}, model_fits, c(2, 2, 4))
model_summ
#                  age      cell interaction
# initial   0.84171865 0.1380891   0.8682455
# post_trim 0.70234593 0.6721195   0.2388473
# mapped    0.78116853 0.7728281   0.3289814
# no_dups   0.08693054 0.3105985   0.2299959


model_summ < 0.05 / 4
#             age  cell interaction
# initial   FALSE FALSE       FALSE
# post_trim FALSE FALSE       FALSE
# mapped    FALSE FALSE       FALSE
# no_dups   FALSE FALSE       FALSE


# > with(subset(cov_summarized, Cell.Type == 'Neuron'), tapply(col, age_group, unique))
#    Infant     Child      Teen     Adult  Prenatal
# "#E31A1C" "#FF7F00" "#33A02C" "#1F78B4"        NA
neun_colors <- with(subset(cov_summarized, Cell.Type == 'Neuron'), tapply(col, age_group, unique))
neun_colors['Prenatal'] <- 'grey50'

age_cell_colors <- with(cov_summarized, tapply(col, age_group_cell, unique))
age_cell_colors['Prenatal'] <- 'grey50'

dir.create('pdf', showWarnings = FALSE)
pdf('pdf/coverage_across_processing_stages_v0.4.0.pdf', width = 11, useDingbats = FALSE)

## Basic
ggplot(cov_summarized, aes(x = type, y = coverage))  + geom_boxplot() + theme_bw(base_size = 20) + xlab('Processing Stage') + ylab('Genome Coverage')

## By cell type
ggplot(cov_summarized, aes(x = type, y = coverage, colour = Cell.Type))  + geom_boxplot() + scale_color_brewer(palette = 'Dark2') + theme_bw(base_size = 20) + xlab('Processing Stage') + ylab('Genome Coverage')

## By age group
ggplot(cov_summarized, aes(x = type, y = coverage, colour = age_group))  + geom_boxplot() + scale_color_manual( values = neun_colors ) + theme_bw(base_size = 20) + xlab('Processing Stage') + ylab('Genome Coverage')

## By age group and cell type
ggplot(cov_summarized, aes(x = type, y = coverage, colour = age_group_cell)) + geom_boxplot() + scale_color_manual( values = age_cell_colors ) + theme_bw(base_size = 20) + xlab('Processing Stage') + ylab('Genome Coverage')

## Using facet by cell type
ggplot(cov_summarized, aes(x = type, y = coverage, colour = age_group_cell)) + geom_boxplot() + scale_color_manual( values = age_cell_colors ) + theme_bw(base_size = 20) + facet_grid(.~Cell.Type) + xlab('Processing Stage') + ylab('Genome Coverage')

## Same but with the simpler colors (no light colors)
ggplot(cov_summarized, aes(x = type, y = coverage, colour = age_group)) + geom_boxplot() + scale_color_manual( values = neun_colors ) + theme_bw(base_size = 20) + facet_grid(.~Cell.Type) + xlab('Processing Stage') + ylab('Genome Coverage')

## Using facet by age group
ggplot(cov_summarized, aes(x = type, y = coverage, colour = age_group_cell)) + geom_boxplot() + scale_color_manual( values = age_cell_colors ) + theme_bw(base_size = 20) + facet_grid(age_group ~ .) + xlab('Processing Stage') + ylab('Genome Coverage')

## Same but with the simpler colors (one by cell type)
ggplot(cov_summarized, aes(x = type, y = coverage, colour = Cell.Type))  + geom_boxplot() + facet_grid(age_group ~ .) + scale_color_brewer(palette = 'Dark2') + theme_bw(base_size = 20) + xlab('Processing Stage') + ylab('Genome Coverage')

## Lines per sample
ggplot(cov_summarized, aes(x = type, y = coverage, group = sample, colour = Cell.Type)) + geom_line() + geom_point() + scale_color_brewer(palette = 'Dark2') + theme_bw(base_size = 20) + xlab('Processing Stage') + ylab('Genome Coverage')

## Same but colored by cell type and split by age group
ggplot(cov_summarized, aes(x = type, y = coverage, group = sample, colour = Cell.Type)) + geom_line() + geom_point() + facet_grid(age_group ~ .) + scale_color_brewer(palette = 'Dark2') + theme_bw(base_size = 20) + xlab('Processing Stage') + ylab('Genome Coverage')
dev.off()

## Resume if needed
if(FALSE) {
    load('rda/cov_bamcount.Rdata', verbose = TRUE)
    load('rda/cov_fastqc.Rdata', verbose = TRUE)
    load('rda/cov_merged.Rdata', verbose = TRUE)
    load('rda/cov_summarized.Rdata', verbose = TRUE)
    load('rda/pd_exp.Rdata', verbose = TRUE)
}

## Check vs prior results
cov_summarized_v0.4.0 <- cov_summarized
load('rda/cov_summarized.Rdata', verbose = TRUE)

## It's identical =)
stopifnot(identical(cov_summarized$auc, cov_summarized_v0.4.0$auc))
stopifnot(identical(cov_summarized$coverage, cov_summarized_v0.4.0$coverage))

## Could have simply checked this earlier :P
cov_bamcount_v0.4.0 <- cov_bamcount
load('rda/cov_bamcount.Rdata', verbose = TRUE)
stopifnot(identical(cov_bamcount_v0.4.0, cov_bamcount))


## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()

# ─ Session info ───────────────────────────────────────────────────────────────────────────────────────────────────────
#  setting  value
#  version  R version 3.5.3 Patched (2019-03-11 r76311)
#  os       CentOS Linux 7 (Core)
#  system   x86_64, linux-gnu
#  ui       X11
#  language (EN)
#  collate  en_US.UTF-8
#  ctype    en_US.UTF-8
#  tz       US/Eastern
#  date     2019-09-04
#
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
#  package              * version   date       lib source
#  assertthat             0.2.1     2019-03-21 [2] CRAN (R 3.5.1)
#  backports              1.1.3     2018-12-14 [2] CRAN (R 3.5.1)
#  Biobase              * 2.42.0    2018-10-30 [2] Bioconductor
#  BiocGenerics         * 0.28.0    2018-10-30 [1] Bioconductor
#  BiocParallel         * 1.16.6    2019-02-10 [1] Bioconductor
#  bitops                 1.0-6     2013-08-17 [2] CRAN (R 3.5.0)
#  broom                * 0.5.2     2019-04-07 [1] CRAN (R 3.5.3)
#  cli                    1.1.0     2019-03-19 [1] CRAN (R 3.5.3)
#  colorout             * 1.2-0     2018-05-02 [1] Github (jalvesaq/colorout@c42088d)
#  colorspace             1.4-1     2019-03-18 [2] CRAN (R 3.5.1)
#  crayon                 1.3.4     2017-09-16 [1] CRAN (R 3.5.0)
#  DelayedArray         * 0.8.0     2018-10-30 [2] Bioconductor
#  digest                 0.6.20    2019-07-04 [1] CRAN (R 3.5.3)
#  dplyr                  0.8.3     2019-07-04 [1] CRAN (R 3.5.3)
#  generics               0.0.2     2018-11-29 [1] CRAN (R 3.5.1)
#  GenomeInfoDb         * 1.18.2    2019-02-12 [1] Bioconductor
#  GenomeInfoDbData       1.2.0     2018-11-02 [2] Bioconductor
#  GenomicRanges        * 1.34.0    2018-10-30 [1] Bioconductor
#  ggplot2              * 3.2.1     2019-08-10 [1] CRAN (R 3.5.3)
#  glue                   1.3.1     2019-03-12 [1] CRAN (R 3.5.1)
#  gtable                 0.3.0     2019-03-25 [2] CRAN (R 3.5.1)
#  hms                    0.4.2     2018-03-10 [2] CRAN (R 3.5.0)
#  htmltools              0.3.6     2017-04-28 [2] CRAN (R 3.5.0)
#  htmlwidgets            1.3       2018-09-30 [1] CRAN (R 3.5.1)
#  httpuv                 1.5.0     2019-03-15 [2] CRAN (R 3.5.1)
#  IRanges              * 2.16.0    2018-10-30 [1] Bioconductor
#  jaffelab             * 0.99.21   2018-05-03 [1] Github (LieberInstitute/jaffelab@7ed0ab7)
#  jsonlite               1.6       2018-12-07 [2] CRAN (R 3.5.1)
#  labeling               0.3       2014-08-23 [2] CRAN (R 3.5.0)
#  later                  0.8.0     2019-02-11 [2] CRAN (R 3.5.1)
#  lattice                0.20-38   2018-11-04 [3] CRAN (R 3.5.3)
#  lazyeval               0.2.2     2019-03-15 [2] CRAN (R 3.5.1)
#  limma                  3.38.3    2018-12-02 [1] Bioconductor
#  magrittr               1.5       2014-11-22 [1] CRAN (R 3.5.0)
#  Matrix                 1.2-15    2018-11-01 [3] CRAN (R 3.5.3)
#  matrixStats          * 0.54.0    2018-07-23 [1] CRAN (R 3.5.1)
#  munsell                0.5.0     2018-06-12 [2] CRAN (R 3.5.1)
#  nlme                   3.1-137   2018-04-07 [3] CRAN (R 3.5.3)
#  pillar                 1.4.2     2019-06-29 [1] CRAN (R 3.5.3)
#  pkgconfig              2.0.2     2018-08-16 [1] CRAN (R 3.5.1)
#  plyr                   1.8.4     2016-06-08 [2] CRAN (R 3.5.0)
#  png                    0.1-7     2013-12-03 [2] CRAN (R 3.5.0)
#  promises               1.0.1     2018-04-13 [2] CRAN (R 3.5.0)
#  purrr                  0.3.2     2019-03-15 [2] CRAN (R 3.5.1)
#  R6                     2.4.0     2019-02-14 [2] CRAN (R 3.5.1)
#  rafalib              * 1.0.0     2015-08-09 [1] CRAN (R 3.5.0)
#  RColorBrewer         * 1.1-2     2014-12-07 [2] CRAN (R 3.5.0)
#  Rcpp                   1.0.2     2019-07-25 [1] CRAN (R 3.5.3)
#  RCurl                  1.95-4.12 2019-03-04 [2] CRAN (R 3.5.1)
#  readr                * 1.3.1     2018-12-21 [1] CRAN (R 3.5.1)
#  reshape2               1.4.3     2017-12-11 [2] CRAN (R 3.5.0)
#  rlang                  0.4.0     2019-06-25 [1] CRAN (R 3.5.3)
#  rmote                * 0.3.4     2018-05-02 [1] deltarho (R 3.5.0)
#  S4Vectors            * 0.20.1    2018-11-09 [1] Bioconductor
#  scales                 1.0.0     2018-08-09 [2] CRAN (R 3.5.1)
#  segmented              0.5-3.0   2017-11-30 [2] CRAN (R 3.5.0)
#  servr                  0.15      2019-08-07 [1] CRAN (R 3.5.3)
#  sessioninfo          * 1.1.1     2018-11-05 [1] CRAN (R 3.5.1)
#  stringi                1.4.3     2019-03-12 [2] CRAN (R 3.5.1)
#  stringr                1.4.0     2019-02-10 [1] CRAN (R 3.5.1)
#  SummarizedExperiment * 1.12.0    2018-10-30 [1] Bioconductor
#  tibble                 2.1.3     2019-06-06 [1] CRAN (R 3.5.3)
#  tidyr                * 0.8.3     2019-03-01 [2] CRAN (R 3.5.1)
#  tidyselect             0.2.5     2018-10-11 [2] CRAN (R 3.5.1)
#  withr                  2.1.2     2018-03-15 [2] CRAN (R 3.5.0)
#  xfun                   0.9       2019-08-21 [1] CRAN (R 3.5.3)
#  XVector                0.22.0    2018-10-30 [1] Bioconductor
#  zlibbioc               1.28.0    2018-10-30 [2] Bioconductor
#
# [1] /users/lcollado/R/x86_64-pc-linux-gnu-library/3.5.x
# [2] /jhpce/shared/jhpce/core/conda/miniconda-3/envs/svnR-3.5.x/R/3.5.x/lib64/R/site-library
# [3] /jhpce/shared/jhpce/core/conda/miniconda-3/envs/svnR-3.5.x/R/3.5.x/lib64/R/library
