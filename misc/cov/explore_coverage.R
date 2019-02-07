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

##### Functions:
## Code adapted from https://github.com/LieberInstitute/RNAseq-pipeline/blob/master/sh/create_count_objects-human.R#L209-L217
splitAt <- function(x, pos) unname(split(x, cumsum(seq_along(x) %in% pos)))

compute_fastqc_seqdist <- function(qcData) {
    R <- lapply(qcData, function(x) scan(x, what = "character", sep= "\n", 
    		quiet = TRUE, strip=TRUE) )
    ## Split list into sublists of metric categories
    zz <- lapply(R, function(x) splitAt(x, which(x==">>END_MODULE")+1))

    ## Extract the sequence length distribution
    seqlen <- lapply(zz, function(x) {
        ## access the sequence distribution part of the report
        y <- x[[8]]
        ## Drop header and module separating lines
        y[-c(1:2, length(y))]
    })

    seqdist <- lapply(seqlen, function(dat) {
    
        ## Process the read length section
        len <- ss(dat, "\t", 1)
    
        ## Sometimes there's a range of read lengths
        res <- do.call(rbind, mapply(function(x, hasdash) {
            if(!hasdash) {
                x <- as.numeric(x)
                result <- data.frame(
                    start = x,
                    end = x,
                    stringsAsFactors = FALSE
                )
            } else {
                ## Get the start/end of the range
                result <- data.frame(
                    start = as.numeric(ss(x, '-', 1)),
                    end = as.numeric(ss(x, '-', 2)),
                    stringsAsFactors = FALSE
                )
            }
            ## Compute the median read length based on the range info
            result$median <- median(c(result$start, result$end))
            return(result)
        }, len, grepl('-', len), SIMPLIFY = FALSE))
    
        res$n <- as.numeric(ss(dat, "\t", 2))
    
        return(res)
    })
    return(seqdist)
}
compute_fastqc_auc <- function(qcData) {
    auc <- sapply(compute_fastqc_seqdist(qcData), function(dat) {
        sum(dat$median * dat$n)
    })
    return(auc)
}

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
get_fastqc_data <- function(id) {
    message(paste(Sys.time(), 'processing sample', id))
    
    ## Build the paths to the FastQC output directories
    id_dirs <- file.path(
        '/dcl02/lieber/WGBS/LIBD_Data/FastQC',
        c('Untrimmed', 'Trimmed'),
    id)
    
    ## Find fastqc data files
    fastqc_files <- dir(id_dirs, pattern = 'data.txt$', recursive = TRUE, full.names = TRUE)
    stopifnot(length(fastqc_files) == 7)
    
    names(fastqc_files) <- gsub(paste0('^', id, '_|_fastqc'), '', ss(fastqc_files, '/', 9))

    auc_dir <- compute_fastqc_auc(fastqc_files)
    
    res <- data.frame(
        sample = id,
        type = names(auc_dir),
        auc = auc_dir,
        stringsAsFactors = FALSE
    )
    rownames(res) <- NULL
    return(res)
}

cov_fastqc <- do.call(rbind, lapply(ids, get_fastqc_data))
# 2019-02-07 11:08:48 processing sample WGC052316L
# 2019-02-07 11:09:09 processing sample WGC059613L
# 2019-02-07 11:09:19 processing sample WGC059614L
# 2019-02-07 11:09:30 processing sample WGC052317L
# 2019-02-07 11:09:50 processing sample WGC055558L
# 2019-02-07 11:10:04 processing sample WGC055559L
# 2019-02-07 11:10:15 processing sample WGC059596L
# 2019-02-07 11:10:28 processing sample WGC055561L
# 2019-02-07 11:10:46 processing sample WGC052318L
# 2019-02-07 11:11:03 processing sample WGC059588L
# 2019-02-07 11:11:27 processing sample WGC059608L
# 2019-02-07 11:11:45 processing sample WGC059594L
# 2019-02-07 11:12:01 processing sample WGC052309L_reseq
# 2019-02-07 11:12:01 processing sample WGC059592L
# 2019-02-07 11:12:19 processing sample WGC052311L
# 2019-02-07 11:12:39 processing sample WGC059612L
# 2019-02-07 11:12:53 processing sample WGC055562L
# 2019-02-07 11:13:07 processing sample WGC055573L
# 2019-02-07 11:13:20 processing sample WGC052324L
# 2019-02-07 11:13:30 processing sample WGC055575L
# 2019-02-07 11:13:51 processing sample WGC059600L
# 2019-02-07 11:14:17 processing sample WGC052328L
# 2019-02-07 11:14:28 processing sample WGC052314L
# 2019-02-07 11:14:42 processing sample WGC059598L
# 2019-02-07 11:14:53 processing sample WGC059603L
# 2019-02-07 11:15:14 processing sample WGC055577L
# 2019-02-07 11:15:31 processing sample WGC059607L
# 2019-02-07 11:15:55 processing sample WGC052312L
# 2019-02-07 11:16:10 processing sample WGC055570L
# 2019-02-07 11:16:22 processing sample WGC055568L
# 2019-02-07 11:16:34 processing sample WGC059601L
# 2019-02-07 11:16:51 processing sample WGC059604L

save(cov_fastqc, file = 'rda/cov_fastqc.Rdata')

## Read in bamcount output data
get_bamcount_data <- function(id) {
    message(paste(Sys.time(), 'processing sample', id))
    
    ## Find the files
    types <- c('duplicatesRemoved', 'Marked_duplicates')
    bamcount_files <- file.path(
        '/dcl01/lieber/ajaffe/lab/brain-epigenomics/misc/cov/auc_files',
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
save(cov_bamcount, file = 'rda/cov_bamcount.Rdata')

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

save(cov_merged, file = 'rda/cov_merged.Rdata')

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


save(pd_exp, file = 'rda/pd_exp.Rdata')

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
    
save(cov_summarized, file = 'rda/cov_summarized.Rdata')

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
pdf('pdf/coverage_across_processing_stages.pdf', width = 11, useDingbats = FALSE)

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
#  date     2019-02-07
#
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
#  package              * version   date       lib source
#  assertthat             0.2.0     2017-04-11 [2] CRAN (R 3.5.0)
#  backports              1.1.3     2018-12-14 [2] CRAN (R 3.5.1)
#  bindr                  0.1.1     2018-03-13 [1] CRAN (R 3.5.0)
#  bindrcpp               0.2.2     2018-03-29 [1] CRAN (R 3.5.0)
#  Biobase              * 2.42.0    2018-10-30 [2] Bioconductor
#  BiocGenerics         * 0.28.0    2018-10-30 [1] Bioconductor
#  BiocParallel         * 1.16.5    2019-01-04 [1] Bioconductor
#  bitops                 1.0-6     2013-08-17 [2] CRAN (R 3.5.0)
#  broom                * 0.5.1     2018-12-05 [1] CRAN (R 3.5.1)
#  cli                    1.0.1     2018-09-25 [1] CRAN (R 3.5.1)
#  colorout             * 1.2-0     2018-05-02 [1] Github (jalvesaq/colorout@c42088d)
#  colorspace             1.4-0     2019-01-13 [2] CRAN (R 3.5.1)
#  crayon                 1.3.4     2017-09-16 [1] CRAN (R 3.5.0)
#  DelayedArray         * 0.8.0     2018-10-30 [2] Bioconductor
#  digest                 0.6.18    2018-10-10 [1] CRAN (R 3.5.1)
#  dplyr                  0.7.8     2018-11-10 [1] CRAN (R 3.5.1)
#  fansi                  0.4.0     2018-10-05 [1] CRAN (R 3.5.1)
#  generics               0.0.2     2018-11-29 [1] CRAN (R 3.5.1)
#  GenomeInfoDb         * 1.18.1    2018-11-12 [1] Bioconductor
#  GenomeInfoDbData       1.2.0     2018-11-02 [2] Bioconductor
#  GenomicRanges        * 1.34.0    2018-10-30 [1] Bioconductor
#  ggplot2              * 3.1.0     2018-10-25 [1] CRAN (R 3.5.1)
#  glue                   1.3.0     2018-07-17 [1] CRAN (R 3.5.1)
#  gtable                 0.2.0     2016-02-26 [2] CRAN (R 3.5.0)
#  hms                    0.4.2     2018-03-10 [2] CRAN (R 3.5.0)
#  htmltools              0.3.6     2017-04-28 [2] CRAN (R 3.5.0)
#  htmlwidgets            1.3       2018-09-30 [1] CRAN (R 3.5.1)
#  httpuv                 1.4.5.1   2018-12-18 [2] CRAN (R 3.5.1)
#  IRanges              * 2.16.0    2018-10-30 [1] Bioconductor
#  jaffelab             * 0.99.21   2018-05-03 [1] Github (LieberInstitute/jaffelab@7ed0ab7)
#  labeling               0.3       2014-08-23 [2] CRAN (R 3.5.0)
#  later                  0.7.5     2018-09-18 [2] CRAN (R 3.5.1)
#  lattice                0.20-38   2018-11-04 [3] CRAN (R 3.5.1)
#  lazyeval               0.2.1     2017-10-29 [2] CRAN (R 3.5.0)
#  limma                  3.38.3    2018-12-02 [1] Bioconductor
#  magrittr               1.5       2014-11-22 [1] CRAN (R 3.5.0)
#  Matrix                 1.2-15    2018-11-01 [3] CRAN (R 3.5.1)
#  matrixStats          * 0.54.0    2018-07-23 [1] CRAN (R 3.5.1)
#  mime                   0.6       2018-10-05 [1] CRAN (R 3.5.1)
#  munsell                0.5.0     2018-06-12 [2] CRAN (R 3.5.0)
#  nlme                   3.1-137   2018-04-07 [3] CRAN (R 3.5.1)
#  pillar                 1.3.1     2018-12-15 [1] CRAN (R 3.5.1)
#  pkgconfig              2.0.2     2018-08-16 [1] CRAN (R 3.5.1)
#  plyr                   1.8.4     2016-06-08 [2] CRAN (R 3.5.0)
#  png                    0.1-7     2013-12-03 [2] CRAN (R 3.5.0)
#  promises               1.0.1     2018-04-13 [2] CRAN (R 3.5.0)
#  purrr                  0.2.5     2018-05-29 [2] CRAN (R 3.5.0)
#  R6                     2.3.0     2018-10-04 [2] CRAN (R 3.5.1)
#  rafalib              * 1.0.0     2015-08-09 [1] CRAN (R 3.5.0)
#  RColorBrewer         * 1.1-2     2014-12-07 [2] CRAN (R 3.5.0)
#  Rcpp                   1.0.0     2018-11-07 [1] CRAN (R 3.5.1)
#  RCurl                  1.95-4.11 2018-07-15 [2] CRAN (R 3.5.1)
#  readr                * 1.3.1     2018-12-21 [1] CRAN (R 3.5.1)
#  reshape2               1.4.3     2017-12-11 [2] CRAN (R 3.5.0)
#  rlang                  0.3.1     2019-01-08 [1] CRAN (R 3.5.1)
#  rmote                * 0.3.4     2018-05-02 [1] deltarho (R 3.5.0)
#  S4Vectors            * 0.20.1    2018-11-09 [1] Bioconductor
#  scales                 1.0.0     2018-08-09 [2] CRAN (R 3.5.1)
#  segmented              0.5-3.0   2017-11-30 [2] CRAN (R 3.5.0)
#  servr                  0.11      2018-10-23 [1] CRAN (R 3.5.1)
#  sessioninfo          * 1.1.1     2018-11-05 [1] CRAN (R 3.5.1)
#  stringi                1.2.4     2018-07-20 [2] CRAN (R 3.5.1)
#  stringr                1.3.1     2018-05-10 [1] CRAN (R 3.5.0)
#  SummarizedExperiment * 1.12.0    2018-10-30 [1] Bioconductor
#  tibble                 2.0.1     2019-01-12 [1] CRAN (R 3.5.1)
#  tidyr                * 0.8.2     2018-10-28 [2] CRAN (R 3.5.1)
#  tidyselect             0.2.5     2018-10-11 [2] CRAN (R 3.5.1)
#  utf8                   1.1.4     2018-05-24 [1] CRAN (R 3.5.0)
#  withr                  2.1.2     2018-03-15 [2] CRAN (R 3.5.0)
#  xfun                   0.4       2018-10-23 [1] CRAN (R 3.5.1)
#  XVector                0.22.0    2018-10-30 [1] Bioconductor
#  zlibbioc               1.28.0    2018-10-30 [2] Bioconductor
#
# [1] /users/lcollado/R/x86_64-pc-linux-gnu-library/3.5.x
# [2] /jhpce/shared/jhpce/core/conda/miniconda-3/envs/svnR-3.5.x/R/3.5.x/lib64/R/site-library
# [3] /jhpce/shared/jhpce/core/conda/miniconda-3/envs/svnR-3.5.x/R/3.5.x/lib64/R/library
