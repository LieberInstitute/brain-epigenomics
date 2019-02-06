# module load bedtools/2.24.0
# wget https://raw.githubusercontent.com/arq5x/bedtools/master/genomes/human.hg19.genome
# bedtools genomecov -ibam /dcl02/lieber/WGBS/LIBD_Data/BAM/WGC052296L.concatenated.sorted.duplicatesRemoved.bam > test_without_genome
# bedtools genomecov -ibam /dcl02/lieber/WGBS/LIBD_Data/BAM/WGC052296L.concatenated.sorted.duplicatesRemoved.bam -g human.hg19.genome > test_with_genome

with <- readr::read_tsv('test_with_genome', col_names = c('chr', 'depth', 'nbases', 'chrsize', 'fraction'))
without <- readr::read_tsv('test_without_genome', col_names = c('chr', 'depth', 'nbases', 'chrsize', 'fraction'))

## Nice, no need for the genome then
stopifnot(identical(with, without))

table( (without$nbases / without$chrsize) - without$fraction < 1e-6)
#
# TRUE
# 5022

## Compute the AUC
n_genomecov <- with(subset(without, chr != 'genome'), sum(depth * nbases))
n_genomecov
## Exactly the same number Chris Wilks calculated with bamcount! ^^
# [1] 4579378967

## Compute the coverage
n_genomecov / subset(without, chr == 'genome')$chrsize[1]
# [1] 1.479274


# module load bamcount/0.2.6
# bamcount /dcl02/lieber/WGBS/LIBD_Data/BAM/WGC052296L.concatenated.sorted.duplicatesRemoved.bam --threads 4 --no-head --auc test_auc
n_bamcount <- readr::read_tsv('test_auc.auc.tsv', col_names = c('cat', 'n'))$n
stopifnot(identical(n_genomecov, n_bamcount))

hg19 <- readr::read_tsv('human.hg19.genome', col_names = c('chr', 'length'))
hg19.size <- sum(hg19$length[hg19$chr %in% paste0('chr', c(1:22, 'X', 'Y', 'M'))])
stopifnot(identical(
    subset(without, chr == 'genome')$chrsize[1],
    hg19.size
))

n_bamcount / hg19.size


qcData <- c('/dcl02/lieber/WGBS/LIBD_Data/FastQC/Untrimmed/WGC052296L/WGC052296L_combined_R1_fastqc/fastqc_data.txt', '/dcl02/lieber/WGBS/LIBD_Data/FastQC/Untrimmed/WGC052296L/WGC052296L_combined_R2_fastqc/fastqc_data.txt', '/dcl02/lieber/WGBS/LIBD_Data/FastQC/Trimmed/WGC052296L/WGC052296L_flashOutput.extendedFrags_fastqc/fastqc_data.txt')


            
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

auc <- compute_fastqc_auc(qcData)

auc / hg19.size

sum(auc[1:2]) / hg19.size

qcData_dir <- dir(c('/dcl02/lieber/WGBS/LIBD_Data/FastQC/Untrimmed/WGC052296L', '/dcl02/lieber/WGBS/LIBD_Data/FastQC/Trimmed/WGC052296L'), pattern = 'data.txt$', recursive = TRUE, full.names = TRUE)
qcData_dir

names(qcData_dir) <- gsub('^[[:alnum:]]*_|_fastqc', '', ss(qcData_dir, '/', 9))

auc_dir <- compute_fastqc_auc(qcData_dir)
auc_dir / hg19.size

dir('/dcl02/lieber/WGBS/LIBD_Data/BAM/', pattern = 'WGC052296L.*bam', full.names = TRUE)
# bamcount /dcl02/lieber/WGBS/LIBD_Data/BAM/WGC052296L.concatenated.sorted.duplicatesRemoved.bam --threads 4 --no-head --auc test_auc
# bamcount /dcl02/lieber/WGBS/LIBD_Data/BAM//WGC052296L_concatenated_Marked_duplicates.bam --threads 4 --no_head --auc test_auc_all

read_bamcount <- function(f) {
    readr::read_tsv(f, col_names = c('cat', 'n'))$n
}
auc_files <- dir(pattern = 'auc.tsv$')
auc_bam <- sapply(auc_files, read_bamcount)
auc_bam

c(auc_dir, auc_bam) / hg19.size
# flashOutput.extendedFrags flashOutput.notCombined_1 flashOutput.notCombined_2
#                 3.4392450                 0.6976701                 0.6253410
#   output_forward_unpaired   output_reverse_unpaired               combined_R1
#                 8.3080693                 0.1219993                19.3254085
#               combined_R2      test_auc_all.auc.tsv          test_auc.auc.tsv
#                19.3254085                 9.7240339                 1.4792738

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()

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
#  date     2019-02-05
#
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
#  package     * version date       lib source
#  assertthat    0.2.0   2017-04-11 [2] CRAN (R 3.5.0)
#  bindr         0.1.1   2018-03-13 [1] CRAN (R 3.5.0)
#  bindrcpp      0.2.2   2018-03-29 [1] CRAN (R 3.5.0)
#  cli           1.0.1   2018-09-25 [1] CRAN (R 3.5.1)
#  colorout    * 1.2-0   2018-05-02 [1] Github (jalvesaq/colorout@c42088d)
#  colorspace    1.4-0   2019-01-13 [2] CRAN (R 3.5.1)
#  crayon        1.3.4   2017-09-16 [1] CRAN (R 3.5.0)
#  digest        0.6.18  2018-10-10 [1] CRAN (R 3.5.1)
#  dplyr         0.7.8   2018-11-10 [1] CRAN (R 3.5.1)
#  ggplot2       3.1.0   2018-10-25 [1] CRAN (R 3.5.1)
#  glue          1.3.0   2018-07-17 [1] CRAN (R 3.5.1)
#  gtable        0.2.0   2016-02-26 [2] CRAN (R 3.5.0)
#  hms           0.4.2   2018-03-10 [2] CRAN (R 3.5.0)
#  htmltools     0.3.6   2017-04-28 [2] CRAN (R 3.5.0)
#  htmlwidgets   1.3     2018-09-30 [1] CRAN (R 3.5.1)
#  httpuv        1.4.5.1 2018-12-18 [2] CRAN (R 3.5.1)
#  later         0.7.5   2018-09-18 [2] CRAN (R 3.5.1)
#  lattice       0.20-38 2018-11-04 [3] CRAN (R 3.5.1)
#  lazyeval      0.2.1   2017-10-29 [2] CRAN (R 3.5.0)
#  magrittr      1.5     2014-11-22 [1] CRAN (R 3.5.0)
#  munsell       0.5.0   2018-06-12 [2] CRAN (R 3.5.0)
#  pillar        1.3.1   2018-12-15 [1] CRAN (R 3.5.1)
#  pkgconfig     2.0.2   2018-08-16 [1] CRAN (R 3.5.1)
#  plyr          1.8.4   2016-06-08 [2] CRAN (R 3.5.0)
#  png           0.1-7   2013-12-03 [2] CRAN (R 3.5.0)
#  promises      1.0.1   2018-04-13 [2] CRAN (R 3.5.0)
#  purrr         0.2.5   2018-05-29 [2] CRAN (R 3.5.0)
#  R6            2.3.0   2018-10-04 [2] CRAN (R 3.5.1)
#  Rcpp          1.0.0   2018-11-07 [1] CRAN (R 3.5.1)
#  readr         1.3.1   2018-12-21 [1] CRAN (R 3.5.1)
#  rlang         0.3.1   2019-01-08 [1] CRAN (R 3.5.1)
#  rmote       * 0.3.4   2018-05-02 [1] deltarho (R 3.5.0)
#  scales        1.0.0   2018-08-09 [2] CRAN (R 3.5.1)
#  servr         0.11    2018-10-23 [1] CRAN (R 3.5.1)
#  sessioninfo   1.1.1   2018-11-05 [1] CRAN (R 3.5.1)
#  tibble        2.0.1   2019-01-12 [1] CRAN (R 3.5.1)
#  tidyselect    0.2.5   2018-10-11 [2] CRAN (R 3.5.1)
#  withr         2.1.2   2018-03-15 [2] CRAN (R 3.5.0)
#  xfun          0.4     2018-10-23 [1] CRAN (R 3.5.1)
#
# [1] /users/lcollado/R/x86_64-pc-linux-gnu-library/3.5.x
# [2] /jhpce/shared/jhpce/core/conda/miniconda-3/envs/svnR-3.5.x/R/3.5.x/lib64/R/site-library
# [3] /jhpce/shared/jhpce/core/conda/miniconda-3/envs/svnR-3.5.x/R/3.5.x/lib64/R/library
