## Based on /users/ajaffe/Lieber/Projects/WGBS/Analysis/test_bsseq.R

###
# qrsh -l bluejay,mem_free=20G,h_vmem=20G,h_fsize=200G -pe local 6
# cd /dcl01/lieber/ajaffe/lab/brain-epigenomics/bumphunting
# module load conda_R/3.4.x

library('bsseq')
library('genefilter')
library('GenomicRanges')
library('devtools')

# read in cleaned phenotype
pd  <- read.csv("/dcl01/lieber/WGBS/LIBD_Data/phenoData/FullMasterList_clean.csv",
    as.is=TRUE)
    
## Keep only prenatal ones
pd <- subset(pd, Age < 0)
dim(pd)

## all report files
pd$reportFiles  <- paste0("/dcl01/lieber/WGBS/LIBD_Data/DNAm_Ratios_duplicates_dropped/Reports/",
	pd$WGC.ID, "/", pd$WGC.ID, ".concatenated.sorted.duplicatesRemoved.CpG_report.txt")
stopifnot(all(file.exists(pd$reportFiles)))


## Read bismark files
bsList <- mclapply(seq_len(nrow(pd)), function(ii) {
    BSobj <- read.bismark(pd$reportFiles[ii], pd$WGC.ID[ii],
        strandCollapse=TRUE, fileType = "cytosineReport")
    cat(".")
	return(BSobj)
}, mc.cores=6, mc.preschedule=FALSE)


BSobj <- combineList(bsList)

## add sample names
stopifnot(identical(sampleNames(BSobj), pd$WGC.ID))
rownames(pd) <-pd$WGC.ID

## append phenotype data
pData(BSobj) <- DataFrame(pd)
dim(BSobj)
# [1] 28217448       20

## Load the filtered data from the postnatal samples
foo <- function() {
    load('BSobj_bsseqSmooth_Neuron_minCov_3.Rdata')
    return(BSobj)
}
postnatal <- foo()
dim(postnatal)
# [1] 18664892       32

## Subset to positions observed in the postnatal samples
ov <- findOverlaps(rowRanges(postnatal), rowRanges(BSobj))
BSobj <- BSobj[subjectHits(ov), ]
dim(BSobj)
# [1] 18664892       20

stopifnot(nrow(BSobj) == nrow(postnatal))

## Smooth
system.time(BSobj <- BSmooth(BSobj, mc.cores = 6, parallelBy = 'sample'))
#      user    system   elapsed
# 59147.641   553.797 12439.142

save(BSobj, file="BSobj_bsseqSmooth_Neuron_minCov_3_prenatal.Rdata")    
    
## Reproducibility info
proc.time()
message(Sys.time())
options(width = 120)
session_info()

# > Sys.time()
# [1] "2017-11-03 15:41:20 EDT"
# > ## Reproducibility info
# > proc.time()
#      user    system   elapsed
# 66662.044  1415.112 17654.373
# > message(Sys.time())
# 2017-11-03 15:41:20
# > options(width = 120)
# > session_info()
# Session info ----------------------------------------------------------------------------------------------------------
#  setting  value
#  version  R version 3.4.2 Patched (2017-10-12 r73550)
#  system   x86_64, linux-gnu
#  ui       X11
#  language (EN)
#  collate  en_US.UTF-8
#  tz       <NA>
#  date     2017-11-03
#
# Packages --------------------------------------------------------------------------------------------------------------
#  package              * version  date       source
#  annotate               1.55.0   2017-08-10 Bioconductor
#  AnnotationDbi          1.39.4   2017-10-16 Bioconductor
#  base                 * 3.4.2    2017-10-13 local
#  Biobase              * 2.37.2   2017-08-06 Bioconductor
#  BiocGenerics         * 0.23.4   2017-10-27 Bioconductor
#  bit                    1.1-12   2014-04-09 CRAN (R 3.4.1)
#  bit64                  0.9-7    2017-05-08 CRAN (R 3.4.1)
#  bitops                 1.0-6    2013-08-17 CRAN (R 3.4.1)
#  blob                   1.1.0    2017-06-17 CRAN (R 3.4.1)
#  bsseq                * 1.13.9   2017-10-20 Bioconductor
#  colorout             * 1.1-2    2017-08-10 Github (jalvesaq/colorout@020a14d)
#  colorspace             1.3-2    2016-12-14 CRAN (R 3.4.1)
#  compiler               3.4.2    2017-10-13 local
#  data.table             1.10.4-3 2017-10-27 CRAN (R 3.4.2)
#  datasets             * 3.4.2    2017-10-13 local
#  DBI                    0.7      2017-06-18 CRAN (R 3.4.1)
#  DelayedArray         * 0.3.21   2017-09-29 Bioconductor
#  devtools             * 1.13.3   2017-08-02 CRAN (R 3.4.1)
#  digest                 0.6.12   2017-01-27 CRAN (R 3.4.1)
#  genefilter           * 1.59.0   2017-08-10 Bioconductor
#  GenomeInfoDb         * 1.13.5   2017-10-16 Bioconductor
#  GenomeInfoDbData       0.99.1   2017-08-06 Bioconductor
#  GenomicRanges        * 1.29.15  2017-10-16 Bioconductor
#  graphics             * 3.4.2    2017-10-13 local
#  grDevices            * 3.4.2    2017-10-13 local
#  grid                   3.4.2    2017-10-13 local
#  gtools                 3.5.0    2015-05-29 CRAN (R 3.4.1)
#  IRanges              * 2.11.19  2017-10-16 Bioconductor
#  lattice                0.20-35  2017-03-25 CRAN (R 3.4.2)
#  limma                  3.33.14  2017-10-16 Bioconductor
#  locfit                 1.5-9.1  2013-04-20 CRAN (R 3.4.1)
#  Matrix                 1.2-11   2017-08-21 CRAN (R 3.4.2)
#  matrixStats          * 0.52.2   2017-04-14 CRAN (R 3.4.1)
#  memoise                1.1.0    2017-04-21 CRAN (R 3.4.1)
#  methods              * 3.4.2    2017-10-13 local
#  munsell                0.4.3    2016-02-13 CRAN (R 3.4.1)
#  parallel             * 3.4.2    2017-10-13 local
#  permute                0.9-4    2016-09-09 CRAN (R 3.4.1)
#  plyr                   1.8.4    2016-06-08 CRAN (R 3.4.1)
#  R.methodsS3            1.7.1    2016-02-16 CRAN (R 3.4.1)
#  R.oo                   1.21.0   2016-11-01 CRAN (R 3.4.1)
#  R.utils                2.5.0    2016-11-07 CRAN (R 3.4.1)
#  Rcpp                   0.12.13  2017-09-28 CRAN (R 3.4.1)
#  RCurl                  1.95-4.8 2016-03-01 CRAN (R 3.4.1)
#  rlang                  0.1.2    2017-08-09 CRAN (R 3.4.1)
#  RSQLite                2.0      2017-06-19 CRAN (R 3.4.1)
#  S4Vectors            * 0.15.14  2017-10-27 Bioconductor
#  scales                 0.5.0    2017-08-24 CRAN (R 3.4.1)
#  splines                3.4.2    2017-10-13 local
#  stats                * 3.4.2    2017-10-13 local
#  stats4               * 3.4.2    2017-10-13 local
#  SummarizedExperiment * 1.7.10   2017-09-29 Bioconductor
#  survival               2.41-3   2017-04-04 CRAN (R 3.4.2)
#  tibble                 1.3.4    2017-08-22 CRAN (R 3.4.1)
#  tools                  3.4.2    2017-10-13 local
#  utils                * 3.4.2    2017-10-13 local
#  withr                  2.0.0    2017-07-28 CRAN (R 3.4.1)
#  XML                    3.98-1.9 2017-06-19 CRAN (R 3.4.1)
#  xtable                 1.8-2    2016-02-05 CRAN (R 3.4.1)
#  XVector                0.17.2   2017-10-27 Bioconductor
#  zlibbioc               1.23.0   2017-08-06 Bioconductor
