## Required libraries
library('bsseq')
library('devtools')
library('data.table')
library('GenomicRanges')
library('BiocParallel')

## Load pheno data
pd <- read.table('/dcl01/lieber/WGBS/ListerData/SraRunTable.txt',
    header = TRUE, sep = '\t')
colnames(pd) <- tolower(gsub('_s$|_l$', '', colnames(pd)))

## Keep only unique entries
pd <- pd[sapply(unique(pd$sample_name), function(x) { which(pd$sample_name == x)[1] }), ]
rownames(pd) <- NULL

## Find bismark report files for CpGs
pd$reportFiles <- file.path(
    '/dcl01/lieber/WGBS/ListerData/DNAm_Ratios_duplicates_dropped/CpGs',
	pd$sample_name, paste0(pd$sample_name,
    '_sorted_and_duplicates_dropped.bismark.cov'))
stopifnot(all(file.exists(pd$reportFiles)))

BSobj_raw <- combineList(bpmapply(function(input_file, sample) {
    library('bsseq')
    read.bismark(input_file, sample, strandCollapse = FALSE,
        fileType = 'cov')
}, pd$reportFiles, pd$sample_name, SIMPLIFY = FALSE,
    BPPARAM = SerialParam()))


## Filter based on CpGs with minCov >=3 
message(paste(Sys.time(), 'filtering CpGs based on those with min coverage 3 in our data'))
load('/dcl01/lieber/ajaffe/lab/brain-epigenomics/bsseq/lister/gr_cps_minCov_3.Rdata')
present <- countOverlaps(rowRanges(BSobj_raw), gr_cpgs) > 0
table(present)
BSobj <- BSobj_raw[present, ]
save(BSobj, file = 'BSobj_lister_minCov_3.Rdata')
rm(BSobj, BSobj_raw, present)


## Now combine
message(paste(Sys.time(), 'Combining nonCpGs results for the chrs'))
files <- dir(pattern = '_BSobj_lister_nonCG_highCov.Rdata')

## Load and filter:
load_filt <- function(f) {
    message(paste(Sys.time(), 'loading', f))
    load(f)
    ## Tweak reportFiles so they can be combined
    colData(BSobj)$reportFiles <- gsub('chr.*', '', colData(BSobj)$reportFiles)
    return(BSobj)
}

BSobj <- do.call(rbind, lapply(files, load_filt))

## Save final combined result
save(BSobj, file = 'allChrs_lister_nonCG_highCov.Rdata')

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
