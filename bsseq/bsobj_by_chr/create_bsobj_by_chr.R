## Required libraries
library('getopt')
library('bsseq')
library('readxl')
library('BiocParallel')
library('devtools')

## Specify parameters
spec <- matrix(c(
    'chr', 'c', 1, 'character', 'chromosome to subset',
    'cores', 't', 1, 'integer', 'Number of cores to use',
	'help' , 'h', 0, 'logical', 'Display help'
), byrow=TRUE, ncol=5)
opt <- getopt(spec)

## if help was asked for print a friendly message
## and exit with a non-zero error code
if (!is.null(opt$help)) {
	cat(getopt(spec, usage=TRUE))
	q(status=1)
}

# read in cleaned phenotype
pd <- read_excel(
    '/users/ajaffe/Lieber/Projects/WGBS/Analysis/WGBS-pd-fixed.xlsx')
pd <- pd[,1:36] # fix excel issues
pd <- pd[!is.na(pd$Data.ID), ]

## keep only postnatal
pd <- pd[pd$Age > 0,]

## Drop Homogenate
pd <- pd[pd$Cell.Type != 'Homogenate', ]

## all report files
pd$reportFiles <- file.path(
    '/dcl01/lieber/WGBS/LIBD_Data/DNAm_Ratios_duplicates_dropped/Reports',
	pd$Data.ID, paste0(pd$Data.ID,
    '.concatenated.sorted.duplicatesRemoved.CpG_report.txt'))

## Load the data
BSobj <- combineList(bplapply(pd$reportFiles, function(input_file) {
    res <- read.bismark(pd$reportFiles[input_file], pd$WGC.ID[input_file],
        strandCollapse=TRUE, fileType = 'cytosineReport')
    res <- res[seqnames(res) == opt$chr, ]
    return(res)
}, BPPARAM = MulticoreParam(opt$cores)))

## add sample names
sampleNames(BSobj) <- pd$Data.ID
rownames(pd) <- pd$Data.ID

## append phenotype data
pData(BSobj) <- DataFrame(pd)
save(BSobj, file = paste0(opt$chr, '_postNatal_cleaned_CpGonly_noHomogenate.Rdata'))

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
