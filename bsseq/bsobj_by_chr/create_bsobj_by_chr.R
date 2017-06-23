## Required libraries
library('getopt')
library('bsseq')
library('readxl')
library('devtools')
library('data.table')
library('GenomicRanges')

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

## For testing
if(FALSE) {
    opt <- list(chr = 'chr21', cores = 1)
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
    '/dcl01/lieber/ajaffe/lab/brain-epigenomics/bsseq/bsobj_by_chr/Reports',
	pd$Data.ID, paste0(pd$Data.ID,
    '.concatenated.sorted.duplicatesRemoved.CX_reportchr', opt$chr, '.txt'))
stopifnot(all(file.exists(pd$reportFiles)))

## Load the data without collapsing by strand
BSobj <- read.bismark(pd$reportFiles, pd$Data.ID, strandCollapse = FALSE,
    fileType = 'cytosineReport', mc.cores = opt$cores)

## Get CX context
context <- fread(pd$reportFiles[1], colClasses = c('factor', 'numeric',
    'factor', 'integer', 'integer', 'factor', 'factor'))
context_gr <- GRanges(seqnames = opt$chr,
    IRanges(start = context[[2]], width = 1), strand = context[[3]],
    c_context = Rle(context[[6]]), trinucleotide_context = Rle(context[[7]]))
    
## Re-order based on strand
context_gr <- c(context_gr[strand(context_gr) == '+'],
    context_gr[strand(context_gr) == '-'])
stopifnot(identical(ranges(granges(BSobj)), ranges(context_gr)))

## Same the context information on the BSobj
rowRanges(BSobj) <- context_gr

## append phenotype data
rownames(pd) <- pd$Data.ID
pData(BSobj) <- DataFrame(pd)
save(BSobj, file = paste0(opt$chr, '_postNatal_cleaned_CX_noHomogenate.Rdata'))

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
