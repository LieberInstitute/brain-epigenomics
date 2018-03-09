## Required libraries
library('getopt')
library('bsseq')
library('readxl')
library('devtools')
library('data.table')
library('GenomicRanges')
library('BiocParallel')

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
pd1 <- read_excel(
    '/users/ajaffe/Lieber/Projects/WGBS/Analysis/WGBS-pd-fixed.xlsx')
pd1 <- pd1[,1:36] # fix excel issues
pd1 <- pd1[!is.na(pd1$Data.ID), ]

pd2 <- read.csv('/dcl01/lieber/WGBS/LIBD_Data/phenoData/FullMasterList_clean.csv')

## Add missing samples
miss1 <- grep('reseq', pd1$Data.ID)
miss2 <- which(pd2$WGC.ID %in% gsub('_reseq', '', pd1$Data.ID[miss1]) )

pd3 <- pd2[miss2, , drop = FALSE]
pd3$WGC.ID <- pd1$Data.ID[miss1]
pd3[, c('avg.Cov', 'cov.sDev', 'Percent.Duplication', 'num.untrimmed.reads', 'num.trimmed.reads', 'alignment.efficiency')] <- NA

## Combine
pd <- rbind(pd2, pd3)
pd$Cell.Type <- ifelse(pd$Cell.Type == 'H', 'Homogenate', ifelse(pd$Cell.Type == 'NeuN-', 'Glia', 'Neuron'))

print('before filtering')
with(pd, table('postnatal' = Age > 0, Cell.Type, useNA = 'ifany'))

## keep only phomogenate
pd <- pd[pd$Cell.Type == 'Homogenate', ]
print('after filtering')
with(pd, table('postnatal' = Age > 0, Cell.Type, useNA = 'ifany'))

## Rename to match previous code
colnames(pd)[colnames(pd) == 'WGC.ID'] <- 'Data.ID'

## all report files
pd$reportFiles <- file.path(
    '/dcl01/lieber/ajaffe/lab/brain-epigenomics/bsseq/bsobj_by_chr/Reports',
	pd$Data.ID, paste0(pd$Data.ID,
    '.concatenated.sorted.duplicatesRemoved.CX_reportchr', opt$chr, '.txt'))
stopifnot(all(file.exists(pd$reportFiles)))

## Load the data without collapsing by strand
BSobj <- combineList(bpmapply(function(input_file, sample) {
    library('bsseq')
    read.bismark(input_file, sample, strandCollapse = FALSE,
        fileType = 'cytosineReport')
}, pd$reportFiles, pd$Data.ID, SIMPLIFY = FALSE,
    BPPARAM = SnowParam(opt$cores, outfile =  Sys.getenv('SGE_STDERR_PATH'))))

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
dir.create('rda', showWarnings = FALSE)
save(BSobj, file = paste0('rda/', opt$chr,
    '_cleaned_CX_Homogenate.Rdata'))

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
