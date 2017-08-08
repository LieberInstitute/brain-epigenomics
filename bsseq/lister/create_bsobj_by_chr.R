## Required libraries
library('getopt')
library('bsseq')
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

## Load pheno data
pd <- read.table('/dcl01/lieber/WGBS/ListerData/SraRunTable.txt',
    header = TRUE, sep = '\t')
colnames(pd) <- tolower(gsub('_s$|_l$', '', colnames(pd)))

## Keep only unique entries
pd <- pd[sapply(unique(pd$sample_name), function(x) { which(pd$sample_name == x)[1] }), ]
rownames(pd) <- NULL

## Find bismark report files for all Cs
pd$reportFiles <- file.path(
    '/dcl01/lieber/WGBS/ListerData/DNAm_Ratios_duplicates_dropped/All_Cs',
	pd$sample_name, paste0(pd$sample_name,
    '_sorted_and_duplicates_dropped.CX_report.chr', opt$chr, '.txt'))
stopifnot(all(file.exists(pd$reportFiles)))

BSobj_raw <- combineList(bpmapply(function(input_file, sample) {
    library('bsseq')
    read.bismark(input_file, sample, strandCollapse = FALSE,
        fileType = 'cytosineReport')
}, pd$reportFiles, pd$sample_name, SIMPLIFY = FALSE,
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
stopifnot(identical(ranges(granges(BSobj_raw)), ranges(context_gr)))

## Same the context information on the BSobj_raw
rowRanges(BSobj_raw) <- context_gr

## append phenotype data
rownames(pd) <- pd$sample_name
pData(BSobj_raw) <- DataFrame(pd)
save(BSobj_raw, file = paste0(opt$chr, '_BSobj_allCs_raw_lister.Rdata'))

## Filter based on our high coverage Cs
message(paste(Sys.time(), 'filtering allCs based on nonCGs with high coverage in our data'))
load('/dcl01/lieber/ajaffe/lab/brain-epigenomics/bsseq/lister/gr_all_highCov.Rdata')
present <- countOverlaps(context_gr, gr_all_highCov) > 0
table(present)
BSobj <- BSobj_raw[present, ]
save(BSobj, file = paste0(opt$chr, '_BSobj_lister_nonCG_highCov.Rdata'))


## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
