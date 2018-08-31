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
    opt <- list(chr = 'chr19', cores = 1)
}

# read in phenotype table

pd = read.table("/dcl02/lieber/WGBS/PublicData/Stroud_Greenberg_Cell_2017/FASTQ/WGBS/SraRunTable.txt", sep="\t", header=T)
pd = pd[,!colnames(pd) %in% c("background_strain","tissue")]
pd$background = "C57BL/6"
rownames(pd) = pd$Run

## all report files

pd$reportFiles = file.path('/dcl02/lieber/WGBS/PublicData/Stroud_Greenberg_Cell_2017/BAMS/WGBS/reports',
                           paste0(pd[pd$LibraryLayout=="PAIRED","Run"],"_1_bismark_bt2_pe.CX_report.txt.chr",opt$chr,".CX_report.txt.gz"))
stopifnot(all(file.exists(pd$reportFiles)))

## Load the data without collapsing by strand
BSobj <- combineList(bpmapply(function(input_file, sample) {
    library('bsseq')
    read.bismark(input_file, sample, strandCollapse = FALSE,
        fileType = 'cytosineReport')
}, pd$reportFiles, pd$Run, SIMPLIFY = FALSE,
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
pData(BSobj) <- DataFrame(pd)
dir.create('rda', showWarnings = FALSE)
save(BSobj, file = paste0('rda/', opt$chr,
    '_cleaned_CX_Stroud.Rdata'))

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
