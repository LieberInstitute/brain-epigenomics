## Required libraries
library('getopt')
library('bsseq')
library('devtools')
library('GenomicRanges')
library('BiocParallel')
library('limma')

## Specify parameters
spec <- matrix(c(
    'chr', 'c', 1, 'character', 'chromosome to subset',
    'cores', 't', 1, 'integer', 'Number of cores to use',
    'context', 'x', 1, 'character', 'Either all, nonCG, CG, CHG, CHH',
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
    opt <- list(chr = 'chr21', cores = 1, context = 'CG')
}

stopifnot(opt$context %in% c('nonCG', 'CG', 'CHG', 'CHH'))

## Load the subset BSobj
load_gr <- function() {
    if(opt$context == 'CG') {
        load('/dcl01/lieber/ajaffe/lab/brain-epigenomics/bsseq/lister/gr_cps_minCov_3.Rdata', verbose = TRUE)
        return(gr_cpgs)
    } else {
        load('/dcl01/lieber/ajaffe/lab/brain-epigenomics/bsseq/lister/gr_all_highCov.Rdata', verbose = TRUE)
        return(gr_all_highCov)
    }
}

load('/dcl01/lieber/ajaffe/lab/brain-epigenomics/lister_phenotype.rda', verbose = TRUE)


load_BSobj <- function() {
    load(paste0('/dcl01/lieber/ajaffe/lab/brain-epigenomics/bsseq/lister/', opt$chr, '_BSobj_allCs_raw_lister.Rdata'), verbose = TRUE)
    
    ## Load the raw data then subset to bases we observed
    gr_ours <- load_gr()
    ov <- findOverlaps(gr_ours, rowRanges(BSobj_raw))
    
    ## Keep only the neurons and non-neurons
    tokeep <- match(pd$geo_accession[which(pd$cellType %in% c('neurons', 'non-neurons'))], colData(BSobj_raw)$sample_name)
    
    res <- BSobj_raw[subjectHits(ov), tokeep]
    colData(res)$Cell.Type <- pd$cellType[which(pd$cellType %in% c('neurons', 'non-neurons'))]
    
    if(opt$context == 'nonCG') {
        res <- res[rowRanges(res)$c_context != 'CG', ]
    } else {
       res <- res[rowRanges(res)$c_context == opt$context, ] 
    }
    return(res)
}

BSobj <- load_BSobj()
gr <- rowRanges(BSobj)

print('Final number of bases of interest and number of samples')
dim(BSobj)

## Genomic position based groups
library('derfinder')
pos <- coverage(gr)[[opt$chr]] > 0
clus <- derfinder:::.clusterMakerRle(pos, maxGap = 1000, ranges = TRUE)

## Size of the groups by number of bases of interest
print('Number of groups')
length(clus)
print('Size of the groups by number of bases of interest')
summary(width(clus))
split.idx <- Rle(seq_len(length(clus)), width(clus))

## Size of the groups by genomic range
gr_list <- split(gr, split.idx)
chunk_width <- width(range(gr_list))
print('Size of groups by genomic range')
summary(unlist(chunk_width))

## Split the methylation data
meth <- DataFrame(getMeth(BSobj, type = 'raw'))

## Replace NaN by 0
for(i in seq_len(ncol(meth))) {
    if(any(!is.finite(meth[[i]]))) {
        meth[[i]][!is.finite(meth[[i]])] <- 0
    }
}
rm(i)
meth_list <- split(meth, split.idx)


print('Keeping only groups with at least 5 observations by group')
table(width(clus) >= 5)
meth_list <- meth_list[width(clus) >= 5]
gr_list <- gr_list[width(clus) >= 5]


if(opt$cores == 1) {
    bpparam <- SerialParam()
} else {
    bpparam <- SnowParam(opt$cores, outfile =  Sys.getenv('SGE_STDERR_PATH'))
}


neurons <- colData(BSobj)$Cell.Type == 'neurons'

auto <- bplapply(meth_list, function(x) {
    ## Drop first row because it's always 1
    auto_res <- apply(x, 2, function(y) { acf(y, plot = FALSE, lag.max = 4)$acf })[-1, ]
    c('neuron' = rowMeans(auto_res[, neurons], na.rm = TRUE),
        'glia' = rowMeans(auto_res[, !neurons], na.rm = TRUE))
}, BPPARAM = bpparam)

## Make it long format
auto_long <- do.call(rbind, lapply(seq_len(length(auto)), function(i) {
    data.frame('acf_neuron' = auto[[i]][1:4],
        'acf_glia' = auto[[i]][5:8],
        'lag' = seq_len(4),
        'clus' = rep(i, 4), 
        'clus_n' = rep(NROW(meth_list[[i]]), 4),
        'clus_range' = rep(width(range(gr_list[[i]], ignore.strand = TRUE)), 4),
    stringsAsFactors = FALSE)
}))
rownames(auto_long) <- NULL
auto_long$context <- opt$context

## Global ACF summary by lag
print('ACF by lag, then abs(ACF) by lag -- neurons')
tapply(auto_long$acf_neuron, auto_long$lag, summary)
tapply(auto_long$acf_neuron, auto_long$lag, function(x) { summary(abs(x))})

print('ACF by lag, then abs(ACF) by lag -- glia')
tapply(auto_long$acf_glia, auto_long$lag, summary)
tapply(auto_long$acf_glia, auto_long$lag, function(x) { summary(abs(x))})

## Save results
dir.create('rda', showWarnings = FALSE)
save(auto_long, file = paste0('rda/auto_long_context_',
    opt$context, '.Rdata'))

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()

