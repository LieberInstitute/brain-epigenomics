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

stopifnot(opt$context %in% c('all', 'nonCG', 'CG', 'CHG', 'CHH'))

message(paste(Sys.time(), 'loading the data'))
load(paste0(opt$chr, '_postNatal_cleaned_CX_noHomogenate.Rdata'))

## Original bases in each context
gr <- rowRanges(BSobj)
context <- c(
    'trinucleotide_context' = sort(table(gr$trinucleotide_context)),
    'c_context' = sort(table(gr$c_context)))

## Filter low coverage bases
cov <- getCoverage(BSobj, type = 'Cov')
cov.ge1 <- cov >= 1
cov.filt <- rowSums(cov.ge1) == ncol(cov)
print("Number of bases with Cov >=1")
table(cov.filt)
#   FALSE    TRUE
# 6801196 7533737
round(table(cov.filt) / length(cov.filt) * 100, 2)
#FALSE  TRUE
#47.44 52.56
BSobj <- BSobj[cov.filt, ]
rm(cov, cov.ge1, cov.filt)


## Subset by context
if(opt$context != 'all') {
    if(opt$context == 'nonCG') {
        BSobj <- BSobj[rowRanges(BSobj)$c_context != 'CG', ]
    } else {
        BSobj <- BSobj[rowRanges(BSobj)$c_context == opt$context, ]
    }
}
print('Final number of bases of interest')
nrow(BSobj)

gr <- rowRanges(BSobj)

context_filt <- c(
    'trinucleotide_context' = sort(table(gr$trinucleotide_context)),
    'c_context' = sort(table(gr$c_context)))
context_filt <- context_filt[names(context)]
    
## Bases in each context
print('Bases in each context')
rbind(context, context_filt,
    'remaining percent' = round(context_filt / context * 100, 2))

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

if(FALSE) {
    ## Groups based on number of observations nearby
    chunksize <- 3000

    ## Determine number of loops
    lastloop <- trunc(nrow(BSobj) /  chunksize)
    ## Fix the lastloop in case that the N is a factor of
    ## chunksize
    if (nrow(BSobj) %% chunksize == 0 & lastloop > 0) {
        lastloop <- lastloop - 1
    }

    ## Build the index for splitting
    split.len <- rep(chunksize, lastloop)
    split.len.sum <- nrow(BSobj) - sum(split.len)
    if (split.len.sum > 0) {
        split.len <- c(split.len, split.len.sum)
    }
    split.idx <- Rle(seq_len(lastloop + 1), split.len)
}

## Size of the groups by genomic range
gr_list <- split(gr, split.idx)
chunk_width <- width(range(gr_list))
print('Size of groups by genomic range')
summary(unlist(chunk_width))

## Split the methylation data
meth <- DataFrame(getMeth(BSobj, type = 'raw'))
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


auto <- bplapply(meth_list, function(x) {
    ## Drop first row because it's always 1
    apply(x, 2, function(y) { acf(y, plot = FALSE, lag.max = 4)$acf })[-1, ]
}, BPPARAM = bpparam)


## Make it long format
auto_long <- do.call(rbind, lapply(seq_len(length(auto)), function(i) {
    data.frame(acf = as.vector(auto[[i]]),
    'lag' = rep(seq_len(4), ncol(auto[[i]])),
    'clus' = rep(i, 4 * ncol(auto[[i]])),
    'clus_n' = rep(nrow(meth_list[[i]]), 4 * ncol(auto[[i]])),
    'clus_range' = rep(width(range(gr_list[[i]])), 4 * ncol(auto[[i]])),
    'sample' = rep(colnames(auto[[i]]), each = 4))
}))
auto_long$context <- opt$context

## Global ACF summary by lag
print('ACF by lag, then abs(ACF) by lag')
tapply(auto_long$acf, auto_long$lag, summary)
tapply(auto_long$acf, auto_long$lag, function(x) { summary(abs(x))})

## Save results
save(auto_long, file = paste0('auto_long_', opt$chr, '_context', opt$context, '.Rdata'))


## extract pheno
pd <- pData(BSobj)

## Explore with limma
models <- list(
    'cell' = with(pd, model.matrix(~ Cell.Type + Age)),
    'age' = with(pd, model.matrix(~ Age + Cell.Type)),
    'interaction' = with(pd, model.matrix(~ Age * Cell.Type)))

fits <- lapply(models, function(mod) {
    lmFit(meth, design = mod)
})
coefs <- c(2, 2, 4)
names(coefs) <- names(fits)


coef_interest <- mapply(function(f, coef) {
    f$coefficients[, coef]
}, fits, coefs)
summary(coef_interest)
summary(abs(coef_interest))

save(fits, models, coef_interest, file = paste0('limma_exploration_', opt$chr, '_context', opt$context, '.Rdata'))



## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
