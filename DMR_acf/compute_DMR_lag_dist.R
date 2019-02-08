## Required libraries
library('getopt')
library('bsseq')
library('sessioninfo')
library('GenomicRanges')
library('BiocParallel')
library('limma')

## Specify parameters
spec <- matrix(c(
    'cores', 't', 1, 'integer', 'Number of cores to use',
    'model', 'm', 1, 'character', 'Either age, cell or interaction',
    'context', 'x', 1, 'character', 'Either nonCG, CG, CHG, CHH',
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
    opt <- list(cores = 1, context = 'CHG', model = 'age')
    opt <- list(cores = 1, context = 'all', model = 'age')
    opt <- list(cores = 1, context = 'all', model = 'cell')
    opt <- list(cores = 1, context = 'all', model = 'interaction')
}

stopifnot(opt$context %in% c('nonCG', 'CG', 'CHG', 'CHH', 'all'))
stopifnot(opt$model %in% c('age', 'cell', 'interaction'))

## Load the subset BSobj
load_DMR <- function(model) {
    if(opt$context == 'CG') {
        load(paste0('/dcl01/lieber/ajaffe/lab/brain-epigenomics/BSobj_subsets/rda/DMR_CpG_', model, '_CpGbased.Rdata'), verbose = TRUE)
        rowRanges(DMR_CpG)$c_context <- 'CG'
        DMR <- DMR_CpG
    } else {
        load(paste0('/dcl01/lieber/ajaffe/lab/brain-epigenomics/BSobj_subsets/rda/DMR_nonCpG_', model, '_CpGbased.Rdata'), verbose = TRUE)
        DMR <- DMR_nonCpG
        if(opt$context != 'nonCG') {
            DMR <- DMR[rowRanges(DMR)$c_context == opt$context, ]
        }
    }
    print('Final number of bases of interest')
    print(nrow(DMR))
    return(DMR)
}
if(opt$context != 'all') {
    DMR <- load_DMR(opt$model)
} else {
    opt$context <- 'CG'
    cg <- load_DMR(opt$model)
    rowRanges(cg)$trinucleotide_context <- Rle(NA)
    opt$context <- 'nonCG'
    ncg <- load_DMR(opt$model)
    colData(cg)$reportFiles <- colData(ncg)$reportFiles
    DMR <- BSseq(
        M = rbind(assays(cg)$M, assays(ncg)$M),
        Cov = rbind(assays(cg)$Cov, assays(ncg)$Cov),
        gr = c(rowRanges(cg), rowRanges(ncg)),
        pData = colData(ncg),
        chr = seqlevels(rowRanges(cg))
    )
    opt$context <- 'all'
    
    ## Andd context info to GR
    gr <- c(rowRanges(cg), rowRanges(ncg))
    ov <- findOverlaps(DMR, gr, type = 'equal')
    
    rowRanges(DMR)$c_context <- gr$c_context[subjectHits(ov)]
    rowRanges(DMR)$trinucleotide_context <- gr$trinucleotide_context[subjectHits(ov)]
    
    ## Sort by chr position
    DMR <- sort(DMR, ignore.strand = TRUE)
    
    print('Final number of bases of interest -- after merging')
    print(nrow(DMR))
    rm(cg, ncg, gr, ov)
}




## Load the bumphunting results
load(paste0('/dcl01/lieber/ajaffe/lab/brain-epigenomics/bumphunting/', 
    'bumps_bsseqSmooth_Neuron_', opt$model, '_250_perm.Rdata'), verbose = TRUE) 

## Build a GRanges object from the table output of bumphunter
bumps_gr <- GRanges(bumps$table[bumps$table$fwer < 0.05, ])

## Group bases by DMR
split.idx <- subjectHits(findOverlaps(DMR, bumps_gr))
gr <- rowRanges(DMR)
gr_list <- split(gr, split.idx)

print(paste(length(unique(split.idx)), 'out of', length(bumps_gr),
    'DMRs have at least one base in the DMR'))

print('Keeping only DMRs with at least 5 observations by group')
n_group <- sapply(split(split.idx, split.idx), length)
table(n_group >= 5)
gr_list <- gr_list[names(gr_list)[n_group >= 5]]


if(opt$cores == 1) {
    bpparam <- SerialParam()
} else {
    bpparam <- MulticoreParam(opt$cores, outfile =  Sys.getenv('SGE_STDERR_PATH'))
}

distance_lag <- function(lag = 1, g) {
    n <- length(g)
    i <- seq_len(n - lag)
    distance(g[i], g[i + lag], ignore.strand = TRUE) + 1
}
compute_d_lags <- function(g, max.lag = 4) {
    res <- lapply(seq_len(max.lag), distance_lag, g = g)
    names(res) <- seq_len(max.lag)
    return(res)
}

## Code that helped me realized I needed ignore.strand = TRUE
# head(which(is.na(compute_d_lags(gg[[1]], 1)[[1]])))
#
# gg[[1]][138:139]

compute_mean_d_lags <- function(g, max.lag = 4) {
    sapply(compute_d_lags(g, max.lag), mean)
}

## Testing code
# gg <- gr_list[1:10]
# res <- do.call(rbind, lapply(gg, compute_mean_d_lags))
# boxplot(res)

start_list <- start(gr_list)


compute_lag_dist <- function(ss, max.lag = 4) {
    i_list <- lapply(seq_len(max.lag), function(lag) {
        IntegerList(lapply(elementNROWS(ss), function(x) {
            seq_len(x - lag)
        }))
    })

    result <- mapply(function(i, lag) {
        mean(ss[i + lag] - ss[i])
    }, i_list, seq_len(max.lag))
    colnames(result) <- seq_len(max.lag)
    return(result)
}


res2 <- compute_lag_dist(start_list[1:10])
stopifnot(identical(res2, res))



lag_dist <- do.call(rbind, bplapply(gr_list, compute_mean_d_lags, BPPARAM = bpparam))

Sys.time()
lag_dist2 <- compute_lag_dist(start_list)
Sys.time()

stopifnot(identical(lag_dist, lag_dist2))

compute_mean_d_lags(g)

lapply(gg, compute_mean_d_lags)

sapply(compute_d_lags(g), mean)

auto <- bplapply(meth_list, function(x) {
    ## Drop first row because it's always 1
    auto_res <- apply(as.matrix(x), 2, function(y) { acf(y, plot = FALSE, lag.max = 4)$acf })[-1, ]
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
auto_long$model <- opt$model

## Global ACF summary by lag
print('ACF by lag, then abs(ACF) by lag -- neurons')
tapply(auto_long$acf_neuron, auto_long$lag, summary)
tapply(auto_long$acf_neuron, auto_long$lag, function(x) { summary(abs(x))})

print('ACF by lag, then abs(ACF) by lag -- glia')
tapply(auto_long$acf_glia, auto_long$lag, summary)
tapply(auto_long$acf_glia, auto_long$lag, function(x) { summary(abs(x))})

## Save results
dir.create('rda', showWarnings = FALSE)
save(auto_long, file = paste0('rda/auto_long_context',
    opt$context, '_', opt$model, '.Rdata'))

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()

