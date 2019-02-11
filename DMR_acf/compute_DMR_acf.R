## Required libraries
library('getopt')
library('bsseq')
library('devtools')
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
    
    print('Final number of bases of interest -- after merging')
    print(nrow(DMR))
    rm(cg, ncg, gr, ov)
}

## Sort by chr position
message(paste(Sys.time(), 'sorting by chr position'))
DMR <- sort(DMR, ignore.strand = TRUE)


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

## Split the methylation data
meth <- DataFrame(getMeth(DMR, type = 'raw'))
meth_list <- split(meth, split.idx)

print('Keeping only DMRs with at least 5 observations by group')
n_group <- sapply(split(split.idx, split.idx), length)
table(n_group >= 5)
meth_list <- meth_list[names(n_group)[n_group >= 5]]
gr_list <- gr_list[names(gr_list)[n_group >= 5]]

print('Summary of number of Cs per DMR')
summary(sapply(meth_list, nrow))

if(opt$context == 'all') {
    byctxt <- do.call(rbind, lapply(gr_list, function(x) { table(x$c_context)}))
    print(nrow(byctxt))
    print(summary(byctxt))
    dir.create('rda', showWarnings = FALSE)
    save(byctxt, file = paste0('rda/byctxt_context', opt$context, '_',
        opt$model, '.Rdata'))
        
        
    panel.hist <- function(x, ...)
    {
        usr <- par("usr"); on.exit(par(usr))
        par(usr = c(usr[1:2], 0, 1.5) )
        h <- hist(x, plot = FALSE, breaks = 50)
        breaks <- h$breaks; nB <- length(breaks)
        y <- h$counts; y <- y/max(y)
        rect(breaks[-nB], 0, breaks[-1], y, col = "light blue", ...)
    }

    dir.create('pdf', showWarnings = FALSE)
    pdf(paste0('pdf/byctxt_context', opt$context, '_', opt$model, '.pdf'),
        width = 10, height = 10)
    pairs(byctxt, panel = panel.smooth, diag.panel = panel.hist, cex.labels = 2, font.labels = 2)
    dev.off()
    
}



if(opt$cores == 1) {
    bpparam <- SerialParam()
} else {
    bpparam <- SnowParam(opt$cores, outfile =  Sys.getenv('SGE_STDERR_PATH'))
}

neurons <- colData(DMR)$Cell.Type == 'Neuron'

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

