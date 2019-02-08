## Required libraries
library('getopt')
library('bsseq')
library('sessioninfo')
library('GenomicRanges')
library('limma')

## Specify parameters
spec <- matrix(c(
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
    opt <- list(context = 'CHG', model = 'age')
    opt <- list(context = 'all', model = 'age')
    opt <- list(context = 'all', model = 'cell')
    opt <- list(context = 'all', model = 'interaction')
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


## Get the positions of the Cs
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

## Compute the mean lag distance for each group/cluster of Cs
lag_dist <- compute_lag_dist(start_list)

## Load the original files and append the info
dir.create('rda/original', showWarnings = FALSE)
r_file <-  paste0('rda/auto_long_context',
    opt$context, '_', opt$model, '.Rdata')
stopifnot(file.exists(r_file))
load(r_file, verbose = TRUE)

## Add the info
auto_long$avg_distance <- as.vector(t(lag_dist))

## Move the original file, then save the new one
system(paste0('mv ', r_file, ' rda/original/', basename(r_file)))
save(auto_long, file = r_file)

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
