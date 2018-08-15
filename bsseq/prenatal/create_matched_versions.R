library('bsseq')
library('GenomicRanges')
library('devtools')
library('getopt')

## Specify parameters
spec <- matrix(c(
    'chr', 'c', 1, 'character', 'chromosome to subset',
    'type', 't', 2, 'character', 'methylation type: either CpG or CpH',
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
    opt <- list(chr = 'chr21', type = 'CpG')
    opt <- list(chr = 'all', type = 'CpG')
}


if(opt$chr != 'all') {
    
    filt_cpg <- paste0('rda/', opt$chr,
        '_cleaned_CX_Prenatal_filtered_CpG.Rdata')
    filt_cph <- paste0('rda/', opt$chr, 
        '_cleaned_CX_Prenatal_filtered_CpH.Rdata')
    
    if(!all(file.exists(c(filt_cpg, filt_cph)))) {
        load('../lister/gr_all_highCov.Rdata', verbose = TRUE)
        load('../lister/gr_cps_minCov_3.Rdata', verbose = TRUE)
        gr_cpgs$c_context <- Rle('CG')
        gr_cpgs$trinucleotide_context <- Rle(NA)
        gr_all <- c(gr_all_highCov, gr_cpgs)
    
        f <- paste0('rda/', opt$chr, '_cleaned_CX_Prenatal.Rdata')
        message(paste(Sys.time(), 'loading', f))
        load(f, verbose = TRUE)
    
        ## Tweak reportFiles so they can be combined
        colData(BSobj)$reportFiles <- gsub('chr.*', '', colData(BSobj)$reportFiles)
        BSobj_all <- BSobj
    
        message(paste(Sys.time(), 'filtering to our CpG bases'))
        BSobj <- BSobj_all[subjectHits(findOverlaps(gr_cpgs, rowRanges(BSobj_all))), ]
        message(paste(Sys.time(), 'saving filtered version'))
        save(BSobj, file = filt_cpg)
    
        message(paste(Sys.time(), 'filtering to our CpH bases'))
        BSobj <- BSobj_all[subjectHits(findOverlaps(gr_all_highCov, rowRanges(BSobj_all))), ]
        message(paste(Sys.time(), 'saving filtered version'))
        save(BSobj, file = filt_cph)
    }
    
    message(paste(Sys.time(), 'loading', filt_cpg))
    load(filt_cpg, verbose = TRUE)
    message(paste(Sys.time(), 'coercing to matrix'))
    assays(BSobj)$M <- as.matrix(assays(BSobj)$M)
    assays(BSobj)$Cov <- as.matrix(assays(BSobj)$Cov)
    message(paste(Sys.time(), 'saving filtered & coerced version'))
    save(BSobj, file = gsub('.Rdata', '_coerced.Rdata', filt_cpg))
        
        
    message(paste(Sys.time(), 'loading', filt_cph))
    load(filt_cph, verbose = TRUE)
    message(paste(Sys.time(), 'coercing to matrix'))
    assays(BSobj)$M <- as.matrix(assays(BSobj)$M)
    assays(BSobj)$Cov <- as.matrix(assays(BSobj)$Cov)
    message(paste(Sys.time(), 'saving filtered & coerced version'))
    save(BSobj, file = gsub('.Rdata', '_coerced.Rdata', filt_cph))

    ## Reproducibility information
    print('Reproducibility information:')
    print(Sys.time())
    print(proc.time())
    options(width = 120)
    print(session_info())
    
    q('no')
}

stopifnot(opt$type %in% c('CpG', 'CpH'))


message(paste(Sys.time(), 'Combining results for the chrs'))
patt <- paste0('_cleaned_CX_Prenatal_filtered_', ifelse(opt$type == 'CpG', 'CpG', 'CpH'), '_coerced.Rdata')

files <- dir('rda', pattern = patt, full.names = TRUE)
message(paste(Sys.time(), 'files to combine:'))
print(files)



## Load and filter:
load_filt <- function(f) {
    message(paste(Sys.time(), 'loading', f))
    load(f, verbose = TRUE)
    
    ## Seems to be needed
    return(updateObject(BSobj))
}

f_filt <- paste0('allChrs_cleaned_CX_Prenatal_', opt$type, '.Rdata')
if(!file.exists(f_filt)) {
    BSobj <- do.call(rbind, lapply(files, load_filt))

    ## Save final combined result
    message(paste(Sys.time(), 'saving', f_filt))
    save(BSobj, file = f_filt)
} else {
    message(paste(Sys.time(), 'loading', f_filt))
    load(f_filt, verbose = TRUE)
}


if(opt$type == 'CpG') {
    ## nonCG next
    message(paste(Sys.time(), 'loading our GR'))
    load('../lister/gr_all_highCov.Rdata', verbose = TRUE)
    
    message(paste(Sys.time(), 'filtering nonCpGs'))
    
    
    Cov <- M <- matrix(0, nrow = length(gr_all_highCov), ncol = ncol(BSobj))
    colnames(Cov) <- colnames(M) <- colnames(BSobj)

    ov <- findOverlaps(gr_all_highCov, rowRanges(BSobj))
    M[queryHits(ov), ] <- assays(BSobj)$M[subjectHits(ov), ]
    Cov[queryHits(ov), ] <- assays(BSobj)$Cov[subjectHits(ov), ]

    BSobj <- BSseq(M = M, Cov = Cov, gr = gr_all_highCov, pData = colData(BSobj))
    message(paste(Sys.time(), 'saving allChrs_matched_cleaned_CX_Prenatal_nonCG_highCov.Rdata'))
    save(BSobj, file = 'allChrs_matched_cleaned_CX_Prenatal_nonCG_highCov.Rdata')
} else {
    ## CG next
    message(paste(Sys.time(), 'loading our GR'))
    load('../lister/gr_cps_minCov_3.Rdata', verbose = TRUE)
    
    message(paste(Sys.time(), 'filtering CpGs'))

    Cov <- M <- matrix(0, nrow = length(gr_cpgs), ncol = ncol(BSobj))
    colnames(Cov) <- colnames(M) <- colnames(BSobj)

    ov <- findOverlaps(gr_cpgs, rowRanges(BSobj))
    M[queryHits(ov), ] <- assays(BSobj)$M[subjectHits(ov), ]
    Cov[queryHits(ov), ] <- assays(BSobj)$Cov[subjectHits(ov), ]

    BSobj <- BSseq(M = M, Cov = Cov, gr = gr_cpgs, pData = colData(BSobj))
    message(paste(Sys.time(), 'saving BSobj_matched_cleaned_CX_Prenatal_minCov_3.Rdata'))
    save(BSobj, file = 'BSobj_matched_cleaned_CX_Prenatal_minCov_3.Rdata')
}

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
