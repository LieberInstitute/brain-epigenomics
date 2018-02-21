library('bsseq')
library('devtools')

files <- dir(pattern = '_postNatal_cleaned_CX_noHomogenate.Rdata')
chrs <- gsub('_.*', '', files)

## Load and filter:
load_filt <- function(f) {
    chr <- gsub('_.*', '', f)
    res_file <- paste0(chr,
        '_postNatal_cleaned_nonCG_noHomogenate_highCov.Rdata')
    
    if(!file.exists(res_file)) {
        message(paste(Sys.time(), 'loading', f))
        load(f)
    
        message(paste(Sys.time(), 'filtering CGs', chr))
    
        ## Drop CG context
        BSobj <- BSobj[rowRanges(BSobj)$c_context != 'CG', ]
        original <- nrow(BSobj)
    
        ## Filter low coverage bases
        message(paste(Sys.time(), 'filtering low coverage', chr))
        cov <- getCoverage(BSobj, type = 'Cov')
        cov.ge5 <- cov >= 5
        cov.filt <- rowSums(cov.ge5) == ncol(cov)
    
        BSobj <- BSobj[cov.filt, ]
    
        filtered <- nrow(BSobj)
        message(paste(Sys.time(), chr, ':', original - filtered, 'bases did not have all samples with coverage >=5 out of', original, 'That is,', round(filtered / original * 100, 2), 'percent is remaining.'))
    
        ## Save results
        message(paste(Sys.time(), 'saving the results in', res_file))
        save(BSobj, file = res_file)
    } else {
        message(paste(Sys.time(), 'loading pre-existing', res_file))
        load(res_file)
    }
    
    ## Tweak reportFiles so they can be combined
    colData(BSobj)$reportFiles <- gsub('chr.*', '', colData(BSobj)$reportFiles)
    
    return(BSobj)
}

## Is this an array job?
array <- Sys.getenv('SGE_TASK_ID')
print(array)

if(array == 'undefined') {
    BSobj <- do.call(rbind, lapply(files, load_filt))

    ## Save final combined result
    save(BSobj, file = 'allChrs_postNatal_cleaned_nonCG_noHomogenate_highCov.Rdata')
} else {
    BSobj <- load_filt(files[as.integer(array)])
}

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
