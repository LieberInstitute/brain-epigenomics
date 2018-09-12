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
    opt <- list(chr = 'chr19', type = 'CpG')
    opt <- list(chr = 'all', type = 'CpG')
}

if(opt$chr != 'all') {
  
  f <- paste0('rda/', opt$chr, '_cleaned_CX_Gabel.Rdata')
  message(paste(Sys.time(), 'loading', f))
  load(f, verbose = TRUE)
  
  BSobj <- updateObject(BSobj, verbose=TRUE)
  colData(BSobj)$reportFiles <- gsub('chr.*', '', colData(BSobj)$reportFiles)
  
  message(paste(Sys.time(), 'coercing to matrix'))
  assays(BSobj)$M <- as.matrix(assays(BSobj)$M)
  assays(BSobj)$Cov <- as.matrix(assays(BSobj)$Cov)
  
  message(paste(Sys.time(), 'Separating CG from CH'))
  BSobj_CG <- BSobj[granges(BSobj)$c_context=="CG",]
  BSobj_CH <- BSobj[granges(BSobj)$c_context!="CG",]
  
  message(paste(Sys.time(), 'saving coerced version'))
  save(BSobj_CG, file = gsub('_cleaned_CX_Gabel.Rdata', '_cleaned_Gabel_CpG_coerced.Rdata', f))
  save(BSobj_CH, file = gsub('_cleaned_CX_Gabel.Rdata', '_cleaned_Gabel_CpH_coerced.Rdata', f))

    
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
patt <- paste0('_cleaned_Gabel_', ifelse(opt$type == 'CpG', 'CpG', 'CpH'), '_coerced.Rdata')

files <- dir('rda', pattern = patt, full.names = TRUE)
message(paste(Sys.time(), 'files to combine:'))
print(files)



## Load and filter:
load_filt <- function(f) {
    message(paste(Sys.time(), 'loading', f))
    load(f, verbose = TRUE)
    if (opt$type == 'CpG') {
      BSobj <- BSobj_CG
      rm(BSobj_CG) } else {
        BSobj <- BSobj_CH
        rm(BSobj_CH)
      }

    ## Seems to be needed
    return(updateObject(BSobj))
}

f_filt <- paste0('rda/allChrs_cleaned_Gabel_', opt$type, '.Rdata')
if(!file.exists(f_filt)) {
    BSobj <- do.call(rbind, lapply(files, load_filt))

    ## Save final combined result
    message(paste(Sys.time(), 'saving', f_filt))
    save(BSobj, file = f_filt)
} else {
    message(paste(Sys.time(), 'loading', f_filt))
    load(f_filt, verbose = TRUE)
}


## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
