# qrsh -l bluejay,mem_free=150G,h_vmem=150G,h_fsize=200G

library('bsseq')
library('GenomicRanges')
library('devtools')
library('getopt')


#### To run interactively 

setwd("/dcl01/lieber/ajaffe/lab/brain-epigenomics/public_WGBS/code")
chrs = c(1:19, "X","Y","M")
types = c('CpG', 'CpH')

## Separate CG from CH and coerce to matrix

for (i in 1:length(chrs)) {
  opt <- list(chr = paste0('chr', chrs[i]))
  
  f <- paste0('rda/', opt$chr, '_cleaned_CX_Stroud.Rdata')
  message(paste(Sys.time(), 'loading', f))
  load(f, verbose = TRUE)
  
  BSobj <- updateObject(BSobj)
  colData(BSobj)$reportFiles <- gsub('chr.*', '', colData(BSobj)$reportFiles)
  
  message(paste(Sys.time(), 'coercing to matrix: chr', chrs[i]))
  assays(BSobj)$M <- as.matrix(assays(BSobj)$M)
  assays(BSobj)$Cov <- as.matrix(assays(BSobj)$Cov)
  
  x = BSobj
  
  message(paste(Sys.time(), 'saving coerced CGs: chr', chrs[i]))
  BSobj = x[granges(x)$c_context=="CG",]
  save(BSobj, file = gsub('_cleaned_CX_Stroud.Rdata', '_cleaned_Stroud_CpG_coerced.Rdata', f))
  message(paste(Sys.time(), 'saving coerced CHs: chr', chrs[i]))
  BSobj = x[granges(x)$c_context!="CG",]
  save(BSobj, file = gsub('_cleaned_CX_Stroud.Rdata', '_cleaned_Stroud_CpH_coerced.Rdata', f))
  
  rm(BSobj,x)
}


## Combine the BSobjs

# Load function
load_filt <- function(f) {
  message(paste(Sys.time(), 'loading', f))
  load(f, verbose = TRUE)
}

for (i in 1:length(types)) {
  opt <- list(type = types[i])
  
  message(paste(Sys.time(), 'Combining results for the chrs'))
  patt <- paste0('_cleaned_Stroud_', ifelse(opt$type == 'CpG', 'CpG', 'CpH'), '_coerced.Rdata')
  
  files <- dir('rda', pattern = patt, full.names = TRUE)
  message(paste(Sys.time(), 'files to combine:'))
  print(files)
  
  f_filt <- paste0('rda/allChrs_cleaned_Stroud_', opt$type, '.Rdata')
  
  if(!file.exists(f_filt)) {
    BSobj <- do.call(rbind, lapply(files, load_filt))
    
    ## Save final combined result
    message(paste(Sys.time(), 'saving', f_filt))
    save(BSobj, file = f_filt)
  } else {
    message(paste(Sys.time(), 'loading', f_filt))
    load(f_filt, verbose = TRUE)
  }
}
