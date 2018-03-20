library('bsseq')
library('GenomicRanges')
library('devtools')

message(paste(Sys.time(), 'Combining results for the chrs'))
files <- dir('rda', pattern = '_cleaned_CX_Homogenate.Rdata', full.names = TRUE)

load('../lister/gr_all_highCov.Rdata', verbose = TRUE)
load('../lister/gr_cps_minCov_3.Rdata', verbose = TRUE)
gr_cpgs$c_context <- Rle('CG')
gr_cpgs$trinucleotide_context <- Rle(NA)
gr_all <- c(gr_all_highCov, gr_cpgs)

## Load and filter:
load_filt <- function(f) {
    f_filtered <- gsub('.Rdata', '_filtered.Rdata', f)
    if(file.exists(f_filtered)) {
        message(paste(Sys.time(), 'loading', f_filtered))
        load(f, verbose = TRUE)
        return(BSobj)
    }
    
    message(paste(Sys.time(), 'loading', f))
    load(f, verbose = TRUE)
    ## Tweak reportFiles so they can be combined
    colData(BSobj)$reportFiles <- gsub('chr.*', '', colData(BSobj)$reportFiles)
    
    message(paste(Sys.time(), 'filtering to our bases'))
    BSobj <- BSobj[subjectHits(findOverlaps(gr_all, rowRanges(BSobj))), ]
    
    message(paste(Sys.time(), 'saving filtered version'))
    save(BSobj, file = f_filtered)
    
    return(BSobj)
}

BSobj <- do.call(rbind, lapply(files, load_filt))
rm(gr_all)

## Save final combined result
message(paste(Sys.time(), 'saving allChrs_cleaned_CX_Homogenate.Rdata'))
save(BSobj, file = 'allChrs_cleaned_CX_Homogenate.Rdata')


## nonCG next
message(paste(Sys.time(), 'filtering nonCpGs'))

gr <- rowRanges(BSobj)

Cov <- M <- matrix(0, nrow = length(gr_all_highCov), ncol = ncol(BSobj))
colnames(Cov) <- colnames(M) <- colnames(BSobj)

ov <- findOverlaps(gr_all_highCov, rowRanges(BSobj))
M[queryHits(ov), ] <- assays(BSobj)$M[subjectHits(ov), ]
Cov[queryHits(ov), ] <- assays(BSobj)$Cov[subjectHits(ov), ]

BSobj <- BSseq(M = M, Cov = Cov, gr = gr_all_highCov, pData = colData(BSobj))
message(paste(Sys.time(), 'saving allChrs_matched_cleaned_CX_Homogenate_nonCG_highCov.Rdata'))
save(BSobj, file = 'allChrs_matched_cleaned_CX_Homogenate_nonCG_highCov.Rdata')

rm(BSobj, gr, Cov, M, ov, gr_all_highCov)

## CG next
message(paste(Sys.time(), 'filtering CpGs'))
load('allChrs_cleaned_CX_Homogenate.Rdata', verbose = TRUE)
Cov <- M <- matrix(0, nrow = length(gr_cpgs), ncol = ncol(BSobj))
colnames(Cov) <- colnames(M) <- colnames(BSobj)

ov <- findOverlaps(gr_cpgs, rowRanges(BSobj))
M[queryHits(ov), ] <- assays(BSobj)$M[subjectHits(ov), ]
Cov[queryHits(ov), ] <- assays(BSobj)$Cov[subjectHits(ov), ]

BSobj <- BSseq(M = M, Cov = Cov, gr = gr_cpgs, pData = colData(BSobj))
message(paste(Sys.time(), 'saving BSobj_matched_cleaned_CX_Homogenate_minCov_3.Rdata'))
save(BSobj, file = 'BSobj_matched_cleaned_CX_Homogenate_minCov_3.Rdata')

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
