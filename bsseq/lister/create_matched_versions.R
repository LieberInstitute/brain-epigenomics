library('bsseq')
library('GenomicRanges')
library('devtools')

## CpGs first
load('gr_cps_minCov_3.Rdata')
load('BSobj_lister_minCov_3.Rdata')

gr <- rowRanges(BSobj)

Cov <- M <- matrix(0, nrow = length(gr_cpgs), ncol = ncol(BSobj))
colnames(Cov) <- colnames(M) <- colnames(BSobj)

ov <- findOverlaps(gr_cpgs, rowRanges(BSobj))
M[queryHits(ov), ] <- assays(BSobj)$M[subjectHits(ov), ]
Cov[queryHits(ov), ] <- assays(BSobj)$Cov[subjectHits(ov), ]

BSobj <- BSseq(M = M, Cov = Cov, gr = gr_cpgs)
save(BSobj, file = 'BSobj_matched_lister_minCov_3.Rdata')

rm(BSobj, gr, Cov, M, ov, gr_cpgs)


## nonCG next
load('gr_all_highCov.Rdata')
load('allChrs_lister_nonCG_highCov.Rdata')

gr <- rowRanges(BSobj)

Cov <- M <- matrix(0, nrow = length(gr_all_highCov), ncol = ncol(BSobj))
colnames(Cov) <- colnames(M) <- colnames(BSobj)

ov <- findOverlaps(gr_cpgs, rowRanges(BSobj))
M[queryHits(ov), ] <- assays(BSobj)$M[subjectHits(ov), ]
Cov[queryHits(ov), ] <- assays(BSobj)$Cov[subjectHits(ov), ]

BSobj <- BSseq(M = M, Cov = Cov, gr = gr_all_highCov)
save(BSobj, file = 'allChrs_matched_lister_nonCG_highCov.Rdata')

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
