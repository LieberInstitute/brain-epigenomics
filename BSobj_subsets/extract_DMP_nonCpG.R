library('bsseq')
library('devtools')
library('limma')

## Load the raw data
load('/dcl01/lieber/ajaffe/lab/brain-epigenomics/bsseq/bsobj_by_chr/allChrs_postNatal_cleaned_nonCG_noHomogenate_highCov.Rdata', verbose = TRUE)

## Load limma results
load('/dcl01/lieber/ajaffe/lab/brain-epigenomics/bsseq/bsobj_by_chr/limma_exploration_nonCG_highCov.Rdata', verbose = TRUE)

## Delete unused parts
rm(models, coef_interest, fits)

## Coefficients of interest
coefs <- c(2, 2, 4)
names(coefs) <- names(ebList)

## Find the DMPs
DMP_index <- mapply(function(f, coef) {
    idx <- which(p.adjust(f$p.value[, coef], method = 'fdr') < 0.05)
}, ebList, coefs, SIMPLIFY = FALSE)

print('Number of DMPs per model (then in millions)')
sapply(DMP_index, length)
sapply(DMP_index, length) / 1e6

## Subset the BSobj according to each model
xx <- mapply(function(idx, model) {
    DMP_nonCpG <- BSobj[idx, ]
    dir.create('rda', showWarnings = FALSE)
    save(DMP_nonCpG, file = paste0('rda/DMP_nonCpG_', model, '.Rdata'))
    return(model)
}, DMP_index, names(DMP_index))

## Reproducibility info
proc.time()
message(Sys.time())
options(width = 120)
session_info()
