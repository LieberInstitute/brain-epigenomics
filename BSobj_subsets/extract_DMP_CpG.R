library('bsseq')
library('devtools')
library('limma')

## Load the raw data
load('/dcl01/lieber/ajaffe/lab/brain-epigenomics/bumphunting/BSobj_bsseqSmooth_Neuron_minCov_3.Rdata', verbose = TRUE)

## Load limma results
load('/dcl01/lieber/ajaffe/lab/brain-epigenomics/bumphunting/limma_Neuron_CpGs_minCov_3.Rdata', verbose = TRUE)

## Delete unused parts
rm(models, coef_interest, t_interest)

## Coefficients of interest
coefs <- c(2, 2, 4)
names(coefs) <- names(fits)

## Find the DMPs
DMP_index <- mapply(function(f, coef) {
    idx <- which(p.adjust(f$p.value[, coef], method = 'fdr') < 0.05)
}, fits, coefs, SIMPLIFY = FALSE)

print('Number of DMPs per model (then in millions)')
sapply(DMP_index, length)
sapply(DMP_index, length) / 1e6

## Subset the BSobj according to each model
xx <- mapply(function(idx, model) {
    DMP_CpG <- BSobj[idx, ]
    dir.create('rda', showWarnings = FALSE)
    save(DMP_CpG, file = paste0('rda/DMP_CpG_', model, '.Rdata'))
    return(model)
}, DMP_index, names(DMP_index))

## Reproducibility info
proc.time()
message(Sys.time())
options(width = 120)
session_info()
