# qrsh -l bluejay,mem_free=150G,h_vmem=150G,h_fsize=100G
# cd /dcl01/lieber/ajaffe/lab/brain-epigenomics/bumphunting
# module load conda_R/3.4.x
library('bsseq')
library('devtools')
library('limma')

## Load the raw data
load('BSobj_bsseqSmooth_Neuron_minCov_3.Rdata')

## extract pheno
pd <- pData(BSobj)

pd$Cell.Type <- factor(pd$Cell.Type)
levels(pd$Cell.Type)

## Get methylation info
meth <- getMeth(BSobj, type = 'raw')
dim(meth)

## Explore with limma
models <- list(
    'cell' = with(pd, model.matrix(~ Cell.Type + Age)),
    'age' = with(pd, model.matrix(~ Age + Cell.Type)),
    'interaction' = with(pd, model.matrix(~ Age * Cell.Type)))

fits <- lapply(models, function(mod) {
    f <- lmFit(meth, design = mod)
    eBayes(f)
})
coefs <- c(2, 2, 4)
names(coefs) <- names(fits)


coef_interest <- mapply(function(f, coef) {
    f$coefficients[, coef]
}, fits, coefs)
summary(coef_interest)
summary(abs(coef_interest))

t_interest <- mapply(function(f, coef) {
    f$t[, coef]
}, fits, coefs)

head(t_interest)

save(fits, models, coef_interest, t_interest, file = 'limma_Neuron_CpGs_minCov_3.Rdata')

## Reproducibility info
proc.time()
message(Sys.time())
options(width = 120)
session_info()

