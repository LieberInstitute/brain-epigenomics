library('bsseq')
library('devtools')
library('limma')
library('jaffelab')
library('shinycsv') # must download
library('RColorBrewer')

## Load data
load('/dcl01/lieber/ajaffe/lab/brain-epigenomics/bsseq/bsobj_by_chr/allChrs_postNatal_cleaned_nonCG_noHomogenate_highCov.Rdata')

## Size of the data
dim(BSobj)

## extract pheno
pd <- pData(BSobj)

## Fix pheno data
pd$Race[pd$Race == "CAUC "] <- 'CAUC'
pd$Sex[pd$Sex == " M"] <- 'M'
pd$RIN <- as.numeric(gsub(" ", "", pd$RIN))
pd$pg.DNA.nuclei.input <- as.numeric(pd$pg.DNA.nuclei.input)
pd$Reads <- as.numeric(pd$Reads)
pd$Percent.GreaterThan.Q30 <- as.numeric(pd$Percent.GreaterThan.Q30)

pd = pd[which(pd$Cell.Type=="Neuron"),]

## Get methylation info
meth <- getMeth(BSobj, type = 'raw')
meth = meth[,which(colnames(meth) %in% rownames(pd))]
rm(BSobj)
meth.g0 <- meth > 0
meth.filt <- rowSums(meth.g0) >= 5
meth.tab <- table(meth.filt)
meth.tab
round(meth.tab / sum(meth.tab) * 100, 2)
meth <- meth[meth.filt, ]
rm(meth.g0, meth.filt)

## PCA
pcs <- prcomp(t(meth))
pcaVars <- getPcaVars(pcs)
names(pcaVars) <- paste0('PC', seq_len(length(pcaVars)))
save(pcs, pcaVars, pd, meth, file="/dcl01/lieber/ajaffe/lab/brain-epigenomics/non-CpG/rdas/nonCG_highCov_neuronsOnly_pca_pd_methTable.Rdata")

## Explore with limma
model <- with(pd, model.matrix(~ Age))

fit <- lmFit(meth, design = model)
pvalue <- fit$p.value[, 2]
qvalue <- qvalue(pvalue)$qvalues

coef_interest <- fit$coefficients[,2]
summary(coef_interest)
summary(abs(coef_interest))

ebResults <- ebayes(fit)
pvalue <- ebResults$p.value[, 2]
padj = p.adjust(pvalue, method = "fdr")
save(fit, model, coef_interest, ebResults, file = '/dcl01/lieber/ajaffe/lab/brain-epigenomics/non-CpG/rdas/limma_exploration_nonCG_highCov_neuronsOnly.Rdata')

