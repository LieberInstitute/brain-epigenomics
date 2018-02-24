library('SGSeq')
library('bsseq')
library('getopt')
library('devtools')
library('limma')
library('gplots')
library('VennDiagram')
library('RColorBrewer')
library('clusterProfiler')
library('org.Hs.eg.db')
library('ggplot2')

opt <- list('feature' = 'gene')

stopifnot(opt$feature %in% c('psi', 'gene', 'exon', 'jx'))
cpgs <- c('CpG', 'nonCpG', 'CpGmarg')

## For re-loading
saved_files <- c(paste0('rda/meqtl_mres_', opt$feature,
    '_using_near_meth11_proteincoding.Rdata'),
    paste0('rda/meqtl_venn_',opt$feature, '_using_near.Rdata'),
    paste0('rda/meqtl_summary_', opt$feature, '_using_near.Rdata'),
    paste0('rda/meqtl_delta_pval_', opt$feature, '_using_near.Rdata'),
    paste0('rda/meqtl_venn_go_', opt$feature, '_using_near.Rdata'),
    paste0('rda/meqtl_age_coef_', opt$feature, '_using_near.Rdata')
)
for(i in saved_files) load(i, verbose = TRUE)
rm(i)

get_age <- function(type) {
    colData(mres[[type]]$meth)$Age
}
get_meth <- function(type, i) {
    as.vector(getMeth(mres[[type]]$meth[i, ], type = 'raw'))
}

get_y <- function(type, i) {
    if(opt$feature == 'psi') {
        res <- variantFreq(mres[[type]]$expr[i, ])
    } else  {
        res <- log2(assays(mres[[type]]$expr)$norm[i, ] + 1)
    }
    return(res)
}

## Check if methylation still explains the expression
## even after adjusting by age
agemeth_coef <- lapply(names(mres), function(type) {
    message(paste(Sys.time(), 'processing', type))
    age <- get_age(type)

    res <- as.data.frame(t(sapply(seq_len(nrow(mres[[type]]$meth)), function(i) {
        fit <- lm(as.vector(get_y(type, i)) ~ get_meth(type, i) + age)
        summary(fit)$coef[2, ]
    })))
    res$ageFDR <- p.adjust(res[, 'Pr(>|t|)'], 'fdr')
    return(res)
})
names(agemeth_coef) <- names(mres)
save(agemeth_coef, file = paste0('rda/meqtl_agemeth_coef_', opt$feature, '_using_near.Rdata'))
