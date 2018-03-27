library("limma")
library("edgeR")


load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/brainseq_pipeline/polyA_unstranded/rse_gene_polyA_dlpfc_n41.Rdata", verbose = T)


# format phenotype table

hompd = colData(rse_gene)

# format counts

geneCounts = assay(rse_gene)
geneCounts = geneCounts[rowSums(geneCounts)>0,]
identical(colnames(geneCounts), rownames(hompd))

## Differential expression: gene level

design <- model.matrix(~ hompd$Age)
dge <- DGEList(counts = geneCounts)
dge <- calcNormFactors(dge)
dat <- voom(dge, design, plot = TRUE)
fit_gene <- lmFit(dat, design)
fit_gene_all <- eBayes(fit_gene)


## now only on postnatal samples

postpd = hompd[hompd$Age>0,]
geneCounts = geneCounts[,colnames(geneCounts) %in% rownames(postpd)]
geneCounts = geneCounts[rowSums(geneCounts)>0,]
identical(colnames(geneCounts), rownames(postpd))

## Differential expression: gene level

design <- model.matrix(~ postpd$Age)
dge <- DGEList(counts = geneCounts)
dge <- calcNormFactors(dge)
dat <- voom(dge, design, plot = TRUE)
fit_gene <- lmFit(dat, design)
fit_gene_post <- eBayes(fit_gene)


## Reformat combined results

homRNAres = as.data.frame(cbind(fit_gene_all$coefficients, fit_gene_all$t, fit_gene_all$p.value))
colnames(homRNAres) = c("Coeff.Intercept", "Coeff", "Tstat.Intercept", "Tstat", "pval.Intercept", "pval")
homRNAres$FDR = p.adjust(homRNAres$pval, method = "fdr")
homRNAres$gencodeID = rownames(homRNAres)

postRNAres = as.data.frame(cbind(fit_gene_post$coefficients, fit_gene_post$t, fit_gene_post$p.value))
colnames(postRNAres) = c("Coeff.Intercept", "Coeff", "Tstat.Intercept", "Tstat", "pval.Intercept", "pval")
postRNAres$FDR = p.adjust(postRNAres$pval, method = "fdr")
postRNAres$gencodeID = rownames(postRNAres)


save(fit_gene_post, fit_gene_all, postRNAres,homRNAres,hompd,postpd, file="/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/homogenate_RNA/DE_limma_results_homogenateRNAseq.rda")




