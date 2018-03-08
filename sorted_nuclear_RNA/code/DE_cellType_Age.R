#results <- do.call(rbind, lapply(regs, function(reg) {
#  message(paste(Sys.time(), 'processing region', reg))
#  design <- with(subset(pd, Region == reg), model.matrix(~ Dx + RIN + age + Sex + Race + mitoRate + gene_Assigned_Percent))
#  
#  counts <- geneCounts[, pd$Region == reg]
#  dge <- DGEList(counts = counts)
#  dge <- calcNormFactors(dge)
#  v <- voom(dge, design, plot = TRUE)
#  fit <- lmFit(v, design)
#  fit <- eBayes(fit)
#  log2FC <- fit$coefficients[, 2]
#  pvalue <- fit$p.value[, 2]
#  qvalue <- qvalue(pvalue)$qvalues
#  gene_res <- cbind(data.frame(log2FC, pvalue, qvalue, Region = rep(reg, nrow(counts))), geneMap)
#  rownames(gene_res) <- NULL
#  return(gene_res)
#}))


library("limma")
library("edgeR")


load("/dcl01/lieber/ajaffe/CellSorting/RNAseq_pipeline/rawCounts_CellSorting_July5_n12.rda")

# format phenotype table
metrics[grep("12", metrics$SampleID),"Age"] = "Neonate"
metrics[grep("13", metrics$SampleID),"Age"] = "Teenager"
metrics[grep("14", metrics$SampleID),"Age"] = "Toddler"
metrics$CellType = ifelse(metrics$NeuN=="NeuN_Minus", "Glia", "Neuron")

polya = metrics[which(metrics$Prep=="PolyA"),]
ribo = metrics[which(metrics$Prep=="Ribo"),]

# format counts
geneCounts = geneCounts[rowSums(geneCounts)>0,]
polyaCounts = geneCounts[,grep("PolyA", colnames(geneCounts))]
riboCounts = geneCounts[,grep("Ribo", colnames(geneCounts))]
exonCounts = exonCounts[rowSums(exonCounts)>0,]
polyaExons = exonCounts[,grep("PolyA", colnames(exonCounts))]
riboExons = exonCounts[,grep("Ribo", colnames(exonCounts))]

match(rownames(polya), colnames(polyaCounts))
match(rownames(ribo), colnames(riboCounts))
match(rownames(polya), colnames(polyaExons))
match(rownames(ribo), colnames(riboExons))

### Differential expression: gene level
# PolyA
design <- model.matrix(~ polya$Age + polya$CellType)
pdge <- DGEList(counts = polyaCounts)
pdge <- calcNormFactors(pdge)
pdat <- voom(pdge, design, plot = TRUE)
fit_gene_polya <- lmFit(pdat, design)
fit_gene_polya <- eBayes(fit_gene_polya)

# Ribo
design <- model.matrix(~ ribo$Age + ribo$CellType)
rdge <- DGEList(counts = riboCounts)
rdge <- calcNormFactors(rdge)
rdat <- voom(rdge, design, plot = TRUE)
fit_gene_ribo <- lmFit(rdat, design)
fit_gene_ribo <- eBayes(fit_gene_ribo)

# Together
combinedCounts = polyaCounts + riboCounts
match(rownames(polya), colnames(combinedCounts))
design <- model.matrix(~ polya$Age + polya$CellType)
rdge <- DGEList(counts = combinedCounts)
rdge <- calcNormFactors(rdge)
rdat <- voom(rdge, design, plot = F)
fit_gene_combined <- lmFit(rdat, design)
fit_gene_combined <- eBayes(fit_gene_combined)


### Differential expression: exon level
# PolyA
design <- model.matrix(~ polya$Age + polya$CellType)
pdee <- DGEList(counts = polyaExons)
pdee <- calcNormFactors(pdee)
pdat.exons <- voom(pdee, design, plot = TRUE)
fit_exon_polya <- lmFit(pdat.exons, design)
fit_exon_polya <- eBayes(fit_exon_polya)

# Ribo
design <- model.matrix(~ ribo$Age + ribo$CellType)
rdee <- DGEList(counts = riboExons)
rdee <- calcNormFactors(rdee)
rdat.exons <- voom(rdee, design, plot = TRUE)
fit_exon_ribo <- lmFit(rdat.exons, design)
fit_exon_ribo <- eBayes(fit_exon_ribo)

# Together
combinedExons = polyaExons + riboExons
match(rownames(polya), colnames(combinedExons))
design <- model.matrix(~ polya$Age + polya$CellType)
rdge <- DGEList(counts = combinedExons)
rdge <- calcNormFactors(rdge)
rdat <- voom(rdge, design, plot = F)
fit_exon_combined <- lmFit(rdat, design)
fit_exon_combined <- eBayes(fit_exon_combined)

## Reformat combined results

nucRNAres = as.data.frame(cbind(fit_gene_combined$coefficients, fit_gene_combined$t, fit_gene_combined$p.value))
colnames(nucRNAres) = c("Coeff.Intercept", "Coeff.AgeTeenager", "Coeff.AgeToddler", "Coeff.CellTypeNeuron",
                        "Tstat.Intercept", "Tstat.AgeTeenager", "Tstat.AgeToddler", "Tstat.CellTypeNeuron",
                        "pval.Intercept", "pval.AgeTeenager", "pval.AgeToddler", "pval.CellTypeNeuron")
nucRNAres$padj.AgeTeenager = p.adjust(nucRNAres$pval.AgeTeenager, method = "fdr")
nucRNAres$padj.AgeToddler = p.adjust(nucRNAres$pval.AgeToddler, method = "fdr")
nucRNAres$padj.CellTypeNeuron = p.adjust(nucRNAres$pval.CellTypeNeuron, method = "fdr")
nucRNAres$gencodeID = rownames(nucRNAres)

nucRNAexonres = as.data.frame(cbind(fit_exon_combined$coefficients, fit_exon_combined$t, fit_exon_combined$p.value))
colnames(nucRNAexonres) = c("Coeff.Intercept", "Coeff.AgeTeenager", "Coeff.AgeToddler", "Coeff.CellTypeNeuron",
                        "Tstat.Intercept", "Tstat.AgeTeenager", "Tstat.AgeToddler", "Tstat.CellTypeNeuron",
                        "pval.Intercept", "pval.AgeTeenager", "pval.AgeToddler", "pval.CellTypeNeuron")
nucRNAexonres$padj.AgeTeenager = p.adjust(nucRNAexonres$pval.AgeTeenager, method = "fdr")
nucRNAexonres$padj.AgeToddler = p.adjust(nucRNAexonres$pval.AgeToddler, method = "fdr")
nucRNAexonres$padj.CellTypeNeuron = p.adjust(nucRNAexonres$pval.CellTypeNeuron, method = "fdr")
nucRNAexonres$gencodeID = rownames(nucRNAexonres)

write.csv(metrics, quote=F, file="/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/sorted.RNAseq.pheno.info.csv")

save(fit_gene_polya,fit_gene_ribo,fit_exon_polya,fit_exon_ribo, geneMap, metrics,
     fit_gene_combined, fit_exon_combined, nucRNAres, nucRNAexonres, 
     file="/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/sorted_nuclear_RNA/DE_limma_results_objects.rda")