library(ggplot2)


### non-CpG: both neurons and glia ###

"/dcl01/lieber/ajaffe/lab/brain-epigenomics/bsseq/bsobj_by_chr/explore_nonCG_highCov.R"
"/dcl01/lieber/ajaffe/lab/brain-epigenomics/bsseq/bsobj_by_chr/pca_nonCG_highCov.pdf"
"/dcl01/lieber/ajaffe/lab/brain-epigenomics/bsseq/bsobj_by_chr/pca_nonCG_highCov_PC1vsAge.pdf"

### non-CpG: only neurons ###

load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/non-CpG/nonCG_highCov_neuronsOnly_pca_pd_methTable.Rdata")

pdf('/dcl01/lieber/ajaffe/lab/brain-epigenomics/non-CpG/figures/pca_nonCG_highCov_neuronsOnly.pdf')
barplot(pcaVars[1:10], col = '#377EB8', ylab = 'Percent of Variance Explained')

plot(pcs$x[, 1] ~ pcs$x[, 2],
     ylab = paste0('PC1: ', pcaVars[1], '% of Var Explained'),
     xlab = paste0('PC2: ', pcaVars[2], '% of Var Explained'), pch = 19)

to_plot <- c('Cell.Type', 'Age', 'Age.Bin', 'PMI', 'Sex', 'Race', 'RIN', 'pH', 'Proportion.of.Neurons', 'yield.nuclei.mg.tissue', 'pg.DNA.nuclei.input', 'X260.280.DNA', 'Library.Technician', 'Flowcell', 'Reads', 'Percent.GreaterThan.Q30', 'avg.Cov', 'cov.sDev', 'Percent.Duplication', 'total_num_trimmed_reads', 'total_num_untrimmed_reads', 'alignment.efficiency')
to_plot <- which(colnames(pd) %in% to_plot)
names(to_plot) <- colnames(pd)[to_plot]

for(pc in 1:4) {
  mapply(function(id, name) {
    plot_twoway(y = pcs$x[, pc], x = pd[, id],
                yvar = paste0('PC', pc, ': ', pcaVars[pc], '% of Var Explained'),
                xvar = name, color = '#377EB8', pal = 'Set1')
    return(invisible(NULL))
  }, to_plot, names(to_plot))
}
dev.off()

pdf('/dcl01/lieber/ajaffe/lab/brain-epigenomics/non-CpG/figures/pca_nonCG_highCov_PC1vsAge_neuronsOnly.pdf')
palette(brewer.pal(8,"Dark2"))
par(mar=c(5,6,2,2), cex.axis=2, cex.lab=2)
plot(pcs$x[, 1] ~ pd$Age, cex=2, pch = 21, bg=factor(pd$Cell.Type),
     ylab = paste0('PC1: ', pcaVars[1], '% of Var Explained'),  xlab = "Age")
legend("right", levels(factor(pd$Cell.Type)), pch = 15, col=1:2, cex=1.4)
dev.off()
