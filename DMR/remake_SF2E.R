library(RColorBrewer)

load('/dcl01/lieber/ajaffe/lab/brain-epigenomics/bumphunting/rda/limma_Neuron_CpGs_minCov_3_ageInfo_dmrs.Rdata')

cols <- brewer.pal(8, 'Dark2')
cutoff <- quantile(abs(dmrs$age_cell_difference_tstat_mean), 0.025)

pdf('/dcl01/lieber/ajaffe/lab/brain-epigenomics/DMR/figures/remade_figS2E.pdf', 
    useDingbats = FALSE, height = 5, width = 5)

plot(x = dmrs$overall_tstat_mean, y = dmrs$age_cell_difference_tstat_mean, 
     col = cols[dmrs$k6cluster], 
     xlab = 'Overall age mean t-statistic', 
     ylab = 'Interaction mean t-statistic', 
     pch = 20, cex = 0.5, 
     main = 't-statistics for global age vs interaction', 
     sub = 'Lines at interaction quantile 2.5% and 2')
abline(v = c(2, -2), col = 'grey80')
abline(v = c(cutoff, - cutoff), col = 'grey80')

dev.off()
