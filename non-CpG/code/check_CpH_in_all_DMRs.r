library(bsseq)
library(pheatmap)
library(RColorBrewer)


## load CpH data
load('/dcl01/lieber/ajaffe/lab/brain-epigenomics/bsseq/bsobj_by_chr/allChrs_postNatal_cleaned_nonCG_noHomogenate_highCov.Rdata')
BSobj_ch = BSobj
meth_ch =getMeth(BSobj_ch, type = 'raw')
methMap_ch = granges(BSobj_ch)

## load dmrs
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/bumphunting/bumps_bsseqSmooth_Neuron_interaction_250_perm.Rdata")
sigInt = bumps$table[bumps$table$fwer < 0.05,]
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/bumphunting/bumps_bsseqSmooth_Neuron_age_250_perm.Rdata")
sigAge = bumps$table[bumps$table$fwer < 0.05,]
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/bumphunting/bumps_bsseqSmooth_Neuron_cell_250_perm.Rdata")
sigCT = bumps$table[bumps$table$fwer < 0.05,]
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/bumphunting/rda/limma_Neuron_CpGs_minCov_3_ageInfo_dmrs.Rdata")
dmrs = split(dmrs, dmrs$k6cluster_label)

## subset to the DMRs
oo_ch = lapply(c(list(CellType = makeGRangesFromDataFrame(sigCT), Age = makeGRangesFromDataFrame(sigAge), Interaction = makeGRangesFromDataFrame(sigInt)), as.list(dmrs)), function(x) findOverlaps(x, methMap_ch))

# first cluster all
meth_DMR_ch = lapply(oo_ch, function(x) meth_ch[queryHits(x),])
dd_ch = lapply(meth_DMR_ch, function(x) dist(t(x)))
dd_ch_mat = lapply(dd_ch, as.matrix)
for (i in 1:length(dd_ch_mat)) { colnames(dd_ch_mat[[i]]) = rownames(dd_ch_mat[[i]]) = paste(pd$Cell.Type, pd$Age.Bin, pd$Working.Num, sep = ":") }


# mean meth in each DMR
meanMeth_ch = lapply(oo_ch, function(oo) do.call("rbind", lapply(split(subjectHits(oo), queryHits(oo)), function(ii) colMeans(t(t(meth_ch[ii,]))))))
dd_ch_mean = lapply(meanMeth_ch, function(x) dist(t(x)))
dd_ch_mean_mat = lapply(dd_ch_mean, as.matrix)
for (i in 1:length(dd_ch_mean_mat)) { colnames(dd_ch_mean_mat[[i]]) = rownames(dd_ch_mean_mat[[i]]) = paste(pd$Cell.Type, pd$Age.Bin, pd$Working.Num, sep = ":") }


pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/non-CpG/figures/nonCG_heatmap_within_DMRs.pdf")
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")) )(255)
for (i in 1:length(dd_ch_mat)) {
  print(pheatmap(dd_ch_mat[[i]], clustering_distance_rows = dd_ch[[i]], clustering_distance_cols = dd_ch[[i]],
                 col = colors, main = paste0("all CpH - Euclidean Distance - ", c("Cell Type","Age","Interaction","Group 1 Interaction","Group 2 Interaction", "Group 3 Interaction", "Group 4 Interaction",
                                                                                      "Group 5 Interaction", "Group 6 Interaction")[i], " DMRs")))
  print(pheatmap(dd_ch_mean_mat[[i]], clustering_distance_rows = dd_ch_mean[[i]], clustering_distance_cols = dd_ch_mean[[i]],
                 col = colors, main = paste0("mean mCpH - Euclidean Distance - ", c("Cell Type", "Age", "Interaction", "Group 1 Interaction", "Group 2 Interaction", "Group 3 Interaction", "Group 4 Interaction",
                                                                                       "Group 5 Interaction", "Group 6 Interaction")[i], " DMRs")))
}
dev.off()

