####
####

library("limma")
library("edgeR")
library(bsseq)
library(recount)
library(pheatmap)
library(RColorBrewer)

## load expression data
load("/dcl01/lieber/ajaffe/CellSorting/RNAseq_pipeline/rse_jx_CellSorting_July5_n12.Rdata")

## load CpG 
load('/dcl01/lieber/ajaffe/lab/brain-epigenomics/bumphunting/BSobj_bsseqSmooth_Neuron_minCov_3.Rdata')
BSobj_cg = BSobj
meth_cg =getMeth(BSobj_cg, type = 'raw')
methMap_cg =granges(BSobj_cg)

## and non-CpG data
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

## subset to the DMRs ###
ooInt_ch = findOverlaps(makeGRangesFromDataFrame(sigInt), methMap_ch)

# first cluster all
meth_Int_ch = meth_ch[queryHits(ooInt_ch),]
dd_Int_ch = dist(t(meth_Int_ch))

dd_Int_ch_mat <- as.matrix(dd_Int_ch)
colnames(dd_Int_ch_mat) = rownames(dd_Int_ch_mat) = paste(
		pd$Cell.Type, pd$Age.Bin, pd$Working.Num, sep = ":")


# width
meanMeth_ch_Int = do.call("rbind", sapply(split(subjectHits(ooInt_ch), 
	factor(queryHits(ooInt_ch), levels=1:nrow(sigInt))), 
		function(ii) colMeans(t(t(meth_ch[ii,])))))
dd_Int_ch_mean = dist(t(meanMeth_ch_Int))
dd_Int_ch_mean_mat <- as.matrix(dd_Int_ch_mean)
colnames(dd_Int_ch_mean_mat) = rownames(dd_Int_ch_mean_mat) = paste(
		pd$Cell.Type, pd$Age.Bin, pd$Working.Num, sep = ":")

save(dd_Int_ch_mat, dd_Int_ch, dd_Int_ch_mean_mat, dd_Int_ch_mean,
     file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/DMR/CH_heatmap_objects.rda")

names = rownames(dd_Int_ch_mat)
names = strsplit(gsub("Young.", "", names), ":", fixed = T)
names = unlist(lapply(names, function(x) paste0(x[1], ": ", x[2])))
names = gsub("Early.", "", names)
names = gsub("Toddler", "Child", names)
names = gsub("Neonate", "Infant", names)
names = data.frame(names, row.names = rownames(dd_Int_ch_mat))

leg = data.frame(labs = c(paste(c('Glia:', 'Neuron:'), 
                                rep(c('Infant', 'Child', 'Teen', 'Adult'), each = 2)),
                          'Homogenate Prenatal'), cols = c(brewer.pal(8, "Paired")[c(5:6, 7:8, 3:4, 1:2)], 
                                                           'grey50'))
colop = list(names = as.character(leg[match(names$names, leg$labs),"cols"]))
names(colop$names) = names$names 

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/non-CpG/figures/updated_euclidean_distance_CpH_cdDMRs.pdf", height = 5, width = 5)
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")) )(255)
pheatmap(dd_Int_ch_mat,
         clustering_distance_rows=dd_Int_ch, clustering_distance_cols=dd_Int_ch,
         col=colors, annotation_legend = F, 
         annotation_row = names, annotation_colors = colop,
         cutree_rows = 2, cutree_cols = 2,
         show_rownames = F, show_colnames = F, annotation_names_row = F,
         main = "2178 Interaction DMRs - all non-CpG")

pheatmap(dd_Int_ch_mean_mat,
         clustering_distance_rows=dd_Int_ch_mean, clustering_distance_cols=dd_Int_ch_mean,
         col=colors, annotation_legend = F, 
         annotation_row = names, annotation_colors = colop,
         cutree_rows = 2, cutree_cols = 2,
         show_rownames = F, show_colnames = F, annotation_names_row = F,
         main = "2178 Interaction DMRs - mean non-CpG")

dev.off()

# (old is called tmp.pdf)

