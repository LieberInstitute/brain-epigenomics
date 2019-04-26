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

## load prenatal CpGs
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/bumphunting/BSobj_bsseqSmooth_Neuron_minCov_3_prenatal.Rdata")
BSobj <- updateObject(BSobj, verbose=TRUE)
meth_cgP = getMeth(BSobj, type = 'raw')
methMap_cgP = granges(BSobj)

meth_cg = cbind(meth_cg, as.matrix(meth_cgP))
pd = rbind(pd[1:32,c("Cell.Type", "Age.Bin", "Working.Num")], data.frame("Cell.Type" = "Homogenate Prenatal",
                                                                     "Age.Bin" = "P", "Working.Num" = 101:120))

## load dmrs
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/bumphunting/bumps_bsseqSmooth_Neuron_interaction_250_perm.Rdata")
sigInt = bumps$table[bumps$table$fwer < 0.05,]

## subset to the DMRs ###
ooInt_cg = findOverlaps(makeGRangesFromDataFrame(sigInt), methMap_cg)

# first cluster all
meth_Int_cg = meth_cg[queryHits(ooInt_cg),]
dd_Int_cg = dist(t(meth_Int_cg))

dd_Int_cg_mat <- as.matrix(dd_Int_cg)
colnames(dd_Int_cg_mat) = rownames(dd_Int_cg_mat) = paste(
  pd$Cell.Type, pd$Age.Bin, pd$Working.Num, sep = ":")


# width
meanMeth_cg_Int = do.call("rbind", sapply(split(subjectHits(ooInt_cg), 
                                                factor(queryHits(ooInt_cg), levels=1:nrow(sigInt))), 
                                          function(ii) colMeans(t(t(meth_cg[ii,])))))
dd_Int_cg_mean = dist(t(meanMeth_cg_Int))
dd_Int_cg_mean_mat <- as.matrix(dd_Int_cg_mean)
colnames(dd_Int_cg_mean_mat) = rownames(dd_Int_cg_mean_mat) = paste(
  pd$Cell.Type, pd$Age.Bin, pd$Working.Num, sep = ":")

save(dd_Int_cg_mat, dd_Int_cg, dd_Int_cg_mean_mat, dd_Int_cg_mean,
     file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/DMR/CG_heatmap_objects.rda")

names = rownames(dd_Int_cg_mat)
names = strsplit(gsub("Young.", "", names), ":", fixed = T)
names = unlist(lapply(names, function(x) paste0(x[1], ": ", x[2])))
names = gsub("Early.", "", names)
names = gsub("Toddler", "Child", names)
names = gsub("Neonate", "Infant", names)
names = gsub(": P", "", names)
names = data.frame(names, row.names = rownames(dd_Int_cg_mat))

leg = data.frame(labs = c(paste(c('Glia:', 'Neuron:'), 
                                rep(c('Infant', 'Child', 'Teen', 'Adult'), each = 2)),
                          'Homogenate Prenatal'), cols = c(brewer.pal(8, "Paired")[c(5:6, 7:8, 3:4, 1:2)], 
                                                           'grey50'))
colop = list(names = as.character(leg[match(names$names, leg$labs),"cols"]))
names(colop$names) = names$names 

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/DMR/figures/updated_euclidean_distance_CpG_cdDMRs.pdf", height = 5, width = 5)
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")) )(255)
pheatmap(dd_Int_cg_mat,
         clustering_distance_rows=dd_Int_cg, clustering_distance_cols=dd_Int_cg,
         col=colors, annotation_legend = F, 
         annotation_row = names, annotation_colors = colop,
         cutree_rows = 2, cutree_cols = 2,
         show_rownames = F, show_colnames = F, annotation_names_row = F,
         main = "cdDMRs - all CpGs")

pheatmap(dd_Int_cg_mean_mat,
         clustering_distance_rows=dd_Int_cg_mean, clustering_distance_cols=dd_Int_cg_mean,
         col=colors, annotation_legend = F, 
         annotation_row = names, annotation_colors = colop,
         cutree_rows = 2, cutree_cols = 2,
         show_rownames = F, show_colnames = F, annotation_names_row = F,
         main = "cdDMRs - mean CpGs")

dev.off()

