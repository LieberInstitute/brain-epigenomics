library('bsseq')
library('bumphunter')
library('devtools')
library('getopt')
library('RColorBrewer')
library('ggplot2')
library('EnsDb.Hsapiens.v75')

load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/DMR/DMR_objects.rda")
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/bumphunting/BSobj_bsseqSmooth_Neuron_minCov_3_prenatal_combined.Rdata")


pd <- pData(BSobj)
dim(BSobj)

## Add colors
colData(BSobj)$col <- brewer.pal(8,"Dark2")[ifelse(colData(BSobj)$Cell.Type == 'Glia', 1, 2)]
colData(BSobj)$col[colData(BSobj)$Age < 0] <- 'grey50'

genes <- genes(EnsDb.Hsapiens.v75)
seqlevels(genes) <- paste0('chr', seqlevels(genes))
regions_gr <- GRanges(seqnames = tab$chr, IRanges(tab$start, tab$end))
regions_gr <- resize(regions_gr, width(regions_gr) + 40000, fix = 'center')
genes <- genes[countOverlaps(genes, regions_gr) > 0]

exons <- exons(EnsDb.Hsapiens.v75)
seqlevels(exons) <- paste0('chr', seqlevels(exons))
exons <- exons[countOverlaps(exons, regions_gr) > 0]



pdf(paste0('pdf/bumps_bsseqSmooth_', opt$subset, '_', opt$model,
           '_', opt$permutations, '_with_bsseq_cell', ifelse(opt$bootstrap, '',
                                                             '_perm'), '.pdf'), width = 14)
palette(brewer.pal(8,"Dark2"))
plotManyRegions(BSobj, regions = tab, extend = 20000, addRegions = subset(bumps$table, fwer < 0.05), annoTrack = list(genes = genes, exons = exons))
dev.off()

colData(BSobj)$age_group <- factor(ifelse(colData(BSobj)$Age < 0, 'Prenatal',
                                          ifelse(colData(BSobj)$Age < 1, 'Infant',
                                                 ifelse(colData(BSobj)$Age <= 12, 'Child',
                                                        ifelse(colData(BSobj)$Age <= 17, 'Teen', 'Adult')))),
                                   levels = c('Infant', 'Child', 'Teen', 'Adult', 'Prenatal'))

colData(BSobj)$age_group_cell <- factor(paste0(colData(BSobj)$age_group, '_',
                                               colData(BSobj)$Cell.Type),
                                        levels = c(paste0(rep(levels(colData(BSobj)$age_group)[1:4], each = 2),
                                                          '_', c('Glia', 'Neuron')), 'Prenatal_H'))
colData(BSobj)$col <- c(brewer.pal(8, "Paired"), 'grey50')[c(5:6, 7:8, 3:4, 1:2, 9)][colData(BSobj)$age_group_cell]

pdf(paste0('pdf/bumps_bsseqSmooth_', opt$subset, '_', opt$model,
           '_', opt$permutations, '_with_bsseq_age_cell', ifelse(opt$bootstrap, '',
                                                                 '_perm'), '.pdf'), width = 14)
palette(c(brewer.pal(8, "Paired")[c(5:6, 7:8, 3:4, 1:2)], 'grey50'))
plot(colData(BSobj)$Age, type = 'p', pch = 21, ylab = 'Age',
     bg = colData(BSobj)$age_group_cell, cex = 3)
legend("topright", levels(colData(BSobj)$age_group_cell), pch = 15, col=1:9, cex=1.4)
plotManyRegions(BSobj, regions = tab, extend = 20000, addRegions = subset(bumps$table, fwer < 0.05), annoTrack = list(genes = genes, exons = exons), regionCol = brewer.pal(8, 'Greys')[2])
dev.off()

## ATAC-seq info: it's the same for all models, so only do this for the age
if(opt$model == 'age') {
  load('/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/ATAC_peaks_methDiffOrdered.rda')
  
  ## Subset annotation to match these regions
  genes <- genes(EnsDb.Hsapiens.v75)
  seqlevels(genes) <- paste0('chr', seqlevels(genes))
  regions_gr <- peaks_methDiffOrdered[seq_len(100)]
  regions_gr <- resize(regions_gr, width(regions_gr) + 4000, fix = 'center')
  genes <- genes[countOverlaps(genes, regions_gr) > 0]
  
  exons <- exons(EnsDb.Hsapiens.v75)
  seqlevels(exons) <- paste0('chr', seqlevels(exons))
  exons <- exons[countOverlaps(exons, regions_gr) > 0]
  
  pdf(paste0('pdf/bumps_bsseqSmooth_', opt$subset, '_', opt$permutations,
             '_ATAC_cell', ifelse(opt$bootstrap, '', '_perm'),
             '.pdf'), width = 14)
  colData(BSobj)$col <- brewer.pal(8,"Dark2")[ifelse(colData(BSobj)$Cell.Type == 'Glia', 1, 2)]
  colData(BSobj)$col[colData(BSobj)$Age < 0] <- 'grey50'
  
  plotManyRegions(BSobj, regions = peaks_methDiffOrdered[seq_len(100)], extend = 2000, addRegions = peaks_methDiffOrdered, annoTrack = list(genes = genes, exons = exons))
  dev.off()
  
  
  pdf(paste0('pdf/bumps_bsseqSmooth_', opt$subset, '_', opt$permutations,
             '_ATAC_age_cell', ifelse(opt$bootstrap, '', '_perm'), '.pdf'),
      width = 14)
  palette(c(brewer.pal(8, "Paired")[c(5:6, 7:8, 3:4, 1:2)], 'grey50'))
  colData(BSobj)$col <- c(brewer.pal(8, "Paired"), 'grey50')[c(5:6, 7:8, 3:4, 1:2, 9)][colData(BSobj)$age_group_cell]
  
  plot(colData(BSobj)$Age, type = 'p', pch = 21, ylab = 'Age',
       bg = colData(BSobj)$age_group_cell, cex = 3)
  legend("topright", levels(colData(BSobj)$age_group_cell), pch = 15, col=1:8, cex=1.4)
  plotManyRegions(BSobj, regions = peaks_methDiffOrdered[seq_len(100)], extend = 2000, addRegions = peaks_methDiffOrdered, annoTrack = list(genes = genes, exons = exons), regionCol = brewer.pal(8, 'Greys')[2])
  dev.off()
}







## Plot interaction PGC2 DMRs

pgc2 = c("GRIN2A","CACNA1C","TCF4","SATB2","D2D2","AS3MT","GRM3","NRGN")

pdf('/dcl01/lieber/ajaffe/lab/brain-epigenomics/DMR/CT_Age_Interaction/figures/PGC2_candidateGenes_plots.pdf', width = 14)

palette(c(brewer.pal(8, "Paired")[c(5:6, 7:8, 3:4, 1:2)], 'grey50'))

plot(colData(BSobj)$Age, type = 'p', pch = 21, ylab = 'Age',
     bg = colData(BSobj)$age_group_cell, cex = 3)
legend("topright", levels(colData(BSobj)$age_group_cell), pch = 15, col=1:9, cex=1.4)
plotManyRegions(BSobj, regions = tab, extend = 20000, addRegions = subset(bumps$table, fwer < 0.05), annoTrack = list(genes = genes, exons = exons), regionCol = brewer.pal(8, 'Greys')[2])
dev.off()