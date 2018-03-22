library('bsseq')
library('bumphunter')
library('devtools')
library('getopt')
library('RColorBrewer')
library('ggplot2')
library('EnsDb.Hsapiens.v75')

## Use interaction DMRs
opt <- list('model' = 'interaction', 'subset' = 'Neuron',
    permutations = 250, 'bootstrap' = FALSE)

## Load DMRs
inputFile <- paste0('bumps_bsseqSmooth_', opt$subset, '_', opt$model,
    '_', opt$permutations, ifelse(opt$bootstrap, '', '_perm'), '.Rdata')
stopifnot(file.exists(inputFile))

load(inputFile, verbose = FALSE)


## Load external data
load('/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/DMR/Birnbaum_PGC2_geneSet_enrichment_interactionClusters_toPLot.rda', verbose = TRUE)

## Load methylation data
load('BSobj_bsseqSmooth_Neuron_minCov_3_prenatal_combined.Rdata', verbose = TRUE)

## Process genes
genes <- genes(EnsDb.Hsapiens.v75)
seqlevels(genes) <- paste0('chr', seqlevels(genes))
regions_gr <- DE_OVERLAP
regions_gr <- resize(regions_gr, width(regions_gr) + 2000, fix = 'center')
genes <- genes[countOverlaps(genes, regions_gr) > 0]

exons <- exons(EnsDb.Hsapiens.v75)
seqlevels(exons) <- paste0('chr', seqlevels(exons))
exons <- exons[countOverlaps(exons, regions_gr) > 0]


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

pdf('pdf/Birnbaum_PGC2_geneSet_enrichment_interactionClusters_toPLot.pdf', width = 14)
palette(c(brewer.pal(8, "Paired")[c(5:6, 7:8, 3:4, 1:2)], 'grey50'))
plot(colData(BSobj)$Age, type = 'p', pch = 21, ylab = 'Age',
    bg = colData(BSobj)$age_group_cell, cex = 3)
legend("topright", levels(colData(BSobj)$age_group_cell), pch = 15, col=1:9, cex=1.4)
plotManyRegions(BSobj, regions = DE_OVERLAP, extend = 1000, addRegions = subset(bumps$table, fwer < 0.05), annoTrack = list(genes = genes, exons = exons), regionCol = brewer.pal(8, 'Greys')[2])
dev.off()
