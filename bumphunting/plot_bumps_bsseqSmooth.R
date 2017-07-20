library('bsseq')
library('bumphunter')
library('devtools')
library('getopt')
library('RColorBrewer')
library('ggplot2')
library('EnsDb.Hsapiens.v75')

## Specify parameters
spec <- matrix(c(
	'model', 'm', 1, 'character', 'Either: cell, age or interaction',
    'subset', 's', 1, 'character', 'Either: Homogenate or Neuron.',
    'cores', 't', 1, 'integer', 'Number of cores to use',
#    'chromosome', 'c', 1, 'character', 'One of chr1 to chr22, chrX, chrY or chrM'
    'permutations', 'p', 1, 'integer', 'Number of permutations to run',
	'help' , 'h', 0, 'logical', 'Display help'
), byrow=TRUE, ncol=5)
opt <- getopt(spec)

## if help was asked for print a friendly message
## and exit with a non-zero error code
if (!is.null(opt$help)) {
	cat(getopt(spec, usage=TRUE))
	q(status=1)
}

## For testing
if(FALSE) {
    opt <- list('model' = 'cell', 'subset' = 'Neuron',
        permutations = 0)
}

## Check inputs
inputFile <- paste0('bumps_bsseqSmooth_', opt$subset, '_', opt$model,
    '_', opt$permutations, '.Rdata')
stopifnot(file.exists(inputFile))

load(inputFile)

load("../processed_beta_values_plusMap.rda")

## tale of top loci
tab = bumps$table[1:100,]
## top CpGs
topInds = mapply(function(s,e) s:e, tab$indexStart, tab$indexEnd)

meanMeth = lapply(topInds, function(ii) colMeans(t(t(meth[ii,]))))
meanMeth = do.call("rbind", meanMeth)

pdfFile <- paste0('bumps_bsseqSmooth_', opt$subset, '_', opt$model,
    '_', opt$permutations, '.pdf')

pdf(pdfFile)
palette(brewer.pal(8,"Dark2"))
for(i in 1:nrow(meanMeth)) {
	plot(meanMeth[i,] ~ pd$Age, pch = 21,
		bg = factor(pd$Cell.Type),ylim = c(0,1),
		ylab = "DNAm Level", xlab = 'Age')
	abline(lm(meanMeth[i,] ~ pd$Age, subset = pd$Cell.Type =="Glia"),
		col = 1)
	abline(lm(meanMeth[i,] ~ pd$Age, subset = pd$Cell.Type =="Neuron"),
		col = 2)
    legend("right", levels(factor(pd$Cell.Type)), pch = 15, col=1:2, cex=1.4)
}
dev.off()


pdf(paste0('bumps_bsseqSmooth_', opt$subset, '_', opt$model,
    '_', opt$permutations, '_CpG_in_DMR.pdf'), width = 14)
for(i in seq_len(100)) {
    matplot(meth[topInds[[i]], ], pch = 20, bg = factor(pd$Cell.Type), ylim = c(0, 1), ylab = 'DNAm Level', xlab = 'CpG in DMR', col = factor(pd$Cell.Type))
}
dev.off()



pdf(paste0('bumps_bsseqSmooth_', opt$subset, '_', opt$model,
    '_', opt$permutations, '_boxplot_by_age.pdf'), width = 14)
for(i in seq_len(100)) {
    df <-  data.frame(
        Meth = as.vector(t(meth[topInds[[i]], ])),
        Age = rep(pd$Age, each = length(topInds[[i]])),
        cell = rep(pd$Cell.Type, each = length(topInds[[i]])),
        dmr = rep(seq_len(length(topInds[[i]])), each = nrow(pd))
    )
    g <- ggplot(data = df, aes(x = factor(Age), y = Meth, fill = cell)) + geom_boxplot() + theme_bw() + scale_fill_brewer(palette = 'Dark2') + xlab('Age') + ylab('DNAm Level')
    print(g)
}
dev.off()

## Use the plotting code from bsseq
load(paste0('BSobj_bsseqSmooth_', opt$subset, '.Rdata'))
## Add colors
colData(BSobj)$col <- brewer.pal(8,"Dark2")[factor(pd$Cell.Type)]

## To get an idea of how much to extend
round(mean(bumps$table$end[1:100] - bumps$table$start[1:100] + 1))

genes <- genes(EnsDb.Hsapiens.v75)
seqlevels(genes) <- paste0('chr', seqlevels(genes))
regions_gr <- GRanges(seqnames = bumps$table$chr[1:100], IRanges(bumps$table$start[1:100], bumps$table$end[1:100]))
regions_gr <- resize(regions_gr, width(regions_gr) + 40000, fix = 'center')
genes <- genes[countOverlaps(genes, regions_gr) > 0]

exons <- exons(EnsDb.Hsapiens.v75)
seqlevels(exons) <- paste0('chr', seqlevels(exons))
exons <- exons[countOverlaps(exons, regions_gr) > 0]



pdf(paste0('bumps_bsseqSmooth_', opt$subset, '_', opt$model,
    '_', opt$permutations, '_with_bsseq_cell.pdf'), width = 14)
palette(brewer.pal(8,"Dark2"))
plotManyRegions(BSobj, regions = bumps$table[1:100, ], extend = 20000, addRegions = bumps$table, annoTrack = list(genes = genes, exons = exons))
dev.off()

colData(BSobj)$age_group <- factor(ifelse(colData(BSobj)$Age < 1, 'Infant',
    ifelse(colData(BSobj)$Age <= 12, 'Child',
    ifelse(colData(BSobj)$Age <= 17, 'Teen', 'Adult'))),
    levels = c('Infant', 'Child', 'Teen', 'Adult'))
    
colData(BSobj)$age_group_cell <- factor(paste0(colData(BSobj)$age_group, '_',
    colData(BSobj)$Cell.Type),
    levels = paste0(rep(levels(colData(BSobj)$age_group), each = 2),
    '_', c('Glia', 'Neuron')))
colData(BSobj)$col <- brewer.pal(8, "Paired")[c(5:6, 7:8, 3:4, 1:2)][colData(BSobj)$age_group_cell]

pdf(paste0('bumps_bsseqSmooth_', opt$subset, '_', opt$model,
    '_', opt$permutations, '_with_bsseq_age_cell.pdf'), width = 14)
palette(brewer.pal(8, "Paired")[c(5:6, 7:8, 3:4, 1:2)])
plot(colData(BSobj)$Age, type = 'p', pch = 21, ylab = 'Age',
    bg = colData(BSobj)$age_group_cell, cex = 3)
legend("bottomright", levels(colData(BSobj)$age_group_cell), pch = 15, col=1:8, cex=1.4)
plotManyRegions(BSobj, regions = bumps$table[1:100, ], extend = 20000, addRegions = bumps$table, annoTrack = list(genes = genes, exons = exons), regionCol = brewer.pal(8, 'Greys')[2])
dev.off()

## ATAC-seq info
load('/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/ATAC_peaks_methDiffOrdered.rda')

## Subset annotation to match these regions
genes <- genes(EnsDb.Hsapiens.v75)
seqlevels(genes) <- paste0('chr', seqlevels(genes))
regions_gr <- peaks_methDiffOrdered[1:100]
regions_gr <- resize(regions_gr, width(regions_gr) + 4000, fix = 'center')
genes <- genes[countOverlaps(genes, regions_gr) > 0]

exons <- exons(EnsDb.Hsapiens.v75)
seqlevels(exons) <- paste0('chr', seqlevels(exons))
exons <- exons[countOverlaps(exons, regions_gr) > 0]

pdf(paste0('bumps_bsseqSmooth_', opt$subset, '_', opt$model,
    '_', opt$permutations, '_ATAC_cell.pdf'), width = 14)
colData(BSobj)$col <- brewer.pal(8,"Dark2")[factor(pd$Cell.Type)]
plotManyRegions(BSobj, regions = peaks_methDiffOrdered[1:100], extend = 2000, addRegions = peaks_methDiffOrdered, annoTrack = list(genes = genes, exons = exons))
dev.off()


pdf(paste0('bumps_bsseqSmooth_', opt$subset, '_', opt$model,
    '_', opt$permutations, '_ATAC_age_cell.pdf'), width = 14)
palette(brewer.pal(8, "Paired")[c(5:6, 7:8, 3:4, 1:2)])
colData(BSobj)$col <- brewer.pal(8, "Paired")[c(5:6, 7:8, 3:4, 1:2)][colData(BSobj)$age_group_cell]
plot(colData(BSobj)$Age, type = 'p', pch = 21, ylab = 'Age',
    bg = colData(BSobj)$age_group_cell, cex = 3)
legend("bottomright", levels(colData(BSobj)$age_group_cell), pch = 15, col=1:8, cex=1.4)
plotManyRegions(BSobj, regions = peaks_methDiffOrdered[1:100], extend = 2000, addRegions = peaks_methDiffOrdered, annoTrack = list(genes = genes, exons = exons), regionCol = brewer.pal(8, 'Greys')[2])
dev.off()



## Reproducibility info
proc.time()
message(Sys.time())
options(width = 120)
session_info()
