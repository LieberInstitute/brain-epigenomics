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
    'bootstrap', 'b', 1, 'logical', 'Whether to use Bootstrap or permutation method',
    'top100', 'i', 1, 'logical', 'Whether to look at the top 100 or another set',
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
        permutations = 0, 'bootstrap' = FALSE, 'top100' = TRUE)
}

dir.create('pdf', showWarnings = FALSE)

## Check inputs
inputFile <- paste0('bumps_bsseqSmooth_', opt$subset, '_', opt$model,
    '_', opt$permutations, ifelse(opt$bootstrap, '', '_perm'), '.Rdata')
stopifnot(file.exists(inputFile))

load(inputFile)

if(opt$permutations != 0) {
    load(paste0('BSobj_bsseqSmooth_', opt$subset, '_minCov_3.Rdata'))
    meth <- getMeth(BSobj, type = 'raw')
    pd <- pData(BSobj)
    rm(BSobj)
} else {
    load('../processed_beta_values_plusMap.rda')
}

## tale of top loci
if(opt$top100) {
    top_set <- seq_len(100)
    top_n <- 100
} else {
    top_set <- c(seq_len(100), unlist(sapply(seq(250, 2500, by = 250),
        function(x) { seq_len(10) + x }, simplify = FALSE)))
    top_n <- length(top_set)
}
tab = bumps$table[top_set,]
## top CpGs
topInds = mapply(function(s,e) s:e, tab$indexStart, tab$indexEnd)

meanMeth = lapply(topInds, function(ii) colMeans(t(t(meth[ii,]))))
meanMeth = do.call("rbind", meanMeth)

pdfFile <- paste0('pdf/bumps_bsseqSmooth_', opt$subset, '_', opt$model,
    '_', opt$permutations, ifelse(opt$bootstrap, '', '_perm'), '.pdf')

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


pdf(paste0('pdf/bumps_bsseqSmooth_', opt$subset, '_', opt$model,
    '_', opt$permutations, '_CpG_in_DMR', ifelse(opt$bootstrap, '', '_perm'),
    '.pdf'), width = 14)
for(i in seq_len(top_n)) {
    matplot(meth[topInds[[i]], ], pch = 20, bg = factor(pd$Cell.Type), ylim = c(0, 1), ylab = 'DNAm Level', xlab = 'CpG in DMR', col = factor(pd$Cell.Type))
}
dev.off()



pdf(paste0('pdf/bumps_bsseqSmooth_', opt$subset, '_', opt$model,
    '_', opt$permutations, '_boxplot_by_age', ifelse(opt$bootstrap, '',
    '_perm'), '.pdf'), width = 14)
for(i in seq_len(top_n)) {
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
load(paste0('BSobj_bsseqSmooth_', opt$subset, '_minCov_3.Rdata'))
#load(paste0('BSobj_bsseqSmooth_', opt$subset, '.Rdata'))
## Add colors
colData(BSobj)$col <- brewer.pal(8,"Dark2")[factor(pd$Cell.Type)]

## To get an idea of how much to extend
round(mean(tab$end - tab$start + 1))

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

colData(BSobj)$age_group <- factor(ifelse(colData(BSobj)$Age < 1, 'Infant',
    ifelse(colData(BSobj)$Age <= 12, 'Child',
    ifelse(colData(BSobj)$Age <= 17, 'Teen', 'Adult'))),
    levels = c('Infant', 'Child', 'Teen', 'Adult'))
    
colData(BSobj)$age_group_cell <- factor(paste0(colData(BSobj)$age_group, '_',
    colData(BSobj)$Cell.Type),
    levels = paste0(rep(levels(colData(BSobj)$age_group), each = 2),
    '_', c('Glia', 'Neuron')))
colData(BSobj)$col <- brewer.pal(8, "Paired")[c(5:6, 7:8, 3:4, 1:2)][colData(BSobj)$age_group_cell]

pdf(paste0('pdf/bumps_bsseqSmooth_', opt$subset, '_', opt$model,
    '_', opt$permutations, '_with_bsseq_age_cell', ifelse(opt$bootstrap, '',
    '_perm'), '.pdf'), width = 14)
palette(brewer.pal(8, "Paired")[c(5:6, 7:8, 3:4, 1:2)])
plot(colData(BSobj)$Age, type = 'p', pch = 21, ylab = 'Age',
    bg = colData(BSobj)$age_group_cell, cex = 3)
legend("bottomright", levels(colData(BSobj)$age_group_cell), pch = 15, col=1:8, cex=1.4)
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

    pdf(paste0('pdf/bumps_bsseqSmooth_', opt$subset, '_', opt$model,
        '_', opt$permutations, '_ATAC_cell', ifelse(opt$bootstrap, '', '_perm'),
        '.pdf'), width = 14)
    colData(BSobj)$col <- brewer.pal(8,"Dark2")[factor(pd$Cell.Type)]
    plotManyRegions(BSobj, regions = peaks_methDiffOrdered[seq_len(100)], extend = 2000, addRegions = peaks_methDiffOrdered, annoTrack = list(genes = genes, exons = exons))
    dev.off()


    pdf(paste0('pdf/bumps_bsseqSmooth_', opt$subset, '_', opt$model,
        '_', opt$permutations, '_ATAC_age_cell', ifelse(opt$bootstrap, '',
        '_perm'), '.pdf'), width = 14)
    palette(brewer.pal(8, "Paired")[c(5:6, 7:8, 3:4, 1:2)])
    colData(BSobj)$col <- brewer.pal(8, "Paired")[c(5:6, 7:8, 3:4, 1:2)][colData(BSobj)$age_group_cell]
    plot(colData(BSobj)$Age, type = 'p', pch = 21, ylab = 'Age',
        bg = colData(BSobj)$age_group_cell, cex = 3)
    legend("bottomright", levels(colData(BSobj)$age_group_cell), pch = 15, col=1:8, cex=1.4)
    plotManyRegions(BSobj, regions = peaks_methDiffOrdered[seq_len(100)], extend = 2000, addRegions = peaks_methDiffOrdered, annoTrack = list(genes = genes, exons = exons), regionCol = brewer.pal(8, 'Greys')[2])
    dev.off()
}




## Reproducibility info
proc.time()
message(Sys.time())
options(width = 120)
session_info()
