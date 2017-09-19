library('bsseq')
library('bumphunter')
library('devtools')
library('RColorBrewer')
library('ggplot2')
library('EnsDb.Hsapiens.v75')

    
dir.create('pdf', showWarnings = FALSE)


## Load CpG data
load(file.path('/dcl01/lieber/ajaffe/lab/brain-epigenomics', 'bumphunting',
   'BSobj_bsseqSmooth_Neuron_minCov_3.Rdata'))

## Regions of interest to plot
regions <- GRanges(c('chr12', 'chr19'),
    IRanges(
        c(2162416, 45392785),
        c(2807115, 45428904)
    )
)
names(regions) <- c('CACNA1C', 'APOE')
regions_long <- resize(regions, width(regions) + 40000, fix = 'center')

## Subset BSobj
ov <- findOverlaps(regions_long, rowRanges(BSobj))
BSobj <- BSobj[subjectHits(ov), ]

## Look only inside the original regions, not outside
gr_meth <- rowRanges(BSobj[subjectHits(findOverlaps(regions, rowRanges(BSobj)))])

## Define new clusters and add them to the regions
clusters <- clusterMaker(seqnames(gr_meth), start(gr_meth), maxGap = 1000)
regs_w_clus <- c(regions, unlist(range(split(gr_meth, clusters))))


## Get annotation
genes <- genes(EnsDb.Hsapiens.v75)
seqlevels(genes) <- paste0('chr', seqlevels(genes))
genes <- genes[countOverlaps(genes, regions_long) > 0]

exons <- exons(EnsDb.Hsapiens.v75)
seqlevels(exons) <- paste0('chr', seqlevels(exons))
exons <- exons[countOverlaps(exons, regions_long) > 0]

## Get the main exons to highlight
exonid <- select(EnsDb.Hsapiens.v75, keys = c('APOE', 'CACNA1C'), keytype = 'GENENAME', columns = c('EXONID', 'GENENAME'))
main_ex <- exons[names(exons) %in% exonid$EXONID]

## Make the plots
colData(BSobj)$col <- brewer.pal(8,"Dark2")[factor(colData(BSobj)$Cell.Type)]
pdf('pdf/CACNA1C_APOE_with_bsseq_cell.pdf', width = 14)
palette(brewer.pal(8,"Dark2"))
plotManyRegions(BSobj, regions = regs_w_clus, extend = 20000, addRegions = main_ex, annoTrack = list(genes = genes, exons = exons), regionCol = brewer.pal(8, 'Greys')[2])
dev.off()

## Add age groups
colData(BSobj)$age_group <- factor(ifelse(colData(BSobj)$Age < 1, 'Infant',
    ifelse(colData(BSobj)$Age <= 12, 'Child',
    ifelse(colData(BSobj)$Age <= 17, 'Teen', 'Adult'))),
    levels = c('Infant', 'Child', 'Teen', 'Adult'))
    
colData(BSobj)$age_group_cell <- factor(paste0(colData(BSobj)$age_group, '_',
    colData(BSobj)$Cell.Type),
    levels = paste0(rep(levels(colData(BSobj)$age_group), each = 2),
    '_', c('Glia', 'Neuron')))
colData(BSobj)$col <- brewer.pal(8, "Paired")[c(5:6, 7:8, 3:4, 1:2)][colData(BSobj)$age_group_cell]

pdf('pdf/CACNA1C_APOE_with_bsseq_age_cell.pdf', width = 14)
palette(brewer.pal(8, "Paired")[c(5:6, 7:8, 3:4, 1:2)])
plot(colData(BSobj)$Age, type = 'p', pch = 21, ylab = 'Age',
    bg = colData(BSobj)$age_group_cell, cex = 3)
legend("bottomright", levels(colData(BSobj)$age_group_cell), pch = 15, col=1:8, cex=1.4)
plotManyRegions(BSobj, regions = regs_w_clus, extend = 20000, addRegions = main_ex, annoTrack = list(genes = genes, exons = exons), regionCol = brewer.pal(8, 'Greys')[2])
dev.off()

## Reproducibility info
proc.time()
message(Sys.time())
options(width = 120)
session_info()
