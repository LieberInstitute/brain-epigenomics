# qrsh -l bluejay,mem_free=40G,h_vmem=40G
library('bsseq')
library('bumphunter')
library('devtools')
library('RColorBrewer')
library('ggplot2')
library('EnsDb.Hsapiens.v75')

    
dir.create('pdf', showWarnings = FALSE)

window_extend <- 1000

## Load CpG data
load(file.path('/dcl01/lieber/ajaffe/lab/brain-epigenomics', 'bumphunting',
   'BSobj_bsseqSmooth_Neuron_minCov_3.Rdata'))

## Regions of interest to plot
regions <- GRanges(c('chr12', 'chr19', 'chr6'),
    IRanges(
        c(2162416, 45392785, 35541362),
        c(2807115, 45428904, 35696397)
    )
)
names(regions) <- c('CACNA1C', 'APOE', 'FKBP5')
regions_long <- resize(regions, width(regions) + window_extend * 2, fix = 'center')

## Subset BSobj
ov <- findOverlaps(regions_long, rowRanges(BSobj))
BSobj <- BSobj[subjectHits(ov), ]

## Look only inside the original regions, not outside
gr_meth <- rowRanges(BSobj[subjectHits(findOverlaps(regions, rowRanges(BSobj)))])

## Define new clusters and add them to the regions
clusters <- clusterMaker(seqnames(gr_meth), start(gr_meth), maxGap = 1000)
regs_w_clus <- unlist(range(split(gr_meth, clusters)))

## APOE only has 1 cluster, so remove it
countOverlaps(regions, regs_w_clus)
regs_w_clus <- regs_w_clus[-queryHits(findOverlaps(regs_w_clus, regions['APOE']))]
regs_w_clus <- c(regions, regs_w_clus)

## Get annotation
genes <- genes(EnsDb.Hsapiens.v75)
seqlevels(genes) <- paste0('chr', seqlevels(genes))
genes <- genes[countOverlaps(genes, regions_long) > 0]

exons <- exons(EnsDb.Hsapiens.v75)
seqlevels(exons) <- paste0('chr', seqlevels(exons))
exons <- exons[countOverlaps(exons, regions_long) > 0]

## Get the main exons to highlight
exonid <- select(EnsDb.Hsapiens.v75, keys = names(regions), keytype = 'GENENAME', columns = c('EXONID', 'GENENAME'))
main_ex <- exons[names(exons) %in% exonid$EXONID]


## Split by gene
for(gene in names(regions)) {
    regs_plot <- regs_w_clus[subjectHits(findOverlaps(regions[gene], regs_w_clus))]

    ## Make the plots
    colData(BSobj)$col <- brewer.pal(8,"Dark2")[factor(colData(BSobj)$Cell.Type)]
    pdf(paste0('pdf/', gene, '_with_bsseq_cell.pdf'), width = 14)
    palette(brewer.pal(8,"Dark2"))
    plotManyRegions(BSobj, regions = regs_plot, extend = window_extend, addRegions = main_ex, annoTrack = list(genes = genes, exons = exons), regionCol = brewer.pal(8, 'Greys')[2])
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

    pdf(paste0('pdf/', gene, '_with_bsseq_age_cell.pdf'), width = 14)
    palette(brewer.pal(8, "Paired")[c(5:6, 7:8, 3:4, 1:2)])
    plot(colData(BSobj)$Age, type = 'p', pch = 21, ylab = 'Age',
        bg = colData(BSobj)$age_group_cell, cex = 3)
    legend("bottomright", levels(colData(BSobj)$age_group_cell), pch = 15, col=1:8, cex=1.4)
    plotManyRegions(BSobj, regions = regs_plot, extend = window_extend, addRegions = main_ex, annoTrack = list(genes = genes, exons = exons), regionCol = brewer.pal(8, 'Greys')[2])
    dev.off()
}

## Reproducibility info
proc.time()
message(Sys.time())
options(width = 120)
session_info()

# Session info ----------------------------------------------------------------------------------------------------------
#  setting  value
#  version  R version 3.3.3 Patched (2017-03-15 r72696)
#  system   x86_64, linux-gnu
#  ui       X11
#  language (EN)
#  collate  en_US.UTF-8
#  tz       <NA>
#  date     2017-10-02
#
# Packages --------------------------------------------------------------------------------------------------------------
#  package                * version  date       source
#  AnnotationDbi          * 1.36.2   2017-02-07 Bioconductor
#  AnnotationHub            2.6.5    2017-05-22 Bioconductor
#  base                   * 3.3.3    2017-05-18 local
#  Biobase                * 2.34.0   2016-10-20 Bioconductor
#  BiocGenerics           * 0.20.0   2017-06-05 Bioconductor
#  BiocInstaller            1.24.0   2016-10-20 Bioconductor
#  BiocParallel             1.8.2    2017-04-28 cran (@1.8.2)
#  biomaRt                  2.30.0   2016-12-06 Bioconductor
#  Biostrings               2.42.1   2016-12-06 Bioconductor
#  bit                      1.1-12   2014-04-09 CRAN (R 3.3.0)
#  bit64                    0.9-7    2017-05-08 CRAN (R 3.3.3)
#  bitops                   1.0-6    2013-08-17 CRAN (R 3.3.0)
#  blob                     1.1.0    2017-06-17 CRAN (R 3.3.3)
#  bsseq                  * 1.10.0   2016-10-20 Bioconductor
#  bumphunter             * 1.14.0   2017-06-05 Bioconductor
#  codetools                0.2-15   2016-10-05 CRAN (R 3.3.3)
#  colorout               * 1.1-2    2016-07-22 Github (jalvesaq/colorout@6d84420)
#  colorspace               1.3-2    2016-12-14 CRAN (R 3.3.1)
#  data.table               1.10.4   2017-02-01 CRAN (R 3.3.1)
#  datasets               * 3.3.3    2017-05-18 local
#  DBI                      0.7      2017-06-18 CRAN (R 3.3.3)
#  devtools               * 1.13.3   2017-08-02 CRAN (R 3.3.3)
#  digest                   0.6.12   2017-01-27 CRAN (R 3.3.1)
#  doRNG                    1.6.6    2017-04-10 CRAN (R 3.3.1)
#  EnsDb.Hsapiens.v75     * 2.1.0    2017-07-18 Bioconductor
#  ensembldb              * 1.6.2    2016-11-18 Bioconductor
#  foreach                * 1.4.3    2015-10-13 CRAN (R 3.3.0)
#  GenomeInfoDb           * 1.10.3   2017-02-14 Bioconductor
#  GenomicAlignments        1.10.1   2017-03-21 Bioconductor
#  GenomicFeatures        * 1.26.4   2017-06-05 Bioconductor
#  GenomicRanges          * 1.26.4   2017-06-05 Bioconductor
#  ggplot2                * 2.2.1    2016-12-30 CRAN (R 3.3.1)
#  graphics               * 3.3.3    2017-05-18 local
#  grDevices              * 3.3.3    2017-05-18 local
#  grid                     3.3.3    2017-05-18 local
#  gtable                   0.2.0    2016-02-26 CRAN (R 3.3.0)
#  gtools                   3.5.0    2015-05-29 CRAN (R 3.3.0)
#  htmltools                0.3.6    2017-04-28 CRAN (R 3.3.1)
#  httpuv                   1.3.5    2017-07-04 CRAN (R 3.3.3)
#  httr                     1.3.1    2017-08-20 CRAN (R 3.3.3)
#  interactiveDisplayBase   1.12.0   2016-10-20 Bioconductor
#  IRanges                * 2.8.2    2017-06-05 Bioconductor
#  iterators              * 1.0.8    2015-10-13 CRAN (R 3.3.0)
#  lattice                  0.20-34  2016-09-06 CRAN (R 3.3.3)
#  lazyeval                 0.2.0    2016-06-12 CRAN (R 3.3.1)
#  limma                  * 3.30.13  2017-03-21 Bioconductor
#  locfit                 * 1.5-9.1  2013-04-20 CRAN (R 3.3.0)
#  magrittr                 1.5      2014-11-22 CRAN (R 3.3.1)
#  Matrix                   1.2-8    2017-01-20 CRAN (R 3.3.3)
#  matrixStats              0.52.2   2017-04-14 CRAN (R 3.3.1)
#  memoise                  1.1.0    2017-04-21 CRAN (R 3.3.1)
#  methods                * 3.3.3    2017-05-18 local
#  mime                     0.5      2016-07-07 cran (@0.5)
#  munsell                  0.4.3    2016-02-13 CRAN (R 3.3.0)
#  parallel               * 3.3.3    2017-05-18 local
#  permute                  0.9-4    2016-09-09 CRAN (R 3.3.1)
#  pkgconfig                2.0.1    2017-03-21 CRAN (R 3.3.3)
#  pkgmaker                 0.22     2014-05-14 CRAN (R 3.3.0)
#  plyr                     1.8.4    2016-06-08 CRAN (R 3.3.1)
#  R.methodsS3              1.7.1    2016-02-16 CRAN (R 3.3.0)
#  R.oo                     1.21.0   2016-11-01 CRAN (R 3.3.1)
#  R.utils                  2.5.0    2016-11-07 CRAN (R 3.3.1)
#  R6                       2.2.2    2017-06-17 CRAN (R 3.3.3)
#  RColorBrewer           * 1.1-2    2014-12-07 CRAN (R 3.3.0)
#  Rcpp                     0.12.12  2017-07-15 CRAN (R 3.3.3)
#  RCurl                    1.95-4.8 2016-03-01 CRAN (R 3.3.0)
#  registry                 0.3      2015-07-08 CRAN (R 3.3.0)
#  rlang                    0.1.2    2017-08-09 CRAN (R 3.3.3)
#  rngtools                 1.2.4    2014-03-06 CRAN (R 3.3.0)
#  Rsamtools                1.26.2   2017-05-10 Bioconductor
#  RSQLite                  2.0      2017-06-19 CRAN (R 3.3.3)
#  rtracklayer              1.34.2   2017-06-05 Bioconductor
#  S4Vectors              * 0.12.2   2017-06-05 Bioconductor
#  scales                   0.5.0    2017-08-24 CRAN (R 3.3.3)
#  shiny                    1.0.5    2017-08-23 CRAN (R 3.3.3)
#  stats                  * 3.3.3    2017-05-18 local
#  stats4                 * 3.3.3    2017-05-18 local
#  stringi                  1.1.5    2017-04-07 CRAN (R 3.3.1)
#  stringr                  1.2.0    2017-02-18 CRAN (R 3.3.1)
#  SummarizedExperiment   * 1.4.0    2017-06-05 Bioconductor
#  tibble                   1.3.3    2017-05-28 CRAN (R 3.3.3)
#  tools                    3.3.3    2017-05-18 local
#  utils                  * 3.3.3    2017-05-18 local
#  withr                    2.0.0    2017-07-28 CRAN (R 3.3.3)
#  XML                      3.98-1.9 2017-06-19 CRAN (R 3.3.3)
#  xtable                   1.8-2    2016-02-05 CRAN (R 3.3.0)
#  XVector                  0.14.1   2017-03-21 Bioconductor
#  yaml                     2.1.14   2016-11-12 CRAN (R 3.3.1)
#  zlibbioc                 1.20.0   2016-10-20 Bioconductor
