library('bsseq')
library('GenomicRanges')
library('bumphunter')

load('/dcl01/lieber/ajaffe/lab/brain-epigenomics/bumphunting/BSobj_bsseqSmooth_Neuron_minCov_3.Rdata')

## Answer "how many Mb of the genome was covered by the range of the clusters for bumphunting?"
clus <- clusterMaker(seqnames(rowRanges(BSobj)), start(rowRanges(BSobj)), maxGap = 1000)
clus_gr <- split(rowRanges(BSobj), clus)
sum(unlist(width(range(clus_gr)))) / 1e6


models <- c('interaction', 'cell', 'age')
dmr_glia_mean_meth <- lapply(models, function(model) {
    load(paste0('/dcl01/lieber/ajaffe/lab/brain-epigenomics/bumphunting/bumps_bsseqSmooth_Neuron_', model, '_250_perm.Rdata'))

    ## Build a GRanges object from the table output of bumphunter
    bumps_gr <- GRanges(bumps$table)

    ## Next find which bases are part of the DMRs
    idx <- countOverlaps(rowRanges(BSobj), bumps_gr) > 0

    ## Get the methylation for the glia samples
    ## only for the bases that are in the DMRs
    meth <- getMeth(BSobj[idx, colData(BSobj)$Cell.Type == 'Glia'], type = 'raw')

    ov <- findOverlaps(bumps_gr, rowRanges(BSobj[idx, colData(BSobj)$Cell.Type == 'Glia']))
    ov_dmr <- split(subjectHits(ov), queryHits(ov))

    ## Finally compute the mean methylation for the bases in the DMR for the glia samples
    dmr_mean_meth <- sapply(ov_dmr, function(dmr_bases) {
        mean(meth[dmr_bases, ])
    })

    stopifnot(length(bumps_gr) == length(dmr_mean_meth))

    ## Separate them into the 3 groups
    dmr_mean_meth_bins <- cut(dmr_mean_meth, c(0, 0.2, 0.8, Inf), include.lowest = TRUE)
    levels(dmr_mean_meth_bins)
    result <- list('mean_meth' = dmr_mean_meth, 'mean_meth_bin' = dmr_mean_meth_bins)
    return(result)
})
names(dmr_glia_mean_meth) <- models
save(dmr_glia_mean_meth, file = 'rdas/dmr_glia_mean_meth.Rdata')


# > sessionInfo()
# R version 3.3.3 Patched (2017-03-15 r72696)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Red Hat Enterprise Linux Server release 6.9 (Santiago)
#
# locale:
#  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C
#  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8
#  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8
#  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C
#  [9] LC_ADDRESS=C               LC_TELEPHONE=C
# [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C
#
# attached base packages:
# [1] stats4    parallel  stats     graphics  grDevices datasets  utils
# [8] methods   base
#
# other attached packages:
#  [1] bumphunter_1.14.0          locfit_1.5-9.1
#  [3] iterators_1.0.8            foreach_1.4.3
#  [5] bsseq_1.10.0               limma_3.30.13
#  [7] SummarizedExperiment_1.4.0 Biobase_2.34.0
#  [9] GenomicRanges_1.26.4       GenomeInfoDb_1.10.3
# [11] IRanges_2.8.2              S4Vectors_0.12.2
# [13] BiocGenerics_0.20.0        colorout_1.1-2
#
# loaded via a namespace (and not attached):
#  [1] Rcpp_0.12.12             plyr_1.8.4               XVector_0.14.1
#  [4] GenomicFeatures_1.26.4   rngtools_1.2.4           R.methodsS3_1.7.1
#  [7] bitops_1.0-6             R.utils_2.5.0            tools_3.3.3
# [10] zlibbioc_1.20.0          biomaRt_2.30.0           bit_1.1-12
# [13] digest_0.6.12            memoise_1.1.0            tibble_1.3.3
# [16] RSQLite_2.0              lattice_0.20-34          rlang_0.1.2
# [19] doRNG_1.6.6              Matrix_1.2-8             DBI_0.7
# [22] registry_0.3             rtracklayer_1.34.2       pkgmaker_0.22
# [25] stringr_1.2.0            Biostrings_2.42.1        gtools_3.5.0
# [28] bit64_0.9-7              grid_3.3.3               data.table_1.10.4
# [31] AnnotationDbi_1.36.2     BiocParallel_1.8.2       XML_3.98-1.9
# [34] blob_1.1.0               magrittr_1.5             GenomicAlignments_1.10.1
# [37] Rsamtools_1.26.2         scales_0.5.0             codetools_0.2-15
# [40] matrixStats_0.52.2       permute_0.9-4            colorspace_1.3-2
# [43] xtable_1.8-2             stringi_1.1.5            RCurl_1.95-4.8
# [46] munsell_0.4.3            R.oo_1.21.0
