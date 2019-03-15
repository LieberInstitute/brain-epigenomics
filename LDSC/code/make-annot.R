# Make .annot file for use with ldsc containing brain-derived genomic features
# Code adapted from Peter Hickey (https://github.com/hansenlab/BrainEpigenomeNN/blob/master/SLDSR/scripts/make-annot.R)

### ============================================================================
### NOTEs
###

# qrsh -l bluejay,mem_free=100G,h_vmem=100G,h_fsize=100G

# - Adapted from tutorial at https://github.com/bulik/ldsc/wiki/LD-Score-Estimation-Tutorial#building-on-top-of-the-finucane-et-al-baseline-model

### ============================================================================
### Setup
###

library(readr)
library(GenomicRanges)
library(R.utils)
library(rtracklayer)
library(AnnotationHub)
library(bumphunter)

options("mc.cores" = 40)

seqlevels <- 1:22

### ============================================================================
### Load and prepare data
###

# ------------------------------------------------------------------------------
# Not separated by sign
#

load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/bumphunting/rda/limma_Neuron_CpGs_minCov_3_ageInfo_dmrs.Rdata")
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/DMR/DMR_objects.rda")
load('/dcl01/lieber/ajaffe/lab/brain-epigenomics/bumphunting/BSobj_bsseqSmooth_Neuron_minCov_3.Rdata')

# Identify all CpG clusters in the genome

gr = granges(BSobj)
cl = clusterMaker( chr = as.character(seqnames(gr)), pos = start(gr),   maxGap = 1000)
gr.clusters = split(gr, cl)
gr.clusters = unlist(range(gr.clusters))
df.clusters = as.data.frame(gr.clusters)
df.clusters$regionID = paste0(df.clusters$seqnames,":",df.clusters$start,"-",df.clusters$end)
df.clusters$rnum = 1:length(gr.clusters)

DMRs <- c(list(
  "Celltype_HypoNeuron" = unique(makeGRangesFromDataFrame(DMR$CellType[which(DMR$CellType$value<0 & DMR$CellType$sig=="FWER < 0.05"),])),
  "Celltype_HypoGlia" = unique(makeGRangesFromDataFrame(DMR$CellType[which(DMR$CellType$value>0 & DMR$CellType$sig=="FWER < 0.05"),])),
  "Age_Decreasing" = unique(makeGRangesFromDataFrame(DMR$Age[which(DMR$Age$value<0 & DMR$Age$sig=="FWER < 0.05"),])), 
  "Age_Increasing" = unique(makeGRangesFromDataFrame(DMR$Age[which(DMR$Age$value>0 & DMR$Age$sig=="FWER < 0.05"),]))),
  as.list(split(dmrs, dmrs$k6cluster_label)))
names(DMRs)[5:10] = c("Gr1_cdDMRs", "Gr2_cdDMRs", "Gr3_cdDMRs", "Gr4_cdDMRs", "Gr5_cdDMRs", "Gr6_cdDMRs")
  
non_DMRs <- subsetByOverlaps(gr.clusters, unlist(as(DMRs, "GRangesList")), invert = TRUE)

set.seed(3000)

nonDMR_subsets <- mclapply(
  list(Celltype_HypoNeuron_nonDMRs = DMRs$Celltype_HypoNeuron, Celltype_HypoGlia_nonDMRs = DMRs$Celltype_HypoGlia, 
       Age_Decreasing_nonDMRs = DMRs$Age_Decreasing, Age_Increasing_nonDMR = DMRs$Age_Increasing,
       Gr1_nonDMRs = DMRs$Gr1_cdDMRs, Gr2_nonDMRs = DMRs$Gr2_cdDMRs, Gr3_nonDMRs = DMRs$Gr3_cdDMRs, 
       Gr4_nonDMRs = DMRs$Gr4_cdDMRs, Gr5_nonDMRs = DMRs$Gr5_cdDMRs, Gr6_nonDMRs = DMRs$Gr6_cdDMRs),
  function(drs) {
    # Keep sampling until get non-DARs with about the same total width as DARs
    tmp <- sample(non_DMRs, length(drs))
    while (sum(width(tmp)) < sum(width(drs))) {
      message(sum(width(drs)) - sum(width(tmp)))
      tmp <- c(tmp, sample(
        subsetByOverlaps(non_DMRs, tmp, invert = TRUE, type = "equal"),
        1000))
    }
    list("sameN" = sample(non_DMRs, length(drs)),
         "sameWidth" = tmp)
  }, mc.cores = getOption("mc.cores"))

our_differential_features <- c(DMRs, unlist(nonDMR_subsets, recursive = FALSE), list(non_DMRs = non_DMRs))
# NOTE: Rename to avoid '.' in names
names(our_differential_features) <- gsub("\\.", "_", names(our_differential_features))


# ------------------------------------------------------------------------------
# Non-differential features
#

load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/UMRs_LMRs_methylSeekR.rda")
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/UMRs_LMRs_methylSeekR_prenatal.rda")

rm("UMRLMRsegments.CG", "UMRLMRsegments.CGpren")


# Non-differential low-methylated regions (LMRs)

LMRs <- list(Prenatal_LMRs = lapply(UMRLMRsegments.CG.100kb.pren, function(y) y[y$type=="LMR"]),
             Neuronal_LMRs = lapply(UMRLMRsegments.CG.100kb[which(names(UMRLMRsegments.CG.100kb) %in% pd[which(pd$Cell.Type=="Neuron"),"Data.ID"])],
                                    function(y) y[y$type=="LMR"]),
             Glial_LMRs = lapply(UMRLMRsegments.CG.100kb[which(names(UMRLMRsegments.CG.100kb) %in% pd[which(pd$Cell.Type=="Glia"),"Data.ID"])],
                                 function(y) y[y$type=="LMR"]))
LMRs <- lapply(LMRs, function(x) reduce(unlist(as(x, "GRangesList"))))

non_dLMRs <- lapply(LMRs, function(gr) subsetByOverlaps(gr, unlist(as(DMRs, "GRangesList")), invert = TRUE))


# Subset of chromHMM track with regulatory region-like states

ah <- AnnotationHub()

# Ten brain samples with chromHMM from RoadMap; identified using
# query(ah, c("chromHMM", "brain"))

list_of_chromHMM <- ah[c("AH46920", "AH46921", "AH46922", "AH46923", "AH46924",
                         "AH46925", "AH46926", "AH46927", "AH46934", "AH46935")]
chromHMM_regulatory_regions <- lapply(list_of_chromHMM, function(x) {
  x <- x[[1]]
  x[x$name %in%
      c("Bivalent Enhancer", "Bivalent/Poised TSS", "Genic enhancers",
        "Flanking Active TSS", "Active TSS", "Strong transcription",
        "Enhancers")]
})

chromHMM_regulatory_regions <- list(
  "chromHMM_union" = reduce(unlist(GRangesList(chromHMM_regulatory_regions))))

non_differential_features <- c(non_dLMRs,
                               chromHMM_regulatory_regions)


# ------------------------------------------------------------------------------
# CNS (central nervous system) cell type group
#

# NOTE: Annotation files for CNS were made available by the LDSC authors.
#       However, the BED files were not, so I construct an approximation to the
#       BED files based on the annotation files
CNS_annot <- mclapply(seqlevels, function(sl) {
  read_tsv(file.path("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/ldsc/Phase1", "cell_type_groups",
                     paste0("CNS.", sl, ".annot.gz")))
}, mc.cores = getOption("mc.cores"))
CNS_BED <- unlist(GRangesList(mclapply(CNS_annot, function(x) {
  y <- Rle(x$CNS)
  s <- cumsum(runLength(y))[runValue(y) == 1] -
    runLength(y)[runValue(y) == 1] + 1
  e <- cumsum(runLength(y))[runValue(y) == 1]
  GRanges(seqnames = paste0("chr", x[["CHR"]][1]),
          ranges = IRanges(x[["BP"]][s],
                           x[["BP"]][e]))
})))

# ------------------------------------------------------------------------------
# List of all brain-derived features
#

categories <- c(our_differential_features,
                non_differential_features,
                list("CNS" = CNS_BED))
stopifnot(all(sapply(categories, isDisjoint)))
stopifnot(all(!grepl("\\.", names(categories))))

saveRDS(categories, "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/ldsc/categories.rds")

### ============================================================================
### Make annotations
###

# NOTE: Don't re-make the CNS annotation
k <- grep("CNS", names(categories), invert = TRUE)
mclapply(seqlevels, function(sl) {
  message(sl)
  cds <- read_tsv(
    paste0("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/ldsc/Phase1/cell_type_groups/CNS.", sl,
           ".annot.gz"))
  cds_gr <- GRanges(paste0("chr", cds$CHR), IRanges(cds$BP, width = 1L))
  annot <- cds[, c("CHR", "BP", "SNP", "CM")]
  annot[names(categories)[k]] <- mclapply(names(categories), function(cn) {
    stopifnot(isDisjoint(categories[[cn]]))
    as.integer(overlapsAny(cds_gr, categories[[cn]]))
  }, mc.cores = 4)
  # 'Marginal' annotation file
  mclapply(names(categories[k]), function(n) {
    fl <- paste0("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/ldsc/annotation/", n, ".Phase1.", sl, ".annot")
    write_tsv(annot[, c("CHR", "BP", "SNP", "CM", n)], fl)
    gzip(fl, overwrite = TRUE)
  }, mc.cores = 4)
}, mc.cores = 10)


# ------------------------------------------------------------------------------
# Annot file for 'base' category
# NOTE: Not added to 'categories' since it is not a brain-derived category
#

mclapply(seqlevels, function(sl) {
  x <- read_tsv(file.path("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/ldsc/", "Phase1", "baseline",
                          paste0("baseline.", sl, ".annot.gz")))
  write_tsv(x[, 1:5], paste0("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/ldsc/annotation/base.Phase1.", sl, ".annot.gz"))
}, mc.cores = getOption("mc.cores"))

### ============================================================================
### Complex set operation features
###

POS_CG_DMRs <- categories[["POS_CG-DMRs"]]
chromHMM <- categories[["chromHMM_union"]]
CNS <- categories[["CNS"]]
H3K27ac <- categories[["H3K27ac"]]

intersect <- GenomicRanges::intersect
union <- GenomicRanges::union
setdiff <- GenomicRanges::setdiff

# Brain only
chromHMM_only <- setdiff(chromHMM, CNS_H3K27ac)
CNS_only <- setdiff(CNS, chromHMM_H3K27ac)
H3K27ac_only <- setdiff(H3K27ac, CNS_chromHMM)

# 2+
two_plus <- setdiff(
  setdiff(
    setdiff(union(CNS, union(H3K27ac, chromHMM)), chromHMM_only),
    CNS_only),
  H3K27ac_only)

# POS_CG-DMRs
CG_only <- setdiff(POS_CG_DMRs, union(chromHMM, union(CNS, H3K27ac)))
CG_shared <- intersect(POS_CG_DMRs, union(chromHMM, union(CNS, H3K27ac)))

# 2+ without CG_shared
two_plus_no_CG <- setdiff(two_plus, CG_shared)

# Brain only without CG_shared
chromHMM_only_no_CG <- setdiff(chromHMM_only, CG_shared)
CNS_only_no_CG <- setdiff(CNS_only, CG_shared)
H3K27ac_only_no_CG <- setdiff(H3K27ac_only, CG_shared)

complex_set_op_features <- list(CG_only = CG_only,
                                CG_shared = CG_shared,
                                chromHMM_only = chromHMM_only,
                                H3K27ac_only = H3K27ac_only,
                                CNS_only = CNS_only,
                                two_plus = two_plus,
                                two_plus_no_CG = two_plus_no_CG,
                                chromHMM_only_no_CG = chromHMM_only_no_CG,
                                CNS_only_no_CG = CNS_only_no_CG,
                                H3K27ac_only_no_CG = H3K27ac_only_no_CG)
saveRDS(complex_set_op_features, "../objects/complex_set_op_features.rds")

# Make annotations
mclapply(seqlevels, function(sl) {
  message(sl)
  cds <- read_tsv(
    paste0("../extdata/Phase1/cell_type_groups/CNS.", sl,
           ".annot.gz"))
  cds_gr <- GRanges(paste0("chr", cds$CHR), IRanges(cds$BP, width = 1L))
  annot <- cds[, c("CHR", "BP", "SNP", "CM")]
  annot[names(complex_set_op_features)] <-
    mclapply(names(complex_set_op_features), function(cn) {
      stopifnot(isDisjoint(complex_set_op_features[[cn]]))
      as.integer(overlapsAny(cds_gr, complex_set_op_features[[cn]]))
    }, mc.cores = 4)
  # 'Marginal' annotation file
  mclapply(names(complex_set_op_features), function(n) {
    fl <- paste0("../output/ldsc/", n, ".Phase1.", sl, ".annot")
    write_tsv(annot[, c("CHR", "BP", "SNP", "CM", n)], fl)
    gzip(fl, overwrite = TRUE)
  }, mc.cores = 4)
}, mc.cores = 10)

### ============================================================================
### Non-differential features without each differential feature
###

differential_features <- categories[c("POS_CG-DMRs", "POSvsNEG_CG-DMRs",
                                      "CH-DMRs",
                                      "POS_DARs", "POSvsNEG_DARs",
                                      "POS_bigDARs", "POSvsNEG_bigDARs")]
non_differential_features <- categories[c("CNS", "H3K27ac", "chromHMM_union")]

n_df <- names(differential_features)
names(n_df) <- n_df
n_ndf <- names(non_differential_features)
names(n_ndf) <- n_ndf
ndf_excluding_df <- unlist(lapply(n_ndf, function(n1) {
  ndf <- non_differential_features[[n1]]
  lapply(n_df, function(n2) {
    df <- differential_features[[n2]]
    setdiff(ndf, df)
  })
}))
names(ndf_excluding_df) <- gsub("\\.", "_excluding_", names(ndf_excluding_df))

saveRDS(
  ndf_excluding_df,
  "../objects/non-differential_features_excluding_differential_features.rds")

# Make annotations
mclapply(seqlevels, function(sl) {
  message(sl)
  cds <- read_tsv(
    paste0("../extdata/Phase1/cell_type_groups/CNS.", sl,
           ".annot.gz"))
  cds_gr <- GRanges(paste0("chr", cds$CHR), IRanges(cds$BP, width = 1L))
  annot <- cds[, c("CHR", "BP", "SNP", "CM")]
  annot[names(ndf_excluding_df)] <-
    mclapply(names(ndf_excluding_df), function(cn) {
      stopifnot(isDisjoint(ndf_excluding_df[[cn]]))
      as.integer(overlapsAny(cds_gr, ndf_excluding_df[[cn]]))
    }, mc.cores = 4)
  # 'Marginal' annotation file
  mclapply(names(ndf_excluding_df), function(n) {
    fl <- paste0("../output/ldsc/", n, ".Phase1.", sl, ".annot")
    write_tsv(annot[, c("CHR", "BP", "SNP", "CM", n)], fl)
    gzip(fl, overwrite = TRUE)
  }, mc.cores = 4)
}, mc.cores = 10)