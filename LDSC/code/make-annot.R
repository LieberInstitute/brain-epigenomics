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
  Celltype = unique(makeGRangesFromDataFrame(DMR$CellType[which(DMR$CellType$sig=="FWER < 0.05"),])),
  Age = unique(makeGRangesFromDataFrame(DMR$Age[which(DMR$Age$sig=="FWER < 0.05"),])),
  cdDMRs = dmrs,
  "Celltype_HypoNeuron" = unique(makeGRangesFromDataFrame(DMR$CellType[which(DMR$CellType$value<0 & DMR$CellType$sig=="FWER < 0.05"),])),
  "Celltype_HypoGlia" = unique(makeGRangesFromDataFrame(DMR$CellType[which(DMR$CellType$value>0 & DMR$CellType$sig=="FWER < 0.05"),])),
  "Age_Decreasing" = unique(makeGRangesFromDataFrame(DMR$Age[which(DMR$Age$value<0 & DMR$Age$sig=="FWER < 0.05"),])), 
  "Age_Increasing" = unique(makeGRangesFromDataFrame(DMR$Age[which(DMR$Age$value>0 & DMR$Age$sig=="FWER < 0.05"),]))),
  as.list(split(dmrs, dmrs$k6cluster_label)))
names(DMRs)[8:13] = c("Gr1_cdDMRs", "Gr2_cdDMRs", "Gr3_cdDMRs", "Gr4_cdDMRs", "Gr5_cdDMRs", "Gr6_cdDMRs")
  
non_DMRs <- subsetByOverlaps(gr.clusters, unlist(as(DMRs, "GRangesList")), invert = TRUE)

set.seed(3000)

nonDMR_subsets <- mclapply(
  list(Celltype_nonDMRs = DMRs$Celltype, Age_nonDMRs = DMRs$Age, cdDMRs_nonDMRs = DMRs$cdDMRs, 
       Celltype_HypoNeuron_nonDMRs = DMRs$Celltype_HypoNeuron, Celltype_HypoGlia_nonDMRs = DMRs$Celltype_HypoGlia, 
       Age_Decreasing_nonDMRs = DMRs$Age_Decreasing, Age_Increasing_nonDMR = DMRs$Age_Increasing,
       Gr1_nonDMRs = DMRs$Gr1_cdDMRs, Gr2_nonDMRs = DMRs$Gr2_cdDMRs, Gr3_nonDMRs = DMRs$Gr3_cdDMRs, 
       Gr4_nonDMRs = DMRs$Gr4_cdDMRs, Gr5_nonDMRs = DMRs$Gr5_cdDMRs, Gr6_nonDMRs = DMRs$Gr6_cdDMRs),
  function(drs) {
    # Keep sampling until get non-DMRs with about the same total width as DMRs
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
names(LMRs) = c("Prenatal_LMRs_all", "Neuronal_LMRs_all", "Glial_LMRs_all")


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

chromHMM_union_reduced = list(chromHMM_union_reduced = subsetByOverlaps(chromHMM_regulatory_regions[[1]], unlist(as(DMRs, "GRangesList")), invert = TRUE))

non_differential_features <- c(LMRs, non_dLMRs,
                               chromHMM_regulatory_regions, chromHMM_union_reduced)


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

CNS_BED_red = list(CNS_reduced = subsetByOverlaps(CNS_BED, unlist(as(DMRs, "GRangesList")), invert = TRUE))

# ------------------------------------------------------------------------------
# List of all brain-derived features
#

categories <- c(our_differential_features,
                non_differential_features,
                list("CNS" = CNS_BED), CNS_BED_red)
stopifnot(all(sapply(categories, isDisjoint)))
stopifnot(all(!grepl("\\.", names(categories))))

saveRDS(categories, "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/ldsc/categories.rds")
write.table(data.frame(names(categories)), quote=F, row.names=F,col.names=F, file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/LDSC/categories.txt")

### ============================================================================
### Make annotations
###

# NOTE: Don't re-make the CNS annotation
k <- which(names(categories)!="CNS")
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
