library(GenomicFeatures)
library(ggplot2)
library(GenomicRanges)
library(bumphunter)
library(data.table)
library(clusterProfiler)
require(org.Hs.eg.db)

load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/DMVs_hypoDMRs_methylSeekR.rda")
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/DMVs_hypoDMRs_methylSeekR_prenatal.rda")
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/DMR/DMR_objects.rda")

rm(DMVs,DMVs_pren)

dmvs = c(DMVs.100kb,DMVs_100kb.pren)
