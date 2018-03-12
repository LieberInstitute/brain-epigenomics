library(GenomicRanges)
library(PWMEnrich)
library(PWMEnrich.Hsapiens.background)
library(ggplot2)
library(SummarizedExperiment)


load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/PWMEnrich/TSS_TFs_bysample.rda")
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/PWMEnrich/filtered_TF_list_expressionInfo.rda")
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/brainseq_pipeline/polyA_unstranded/rse_gene_polyA_dlpfc_n41.Rdata")
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/UMRs_LMRs_annotated.rda")

# get TPM

names(assays(rse_gene))


## construct matrix of TF-TF regulatory networks for each sample

pmat = vector("list",length(TSS_TFs_bysample))
for (i in 1:length(TSS_TFs_bysample)) {
  pmat[[i]] = list(vector("list", length(TSS_TFs_bysample[[i]]$sequences)))
  for (j in 1:length(TSS_TFs_bysample[[i]]$sequences)) {
    pmat[[i]][[j]] = sequenceReport(TSS_TFs_bysample[[i]], seq.id=j)
  }
}
pmat = lapply(pmat, function(x) lapply(x, function(y) as.data.frame(y)[order(y$id),]))
pmat = lapply(pmat, function(t) lapply(t, function(x) x[which(x$target %in% names(targettogeneID)),]))

pvalMat = lapply(pmat, function(y) lapply(y, do.call(cbind, lapply(y, p.adjust(y$p.value, method = "fdr")))))
rownames(pvalMat) = pmat[[1]][[1]]$target

for (i in 1:length(TSS_TFs_bysample)) { colnames(pvalMat[[i]]) = names(TSS_TFs_bysample[[i]]$sequences) }
gr = lapply(TSS_TFs_bysample, function(x) GRanges(names(x$sequences)))
elementNROWS(gr)

TSSs = makeGRangesFromDataFrame(data.frame(seqnames = geneMap$Chr, start = geneMap$Start-5000, end = geneMap$Start+5000, strand = geneMap$Strand))
names(TSSs) = geneMap$gencodeID
tfTSS = TSSs[which(names(TSSs) %in% targettogeneID)]

oo = lapply(gr, function(x) findOverlaps(tfTSS, x)) 
splGR = lapply(gr, function(x) split(x[subjectHits(oo)], queryHits(oo)))



network = TSS_TFs_bysample



all_int_seqRep = list()
for (i in 1:length(all_int$sequences)) {
  all_int_seqRep[[i]] = sequenceReport(all_int, seq.id=i)
}
targ = all_int_seqRep[[1]]$target    
pvalMat = do.call(cbind, lapply(all_int_seqRep, function(x) x$p.value[match(targ, x$target)]))
rownames(pvalMat) = targ
colnames(pvalMat) = names(all_int$sequences)

lMat = -log10(pvalMat)


# Figures 5D (construction of TF-TF regulatory networks) TF A was predicted to regulate TF B when: 
# (1) TF A was expressed in a cell type (≥30 TPM), 
# (2) TF A had a predicted footprint (FP A) in a cell type-specific ATAC-seq peak, 
# (3) the ATAC-seq peak was within 10 kb of the TSS for TF B, and 
# (4) TF B was expressed in that cell type (≥30 TPM). 
# The resulting set of predicted regulatory interactions was visualized as a network (igraph package in R), 
# omitting TFs with more than 20 connections to ease visualization. To define a pan-neuronal regulatory network, 
# we identified footprints common to all three cell types that occurred in shared ATAC-seq peaks and did not overlap ubiquitous DNaseI peaks 
# (peaks occurring in at least 40 out of 53 processed DNaseI-seq samples). The full networks are listed in Table S4.



# other tools
# GeneNetworkBuilder: Appliation for discovering direct or indirect targets of transcription factors using ChIP-chip or ChIP-seq, and microarray or RNA-seq gene expression data. Inputting a list of genes of potential targets of one TF from ChIP-chip or ChIP-seq, and the gene expression results, GeneNetworkBuilder generates a regulatory network of the TF.
# http://bioconductor.org/packages/release/bioc/html/GeneNetworkBuilder.html

# RTN: reconstruction of transcriptional networks and analysis of master regulators.
# Circos plots
# HiveR and Gephi plots, http://www.vesnam.com/Rblog/viznets3/
# mfinder: Network motifs detection tool