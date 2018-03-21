library(GenomicRanges)
library(PWMEnrich)
library(PWMEnrich.Hsapiens.background)
library(ggplot2)
library(SummarizedExperiment)
library(jaffelab)

load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/PWMEnrich/TSS_TFs_bysample.rda")
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/PWMEnrich/filtered_TF_list_expressionInfo.rda", verbose = T)
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/brainseq_pipeline/polyA_unstranded/rse_gene_polyA_dlpfc_n41.Rdata", verbose = T)
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/UMRs_LMRs_annotated.rda", verbose = T)


## construct matrix of TF-TF regulatory networks for each sample

group = groupReport(TSS_TFs_bysample[[1]])
group = as.data.frame(group)[order(group$id),]
group = group[which(group$target %in% names(prepostTargettoGeneID)),]

pmat1.10 = pmat11.20 = pmat21.30 =pmat31.40 = pmat41.52 = list()
# 1-10
for (i in 1:10) {
  pmat1.10[[i]] = list(vector("list", length(TSS_TFs_bysample[[i]]$sequences)))
  for (j in 1:length(TSS_TFs_bysample[[i]]$sequences)) {
    pmat1.10[[i]][[j]] = sequenceReport(TSS_TFs_bysample[[i]], seq.id=j)
  }
}
pmat1.10 = pmat1.10[which(elementNROWS(pmat1.10)>0)]
names(pmat1.10) = names(TSS_TFs_bysample)[1:10]
pmat1.10 = lapply(pmat1.10, function(x) lapply(x, function(y) as.data.frame(y)[order(y$id),]))
pmat1.10 = lapply(pmat1.10, function(y) do.call(cbind, lapply(y, function(z) p.adjust(z$p.value, method = "fdr"))))
for (i in 1:length(pmat1.10)) { rownames(pmat1.10[[i]]) = group$target }
cols = list()
for (i in 1:10) { cols[[i]] = names(TSS_TFs_bysample[[i]]$sequences) }
cols = cols[which(elementNROWS(cols)>0)]
for (i in 1:length(pmat1.10)) { colnames(pmat1.10[[i]]) = cols[[i]] }
save(pmat1.10, file="/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/PWMEnrich/pmat1.10.rda")
# 11-20
for (i in 11:20) {
  pmat11.20[[i]] = list(vector("list", length(TSS_TFs_bysample[[i]]$sequences)))
  for (j in 1:length(TSS_TFs_bysample[[i]]$sequences)) {
    pmat11.20[[i]][[j]] = sequenceReport(TSS_TFs_bysample[[i]], seq.id=j)
  }
}
names(pmat11.20) = names(TSS_TFs_bysample)[11:20]
pmat11.20 = lapply(pmat11.20, function(x) lapply(x, function(y) as.data.frame(y)[order(y$id),]))
pmat11.20 = lapply(pmat11.20, function(y) do.call(cbind, lapply(y, function(z) p.adjust(z$p.value, method = "fdr"))))
for (i in 1:length(pmat11.20)) { rownames(pmat11.20[[i]]) = group$target }
for (i in 11:20) { cols[[i]] = names(TSS_TFs_bysample[[i]]$sequences) }
cols = cols[which(elementNROWS(cols)>0)]
for (i in 1:length(pmat11.20)) { colnames(pmat11.20[[i]]) = cols[[i]] }
save(pmat11.20, file="/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/PWMEnrich/pmat11.20.rda")
# 21-30
for (i in 21:30) {
  pmat21.30[[i]] = list(vector("list", length(TSS_TFs_bysample[[i]]$sequences)))
  for (j in 1:length(TSS_TFs_bysample[[i]]$sequences)) {
    pmat21.30[[i]][[j]] = sequenceReport(TSS_TFs_bysample[[i]], seq.id=j)
  }
}
pmat21.30 = pmat21.30[which(elementNROWS(pmat21.30)>0)]
names(pmat21.30) = names(TSS_TFs_bysample)[21:30]
pmat21.30 = lapply(pmat21.30, function(x) lapply(x, function(y) as.data.frame(y)[order(y$id),]))
pmat21.30 = lapply(pmat21.30, function(y) do.call(cbind, lapply(y, function(z) p.adjust(z$p.value, method = "fdr"))))
for (i in 1:length(pmat21.30)) { rownames(pmat21.30[[i]]) = group$target }
cols = list()
for (i in 21:30) { cols[[i]] = names(TSS_TFs_bysample[[i]]$sequences) }
cols = cols[which(elementNROWS(cols)>0)]
for (i in 1:length(pmat21.30)) { colnames(pmat21.30[[i]]) = cols[[i]] }
save(pmat21.30, file="/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/PWMEnrich/pmat21.30.rda")
# 31-40
for (i in 31:40) {
  pmat31.40[[i]] = list(vector("list", length(TSS_TFs_bysample[[i]]$sequences)))
  for (j in 1:length(TSS_TFs_bysample[[i]]$sequences)) {
    pmat31.40[[i]][[j]] = sequenceReport(TSS_TFs_bysample[[i]], seq.id=j)
  }
}
pmat31.40 = pmat31.40[which(elementNROWS(pmat31.40)>0)]
names(pmat31.40) = names(TSS_TFs_bysample)[31:40]
pmat31.40 = lapply(pmat31.40, function(x) lapply(x, function(y) as.data.frame(y)[order(y$id),]))
pmat31.40 = lapply(pmat31.40, function(y) do.call(cbind, lapply(y, function(z) p.adjust(z$p.value, method = "fdr"))))
for (i in 1:length(pmat31.40)) { rownames(pmat31.40[[i]]) = group$target }
cols = list()
for (i in 31:40) { cols[[i]] = names(TSS_TFs_bysample[[i]]$sequences) }
cols = cols[which(elementNROWS(cols)>0)]
for (i in 1:length(pmat31.40)) { colnames(pmat31.40[[i]]) = cols[[i]] }
save(pmat31.40, file="/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/PWMEnrich/pmat31.40.rda")
# 41-52
for (i in 41:52) {
  pmat41.52[[i]] = list(vector("list", length(TSS_TFs_bysample[[i]]$sequences)))
  for (j in 1:length(TSS_TFs_bysample[[i]]$sequences)) {
    pmat41.52[[i]][[j]] = sequenceReport(TSS_TFs_bysample[[i]], seq.id=j)
  }
}
pmat41.52 = pmat41.52[which(elementNROWS(pmat41.52)>0)]
names(pmat41.52) = names(TSS_TFs_bysample)[41:52]
pmat41.52 = lapply(pmat41.52, function(x) lapply(x, function(y) as.data.frame(y)[order(y$id),]))
pmat41.52 = lapply(pmat41.52, function(y) do.call(cbind, lapply(y, function(z) p.adjust(z$p.value, method = "fdr"))))
for (i in 1:length(pmat41.52)) { rownames(pmat41.52[[i]]) = group$target }
cols = list()
for (i in 41:52) { cols[[i]] = names(TSS_TFs_bysample[[i]]$sequences) }
cols = cols[which(elementNROWS(cols)>0)]
for (i in 1:length(pmat41.52)) { colnames(pmat41.52[[i]]) = cols[[i]] }
save(pmat41.52, file="/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/PWMEnrich/pmat41.52.rda")

pmat = c(pmat1.10, pmat11.20, pmat21.30, pmat31.40, pmat41.52)

gr = lapply(TSS_TFs_bysample, function(x) GRanges(names(x$sequences)))
elementNROWS(gr)

save(pmat,gr, file="/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/PWMEnrich/pmat.rda")

#####
## label regions with their appropriate gene TSS
gr = GRangesList(gr)
grFull = unlist(gr)
grFull$SampleID = rep(names(gr), times = lengths(gr))
grFull$chrpos = paste0(seqnames(grFull), ":", start(grFull), "-", end(grFull))
grFull$FullID = paste0(grFull$SampleID, "_", grFull$chrpos)

for(i in seq(along=pmat)) {
	colnames(pmat[[i]]) = paste0(names(pmat)[i], "_", colnames(pmat[[i]]))
}

## filter to TFs
geneMap= rowRanges(rse_gene)

# promoters, -5kb to +5kb
genePromoters = GRanges(seqnames(geneMap),
	IRanges(start = ifelse(strand(geneMap) == "+",
		start(geneMap)-5000, end(geneMap)-5000),
	end = ifelse(strand(geneMap) == "+",
		start(geneMap)+5000, end(geneMap)+5000)),
		strand = strand(geneMap))
mcols(genePromoters) = mcols(geneMap)
names(genePromoters) = names(geneMap)
tfTSS = genePromoters[which(names(genePromoters) %in% prepostTargettoGeneID)]
tfTSS$otherSymbol = names(prepostTargettoGeneID)[match(tfTSS$gencodeID, prepostTargettoGeneID)]

## get symbol
oo1 = findOverlaps(grFull, tfTSS)
grFull$Symbol = NA
grFull$Symbol[queryHits(oo1)] = tfTSS$otherSymbol[subjectHits(oo1)]
grFull$gencodeID = NA
grFull$gencodeID[queryHits(oo1)] = tfTSS$gencodeID[subjectHits(oo1)]

## find overlaps
oo = findOverlaps(tfTSS, grFull)

## filter
tfTSS_sub = tfTSS[unique(queryHits(oo)),]
idUniq = unique(grFull$FullID[subjectHits(oo)])
pmat_sub = lapply(pmat, function(x) x[,colnames(x) %in% idUniq])
pmat_sub = do.call("rbind", lapply(pmat_sub,t))
grFull_sub = grFull[grFull$FullID %in% idUniq,]
pmat_sub = pmat_sub[grFull_sub$FullID,]

#### 
sIndexes = splitit(grFull_sub$SampleID)
adjMatList = vector("list", length(sIndexes))
names(adjMatList) = names(sIndexes)
for(i in seq(along=sIndexes)) {
	adjMat = matrix(FALSE, nr = length(tfTSS), ncol = length(tfTSS),
		dimnames = list(tfTSS$otherSymbol, tfTSS$otherSymbol))
	ii = sIndexes[[i]]
	g = grFull_sub[ii]
	p = pmat_sub[ii,] < 0.01
	pp = sapply(splitit(colnames(p)), function(jj) rowSums(t(t(p[,jj]))) > 0)
	pp = pp[,colnames(adjMat)]
	pp = t(sapply(splitit(g$Symbol), function(jj) {
		if(length(jj) > 1) colSums(pp[jj,]) > 0 else pp[jj,] > 0
	}))
	adjMat[rownames(pp), colnames(pp)] = pp
	adjMatList[[i]] = adjMat
}
save(adjMatList, file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/adjMat_List_TFs_LMRs.rda")


### Construct networks for DMRs

load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/DMR/DMR_objects.rda")
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/PWMEnrich/DMR_PWMEnrich_objects.rda")


DMRpmat = vector("list", 5)
names(DMRpmat) = c("")
for (i in 1:5) {
  DMRpmat[[i]] = list(vector("list", length(TSS_TFs_bysample[[i]]$sequences)))
  for (j in 1:length(TSS_TFs_bysample[[i]]$sequences)) {
    pmat41.52[[i]][[j]] = sequenceReport(TSS_TFs_bysample[[i]], seq.id=j)
  }
}
pmat41.52 = pmat41.52[which(elementNROWS(pmat41.52)>0)]
names(pmat41.52) = names(TSS_TFs_bysample)[41:52]
pmat41.52 = lapply(pmat41.52, function(x) lapply(x, function(y) as.data.frame(y)[order(y$id),]))
pmat41.52 = lapply(pmat41.52, function(t) lapply(t, function(x) x[which(x$target %in% names(prepostTargettoGeneID)),]))
pmat41.52 = lapply(pmat41.52, function(y) do.call(cbind, lapply(y, function(z) p.adjust(z$p.value, method = "fdr"))))
for (i in 1:length(pmat41.52)) { rownames(pmat41.52[[i]]) = group$target }
cols = list()
for (i in 41:52) { cols[[i]] = names(TSS_TFs_bysample[[i]]$sequences) }
cols = cols[which(elementNROWS(cols)>0)]
for (i in 1:length(pmat41.52)) { colnames(pmat41.52[[i]]) = cols[[i]] }
save(pmat41.52, file="/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/PWMEnrich/pmat41.52.rda")

ct.age = list()
for (i in 1:length(all_split[-grep("Interaction", names(all_split))])) {
  ct.age[[i]] = list(vector("list", length(all_split[-grep("Interaction", names(all_split))][[i]]$sequences)))
  for (j in 1:length(all_split[-grep("Interaction", names(all_split))][[i]]$sequences)) {
    ct.age[[i]][[j]] = sequenceReport(all_split[-grep("Interaction", names(all_split))][[i]], seq.id=j)
  }
}
names(ct.age) = names(all_split[-grep("Interaction", names(all_split))])
ct.age = lapply(ct.age, function(x) lapply(x, function(y) as.data.frame(y)[order(y$id),]))
ct.age = lapply(ct.age, function(t) lapply(t, function(x) x[which(x$target %in% names(targettogeneID)),]))
ct.age.Mat = lapply(ct.age, function(y) do.call(cbind, lapply(y, function(z) p.adjust(z$p.value, method = "fdr"))))
for (i in 1:length(ct.age.Mat)) { rownames(ct.age.Mat[[i]]) = ct.age[[1]][[1]]$target }
for (i in 1:length(ct.age.Mat)) { 
  colnames(ct.age.Mat[[i]]) = names(all_split[-grep("Interaction", names(all_split))][[i]]$sequences) }
ct.age.gr = lapply(all_split[-grep("Interaction", names(all_split))], function(x) GRanges(names(x$sequences)))
elementNROWS(ct.age.gr)


gr = lapply(TSS_TFs_bysample, function(x) GRanges(names(x$sequences)))
elementNROWS(gr)

save(pmat,gr, file="/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/PWMEnrich/pmat.rda")















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