####

library("limma")
library("edgeR")
library(bsseq)
library(recount)
library(jaffelab)

## load CpG data
load('/dcl01/lieber/ajaffe/lab/brain-epigenomics/bumphunting/BSobj_bsseqSmooth_Neuron_minCov_3.Rdata')
BSobj_cpg = BSobj
pd <- pData(BSobj_cpg)
gr_cpg <- granges(BSobj_cpg)
meth_cpg <- getMeth(BSobj_cpg, type = 'smooth')
rm(BSobj)

## load non-CpG data
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/bsseq/bsobj_by_chr/allChrs_postNatal_cleaned_nonCG_noHomogenate_highCov.Rdata")
BSobj_non = BSobj
gr_non <- granges(BSobj_non)
meth_non <- getMeth(BSobj_non, type = 'raw')
rm(BSobj)

### get mean summaries
pd$newAgeGroup = cut(pd$Age, c(0,1,10,30),
	label = c("Infant", "Child", "Teen"))
gIndexes = splitit(paste0(pd$Cell.Type, ":", pd$newAgeGroup))

mean_meth_cpg = sapply(gIndexes, function(ii) {
	if(length(ii) > 1) rowMeans(meth_cpg[,ii]) else meth_cpg[,ii]
})
mean_meth_non = sapply(gIndexes, function(ii) {
	if(length(ii) > 1) rowMeans(meth_non[,ii]) else meth_non[,ii]
})

## load expression data
library(SummarizedExperiment)
load("/dcl01/lieber/ajaffe/CellSorting/RNAseq_pipeline/rse_gene_CellSorting_July5_n12.Rdata")
load("/dcl01/lieber/ajaffe/CellSorting/RNAseq_pipeline/rse_exon_CellSorting_July5_n12.Rdata")
load("/dcl01/lieber/ajaffe/CellSorting/RNAseq_pipeline/rse_jx_CellSorting_July5_n12.Rdata")

# and SGESeq
xx=load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/sorted_nuclear_RNA/SGSeq_objects_sortedRNA.rda")

# make RSE
exprpos <- data.frame(
		spliceid = mcols(sgvc)$variantName,
		chr = as.character(seqnames(unlist(range(rowRanges(sgvc))))),
		start = min(start(rowRanges(sgvc))),
		end = max(end(rowRanges(sgvc))),
		stringsAsFactors = FALSE
	)
	
psi = variantFreq(sgvc)
rownames(psi) = exprpos$spliceid
gr = makeGRangesFromDataFrame(exprpos)
mcols(gr) = mcols(sgvc)
rse_psi = SummarizedExperiment(assays = list(psi = psi),
	rowRanges = gr, colData = si)

		
## keep only main chrs
mainChrs = paste0("chr", c(1:22, "X","Y","M"))
rse_gene = rse_gene[seqnames(rse_gene) %in% mainChrs,]
rse_exon = rse_exon[seqnames(rse_exon) %in% mainChrs,]
rse_jx = rse_jx[seqnames(rse_jx) %in% mainChrs,]
rse_psi = rse_psi[seqnames(rse_psi) %in% mainChrs,]

# get expression, 12:infant, 13: teen, 14: child
rse_gene$newAgeGroup = "Infant"
rse_gene$newAgeGroup[rse_gene$SampleID == 13] = "Teen"
rse_gene$newAgeGroup[rse_gene$SampleID == 14] = "Child"
rse_gene$Cell.Type = ifelse(rse_gene$NeuN == "NeuN_Plus", 
	"Neuron", "Glia")

## get exprs itself
rowData(rse_jx)$bp_length = 100
jRp10m = getRPKM(rse_jx)
jMap = rowRanges(rse_jx)
geneRpkm = getRPKM(rse_gene, "Length")
geneMap = rowRanges(rse_gene)
exonRpkm = getRPKM(rse_exon, "Length")
exonMap = rowRanges(rse_exon)
psiMat = assays(rse_psi)$psi
psiMap = rowRanges(rse_psi)

## split by group	
gIndexesExprs = splitit(paste0(rse_gene$Cell.Type, 
	":", rse_gene$newAgeGroup))
mean_jRp10m = sapply(gIndexesExprs, function(ii) {
	if(length(ii) > 1) rowMeans(jRp10m[,ii]) else jRp10m[,ii]
})
mean_geneRpkm = sapply(gIndexesExprs, function(ii) {
	if(length(ii) > 1) rowMeans(geneRpkm[,ii]) else geneRpkm[,ii]
})
mean_exonRpkm = sapply(gIndexesExprs, function(ii) {
	if(length(ii) > 1) rowMeans(exonRpkm[,ii]) else exonRpkm[,ii]
})
mean_psiMat = psiMat[,c(3,1,2,6,4,5)]
colnames(mean_psiMat) = colnames(mean_exonRpkm)

## save
save(mean_jRp10m, mean_geneRpkm, mean_exonRpkm,
	mean_psiMat, psiMap,
	jMap, exonMap, geneMap, mean_meth_cpg, mean_meth_non,
	gr_cpg, gr_non, 	compress=TRUE,
	file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/sorted_nuclear_RNA/matched_up_genomicData_meaned.rda")

###############################
### link DNAm to expression ###
###############################

load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/sorted_nuclear_RNA/matched_up_genomicData_meaned.rda")

## cuts
dnamCuts = c(0,0.2,0.8,Inf)
rnaseqLogCuts = c(0,1,2,5,Inf)

#######
# Gene Level ##
###########

# promoters, -2000bp to +200bp
genePromoters = GRanges(seqnames(geneMap),
	IRanges(start = ifelse(strand(geneMap) == "+",
		start(geneMap)-2000, end(geneMap)-200),
	end = ifelse(strand(geneMap) == "+",
		start(geneMap)+200, end(geneMap)+2000)),
		strand = strand(geneMap))
mcols(genePromoters) = mcols(geneMap)

## cpg
ooGene_promoter_cpg = findOverlaps(genePromoters, gr_cpg)
mean_meth_cpg_genePromoter = t(sapply(split(subjectHits(ooGene_promoter_cpg),
	queryHits(ooGene_promoter_cpg)), function(ii) {
		if(length(ii) > 1) colMeans(mean_meth_cpg[ii,]) else mean_meth_cpg[ii,]
}))
mean_geneRpkm_cpg_genePromoter = mean_geneRpkm[
	as.numeric(rownames(mean_meth_cpg_genePromoter)),]
	
## check correlation
signif(diag(cor(mean_meth_cpg_genePromoter,
	log2(mean_geneRpkm_cpg_genePromoter+1))),3)
smoothScatter(	mean_meth_cpg_genePromoter[,6],
	log2(mean_geneRpkm_cpg_genePromoter[,6]+1))
	
## binned
mean_meth_cpg_genePromoter_cut = apply(mean_meth_cpg_genePromoter, 2, 
	cut, breaks = dnamCuts, labels = FALSE,include.lowest=TRUE)
mean_exprs_cut_cpg_genePromoter	= apply(log2(mean_geneRpkm_cpg_genePromoter+1),2,
	cut, breaks = rnaseqLogCuts, labels = FALSE,include.lowest=TRUE)
signif(diag(cor(mean_meth_cpg_genePromoter_cut,
	mean_exprs_cut_cpg_genePromoter+1)),2)
table(mean_meth_cpg_genePromoter_cut[,1] , mean_exprs_cut_cpg_genePromoter[,1])
	
## non-cpg
ooGene_promoter_non = findOverlaps(genePromoters, gr_non)
mean_meth_non_genePromoter = t(sapply(split(subjectHits(ooGene_promoter_non),
	queryHits(ooGene_promoter_non)), function(ii) {
		if(length(ii) > 1) colMeans(mean_meth_non[ii,]) else mean_meth_non[ii,]
}))
mean_geneRpkm_non_genePromoter = mean_geneRpkm[
	as.numeric(rownames(mean_meth_non_genePromoter)),]
## check correlation
signif(diag(cor(mean_meth_non_genePromoter,
	log2(mean_geneRpkm_non_genePromoter+1))),3)
smoothScatter(	mean_meth_non_genePromoter[,6],
	log2(mean_geneRpkm_non_genePromoter[,6]+1))

########################
## check CpH context ###
cIndexes = splitit(gr_non$trinucleotide_context)
	
corByContext_promoters = t(sapply(cIndexes, function(ii) {
	gr_non_sub = gr_non[ii]
	mean_meth_non_sub = mean_meth_non[ii,]

	## non-cpg
	ooGene_promoter_non = findOverlaps(genePromoters, gr_non_sub)
	mean_meth_non_genePromoter = t(sapply(split(subjectHits(ooGene_promoter_non),
		queryHits(ooGene_promoter_non)), function(ii) {
			if(length(ii) > 1) colMeans(mean_meth_non_sub[ii,]) else mean_meth_non_sub[ii,]
	}))
	mean_geneRpkm_non_genePromoter = mean_geneRpkm[
		as.numeric(rownames(mean_meth_non_genePromoter)),]
	## check correlation
	signif(diag(cor(mean_meth_non_genePromoter,
		log2(mean_geneRpkm_non_genePromoter+1))),2)
}))
		
#######
### gene body, rest of gene less promoter
geneBody = geneMap
longIndex = which(width(geneBody) > 250)
start(geneBody)[longIndex] = ifelse(strand(geneBody) == "+",
		start(geneBody)+200, start(geneBody))[longIndex]
end(geneBody)[longIndex] = ifelse(strand(geneBody) == "+",
		end(geneBody), end(geneBody)-200)[longIndex]
		
## cpg
ooGene_body_cpg = findOverlaps(geneBody, gr_cpg)
mean_meth_cpg_geneBody = t(sapply(split(subjectHits(ooGene_body_cpg),
	queryHits(ooGene_body_cpg)), function(ii) {
		if(length(ii) > 1) colMeans(mean_meth_cpg[ii,]) else mean_meth_cpg[ii,]
}))
mean_geneRpkm_cpg_geneBody = mean_geneRpkm[
	as.numeric(rownames(mean_meth_cpg_geneBody)),]
## check correlation
signif(diag(cor(mean_meth_cpg_geneBody,log2(mean_geneRpkm_cpg_geneBody+1))),3)

## non-non
ooGene_body_non = findOverlaps(geneBody, gr_non)
mean_meth_non_geneBody = t(sapply(split(subjectHits(ooGene_body_non),
	queryHits(ooGene_body_non)), function(ii) {
		if(length(ii) > 1) colMeans(mean_meth_non[ii,]) else mean_meth_non[ii,]
}))
mean_geneRpkm_non_geneBody = mean_geneRpkm[
	as.numeric(rownames(mean_meth_non_geneBody)),]
## check correlation
signif(diag(cor(mean_meth_non_geneBody,
	log2(mean_geneRpkm_non_geneBody+1))),3)
smoothScatter(	mean_meth_non_geneBody[,6],
	log2(mean_geneRpkm_non_geneBody[,6]+1))
	
# by context	
corByContext_body = t(sapply(cIndexes, function(ii) {
	gr_non_sub = gr_non[ii]
	mean_meth_non_sub = mean_meth_non[ii,]

	## non-cpg
	ooGene_body_non = findOverlaps(geneBody, gr_non_sub)
	mean_meth_non_geneBody = t(sapply(split(subjectHits(ooGene_body_non),
	queryHits(ooGene_body_non)), function(ii) {
		if(length(ii) > 1) colMeans(mean_meth_non_sub[ii,]) else mean_meth_non_sub[ii,]
	}))
	mean_geneRpkm_non_geneBody = mean_geneRpkm[
		as.numeric(rownames(mean_meth_non_geneBody)),]
	## check correlation
	signif(diag(cor(mean_meth_non_geneBody,
		log2(mean_geneRpkm_non_geneBody+1))),2)
}))

## flanking
geneBodyLeft = geneBodyRight = geneBody
start(geneBodyLeft) = start(geneBody)-100000
end(geneBodyLeft) = start(geneBody)
start(geneBodyRight) = end(geneBody)
end(geneBodyRight) = end(geneBody) + 100000

ooGene_cpg_left = findOverlaps(geneBodyLeft, gr_cpg)
ooGene_cpg_right = findOverlaps(geneBodyRight, gr_cpg)
ooGene_cpg = as.data.frame(rbind(as.matrix(ooGene_cpg_left),
	as.matrix(ooGene_cpg_right)))

mean_meth_cpg_gene_flanking = t(sapply(split(ooGene_cpg$subjectHits,
	ooGene_cpg$queryHits), function(ii) {
		if(length(ii) > 1) colMeans(mean_meth_cpg[ii,]) else mean_meth_cpg[ii,]
}))
n = intersect(rownames(mean_meth_cpg_geneBody), rownames(mean_meth_cpg_gene_flanking))
mean_meth_cpg_gene_flanking = mean_meth_cpg_gene_flanking[n,]
mean_meth_cpg_geneBody = mean_meth_cpg_geneBody[n,]

boxplot(mean_meth_cpg_geneBody)
boxplot(mean_meth_cpg_geneBody/mean_meth_cpg_gene_flanking)

mean_meth_ratioToFlanking = mean_meth_cpg_geneBody/mean_meth_cpg_gene_flanking
mean_geneRpkm_geneBody_ratio = mean_geneRpkm[
	as.numeric(rownames(mean_meth_ratioToFlanking)),]
	
## check correlation
signif(diag(cor(mean_meth_ratioToFlanking,
	log2(mean_geneRpkm_geneBody_ratio+1))),3)
	
####################
### Exon level #####
####################	
	
## cpg
ooExon_cpg = findOverlaps(exonMap, gr_cpg, maxgap=500)
mean_meth_cpg_exon = t(sapply(split(subjectHits(ooExon_cpg),
	queryHits(ooExon_cpg)), function(ii) {
		if(length(ii) > 1) colMeans(mean_meth_cpg[ii,]) else mean_meth_cpg[ii,]
}))
mean_exonRpkm_cpg = mean_exonRpkm[
	as.numeric(rownames(mean_meth_cpg_exon)),]
## check correlation
signif(diag(cor(mean_meth_cpg_exon,
	log2(mean_exonRpkm_cpg+1))),3)

## non-non
ooExon_non = findOverlaps(exonMap, gr_non, maxgap=500)
mean_meth_non_exon = t(sapply(split(subjectHits(ooExon_non),
	queryHits(ooExon_non)), function(ii) {
		if(length(ii) > 1) colMeans(mean_meth_non[ii,]) else mean_meth_non[ii,]
}))
mean_exonRpkm_non = mean_exonRpkm[
	as.numeric(rownames(mean_meth_non_exon)),]
## check correlation
signif(diag(cor(mean_meth_non_exon,
	log2(mean_exonRpkm_non+1))),3)

## non-non
corByContext_exon = t(sapply(cIndexes, function(ii) {
	gr_non_sub = gr_non[ii]
	mean_meth_non_sub = mean_meth_non[ii,]

	ooExon_non = findOverlaps(exonMap, gr_non_sub, maxgap=500)
	mean_meth_non_exon = t(sapply(split(subjectHits(ooExon_non),
		queryHits(ooExon_non)), function(ii) {
			if(length(ii) > 1) colMeans(mean_meth_non_sub[ii,]) else mean_meth_non_sub[ii,]
	}))
	mean_exonRpkm_non = mean_exonRpkm[
		as.numeric(rownames(mean_meth_non_exon)),]
	## check correlation
	signif(diag(cor(mean_meth_non_exon,
		log2(mean_exonRpkm_non+1))),3)
}))
corByContext_exon
####################
### Junction Level #
####################	

## just do the intron sequence first

## cpg
ooJxn_cpg = findOverlaps(jMap, gr_cpg)
mean_meth_cpg_jxn = t(sapply(split(subjectHits(ooJxn_cpg),
	queryHits(ooJxn_cpg)), function(ii) {
		if(length(ii) > 1) colMeans(mean_meth_cpg[ii,]) else mean_meth_cpg[ii,]
}))
mean_jRp10m_cpg = mean_jRp10m[
	as.numeric(rownames(mean_meth_cpg_jxn)),]
## check correlation
signif(diag(cor(mean_meth_cpg_jxn,
	log2(mean_jRp10m_cpg+1))),3)

## non-non
ooJxn_non = findOverlaps(jMap, gr_non)
mean_meth_non_jxn = t(sapply(split(subjectHits(ooJxn_non),
	queryHits(ooJxn_non)), function(ii) {
		if(length(ii) > 1) colMeans(mean_meth_non[ii,]) else mean_meth_non[ii,]
}))
mean_jxnRp10m_non = mean_jRp10m[
	as.numeric(rownames(mean_meth_non_jxn)),]
## check correlation
signif(diag(cor(mean_meth_non_jxn,
	log2(mean_jxnRp10m_non+1))),3)

## just do the 500bp flanking each intron end

## cpg
jMapLeft = jMap
end(jMapLeft) = start(jMapLeft)+50
jMapRight = jMap
start(jMapRight) = end(jMapRight)-50

ooJxn_cpg_left = findOverlaps(jMapLeft, gr_cpg)
ooJxn_cpg_right = findOverlaps(jMapRight, gr_cpg)
ooJxn_cpg = as.data.frame(rbind(as.matrix(ooJxn_cpg_left),
	as.matrix(ooJxn_cpg_right)))

mean_meth_cpg_jxn = t(sapply(split(ooJxn_cpg$subjectHits,
	ooJxn_cpg$queryHits), function(ii) {
		if(length(ii) > 1) colMeans(mean_meth_cpg[ii,]) else mean_meth_cpg[ii,]
}))
mean_jRp10m_cpg = mean_jRp10m[
	as.numeric(rownames(mean_meth_cpg_jxn)),]
## check correlation
signif(diag(cor(mean_meth_cpg_jxn,
	log2(mean_jRp10m_cpg+1))),3)

## non-non
ooJxn_non_left = findOverlaps(jMapLeft, gr_non)
ooJxn_non_right = findOverlaps(jMapRight, gr_non)
ooJxn_non = as.data.frame(rbind(as.matrix(ooJxn_non_left),
	as.matrix(ooJxn_non_right)))
	
mean_meth_non_jxn = t(sapply(split(ooJxn_non$subjectHits,
	ooJxn_non$queryHits), function(ii) {
		if(length(ii) > 1) colMeans(mean_meth_non[ii,]) else mean_meth_non[ii,]
}))
mean_jxnRp10m_non = mean_jRp10m[
	as.numeric(rownames(mean_meth_non_jxn)),]
## check correlation
signif(diag(cor(mean_meth_non_jxn,
	log2(mean_jxnRp10m_non+1))),3)

## by context
## non-non
corByContext_jxn = t(sapply(cIndexes, function(ii) {
	gr_non_sub = gr_non[ii]
	mean_meth_non_sub = mean_meth_non[ii,]
	
	ooJxn_non_left = findOverlaps(jMapLeft, gr_non_sub)
	ooJxn_non_right = findOverlaps(jMapRight, gr_non_sub)
	ooJxn_non = as.data.frame(rbind(as.matrix(ooJxn_non_left),
		as.matrix(ooJxn_non_right)))
		
	mean_meth_non_jxn = t(sapply(split(ooJxn_non$subjectHits,
		ooJxn_non$queryHits), function(ii) {
			if(length(ii) > 1) colMeans(mean_meth_non_sub[ii,]) else mean_meth_non_sub[ii,]
	}))
	mean_jxnRp10m_non = mean_jRp10m[
		as.numeric(rownames(mean_meth_non_jxn)),]
	## check correlation
	signif(diag(cor(mean_meth_non_jxn,
		log2(mean_jxnRp10m_non+1))),3)
}))
corByContext_jxn