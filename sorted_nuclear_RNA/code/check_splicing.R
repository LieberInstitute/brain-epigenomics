####

library(limma)
library(edgeR)
library(bsseq)
library(recount)
library(jaffelab)
library(SummarizedExperiment)
library(SGSeq)


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
load("../Downloads/matched_up_genomicData_meaned.rda", verbose=T)

## cuts
dnamCuts = c(0,0.2,0.8,Inf)
rnaseqLogCuts = c(0,1,2,5,Inf)

## set up genomic features granges

# promoters, -2000bp to +200bp
genePromoters = GRanges(seqnames(geneMap), 
                        IRanges(start = ifelse(strand(geneMap) == "+",
                                               start(geneMap)-2000, end(geneMap)-200),
                                end = ifelse(strand(geneMap) == "+",
                                             start(geneMap)+200, end(geneMap)+2000)),
                        strand = strand(geneMap))
mcols(genePromoters) = mcols(geneMap)

# Gene bodies
geneBody = geneMap
longIndex = which(width(geneBody) > 250)
start(geneBody)[longIndex] = ifelse(strand(geneBody) == "+",
                                    start(geneBody)+200, start(geneBody))[longIndex]
end(geneBody)[longIndex] = ifelse(strand(geneBody) == "+",
                                  end(geneBody), end(geneBody)-200)[longIndex]

# Junctions (the 500bp flanking each intron end)
jMapLeft = jMap
end(jMapLeft) = start(jMapLeft)+50
jMapRight = jMap
start(jMapRight) = end(jMapRight)-50

# context indices
cIndexes = splitit(gr_non$trinucleotide_context)
cIndexes = cIndexes[which(names(cIndexes) %in% c("CAC", "CAG"))]



######## Promoters ########

## CG 

ooGene_promoter_cpg = findOverlaps(genePromoters, gr_cpg)

mean_meth_cpg_genePromoter = t(sapply(split(subjectHits(ooGene_promoter_cpg),
                                            queryHits(ooGene_promoter_cpg)), function(ii) {
                                              if(length(ii) > 1) colMeans(mean_meth_cpg[ii,]) 
                                              else mean_meth_cpg[ii,]
                                              }))
mean_geneRpkm_cpg_genePromoter = mean_geneRpkm[as.numeric(rownames(mean_meth_cpg_genePromoter)),]
rownames(mean_meth_cpg_genePromoter) = rownames(mean_geneRpkm_cpg_genePromoter)

signif(diag(cor(mean_meth_cpg_genePromoter, log2(mean_geneRpkm_cpg_genePromoter+1))),3)
	
	
## CH

ooGene_promoter_non = findOverlaps(genePromoters, gr_non)

mean_meth_non_genePromoter = t(sapply(split(subjectHits(ooGene_promoter_non),
                                            queryHits(ooGene_promoter_non)), function(ii) {
                                              if(length(ii) > 1) colMeans(mean_meth_non[ii,]) 
                                              else mean_meth_non[ii,]
                                              }))
mean_geneRpkm_non_genePromoter = mean_geneRpkm[as.numeric(rownames(mean_meth_non_genePromoter)),]
rownames(mean_meth_non_genePromoter) = rownames(mean_geneRpkm_non_genePromoter)

signif(diag(cor(mean_meth_non_genePromoter, log2(mean_geneRpkm_non_genePromoter+1))),3)


## CAG and CAC

corByContext_promoters = t(sapply(cIndexes, function(ii) {
  gr_non_sub = gr_non[ii]
	mean_meth_non_sub = mean_meth_non[ii,]
	ooGene_promoter_non = findOverlaps(genePromoters, gr_non_sub)
	mean_meth_context_genePromoter = t(sapply(split(subjectHits(ooGene_promoter_non),
	                                            queryHits(ooGene_promoter_non)), function(ii) {
	                                              if(length(ii) > 1) colMeans(mean_meth_non_sub[ii,]) 
	                                              else mean_meth_non_sub[ii,]
	                                              }))
	mean_geneRpkm_context_genePromoter = mean_geneRpkm[as.numeric(rownames(mean_meth_context_genePromoter)),]
	signif(diag(cor(mean_meth_context_genePromoter, log2(mean_geneRpkm_context_genePromoter+1))),3)
	}))
		
mean_meth_context_genePromoter = t(sapply(cIndexes, function(ii) {
  gr_non_sub = gr_non[ii]
  mean_meth_non_sub = mean_meth_non[ii,]
  ooGene_promoter_non = findOverlaps(genePromoters, gr_non_sub)
  t(sapply(split(subjectHits(ooGene_promoter_non),
                 queryHits(ooGene_promoter_non)), function(ii) {
                   if(length(ii) > 1) colMeans(mean_meth_non_sub[ii,]) 
                   else mean_meth_non_sub[ii,]
                   }))
  }))

mean_geneRpkm_context_genePromoter = lapply(mean_meth_context_genePromoter, 
                                            function(ii) mean_geneRpkm[as.numeric(rownames(ii)),])
names(mean_meth_context_genePromoter) = names(mean_geneRpkm_context_genePromoter) = names(cIndexes)
rownames(mean_meth_context_genePromoter$CAC) = rownames(mean_geneRpkm_context_genePromoter$CAC)
rownames(mean_meth_context_genePromoter$CAG) = rownames(mean_geneRpkm_context_genePromoter$CAG)


######## gene body (rest of gene less promoter) ########

## CG

ooGene_body_cpg = findOverlaps(geneBody, gr_cpg)
mean_meth_cpg_geneBody = t(sapply(split(subjectHits(ooGene_body_cpg),
                                        queryHits(ooGene_body_cpg)), function(ii) {
                                          if(length(ii) > 1) colMeans(mean_meth_cpg[ii,]) 
                                          else mean_meth_cpg[ii,]
                                          }))
mean_geneRpkm_cpg_geneBody = mean_geneRpkm[as.numeric(rownames(mean_meth_cpg_geneBody)),]
rownames(mean_meth_cpg_geneBody) = rownames(mean_geneRpkm_cpg_geneBody)

signif(diag(cor(mean_meth_cpg_geneBody,log2(mean_geneRpkm_cpg_geneBody+1))),3)


## CH

ooGene_body_non = findOverlaps(geneBody, gr_non)
mean_meth_non_geneBody = t(sapply(split(subjectHits(ooGene_body_non),
                                        queryHits(ooGene_body_non)), function(ii) {
                                          if(length(ii) > 1) colMeans(mean_meth_non[ii,]) 
                                          else mean_meth_non[ii,]
                                          }))
mean_geneRpkm_non_geneBody = mean_geneRpkm[as.numeric(rownames(mean_meth_non_geneBody)),]
rownames(mean_meth_non_geneBody) = rownames(mean_geneRpkm_non_geneBody)

signif(diag(cor(mean_meth_non_geneBody, log2(mean_geneRpkm_non_geneBody+1))),3)


## CAG and CAC

corByContext_body = t(sapply(cIndexes, function(ii) {
  gr_non_sub = gr_non[ii]
  mean_meth_non_sub = mean_meth_non[ii,]
  ooGene_body_non = findOverlaps(geneBody, gr_non_sub)
  mean_meth_context_geneBody = t(sapply(split(subjectHits(ooGene_body_non),
                                              queryHits(ooGene_body_non)), function(ii) {
                                                if(length(ii) > 1) colMeans(mean_meth_non_sub[ii,]) 
                                                else mean_meth_non_sub[ii,]
                                              }))
  mean_geneRpkm_context_geneBody = mean_geneRpkm[as.numeric(rownames(mean_meth_context_geneBody)),]
  signif(diag(cor(mean_meth_context_geneBody,
                  log2(mean_geneRpkm_context_geneBody+1))),3)
}))

mean_meth_context_geneBody = t(sapply(cIndexes, function(ii) {
  gr_non_sub = gr_non[ii]
  mean_meth_non_sub = mean_meth_non[ii,]
  ooGene_promoter_non = findOverlaps(geneBody, gr_non_sub)
  t(sapply(split(subjectHits(ooGene_promoter_non),
                 queryHits(ooGene_promoter_non)), function(ii) {
                   if(length(ii) > 1) colMeans(mean_meth_non_sub[ii,]) 
                   else mean_meth_non_sub[ii,]
                 }))
}))

mean_geneRpkm_context_geneBody = lapply(mean_meth_context_geneBody, 
                                            function(ii) mean_geneRpkm[as.numeric(rownames(ii)),])
names(mean_meth_context_geneBody) = names(mean_geneRpkm_context_geneBody) = names(cIndexes)
rownames(mean_meth_context_geneBody$CAG) = rownames(mean_geneRpkm_context_geneBody$CAG)
rownames(mean_meth_context_geneBody$CAC) = rownames(mean_geneRpkm_context_geneBody$CAC)



######## Exons ########

## CG

ooExon_cpg = findOverlaps(exonMap, gr_cpg, maxgap=500)
mean_meth_cpg_exon = t(sapply(split(subjectHits(ooExon_cpg),
                                    queryHits(ooExon_cpg)), function(ii) {
                                      if(length(ii) > 1) colMeans(mean_meth_cpg[ii,]) 
                                      else mean_meth_cpg[ii,]
                                      }))
mean_exonRpkm_cpg = mean_exonRpkm[as.numeric(rownames(mean_meth_cpg_exon)),]
rownames(mean_meth_cpg_exon) = rownames(mean_exonRpkm_cpg)

signif(diag(cor(mean_meth_cpg_exon, log2(mean_exonRpkm_cpg+1))),3)


## CH

ooExon_non = findOverlaps(exonMap, gr_non, maxgap=500)
mean_meth_non_exon = t(sapply(split(subjectHits(ooExon_non),
                                    queryHits(ooExon_non)), function(ii) {
                                      if(length(ii) > 1) colMeans(mean_meth_non[ii,]) 
                                      else mean_meth_non[ii,]
                                      }))
mean_exonRpkm_non = mean_exonRpkm[as.numeric(rownames(mean_meth_non_exon)),]
rownames(mean_meth_non_exon) = rownames(mean_exonRpkm_non)

signif(diag(cor(mean_meth_non_exon, log2(mean_exonRpkm_non+1))),3)


## CAG and CAC

corByContext_exon = t(sapply(cIndexes, function(ii) {
  gr_non_sub = gr_non[ii]
  mean_meth_non_sub = mean_meth_non[ii,]
  ooGene_exon_non = findOverlaps(exonMap, gr_non_sub, maxgap=500)
  mean_meth_context_exon = t(sapply(split(subjectHits(ooGene_exon_non),
                                              queryHits(ooGene_exon_non)), function(ii) {
                                                if(length(ii) > 1) colMeans(mean_meth_non_sub[ii,]) 
                                                else mean_meth_non_sub[ii,]
                                              }))
  mean_exonRpkm_context_exon = mean_exonRpkm[as.numeric(rownames(mean_meth_context_exon)),]
  signif(diag(cor(mean_meth_context_exon,
                  log2(mean_exonRpkm_context_exon+1))),3)
}))

mean_meth_context_exon = t(sapply(cIndexes, function(ii) {
  gr_non_sub = gr_non[ii]
  mean_meth_non_sub = mean_meth_non[ii,]
  ooGene_promoter_non = findOverlaps(geneBody, gr_non_sub)
  t(sapply(split(subjectHits(ooGene_promoter_non),
                 queryHits(ooGene_promoter_non)), function(ii) {
                   if(length(ii) > 1) colMeans(mean_meth_non_sub[ii,]) 
                   else mean_meth_non_sub[ii,]
                 }))
}))

mean_geneRpkm_context_exon = lapply(mean_meth_context_exon, 
                                        function(ii) mean_geneRpkm[as.numeric(rownames(ii)),])
names(mean_meth_context_exon) = names(mean_geneRpkm_context_exon) = names(cIndexes)
rownames(mean_meth_context_exon$CAG) = rownames(mean_geneRpkm_context_exon$CAG)
rownames(mean_meth_context_exon$CAC) = rownames(mean_geneRpkm_context_exon$CAC)


######## Junctions ########

## CG

ooJxn_cpg_left = findOverlaps(jMapLeft, gr_cpg)
ooJxn_cpg_right = findOverlaps(jMapRight, gr_cpg)
ooJxn_cpg = as.data.frame(rbind(as.matrix(ooJxn_cpg_left), as.matrix(ooJxn_cpg_right)))

mean_meth_cpg_jxn = t(sapply(split(ooJxn_cpg$subjectHits, 
                                   ooJxn_cpg$queryHits), function(ii) {
                                     if(length(ii) > 1) colMeans(mean_meth_cpg[ii,]) 
                                     else mean_meth_cpg[ii,]
                                     }))
mean_jRp10m_cpg = mean_jRp10m[as.numeric(rownames(mean_meth_cpg_jxn)),]
rownames(mean_meth_cpg_jxn) = rownames(mean_jRp10m_cpg)

signif(diag(cor(mean_meth_cpg_jxn, log2(mean_jRp10m_cpg+1))),3)

## CH

ooJxn_non_left = findOverlaps(jMapLeft, gr_non)
ooJxn_non_right = findOverlaps(jMapRight, gr_non)
ooJxn_non = as.data.frame(rbind(as.matrix(ooJxn_non_left), as.matrix(ooJxn_non_right)))
	
mean_meth_non_jxn = t(sapply(split(ooJxn_non$subjectHits,
                                   ooJxn_non$queryHits), function(ii) {
                                     if(length(ii) > 1) colMeans(mean_meth_non[ii,]) 
                                     else mean_meth_non[ii,]
                                     }))
mean_jxnRp10m_non = mean_jRp10m[as.numeric(rownames(mean_meth_non_jxn)),]
rownames(mean_meth_non_jxn) = rownames(mean_jxnRp10m_non)

signif(diag(cor(mean_meth_non_jxn, log2(mean_jxnRp10m_non+1))),3)


## CAG and CAC

corByContext_junction = t(sapply(cIndexes, function(ii) {
  gr_non_sub = gr_non[ii]
  mean_meth_non_sub = mean_meth_non[ii,]
  ooJxn_non_left = findOverlaps(jMapLeft, gr_non_sub)
  ooJxn_non_right = findOverlaps(jMapRight, gr_non_sub)
  ooJxn_non = as.data.frame(rbind(as.matrix(ooJxn_non_left),
                                  as.matrix(ooJxn_non_right)))
  mean_meth_context_junction = t(sapply(split(ooJxn_non$subjectHits,
                                              ooJxn_non$queryHits), function(ii) {
                                            if(length(ii) > 1) colMeans(mean_meth_non_sub[ii,]) 
                                            else mean_meth_non_sub[ii,]
                                          }))
  mean_junctionRpkm_context_junction = mean_jRp10m[as.numeric(rownames(mean_meth_context_junction)),]
  signif(diag(cor(mean_meth_context_junction,
                  log2(mean_junctionRpkm_context_junction+1))),3)
}))

mean_meth_context_junction = t(sapply(cIndexes, function(ii) {
  gr_non_sub = gr_non[ii]
  mean_meth_non_sub = mean_meth_non[ii,]
  ooJxn_non_left = findOverlaps(jMapLeft, gr_non_sub)
  ooJxn_non_right = findOverlaps(jMapRight, gr_non_sub)
  ooJxn_non = as.data.frame(rbind(as.matrix(ooJxn_non_left),
                                  as.matrix(ooJxn_non_right)))
  mean_meth_context_junction = t(sapply(split(ooJxn_non$subjectHits,
                                              ooJxn_non$queryHits), function(ii) {
                                                if(length(ii) > 1) colMeans(mean_meth_non_sub[ii,]) 
                                                else mean_meth_non_sub[ii,]
                                              }))
  }))

mean_Rpkm_context_junction = lapply(mean_meth_context_junction, 
                                    function(ii) mean_jRp10m[as.numeric(rownames(ii)),])
names(mean_meth_context_junction) = names(mean_Rpkm_context_junction) = names(cIndexes)
rownames(mean_meth_context_junction$CAG) = rownames(mean_Rpkm_context_junction$CAG)
rownames(mean_meth_context_junction$CAC) = rownames(mean_Rpkm_context_junction$CAC)


#### Make plots of all correlation (subsampled to 10,000 each)

## Organize CG object

CG = list(Promoter = cbind(reshape2::melt(mean_meth_cpg_genePromoter), 
                           Expression = reshape2::melt(log2(mean_geneRpkm_cpg_genePromoter+1))$value),
          "Gene Body" = cbind(reshape2::melt(mean_meth_cpg_geneBody),
                              Expression = reshape2::melt(log2(mean_geneRpkm_cpg_geneBody+1))$value),
          Exon = cbind(reshape2::melt(mean_meth_cpg_exon), Expression = reshape2::melt(log2(mean_exonRpkm_cpg+1))$value),
          Junction = cbind(reshape2::melt(mean_meth_cpg_jxn), Expression = reshape2::melt(log2(mean_jRp10m_cpg+1))$value))
range(elementNROWS(CG)) # 301230-3204642
CGsamp = lapply(CG, function(x) sample(1:nrow(x), 10000))
CG = mapply(function(r, s) r[s,], CG, CGsamp, SIMPLIFY = F)

CG = do.call(rbind, Map(cbind, CG, Location = as.list(names(CG))))
CG$Var2 = factor(CG$Var2, levels = c("Glia:Infant","Glia:Child","Glia:Teen","Neuron:Infant","Neuron:Child","Neuron:Teen"))

CGcor = reshape2::melt(data.frame(rbind(signif(diag(cor(mean_meth_cpg_genePromoter, log2(mean_geneRpkm_cpg_genePromoter+1))),3),
                                        signif(diag(cor(mean_meth_cpg_geneBody, log2(mean_geneRpkm_cpg_geneBody+1))),3),
                                        signif(diag(cor(mean_meth_cpg_exon, log2(mean_exonRpkm_cpg+1))),3),
                                        signif(diag(cor(mean_meth_cpg_jxn, log2(mean_jRp10m_cpg+1))),3)), 
                                  Location = c("Promoter", "Gene Body", "Exon", "Junction")))
colnames(CGcor) = c("Location","Var2","cor")
CGcor$Var2 = gsub(".",":", CGcor$Var2, fixed=T)
CGcor$Var2 = factor(CGcor$Var2, levels = c("Glia:Infant","Glia:Child","Glia:Teen","Neuron:Infant","Neuron:Child","Neuron:Teen"))


## Organize CH object

CH = list(Promoter = cbind(reshape2::melt(mean_meth_non_genePromoter), 
                           Expression = reshape2::melt(log2(mean_geneRpkm_non_genePromoter+1))$value),
          "Gene Body" = cbind(reshape2::melt(mean_meth_non_geneBody),
                              Expression = reshape2::melt(log2(mean_geneRpkm_non_geneBody+1))$value),
          Exon = cbind(reshape2::melt(mean_meth_non_exon), Expression = reshape2::melt(log2(mean_exonRpkm_non+1))$value),
          Junction = cbind(reshape2::melt(mean_meth_non_jxn), Expression = reshape2::melt(log2(mean_jxnRp10m_non+1))$value))

range(elementNROWS(CH)) # 232668-2151780
CHsamp = lapply(CH, function(x) sample(1:nrow(x), 10000))
CH = mapply(function(r, s) r[s,], CH, CHsamp, SIMPLIFY = F)
CH = do.call(rbind, Map(cbind, CH, Location = as.list(names(CH))))
CH$Var2 = factor(CH$Var2, levels = c("Glia:Infant","Glia:Child","Glia:Teen","Neuron:Infant","Neuron:Child","Neuron:Teen"))

CHcor = reshape2::melt(data.frame(rbind(signif(diag(cor(mean_meth_non_genePromoter, log2(mean_geneRpkm_non_genePromoter+1))),3),
                                        signif(diag(cor(mean_meth_non_geneBody, log2(mean_geneRpkm_non_geneBody+1))),3),
                                        signif(diag(cor(mean_meth_non_exon, log2(mean_exonRpkm_non+1))),3),
                                        signif(diag(cor(mean_meth_non_jxn, log2(mean_jxnRp10m_non+1))),3)), 
                                  Location = c("Promoter", "Gene Body", "Exon", "Junction")))
colnames(CHcor) = c("Location","Var2","cor")
CHcor$Var2 = gsub(".",":", CHcor$Var2, fixed=T)
CHcor$Var2 = factor(CHcor$Var2, levels = c("Glia:Infant","Glia:Child","Glia:Teen","Neuron:Infant","Neuron:Child","Neuron:Teen"))


## Organize CAC object

CAC = list(Promoter = cbind(reshape2::melt(mean_meth_context_genePromoter$CAC), 
                            Expression = reshape2::melt(log2(mean_geneRpkm_context_genePromoter$CAC+1))$value),
           "Gene Body" = cbind(reshape2::melt(mean_meth_context_geneBody$CAC),
                               Expression = reshape2::melt(log2(mean_geneRpkm_context_geneBody$CAC+1))$value),
           Exon = cbind(reshape2::melt(mean_meth_context_exon$CAC), 
                        Expression = reshape2::melt(log2(mean_geneRpkm_context_exon$CAC+1))$value),
           Junction = cbind(reshape2::melt(mean_meth_context_junction$CAC), 
                            Expression = reshape2::melt(log2(mean_Rpkm_context_junction$CAC+1))$value))
range(elementNROWS(CAC)) # 190428-449214
CACsamp = lapply(CAC, function(x) sample(1:nrow(x), 10000))
CAC = mapply(function(r, s) r[s,], CAC, CACsamp, SIMPLIFY = F)
CAC = do.call(rbind, Map(cbind, CAC, Location = as.list(names(CAC))))
CAC$Var2 = factor(CAC$Var2, levels = c("Glia:Infant","Glia:Child","Glia:Teen","Neuron:Infant","Neuron:Child","Neuron:Teen"))
head(CAC)


## Organize CAG object

CAG = list(Promoter = cbind(reshape2::melt(mean_meth_context_genePromoter$CAG), 
                            Expression = reshape2::melt(log2(mean_geneRpkm_context_genePromoter$CAG+1))$value),
           "Gene Body" = cbind(reshape2::melt(mean_meth_context_geneBody$CAG),
                               Expression = reshape2::melt(log2(mean_geneRpkm_context_geneBody$CAG+1))$value),
           Exon = cbind(reshape2::melt(mean_meth_context_exon$CAG), 
                        Expression = reshape2::melt(log2(mean_geneRpkm_context_exon$CAG+1))$value),
           Junction = cbind(reshape2::melt(mean_meth_context_junction$CAG), 
                            Expression = reshape2::melt(log2(mean_Rpkm_context_junction$CAG+1))$value))
range(elementNROWS(CAG)) # 218046-735090
CAGsamp = lapply(CAG, function(x) sample(1:nrow(x), 10000))
CAG = mapply(function(r, s) r[s,], CAG, CAGsamp, SIMPLIFY = F)

CAG = do.call(rbind, Map(cbind, CAG, Location = as.list(names(CAG))))
CAG$Var2 = factor(CAG$Var2, levels = c("Glia:Infant","Glia:Child","Glia:Teen","Neuron:Infant","Neuron:Child","Neuron:Teen"))
head(CAG)

CAHcor = reshape2::melt(data.frame(Context = rep.int(c("CAC","CAG"), 4), 
                                   rbind(corByContext_promoters, corByContext_body,
                                         corByContext_exon, corByContext_junction),
                                   Location = c("Promoter","Promoter", "Gene Body","Gene Body",
                                                "Exon","Exon", "Junction","Junction")))
CAHcor$variable = gsub(".",":", CAHcor$variable, fixed=T)
CAHcor$variable = factor(CAHcor$variable, levels = c("Glia:Infant","Glia:Child","Glia:Teen","Neuron:Infant","Neuron:Child","Neuron:Teen"))
colnames(CAHcor)[3] = "Var2"


pdf("./brain-epigenomics/sorted_nuclear_RNA/figures/correlation_density_plot.pdf",
    height = 6.75, width = 15.5)
ggplot(data = CG, aes(x = value, y = Expression)) +
  geom_point(alpha=0.05) +
  facet_grid(Location ~ Var2) +
  geom_text(data = CGcor, aes(x=0.8, y=14, label=paste0("rho==",cor)), parse = T) +
  theme_bw() + xlab("Methylation") + ylab("Expression") +
  theme(title = element_text(size = 18), 
        text = element_text(size = 18),
        legend.title = element_blank()) +
  geom_smooth(method=lm) +
  ggtitle("CG Methylation vs. Expression")

ggplot(data = CH, aes(x = value, y = Expression)) +
  geom_point(alpha=0.05) + ylim(0,16) +
  facet_grid(Location ~ Var2) +
  geom_text(data = CHcor, aes(x=0.7, y=14, label=paste0("rho==",cor)), parse = T) +
  theme_bw() + xlab("Methylation") + ylab("Expression") +
  theme(title = element_text(size = 18), 
        text = element_text(size = 18),
        legend.title = element_blank()) +
  geom_smooth(method=lm) +
  ggtitle("CH Methylation vs. Expression")

ggplot(data = CAC, aes(x = value, y = Expression)) +
  geom_point(alpha=0.05) + ylim(0,16) +
  facet_grid(Location ~ Var2) +
  geom_text(data = CAHcor[which(CAHcor$Context=="CAC"),], aes(x=0.7, y=14, label=paste0("rho==",value)), parse = T) +
  theme_bw() + xlab("Methylation") + ylab("Expression") +
  theme(title = element_text(size = 18), 
        text = element_text(size = 18),
        legend.title = element_blank()) +
  geom_smooth(method=lm) +
  ggtitle("CAC Methylation vs. Expression")

ggplot(data = CAG, aes(x = value, y = Expression)) +
  geom_point(alpha=0.05) + ylim(0,16) +
  facet_grid(Location ~ Var2) +
  geom_text(data = CAHcor[which(CAHcor$Context=="CAG"),], aes(x=0.7, y=14, label=paste0("rho==",value)), parse = T) +
  theme_bw() + xlab("Methylation") + ylab("Expression") +
  theme(title = element_text(size = 18), 
        text = element_text(size = 18),
        legend.title = element_blank()) +
  geom_smooth(method=lm) +
  ggtitle("CAG Methylation vs. Expression")

dev.off()



##### Extra code #####

### Promoters

smoothScatter(	mean_meth_cpg_genePromoter[,6],
               log2(mean_geneRpkm_cpg_genePromoter[,6]+1),ylim = c(0,5))
## binned
mean_meth_cpg_genePromoter_cut = apply(mean_meth_cpg_genePromoter, 2, 
                                       cut, breaks = dnamCuts, labels = FALSE,include.lowest=TRUE)
mean_exprs_cut_cpg_genePromoter	= apply(log2(mean_geneRpkm_cpg_genePromoter+1),2,
                                        cut, breaks = rnaseqLogCuts, labels = FALSE,include.lowest=TRUE)
signif(diag(cor(mean_meth_cpg_genePromoter_cut,
                mean_exprs_cut_cpg_genePromoter+1)),2)
table(mean_meth_cpg_genePromoter_cut[,1] , mean_exprs_cut_cpg_genePromoter[,1])
prop.table(table(mean_meth_cpg_genePromoter_cut[,1] , mean_exprs_cut_cpg_genePromoter[,1]),2)

smoothScatter(	mean_meth_non_genePromoter[,6],
               log2(mean_geneRpkm_non_genePromoter[,6]+1))


## Gene body 
smoothScatter(	mean_meth_non_geneBody[,6],
               log2(mean_geneRpkm_non_geneBody[,6]+1))

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

###################
## constitutive exons ##
exonMap$numTx_gene = geneMap$NumTx[match(exonMap$gencodeID, names(geneMap))]
exonMap$constit = exonMap$numTx_gene == exonMap$NumTx

exonMapMain = exonMap[exonMap$constit,]
exonMapAlt = exonMap[!exonMap$constit,]

ooExon_cpg_main = findOverlaps(exonMapMain, gr_cpg, maxgap=500)
mean_meth_cpg_exon_main = t(sapply(split(subjectHits(ooExon_cpg_main),
                                         queryHits(ooExon_cpg_main)), function(ii) {
                                           if(length(ii) > 1) colMeans(mean_meth_cpg[ii,]) else mean_meth_cpg[ii,]
                                         }))
mean_exonRpkm_cpg_main = mean_exonRpkm[
  as.numeric(rownames(mean_meth_cpg_exon_main)),]
## check correlation
signif(diag(cor(mean_meth_cpg_exon_main,
                log2(mean_exonRpkm_cpg_main+1))),3)

## non: main
ooExon_non_main = findOverlaps(exonMapMain, gr_non, maxgap=500)
mean_meth_non_exon_main = t(sapply(split(subjectHits(ooExon_non_main),
                                         queryHits(ooExon_non_main)), function(ii) {
                                           if(length(ii) > 1) colMeans(mean_meth_non[ii,]) else mean_meth_non[ii,]
                                         }))
mean_exonRpkm_non_main = mean_exonRpkm[
  as.numeric(rownames(mean_meth_non_exon_main)),]
## check correlation
signif(diag(cor(mean_meth_non_exon_main,
                log2(mean_exonRpkm_non_main+1))),3)

## alt	

ooExon_cpg_alt = findOverlaps(exonMapAlt, gr_cpg, maxgap=500)
mean_meth_cpg_exon_alt = t(sapply(split(subjectHits(ooExon_cpg_alt),
                                        queryHits(ooExon_cpg_alt)), function(ii) {
                                          if(length(ii) > 1) colMeans(mean_meth_cpg[ii,]) else mean_meth_cpg[ii,]
                                        }))
mean_exonRpkm_cpg_alt = mean_exonRpkm[
  as.numeric(rownames(mean_meth_cpg_exon_alt)),]
## check correlation
signif(diag(cor(mean_meth_cpg_exon_alt,
                log2(mean_exonRpkm_cpg_alt+1))),3)

## non: alt
ooExon_non_alt = findOverlaps(exonMapAlt, gr_non, maxgap=500)
mean_meth_non_exon_alt = t(sapply(split(subjectHits(ooExon_non_alt),
                                        queryHits(ooExon_non_alt)), function(ii) {
                                          if(length(ii) > 1) colMeans(mean_meth_non[ii,]) else mean_meth_non[ii,]
                                        }))
mean_exonRpkm_non_alt = mean_exonRpkm[
  as.numeric(rownames(mean_meth_non_exon_alt)),]
## check correlation
signif(diag(cor(mean_meth_non_exon_alt,
                log2(mean_exonRpkm_non_alt+1))),3)


## Junctions
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


## Do more things

library(ggplot2)

df = openxlsx::read.xlsx("./Dropbox/wgbs_development/Third\ submission/Figures/Supplementary\ Tables/TableS6old.xlsx")
df = reshape2::melt(df, id.vars = c("Type","Location"))
df$value = as.numeric(as.character(df$value))
df$CT = unlist(lapply(strsplit(as.character(df$variable), ":", fixed=T), function(x) x[1]))
df$Age = unlist(lapply(strsplit(as.character(df$variable), ":", fixed=T), function(x) x[2]))
df$Age = factor(df$Age, levels = c("Infant", "Child", "Teen"))
df$Location = gsub("Body", "Gene Body", df$Location)
df$Location = factor(df$Location, levels = c("Promoter","Gene Body","Exon","Junction"))
df$Type = factor(df$Type, levels = c("CpG","CpH","CAC","CAG"))

pdf("./Dropbox/brain-epigenomics/sorted_nuclear_RNA/figures/correlation_plot.pdf", width=5.5, height=7)
ggplot(df, aes(x = Age, y = value)) + 
  geom_point(aes(x = Age, y = value, colour = CT)) + 
  geom_line(aes(x = as.numeric(Age), y = value, colour = CT)) +
  geom_line(aes(x = as.numeric(Age), y = 0), linetype="dashed") +
  facet_grid(Location ~ Type) +
  scale_colour_brewer(8, palette="Dark2") +
  ylab("Correlation") + xlab("") +
  theme(title = element_text(size = 18),
        text = element_text(size = 14),
        legend.text = element_text(size=14),
        legend.title=element_blank(),
        legend.position = "bottom",
        axis.text.x = element_text(angle = 30, hjust = 1))
dev.off()


range(df[which(df$Location!="Junction" & df$Type=="CpH" & df$CT=="Neuron"),"value"])
# -0.322 -0.103

range(df[which(df$Location %in% c("Promoter", "Gene Body") & df$Type=="CpG" & df$CT=="Neuron"),"value"])
# -0.418 -0.223




