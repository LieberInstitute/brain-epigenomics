###
library(GenomicRanges)
library(jaffelab)


### load regions to check
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/DMR/DMR_objects.rda", verbose = T)
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/PMDs_100kb_noGaps.rda", verbose = T)
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/UMRs_LMRs_annotated.rda", verbose = T)
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/dmvs_annotated.rda", verbose = T)
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/bumphunting/rda/limma_Neuron_CpGs_minCov_3_ageInfo_dmrs.Rdata")
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/DMV_gene_comps.rda")
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/PMD_gene_comps.rda")


## Prepare DMRs

cdDmrs = split(dmrs, dmrs$k6cluster_label)
cdDmrs$All = dmrs
dmrsList = endoapply(cdDmrs, granges)

ageDMRs = makeGRangesFromDataFrame(DMR$Age[DMR$Age$sig=="FWER < 0.05",])
dmrsList$Age = ageDMRs[!duplicated(ageDMRs)]
cellDMRs = makeGRangesFromDataFrame(DMR$CellType[DMR$CellType$sig=="FWER < 0.05",])
dmrsList$CellType = cellDMRs[!duplicated(cellDMRs)]

## Prepare methylation features

methfeatures = list(UMR = lapply(uDMR, makeGRangesFromDataFrame, keep=T), 
				LMR = lapply(lDMR, makeGRangesFromDataFrame, keep=T), 
                PMD = lapply(total.nogaps[-grep("GSM", names(total.nogaps))], function(x) x[which(x$type=="PMD")]), 
                DMV = lapply(dmvDMR, makeGRangesFromDataFrame, keep=T))
methfeatures$UMR = mapply(function(u,d) u[!u %in% d],  methfeatures$UMR, 
						methfeatures$DMV, SIMPLIFY = F) # because DMV is a special case of UMR, remove DMVs from UMR list
rmethfeatures = lapply(methfeatures, function(m) lapply(m, reduce))
methFeat = unlist(rmethfeatures)

## Object lists to check for chromatin state
allList = c(dmrsList, GRangesList(methFeat))
length(allList)

#########################################
## load roadmap genomic state objects
load("/dcl01/lieber/ajaffe/PublicData/EpigenomeRoadmap/ChromHMM/chromHMM_15state_coverageDF.rda")

dat = read.delim("/dcl01/lieber/ajaffe/PublicData/EpigenomeRoadmap/ChromHMM/jul2013.roadmapData.qc - Consolidated_EpigenomeIDs_summary_Table.tsv",
	as.is=TRUE)
colnames(dat)[2] = "EID"
dat = dat[dat$EID != "",]

brainIDs = which(dat$ANATOMY == "BRAIN")

#####################
### OVERLAPS ########
#####################

## only autosomal
allList= endoapply(allList, keepSeqlevels, names(stateCovDf), pruning.mode="coarse")
length(allList)
sum(width(allList))

# split by chr
enrichmentsList = mclapply(allList, function(x) {
	
	byChr = split(x, seqnames(x))
	byChr = byChr[names(stateCovDf)]

	bpOverlapByChr = mapply(function(regByChr, stateByChr) {
		cat(".")	
		enr = stateByChr[ranges(regByChr),]
		bg = stateByChr

		list(enrTab = sapply(enr, table),
			bgTab  = sapply(bg, table))
	}, byChr, stateCovDf,SIMPLIFY=FALSE)

	bpOverlapEnr = sapply(bpOverlapByChr, "[", "enrTab")
	bpOverlapArray = array(unlist(bpOverlapEnr), 
		dim = c(nrow(bpOverlapEnr[[1]]), 
		ncol(bpOverlapEnr[[1]]), length(bpOverlapEnr)))
	enrTab = apply(bpOverlapArray, 1:2, sum)
	dimnames(enrTab) = list(rownames(bpOverlapEnr[[1]]),
		colnames(bpOverlapEnr[[1]]))
	
		
	# states in rest of genome	
	bpOverlapBg = sapply(bpOverlapByChr, "[", "bgTab")
	bpOverlapArray = array(unlist(bpOverlapBg), 
		dim = c(nrow(bpOverlapBg[[1]]), 
		ncol(bpOverlapBg[[1]]), length(bpOverlapBg)))
	bgTab = apply(bpOverlapArray, 1:2, function(x) sum(as.numeric(x)))
	dimnames(bgTab) = list(rownames(bpOverlapBg[[1]]),
		colnames(bpOverlapBg[[1]]))

	# take ratio
	statTab = as.data.frame(t(prop.table(enrTab,2) / 
		prop.table(bgTab,2)))
	statTab$Sample = dat$Standardized.Epigenome.name
	statTab = statTab[,c(16,1:15)]
	return(statTab)
}, mc.cores=5)
save(enrichmentsList, allList, dat, compress=TRUE,
	file = "../rdas/enrichmentState_byObj_allRoadmap.rda")

###############
# make plots ##
###############

# for the roadmap state analysis, I'd like 
# 1) comparison of cell type, overall age, and interaction DMRs in brain (adult DLPFC?) to other tissues
# 2) in adult DLPFC and fetal cortex samples, enrichment of UMR, LMR, PMD, and DMV states

# long stats
enrichListLong = lapply(enrichmentsList, function(x) {
	y = data.frame(State = rep(colnames(x)[-1],each=nrow(x)),
			Enrichment = as.numeric(as.matrix(x[,-1])),
			Sample = rep(rownames(x), times= ncol(x)-1),
			stringsAsFactors=FALSE)
			
	y$log2Enrich = log2(y$Enrichment+0.001) # use offset
	y$State = factor(y$State, levels=unique(y$State))
	y$isBrain = y$Sample %in% dat$EID[brainIDs]
	y$isDLPFC = y$Sample == "E073"
	return(y)
})

## comparison of cell type, overall age, and interaction DMRs in brain (adult DLPFC?) to other tissues
enrichDMRs = enrichListLong[1:9]
names(enrichDMRs)[7] = "All_cdDMRs")
# plot
pdf("../plots/DMRs_roadmap_enrichments.pdf",h=5,w=10,useDingbats=FALSE)
for(i in seq(along=enrichDMRs)) {
	x = enrichDMRs[[i]]
	par(mar=c(10,6,2,2),cex.axis=1.6,cex.lab=1.6,cex.main=1.6)
	boxplot(log2Enrich ~ State, data = x, las=3, outline=FALSE,
		main = names(enrichDMRs)[i],ylim=c(-10,10),
		ylab = "Coverage Ratio (log2)")
	points(log2Enrich ~ jitter(as.numeric(State), amount=0.15), 
		data = x, pch = ifelse(isDLPFC, 22, 21), 
		bg = ifelse(isBrain, "red", "grey"),
		cex = ifelse(!isBrain, 0.5, ifelse(!isDLPFC, 1,1.5)))
	abline(h=0,lty=2,col="blue")
}
dev.off()





# long, all
allLong = data.frame(State = rep(colnames(statTab_all)[-1],each=nrow(statTab_all)),
	Enrichment = as.numeric(as.matrix(statTab_all[,-1])),
	Sample = rep(rownames(statTab_all), times= ncol(statTab_all)-1))
	
allLong$log2Enrich = log2(allLong$Enrichment)
allLong$State = factor(allLong$State, levels=unique(allLong$State))
allLong$inBrain = as.numeric(allLong$Sample %in% dat$EID[brainIDs])+1

abline(h=0,lty=2,col="blue")
legend("topright", "Brain", col = "red", pch =15,cex=1.6)
boxplot(log2Enrich ~ State, data = allLong, las=3, outline=FALSE,
	main = "NHGRI GWAS",ylim=c(-4,4),
	ylab = "Coverage Ratio (log2)")
points(log2Enrich ~ jitter(as.numeric(State), amount=0.15), 
	data = allLong, pch = 21, bg = ifelse(inBrain == 2, "red", "grey"),
	cex = ifelse(inBrain == 2, 1, 0.5))
abline(h=0,lty=2,col="blue")
dev.off()
	
#### extract numbers for SCZD
round(2^tapply(pgcLong$log2Enrich, pgcLong$State, 
	function(x) mean(x[is.finite(x)])),3)
round(2^tapply(pgcLong$log2Enrich[pgcLong$inBrain==1], 
	pgcLong$State[pgcLong$inBrain==1], 
	function(x) mean(x[is.finite(x)])),3)
	
#### extract numbers for NHGRI
round(2^tapply(allLong$log2Enrich, allLong$State, 
	function(x) mean(x[is.finite(x)])),3)
round(2^tapply(allLong$log2Enrich[allLong$inBrain==1], 
	allLong$State[allLong$inBrain==1], 
	function(x) mean(x[is.finite(x)])),3)
	