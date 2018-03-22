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

dmrs = split(dmrs, dmrs$k6cluster_label)
oo = lapply(dmrs, function(x) findOverlaps(x, makeGRangesFromDataFrame(DMR$Interaction)))
dmrs = c(lapply(DMR[which(names(DMR) %in% c("CellType","Age"))], function(x) x[which(x$sig=="FWER < 0.05"),]), lapply(oo, function(x) DMR$Interaction[subjectHits(x),]))
names(dmrs) = c("CellType","Age","Gr1","Gr2","Gr3","Gr4","Gr5","Gr6")
rdmrs = lapply(dmrs, function(x) reduce(makeGRangesFromDataFrame(x)))

## Prepare methylation features

methfeatures = list(UMR = lapply(uDMR, makeGRangesFromDataFrame, keep=T), LMR = lapply(lDMR, makeGRangesFromDataFrame, keep=T), 
                    PMD = lapply(total.nogaps[-grep("GSM", names(total.nogaps))], function(x) x[which(x$type=="PMD")]), 
                    DMV = lapply(dmvDMR, makeGRangesFromDataFrame, keep=T))
methfeatures$UMR = mapply(function(u,d) u[!u %in% d],  methfeatures$UMR, methfeatures$DMV, SIMPLIFY = F) # because DMV is a special case of UMR, remove DMVs from UMR list
rmethfeatures = lapply(methfeatures, function(m) lapply(m, reduce))
do.call(cbind, lapply(rmethfeatures, elementNROWS))


## Object lists to check for chromatin state

rdmrs # Cell type and overall age DMRs, and the 6 kmeans cluster groups of interaction DMRs
rmethfeatures # a list of all UMRs, LMRs, DMVs and PMDs in each samples (52 total)
DMV.CTcomps # genes differentially included in DMVs by cell type
DMV.Agecomps # genes differentially included in DMVs by age in neurons
PMD.CTcomps # genes differentially included in PMDs by cell type
PMD.Agecomps # genes differentially included in PMDs by age in neurons




## load eqtls
load("/users/ajaffe/Lieber/Projects/RNAseq/DLPFC_eQTL_paper/joint/rdas/gwas_hits_allEqtl_subset_fdr.rda")
load("/users/ajaffe/Lieber/Projects/RNAseq/DLPFC_eQTL_paper/joint/rdas/PGC_SZ_hits_allEqtl_subset_fdr.rda")

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

## check PGC 
pgcEqtl_novel = pgcEqtls[which(pgcEqtls$Class == "Novel" & pgcEqtls$Type=="ER"),]
pgcEqtl_novel = pgcEqtl_novel[!duplicated(pgcEqtl_novel$Feature),]
gr_pgc_novel = GRanges(ss(pgcEqtl_novel$Coordinates,":"),
	IRanges(as.numeric(ss(ss(pgcEqtl_novel$Coordinates,":",2), "-")),
			as.numeric(ss(ss(ss(pgcEqtl_novel$Coordinates,":",2), "-",2),"\\("))))
seqlevels(gr_pgc_novel) = names(stateCovDf	)

length(gr_pgc_novel)
sum(width(gr_pgc_novel))

# split by chr
pgcByChr = split(gr_pgc_novel, seqnames(gr_pgc_novel))
pgcByChr = pgcByChr[names(stateCovDf)]

bpOverlapByChr_pgc = mapply(function(regByChr, stateByChr) {
	cat(".")	
	enr = stateByChr[ranges(regByChr),]
	bg = stateByChr

	list(enrTab = sapply(enr, table),
		bgTab  = sapply(bg, table))
}, pgcByChr, stateCovDf,SIMPLIFY=FALSE)


# states in s1+s2
bpOverlapEnr = sapply(bpOverlapByChr_pgc, "[", "enrTab")
bpOverlapArray = array(unlist(bpOverlapEnr), 
	dim = c(nrow(bpOverlapEnr[[1]]), 
	ncol(bpOverlapEnr[[1]]), length(bpOverlapEnr)))
enrTab = apply(bpOverlapArray, 1:2, sum)
dimnames(enrTab) = list(rownames(bpOverlapEnr[[1]]),
	colnames(bpOverlapEnr[[1]]))

# states in rest of genome	
bpOverlapBg = sapply(bpOverlapByChr_pgc, "[", "bgTab")
bpOverlapArray = array(unlist(bpOverlapBg), 
	dim = c(nrow(bpOverlapBg[[1]]), 
	ncol(bpOverlapBg[[1]]), length(bpOverlapBg)))
bgTab = apply(bpOverlapArray, 1:2, function(x) sum(as.numeric(x)))
dimnames(bgTab) = list(rownames(bpOverlapBg[[1]]),
	colnames(bpOverlapBg[[1]]))

# take ratio
statTab_pgc = as.data.frame(t(prop.table(enrTab,2) / 
	prop.table(bgTab,2)))
statTab_pgc$Sample = dat$Standardized.Epigenome.name
statTab_pgc = statTab_pgc[,c(16,1:15)]

###############
## check other GWAS
allEqtl_novel = allEqtlGwas[which(allEqtlGwas$Class == "Novel" & allEqtlGwas$Type=="ER"),]
allEqtl_novel = allEqtl_novel[!duplicated(allEqtl_novel$Feature),]
gr_all_novel = GRanges(ss(allEqtl_novel$Coordinates,":"),
	IRanges(as.numeric(ss(ss(allEqtl_novel$Coordinates,":",2), "-")),
			as.numeric(ss(ss(ss(allEqtl_novel$Coordinates,":",2), "-",2),"\\("))))
seqlevels(gr_all_novel) = names(stateCovDf	)

length(gr_all_novel)
sum(width(gr_all_novel))

# split by chr
allByChr = split(gr_all_novel, seqnames(gr_all_novel))
allByChr = allByChr[names(stateCovDf)]

bpOverlapByChr_all = mapply(function(regByChr, stateByChr) {
	cat(".")	
	enr = stateByChr[ranges(regByChr),]
	bg = stateByChr

	list(enrTab = sapply(enr, table),
		bgTab  = sapply(bg, table))
}, allByChr, stateCovDf,SIMPLIFY=FALSE)


# states in s1+s2
bpOverlapEnr = sapply(bpOverlapByChr_all, "[", "enrTab")
bpOverlapArray = array(unlist(bpOverlapEnr), 
	dim = c(nrow(bpOverlapEnr[[1]]), 
	ncol(bpOverlapEnr[[1]]), length(bpOverlapEnr)))
enrTab = apply(bpOverlapArray, 1:2, sum)
dimnames(enrTab) = list(rownames(bpOverlapEnr[[1]]),
	colnames(bpOverlapEnr[[1]]))

# states in rest of genome	
bpOverlapBg = sapply(bpOverlapByChr_all, "[", "bgTab")
bpOverlapArray = array(unlist(bpOverlapBg), 
	dim = c(nrow(bpOverlapBg[[1]]), 
	ncol(bpOverlapBg[[1]]), length(bpOverlapBg)))
bgTab = apply(bpOverlapArray, 1:2, function(x) sum(as.numeric(x)))
dimnames(bgTab) = list(rownames(bpOverlapBg[[1]]),
	colnames(bpOverlapBg[[1]]))

# take ratio
statTab_all = as.data.frame(t(prop.table(enrTab,2) / 
	prop.table(bgTab,2)))
statTab_all$Sample = dat$Standardized.Epigenome.name
statTab_all = statTab_all[,c(16,1:15)]

boxplot(log2(statTab_all[,-1]))

###############
# make plots ##
###############

# long, pgc
pgcLong = data.frame(State = rep(colnames(statTab_pgc)[-1],each=nrow(statTab_pgc)),
	Enrichment = as.numeric(as.matrix(statTab_pgc[,-1])),
	Sample = rep(rownames(statTab_pgc), times= ncol(statTab_pgc)-1))
	
pgcLong$log2Enrich = log2(pgcLong$Enrichment)
pgcLong$State = factor(pgcLong$State, levels=unique(pgcLong$State))
pgcLong$inBrain = as.numeric(pgcLong$Sample %in% dat$EID[brainIDs])+1

# long, all
allLong = data.frame(State = rep(colnames(statTab_all)[-1],each=nrow(statTab_all)),
	Enrichment = as.numeric(as.matrix(statTab_all[,-1])),
	Sample = rep(rownames(statTab_all), times= ncol(statTab_all)-1))
	
allLong$log2Enrich = log2(allLong$Enrichment)
allLong$State = factor(allLong$State, levels=unique(allLong$State))
allLong$inBrain = as.numeric(allLong$Sample %in% dat$EID[brainIDs])+1

# plot
pdf("GWAS_unannotated_roadmap_enrichments.pdf",h=5,w=10,useDingbats=FALSE)
par(mar=c(10,6,2,2),cex.axis=1.6,cex.lab=1.6,cex.main=1.6)
boxplot(log2Enrich ~ State, data = pgcLong, las=3, outline=FALSE,
	main = "SCZD GWAS",ylim=c(-4,4),
	ylab = "Coverage Ratio (log2)")
points(log2Enrich ~ jitter(as.numeric(State), amount=0.15), 
	data = pgcLong, 
	pch = 21, bg = ifelse(inBrain == 2, "red", "grey"),
	cex = ifelse(inBrain == 2, 1, 0.5))
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
	