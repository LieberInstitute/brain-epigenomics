###	
library(jaffelab)
library(GenomicRanges)
library(bsseq)
library(bumphunter)

## load DMRs
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/DMR/DMR_objects.rda")

# load roadmap states
load("/dcl01/lieber/ajaffe/PublicData/EpigenomeRoadmap/ChromHMM/chromHMM_15state_coverageDF.rda")
dat = read.delim("/dcl01/lieber/ajaffe/PublicData/EpigenomeRoadmap/ChromHMM/jul2013.roadmapData.qc - Consolidated_EpigenomeIDs_summary_Table.tsv",
	as.is=TRUE)
colnames(dat)[2] = "EID"
dat = dat[dat$EID != "",]


## load data
load('/dcl01/lieber/ajaffe/lab/brain-epigenomics/bumphunting/BSobj_bsseqSmooth_Neuron_minCov_3.Rdata')

# Identify all CpG clusters in the genome
gr <- granges(BSobj)
cl = clusterMaker( chr = as.character(seqnames(gr)), 
                   pos = start(gr),   maxGap = 1000)
gr.clusters = split(gr, cl)
gr.clusters = unlist(range(gr.clusters))

# split out DMRs per background
dmrList = GRangesList(lapply(DMR, makeGRangesFromDataFrame, keep=TRUE))
bgWithDmr = endoapply(dmrList, function(x) {
	x = granges(x)
	m = c(x, gr.clusters)
	d = disjoin(m)
	d$inDMR = countOverlaps(d,x) > 0
	return(d)
})

################
### cell type DMRs 
cellDMRsByChr = split(bgWithDmr$CellType, seqnames(bgWithDmr$CellType))

## autosomal
cellDMRs_bpOverlapByChr = mapply(function(regByChr, stateByChr) {
	cat(".")
	inDMR = regByChr[regByChr$inDMR]
	outDMR = regByChr[!regByChr$inDMR]
	
	enr = stateByChr[ranges(inDMR),]
	bg = stateByChr[ranges(outDMR),]

	list(enrTab = sapply(enr, table),
		bgTab  = sapply(bg, table))
}, cellDMRsByChr[paste0("chr",1:22)], stateCovDf,SIMPLIFY=FALSE)

# states in cell DMRs
cellDMRs_bpOverlapEnr = sapply(cellDMRs_bpOverlapByChr, "[", "enrTab")
bpOverlapArray = array(unlist(cellDMRs_bpOverlapEnr), 
	dim = c(nrow(cellDMRs_bpOverlapEnr[[1]]), 
	ncol(cellDMRs_bpOverlapEnr[[1]]), length(cellDMRs_bpOverlapEnr)))
cellDMRs_enrTab = apply(bpOverlapArray, 1:2, sum)
dimnames(cellDMRs_enrTab) = list(rownames(cellDMRs_bpOverlapEnr[[1]]),
	colnames(cellDMRs_bpOverlapEnr[[1]]))

# states in rest of genome	
cellDMRs_bpOverlapBg = sapply(cellDMRs_bpOverlapByChr, "[", "bgTab")
bpOverlapArray = array(unlist(cellDMRs_bpOverlapBg), 
	dim = c(nrow(cellDMRs_bpOverlapBg[[1]]), 
	ncol(cellDMRs_bpOverlapBg[[1]]), length(cellDMRs_bpOverlapBg)))
cellDMRs_bgTab = apply(bpOverlapArray, 1:2, function(x) sum(as.numeric(x)))
dimnames(cellDMRs_bgTab) = list(rownames(cellDMRs_bpOverlapBg[[1]]),
	colnames(cellDMRs_bpOverlapBg[[1]]))

# take ratio
cellDMRs_statTab = as.data.frame(t(prop.table(cellDMRs_enrTab,2) / 
	prop.table(cellDMRs_bgTab,2)))
cellDMRs_statTab$Sample = dat$Standardized.Epigenome.name[match(rownames(cellDMRs_statTab), dat$EID)]
cellDMRs_statTab = cellDMRs_statTab[,c(16,1:15)]

cellDMRs_statTab["E073",]


##################
### age DMRs

ageDMRsByChr = split(bgWithDmr$Age, seqnames(bgWithDmr$Age))

## autosomal
ageDMRs_bpOverlapByChr = mapply(function(regByChr, stateByChr) {
	cat(".")
	inDMR = regByChr[regByChr$inDMR]
	outDMR = regByChr[!regByChr$inDMR]
	
	enr = stateByChr[ranges(inDMR),]
	bg = stateByChr[ranges(outDMR),]

	list(enrTab = sapply(enr, table),
		bgTab  = sapply(bg, table))
}, ageDMRsByChr[paste0("chr",1:22)], stateCovDf,SIMPLIFY=FALSE)

# states in cell DMRs
ageDMRs_bpOverlapEnr = sapply(ageDMRs_bpOverlapByChr, "[", "enrTab")
bpOverlapArray = array(unlist(ageDMRs_bpOverlapEnr), 
	dim = c(nrow(ageDMRs_bpOverlapEnr[[1]]), 
	ncol(ageDMRs_bpOverlapEnr[[1]]), length(ageDMRs_bpOverlapEnr)))
ageDMRs_enrTab = apply(bpOverlapArray, 1:2, sum)
dimnames(ageDMRs_enrTab) = list(rownames(ageDMRs_bpOverlapEnr[[1]]),
	colnames(ageDMRs_bpOverlapEnr[[1]]))

# states in rest of genome	
ageDMRs_bpOverlapBg = sapply(ageDMRs_bpOverlapByChr, "[", "bgTab")
bpOverlapArray = array(unlist(ageDMRs_bpOverlapBg), 
	dim = c(nrow(ageDMRs_bpOverlapBg[[1]]), 
	ncol(ageDMRs_bpOverlapBg[[1]]), length(ageDMRs_bpOverlapBg)))
ageDMRs_bgTab = apply(bpOverlapArray, 1:2, function(x) sum(as.numeric(x)))
dimnames(ageDMRs_bgTab) = list(rownames(ageDMRs_bpOverlapBg[[1]]),
	colnames(ageDMRs_bpOverlapBg[[1]]))

# take ratio
ageDMRs_statTab = as.data.frame(t(prop.table(ageDMRs_enrTab,2) / 
	prop.table(ageDMRs_bgTab,2)))
ageDMRs_statTab$Sample = dat$Standardized.Epigenome.name[match(rownames(ageDMRs_statTab), dat$EID)]
ageDMRs_statTab = ageDMRs_statTab[,c(16,1:15)]

ageDMRs_statTab["E073",]

##################
### interaction DMRs

intDMRsByChr = split(bgWithDmr$Interaction, seqnames(bgWithDmr$Interaction))

## autosomal
intDMRs_bpOverlapByChr = mapply(function(regByChr, stateByChr) {
	cat(".")
	inDMR = regByChr[regByChr$inDMR]
	outDMR = regByChr[!regByChr$inDMR]
	
	enr = stateByChr[ranges(inDMR),]
	bg = stateByChr[ranges(outDMR),]

	list(enrTab = sapply(enr, table),
		bgTab  = sapply(bg, table))
}, intDMRsByChr[paste0("chr",1:22)], stateCovDf,SIMPLIFY=FALSE)

# states in cell DMRs
intDMRs_bpOverlapEnr = sapply(intDMRs_bpOverlapByChr, "[", "enrTab")
bpOverlapArray = array(unlist(intDMRs_bpOverlapEnr), 
	dim = c(nrow(intDMRs_bpOverlapEnr[[1]]), 
	ncol(intDMRs_bpOverlapEnr[[1]]), length(intDMRs_bpOverlapEnr)))
intDMRs_enrTab = apply(bpOverlapArray, 1:2, sum)
dimnames(intDMRs_enrTab) = list(rownames(intDMRs_bpOverlapEnr[[1]]),
	colnames(intDMRs_bpOverlapEnr[[1]]))

# states in rest of genome	
intDMRs_bpOverlapBg = sapply(intDMRs_bpOverlapByChr, "[", "bgTab")
bpOverlapArray = array(unlist(intDMRs_bpOverlapBg), 
	dim = c(nrow(intDMRs_bpOverlapBg[[1]]), 
	ncol(intDMRs_bpOverlapBg[[1]]), length(intDMRs_bpOverlapBg)))
intDMRs_bgTab = apply(bpOverlapArray, 1:2, function(x) sum(as.numeric(x)))
dimnames(intDMRs_bgTab) = list(rownames(intDMRs_bpOverlapBg[[1]]),
	colnames(intDMRs_bpOverlapBg[[1]]))

# take ratio
intDMRs_statTab = as.data.frame(t(prop.table(intDMRs_enrTab,2) / 
	prop.table(intDMRs_bgTab,2)))
intDMRs_statTab$Sample = dat$Standardized.Epigenome.name[match(rownames(intDMRs_statTab), dat$EID)]
intDMRs_statTab = intDMRs_statTab[,c(16,1:15)]

intDMRs_statTab["E073",]

dlpfcStats = rbind(cellDMRs_statTab["E073",-1], ageDMRs_statTab["E073",-1],
	intDMRs_statTab["E073",-1])
rownames(dlpfcStats) = names(dmrList)
dlpfcStats = t(signif(dlpfcStats,3))

tmpList = vector("list", nrow(cellDMRs_statTab))
for(i in 1:nrow(cellDMRs_statTab)) {
	tmp =  rbind(cellDMRs_statTab[i,-1], ageDMRs_statTab[i,-1],
				intDMRs_statTab[i,-1])
	rownames(tmp) = names(dmrList)
	tmpList[[i]] = t(signif(tmp,3)) / dlpfcStats
}
names(tmpList) = rownames(cellDMRs_statTab)

########
# relative enrichments ###
############

###### cell type
cellVsDlpfc = sapply(tmpList, function(x) x[,1])

## general
rowSums(cellVsDlpfc > 2)
rowSums(cellVsDlpfc < 0.5)

## brain
rowMeans(cellVsDlpfc[,grep("Brain", dat$Standardized.Epigenome.name)] > 2)
rowMeans(cellVsDlpfc[,grep("Brain", dat$Standardized.Epigenome.name)] < 0.5)

## not brain
rowMeans(cellVsDlpfc[,-grep("Brain", dat$Standardized.Epigenome.name)] > 2)
rowMeans(cellVsDlpfc[,-grep("Brain", dat$Standardized.Epigenome.name)] < 0.5)

#### age
ageVsDlpfc = sapply(tmpList, function(x) x[,2])
rowMeans(ageVsDlpfc > 2)
rowMeans(ageVsDlpfc < 0.5)

## brain
rowMeans(ageVsDlpfc[,grep("Brain", dat$Standardized.Epigenome.name)] > 2)
rowMeans(ageVsDlpfc[,grep("Brain", dat$Standardized.Epigenome.name)] < 0.5)

## not brain
rowMeans(ageVsDlpfc[,-grep("Brain", dat$Standardized.Epigenome.name)] > 2)
rowMeans(ageVsDlpfc[,-grep("Brain", dat$Standardized.Epigenome.name)] < 0.5)

####### interaction
interVsDlpfc = sapply(tmpList, function(x) x[,3])
rowMeans(interVsDlpfc > 2)
rowMeans(interVsDlpfc < 0.5)

## brain
rowMeans(interVsDlpfc[,grep("Brain", dat$Standardized.Epigenome.name)] > 2)
rowMeans(interVsDlpfc[,grep("Brain", dat$Standardized.Epigenome.name)] < 0.5)

## not brain
rowMeans(interVsDlpfc[,-grep("Brain", dat$Standardized.Epigenome.name)] > 2)
rowMeans(interVsDlpfc[,-grep("Brain", dat$Standardized.Epigenome.name)] < 0.5)


atriumStats = rbind(cellDMRs_statTab["E104",-1], ageDMRs_statTab["E104",-1],
	intDMRs_statTab["E104",-1])
rownames(atriumStats) = names(dmrList)
atriumStats = t(signif(atriumStats,3))

### other tissues
otherList = list(cellDMRs_statTab, ageDMRs_statTab, intDMRs_statTab)

ageVsCell = ageDMRs_statTab[,-1] / cellDMRs_statTab[,-1]
cor(t(ageVsCell),t(ageVsCell["E073",]))