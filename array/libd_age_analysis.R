######

library(minfi)
library(limma)
library(RColorBrewer)
library(jaffelab)

## read in data
pd = read.csv("/users/ajaffe/Lieber/Projects/450k/full_sample_sheet_450k_brain.csv",
	header=TRUE,as.is=TRUE)

## keep only SZ and Control, CAUC+AA
pd = pd[pd$Race %in% c("CAUC", "AA") & 
	pd$Dx %in% c("Schizo","Control"),]
	
## subset from 0 to 24
pd = pd[pd$Age > 0 & pd$Age < 24,]

# update path as theyre on nexsan now
pd$BasePath=gsub("DATA", "/legacy/nexsan2/disk3/ajaffe/450k/BrainData", pd$BasePath)

# read in data	
RGset = read.metharray(pd$BasePath)

# add control genes
controlProbes = minfi:::.extractFromRGSet450k(RGset)
negControlPCs = prcomp(t(log2(rbind(controlProbes$greenControls$NEGATIVE,
	controlProbes$redControls$NEGATIVE)+1)))$x[,1:4]
colnames(negControlPCs) = paste0("negControl_", 	
	colnames(negControlPCs))
pd = cbind(pd, negControlPCs)


## qc, to drop replicates
mset = mapToGenome(RGset)
qc = minfiQC(mset)
qcMat = as.data.frame(qc$qc)
pd = cbind(pd,qcMat)

### KEEP DUPLICATES IN FOR NOW BUT FLAG
# distance to centroid to flag duplicates
bIndexes = split0(pd$BrNum)
means = colMeans(pd[,c("mMed","uMed")])
bestQC= sapply(bIndexes, function(x) {
	if(length(x) > 1) {
		tmp = rowSums((pd[x,c("mMed","uMed")] - means)^2) 
	}	else tmp = 1
	x[which.min(tmp)]
})
pd$bestQC = rep(FALSE)
pd$bestQC[bestQC] = TRUE

## normalize
Mset = preprocessQuantile(RGset, merge=TRUE)
pData(Mset) = DataFrame(pd)

#######################
#### age modeling #####
#######################

library(limma)

## filter people
keepIndex = which(pData(Mset)$bestQC & 
	pData(Mset)$Gender == pData(Mset)$predictedSex & 
	pData(Mset)$Dx == "Control")
Mset = Mset[,keepIndex]
	
## extract data on remaining people
pd = as.data.frame(pData(Mset))
p = getBeta(Mset)
map = as.data.frame(rowData(Mset))
colnames(map)[1] = "chr"

## age modeling
mod = model.matrix(~Age, data=pd)
fit = lmFit(p,mod)
eb = ebayes(fit)
outLibd = rowRanges(Mset)
mcols(outLibd)$ageChange = fit$coef[,2]
mcols(outLibd)$tstat = eb$t[,2]
mcols(outLibd)$pval = eb$p[,2]

save(outLibd, file="LIBD_DNAm_age_stats_postnatal_GRanges.rda")