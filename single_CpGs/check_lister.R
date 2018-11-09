###
library(bsseq)
library(GEOquery)
library(recount)
library(jaffelab)

## load data
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/bsseq/lister/BSobj_matched_lister_minCov_3.Rdata")

pheno = pData(getGEO("GSE47966")[[1]])
pd = pheno[colnames(BSobj),c(1,2,8:13)]

for(i in 1:ncol(pd)) pd[,i] = as.character(pd[,i])
colnames(pd)[5:8] = ss(as.character(pd[1,5:8]), ": ")
for(i in 5:8) pd[,i] = ss(pd[,i], ": ", 2)
## fix last row
pd$gender[13] = "female"
pd$age[13] = NA
names(pd)[6] = "cellType"

## redo age
pd$ageNumeric = ss(pd$age, " ")
pd$ageNumeric = as.numeric(pd$ageNumeric)
pd$ageNumeric[grep("week", pd$age)] = pd$ageNumeric[grep("week", pd$age)]/52
pd$ageNumeric[grep("day", pd$age)] = pd$ageNumeric[grep("week", pd$age)]/365

pd$cellType = ss(pd$cellType, " ")
save(pd, file = "lister_phenotype.rda")

meth = getMeth(BSobj, type = 'raw')

#######################
## cell type dmrs ####
########################

library(genefilter)

dmrs_cell = read.csv("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/DMR/CellType/CellType_DMRs_fwer_0.05.csv",
	as.is=TRUE)
cellInds = mapply(function(s,e) s:e, dmrs_cell$indexStart, 
	dmrs_cell$indexEnd)

pdCell = pd[pd$cellType %in% c("neurons","non-neurons"),]
methCell = meth[,pd$cellType %in% c("neurons","non-neurons")]
ttCell = rowttests(methCell, factor(pdCell$cellType))

meanStats_cell = as.data.frame(t(sapply(cellInds, function(ii) {
	c(value = mean(ttCell$dm[ii], na.rm=TRUE), 
		area = sum(ttCell$dm[ii], na.rm=TRUE), 
		numCpGs = sum(!is.na(ttCell$stat[ii])),
		stat = mean(ttCell$stat[ii], na.rm=TRUE), 
		pvalue = mean(ttCell$p.value[ii], na.rm=TRUE))
})))

## validate
pdf("plots/libd_vs_lister_age.pdf")
par(mar=c(5,6,4,2),cex.axis=2,cex.lab=2,cex.main=1.7)
plot(dmrs_cell$value, meanStats_cell$value,
	xlab = "Discovery", ylab="Replication",
	pch = 21, bg="grey",cex=0.8,
	main = "Mean DNAm Difference\nNeuron vs Glia")
abline(h=0,v=0,col="blue")

plot(sign(dmrs_cell$value)*dmrs_cell$area, meanStats_cell$area,
	xlab = "Discovery", ylab="Replication",
	pch = 21, bg="grey",cex=0.8,
	main = "Area DNAm Difference\nNeuron vs Glia")
abline(0,1,col="blue")
dev.off()

cor(dmrs_cell$value, meanStats_cell$value,use="comp")
table(dmrs_cell$value > 0, meanStats_cell$value > 0)

#######################
## age dmrs ####
########################

library(limma)

dmrs_age = read.csv("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/DMR/Age/Age_DMRs_fwer_0.05.csv",
	as.is=TRUE)
ageInds = mapply(function(s,e) s:e, dmrs_age$indexStart, 
	dmrs_age$indexEnd)

pdAge = pd[!is.na(pd$ageNumeric),]
methAge = meth[,!is.na(pd$ageNumeric)]

meanStats_age = as.data.frame(t(sapply(ageInds, function(ii) {
	fitAge = lmFit(methAge[ii,], model.matrix(~pdAge$ageNumeric))

	c(value = mean(fitAge$coef[,2], na.rm=TRUE),
		area = sum(fitAge$coef[,2], na.rm=TRUE), 
		numCpGs = sum(!is.na(fitAge$coef[,2])))
})))

## validate
pdf("plots/libd_vs_lister_age.pdf")
par(mar=c(5,6,4,2),cex.axis=2,cex.lab=2,cex.main=1.7)
plot(10*dmrs_age$value,10*meanStats_age$value ,
	xlab = "Discovery", ylab="Replication",
	pch = 21, bg="grey",cex=0.8,
	main = "Mean DNAm Change/Decade")
abline(h=0,v=0,col="blue")
dev.off()