#### get data
library(jaffelab)

#### read in brainspan data
p1 = read.delim("/legacy/nexsan2/disk3/ajaffe/BrainSpan/DNAm/1109_methylation_beta_values.txt",
	as.is=TRUE, skip=14, header=TRUE,row.names=1)
pd1 = read.delim("/legacy/nexsan2/disk3/ajaffe/BrainSpan/DNAm/1109_methylation_beta_values.txt",
	as.is=TRUE, header=FALSE,nrow = 14)

p2 = read.delim("/legacy/nexsan2/disk3/ajaffe/BrainSpan/DNAm/1110_methylation_beta_values.txt",
	as.is=TRUE, skip=14, header=TRUE,row.names=1)
pd2 = read.delim("/legacy/nexsan2/disk3/ajaffe/BrainSpan/DNAm/1110_methylation_beta_values.txt",
	as.is=TRUE, header=FALSE,nrow = 14)
	
p= cbind(p1,p2)
colnames(p) = ss(colnames(p), "\\.")

## pheno
pd = cbind(pd1, pd2[,-1])
pd = pd[-4,] # duplicate well
rownames(pd) = pd[,1]
pd = pd[,-1]
pd = t(pd)
pd = as.data.frame(pd, stringsAsFactors=FALSE)
rownames(pd) = NULL
pd$Plate = ss(pd$Plate, "\\.")
names(pd) = gsub(" ", "", names(pd))
pd$SampleID = gsub(" ", "_", pd$SampleID)

identical(gsub("X","", colnames(p)), pd$CompleteBarcode) # TRUE
colnames(p) = rownames(pd) = pd$SampleID

## fix up region
pd$Regions = ss(rownames(pd),"_",2)
pd$Regions[pd$Regions=="AIC"] = "A1C"
pd$Regions[pd$Regions=="MIC"] = "M1C"
pd$Regions[pd$Regions=="VIC"] = "V1C"
pd$Regions[pd$Regions=="SIC"] = "S1C"
pd$Regions[pd$Regions=="STS"] = "STR"
pd$Regions = factor(pd$Regions, levels = c("DFC","VFC","MFC",
	"OFC","M1C","S1C", "IPC", "A1C", "STC", "ITC", "V1C", "HIP",
	"AMY", "STR", "MD", "CBC"))
pd$Specimen.Code = ss(pd$SampleID, "_")

pheno = read.csv("/users/ajaffe/Lieber/Projects/450k/ECD2014/brainspan_pheno_meth.csv",as.is=TRUE)
tmp = as.numeric(ss(pheno$Age," ", 1))
tmp[ss(pheno$Age," ", 2)=="M"] = tmp[ss(pheno$Age," ", 2)=="M"]/12
pheno$Age= tmp

pd$Age = pheno$Age[match(pd$Specimen.Code, pheno$Specimen.Code)]
pd$Specimen.ID = pheno$Specimen.ID[match(pd$Specimen.Code, pheno$Specimen.Code)]

save(p,pd,file="brainspan_dnam_data_postnatal.rda",compress=TRUE)

## modeling age
library(limma)

## overall
mod = model.matrix(~Age + factor(Regions),data=pd)
fit = lmFit(p[,!is.na(pd$Age)], mod)
eb = ebayes(fit)

## just DFC
mod1 = model.matrix(~Age,data=pd[pd$Regions == "DFC",])
fit1 = lmFit(p[,rownames(mod1)], mod1)
eb1 = ebayes(fit1)

## output
outSpan = data.frame(ageChange_all = fit$coef[,2], 
	tstat_all = eb$t[,2], pval_all = eb$p[,2],
	ageChange_dfc = fit1$coef[,2], tstat_dfc = eb1$t[,2],
	pval_dfc = eb1$p[,2])
save(outSpan, file="brainspan_DNAm_age_stats_postnatal.rda")