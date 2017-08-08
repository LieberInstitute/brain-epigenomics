##

library(jaffelab)

## old
pdOld= read.csv("/dcl01/lieber/WGBS/LIBD_Data/phenoData/FullMasterList_clean.csv")
pdOld = pdOld[pdOld$Age < 0,]
pdOld$FastQC_path = paste0("/dcl01/lieber/WGBS/LIBD_Data/FastQC/Untrimmed/",
	pdOld$WGC.ID, "/", pdOld$WGC.ID, "_combined")
pdOld$Batch = "Old"
colnames(pdOld)[1:2] = c("SampleID", "BrNum")

## new
pdNew = read.csv("/dcl01/lieber/WGBS/Macrogen_RawData/Macrogen_rawData_samples.csv",as.is=TRUE)
colnames(pdNew)[1] = "SampleID"
pdNew = pdNew[match(pdOld$BrNum, pdNew$BrNum),]
pdNew$FastQC_path = paste0("/dcl01/lieber/WGBS/psychENCODE_FetalSamples/FastQC/Untrimmed/",
	pdNew$BrNum, "_DLPFC/", pdNew$SampleID)
pdNew$Batch = "New"

n = intersect(names(pdOld), names(pdNew))
pd = rbind(pdOld[,n], pdNew[,n])
N = nrow(pd)

#####################
## read in FastQC ###
#####################
flags = c("FQCbasicStats","perBaseQual","perTileQual","perSeqQual",
			"perBaseContent","GCcontent","Ncontent","SeqLengthDist",
			"SeqDuplication","OverrepSeqs","AdapterContent","KmerContent")
fastqcdata = c("SeqLength","percentGC","phred1","phred2","phred3","phred4",
				"phredGT30","phredGT35","Adapter1","Adapter2","Adapter3")
splitAt = function(x, pos) unname(split(x, cumsum(seq_along(x) %in% pos)))

#### Summary flags (PASS/WARN/FAIL) ####
qcFlagsR1 = paste0(pd$FastQC_path,"_R1_fastqc/summary.txt")
qcFlagsR2 = paste0(pd$FastQC_path,"_R2_fastqc/summary.txt")
		
R1 = lapply(qcFlagsR1, function(x) scan(x, what="character", sep="\n", quiet=TRUE, strip=TRUE) )
R2 = lapply(qcFlagsR2, function(x) scan(x, what="character", sep="\n", quiet=TRUE, strip=TRUE) )	
o1 = lapply(R1, function(x) ss(x, "\t"))	
o1 = matrix(unlist(o1), ncol = 12, byrow = TRUE)
o2 = lapply(R2, function(x) ss(x, "\t"))	
o2 = matrix(unlist(o2), ncol = 12, byrow = TRUE)

## combine
o3 = paste0(o1,"/",o2)
dim(o3) = c(N,12)

o3[o3=="PASS/PASS"] = "PASS"
o3[o3=="WARN/WARN"] = "WARN"
o3[o3=="FAIL/FAIL"] = "FAIL"
colnames(o3)=paste0(flags)
pd = cbind(pd,o3)

#### FastQC pd/data ####
rr = c("R1","R2")
for (i in c(1:2)) {
	qcData = paste0(pd$FastQC_path,"_", rr[i], "_fastqc/fastqc_data.txt")
			
	R = sapply(qcData, function(x) scan(x, what="character", sep="\n", quiet=TRUE, strip=TRUE) )	
	names(R) = pd$SAMPLE_ID
	## Split list into sublists of metric categories
	zz = lapply(R, function(x) splitAt(x, which(x==">>END_MODULE")+1))
	
	# sequence length
	seqlen = lapply(zz, function(x) x[[1]][9])
	seqlen = sapply(seqlen, function(x) ss(x, "\t", 2))
	# percent GC
	gcp = lapply(zz, function(x) x[[1]][10])
	gcp = sapply(gcp, function(x) ss(x, "\t", 2))
	
	# median phred scores (at roughly 1/4, 1/2, 3/4, and end of seq length)
	# get positions 
	len = round((length(zz[[1]][[2]])-3) / 4)
	pos = c(len+3, 2*len+3, 3*len+3, length(zz[[1]][[2]])-1)
	nameSuf = ss(zz[[1]][[2]][pos], "\t", 1)
	fastqcdata[3:6] = paste0("phred", nameSuf)
	phred = lapply(zz, function(x) x[[2]][pos])
	phred = lapply(phred, function(x) ss(x, "\t", 3))
	phred = matrix(unlist(phred), ncol=4,byrow=T)
	
	# proportion of reads above phred 30 and 35
	phred2 = lapply(zz, function(x) x[[4]][3:(length(x[[4]])-1)])
	phred2 = lapply(phred2, function(x) 
			data.frame(score=ss(x, "\t", 1),count=ss(x, "\t", 2)))
	phred2 = lapply(phred2, function(x) 
			data.frame(x, cumulRev = rev(cumsum(rev(as.numeric(levels(x$count))[x$count]))) ))
	phred2 = lapply(phred2, function(x) 
			data.frame(x, prop = x$cumulRev/x$cumulRev[1] ))
	phred2 = lapply(phred2, function(x) x[which(x$score%in%c(30,35)),4] )
	phred2 = matrix(unlist(phred2), ncol=2, byrow=T)
	
	# Illumina adapter content
	# get positions 
	len = round((length(zz[[1]][[11]])-3) / 5)
	pos = c(3*len+2, 4*len+2, length(zz[[1]][[11]])-1)
	nameSuf = ss(zz[[1]][[11]][pos], "\t", 1)
	fastqcdata[9:11] = paste0("Adapter", nameSuf)
	adap = lapply(zz, function(x) x[[11]][pos])
	adap = lapply(adap, function(x) ss(x, "\t", 2))
	adap = matrix(unlist(adap), ncol=3, byrow=T)
	adap = matrix(as.numeric(adap), ncol=3, byrow=F)
	
	combined = data.frame(SeqLen=unlist(seqlen), GCprec=unlist(gcp), phred, phred2, adap)
	rownames(combined)=NULL
	names(combined) = paste0(fastqcdata,"_R",i)
	pd = cbind(pd,combined)
}

pdFetal = pd
facIndex = which(sapply(pdFetal,class)=="factor")
for(i in facIndex) pdFetal[,i] = as.character(pdFetal[,i])
for(i in 17:ncol(pdFetal)) pdFetal[,i] = as.numeric(pdFetal[,i])
save(pdFetal, file="rdas/fetal_qc_metrics_oldAndNew.rda")

boxplot(phred150_R2~factor(Batch),data=pdFetal)
boxplot(Adapter138_R2~factor(Batch),data=pdFetal)

# split back
pdOld = pdFetal[pdFetal$Batch == "Old",]
pdNew = pdFetal[pdFetal$Batch == "New",]
identical(pdOld$BrNum, pdNew$BrNum)

plot(pdOld$phred150_R2, pdNew$phred150_R2)
plot(pdOld$Adapter138_R2, pdNew$Adapter138_R2)
plot(pdOld$phredGT35_R2, pdNew$phredGT35_R2)
plot(pdOld$phredGT35_R2, pdNew$phredGT35_R2)
plot(pdOld$phredGT35_R1, pdNew$phredGT35_R1)
plot(pdOld$"Adapter110-111_R2", pdNew$"Adapter110-111_R2")