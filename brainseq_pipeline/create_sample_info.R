###
library(jaffelab)
library(SummarizedExperiment)

## directory: /dcl01/lieber/ajaffe/lab/brain-epigenomics/brainseq_pipeline

################
#### phenoData
master = read.csv("/dcl01/lieber/WGBS/LIBD_Data/phenoData/FullMasterList_clean.csv")
colnames(master) = gsub("\\.", "\\_", colnames(master))

################
#### PolyA
load("/dcl01/ajaffe/data/lab/brainseq_phase1/count_data/dlpfc_polyA_brainseq_phase1_hg38_rseGene_merged_n732.rda")
polya = colData(rse_gene)[,1:9]
polya = polya[which(polya$BrNum %in% master$Brain_Num),]
polya$Experiment = "PolyA"
rownames(polya) = paste0(polya$RNum, "_", polya$Experiment)

man_polya = read.table("/dcl01/ajaffe/data/lab/brainseq_phase1/preprocessed_data/.samples_unmerged.manifest")

rnums = polya$RNum
matches = unique(grep(paste(rnums,collapse="|"), man_polya$V5, value=TRUE))
man_polya = man_polya[which(man_polya$V5 %in% matches),]
man_polya = man_polya[order(man_polya$V5),]
man_polya$V5 = ss(as.character(man_polya$V5),"_",1)

# check
nrow(polya)  # 41
length(unique(man_polya$V5))	 # 41	

write.table(man_polya, file="polyA_unstranded/samples.manifest", quote=FALSE, row.names=FALSE, sep="\t")

	

################					
#### RiboZero
load("/dcl01/ajaffe/data/lab/brainseq_phase2/count_data/dlpfc_ribozero_brainseq_phase2_hg38_rseGene_merged_n449.rda")
ribo = colData(rse_gene)[,c(60:67,1)]
ribo = ribo[which(ribo$BrNum %in% master$Brain_Num),]
ribo$Experiment = "RiboZero"
rownames(ribo) = paste0(ribo$RNum, "_", ribo$Experiment)

man_ribo = read.table("/dcl01/ajaffe/data/lab/brainseq_phase2/preprocessed_data/DLPFC_RiboZero/samples.manifest")

rnums = ribo$RNum
matches = unique(grep(paste(rnums,collapse="|"), man_ribo$V5, value=TRUE))
man_ribo = man_ribo[which(man_ribo$V5 %in% matches),]
man_ribo = man_ribo[order(man_ribo$V5),]
man_ribo$V5 = ss(as.character(man_ribo$V5),"_",1)
						
# check
nrow(ribo)  # 31
length(unique(man_ribo$V5))	# 31	

write.table(man_ribo, file="riboZero_stranded/samples.manifest", quote=FALSE, row.names=FALSE, sep="\t")
	


################
#### combine

sample_information = rbind(as.data.frame(polya),as.data.frame(ribo))
sample_information$RIN = as.numeric(sample_information$RIN)
sample_information$SAMPLE_ID = sapply(sample_information$SAMPLE_ID, paste, collapse=", ")

write.csv(sample_information, file="sample_information.csv")

# sample_manifest = rbind(man_polya, man_ribo)
# write.table(sample_manifest, file="samples.manifest", quote=FALSE, row.names=FALSE, sep="\t")









