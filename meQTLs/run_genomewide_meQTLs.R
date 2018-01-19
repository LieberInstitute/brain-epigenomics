###
library(bsseq)
library(MatrixEQTL)
library(readr)
library(stringr)

## load DNAm data
load('/dcl01/lieber/ajaffe/lab/brain-epigenomics/bumphunting/BSobj_bsseqSmooth_Neuron_minCov_3.Rdata')
pd = pData(BSobj)
## get mean meth per DMR
meth <- getMeth(BSobj, type = 'raw')

## subset to just neurons
meth = meth[,pd$Cell.Type == "Neuron"]
methMap = granges(BSobj)
pd = pd[pd$Cell.Type == "Neuron",]

## load genotype data
load("/dcl01/lieber/WGBS/Genotypes/Merged/Jaffe_DNAmDevelR21_Genotype_Imputed_1000G_Phase3_Merged.rda")
snp = snp[,pd$Brain.ID]
mds = mds[pd$Brain.ID,]

### do PCA
oo = order(rowSds(meth),decreasing=TRUE)[1:1e6]
pca = prcomp(t(meth[oo,]))

pd$Sex = str_trim(pd$Sex)
mod = model.matrix(~mds$snpPC1 + mds$snpPC2 + mds$snpPC3 + pca$x[,1:5] + pd$Age + pd$Sex)
# covs = SlicedData$new(t(mod[,-1]))

## set up matrix eqtl run
# theSnps = SlicedData$new(as.matrix(snp))
# theSnps$ResliceCombined(sliceSize = 10000)

snpspos = snpMap[,c("SNP","CHR","POS")]
snpspos$CHR = paste0("chr",snpspos$CHR)
colnames(snpspos) = c("name","chr","pos")

### meth position
posMeth = as.data.frame(methMap)
posMeth$Name = paste0(posMeth$seqnames, ":", posMeth$start)
posMeth = posMeth[,c(6,1:3)]
names(posMeth)[2] = "chr"

rownames(meth) = posMeth$Name
# theMeth = SlicedData$new(as.matrix(meth))
# theMeth$ResliceCombined(sliceSize = 10000)

# meMeth = Matrix_eQTL_main(snps=theSnps, gene = theMeth, 
	# cvrt = covs, output_file_name.cis =  ".txt" ,
	# pvOutputThreshold.cis = 0.001, pvOutputThreshold=0,
	# snpspos = snpspos, genepos = posMeth, 
	# useModel = modelLINEAR,	cisDist=10000,
	# pvalue.hist = 100,min.pv.by.genesnp = TRUE)
# save(meMeth, file="meQTLs_neurons_cis.rda")

## check out meQTLs
meqtls = meMeth$cis$eqtls
meqtls$snps = as.character(meqtls$snps)
meqtls$gene = as.character(meqtls$gene)
colnames(meqtls)[2] = "CpG"

meqtls = meqtls[meqtls$FDR < 0.05,]

length(unique(meqtls$CpG))
length(unique(meqtls$snps))

meqtls$methChr = posMeth$chr[match(meqtls$CpG, posMeth$Name)]
meqtls$methPos = posMeth$start[match(meqtls$CpG, posMeth$Name)]
meqtls$snpPos = snpspos$pos[match(meqtls$snps, snpspos$name)]
meqtls$snpMinusMeth = meqtls$methPos - meqtls$snpPos