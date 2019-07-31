library('bsseq')
library('devtools')
library('bumphunter')
library('GenomicRanges')
library('jaffelab')
library('SNPlocs.Hsapiens.dbSNP144.GRCh37')
library(matrixStats)

###################
# first overall ###
###################
load("../bumphunting/BSobj_bsseqSmooth_Neuron_minCov_3.Rdata")
BSobj_CpG = BSobj

## get mean cov
cpg_cov = getCoverage(BSobj_CpG)
mean_cpg_cov = rowMeans(cpg_cov)
quantile(mean_cpg_cov)

## get SNP data to get a sense of MAF
load("/dcl02/lieber/WGBS/Genotypes/Merged/Jaffe_DNAmDevelR21_Genotype_Imputed_1000G_Phase3_Merged.rda")
br = unique(colData(BSobj_CpG)$Brain.ID)
snp = snp[,br]
snpMap$MAF = rowSums(snp, na.rm=TRUE) / (2*rowSums(!is.na(snp)))

## change seqnames
r = rowRanges(BSobj_CpG)
seqlevels(r) = gsub("chr", "", seqlevels(r))
seqlevels(r)[which(seqlevels(r) == "M")] = "MT"
s_cpg = snpsByOverlaps(SNPlocs.Hsapiens.dbSNP144.GRCh37, r)
	
## what fraction of CpGs have snps?
length(s_cpg) / nrow(BSobj_CpG )
mean(s_cpg$alleles_as_ambig == "Y")

##check against our data
s_cpg$inSnpData = s_cpg$RefSNP_id %in% snpMap$name
table(s_cpg$inSnpData )
tt_cpg = table(s_cpg$alleles_as_ambig, s_cpg$inSnpData)
tt_cpg
prop.table(tt_cpg,1)

############################
###### DMPs and DMRs #######

## get CpG lists
cpgFiles = list.files("../BSobj_subsets/rda",
	pattern = "_CpG_",full=TRUE)
names(cpgFiles) = gsub(".Rdata", "", ss(cpgFiles, "/",4),fixed=TRUE)
## read in
cpgList = lapply(cpgFiles, function(x) {
	z = load(x)
	return(get(z))
})
sapply(cpgList, nrow)

## change seqnames
snpOverlaps = lapply(cpgList, function(x) {
	r = rowRanges(x)
	seqlevels(r) = gsub("chr", "", seqlevels(r))
	seqlevels(r)[1] = "MT"
	snpsByOverlaps(SNPlocs.Hsapiens.dbSNP144.GRCh37, r)
})
	
## what fraction of CpGs have snps?
lengths(snpOverlaps) / sapply(cpgList, nrow)

snpOverlaps = lapply(snpOverlaps, function(x) {
	x$inSnpData = x$RefSNP_id %in% snpMap$name
	x
})

sapply(snpOverlaps, function(x) mean(x$inSnpData))

## make a little table
stats = data.frame(numCpGs =sapply(cpgList, nrow))
stats$numSnps = lengths(snpOverlaps)
stats$numCTSnps = sapply(snpOverlaps, function(x) sum(x$alleles_as_ambig == "Y"))
stats$numWithMAF05 = sapply(snpOverlaps, function(x) sum(x$variable))
stats$prop = signif(stats$numWithMAF05 / stats$numCpGs,4)

####################
##### non CpG  #####
####################

## just need gr
load('/dcl01/lieber/ajaffe/lab/brain-epigenomics/bsseq/bsobj_by_chr/allChrs_postNatal_cleaned_nonCG_noHomogenate_highCov.Rdata',verbose=TRUE)
BSobj_nonCpG = BSobj
dim(BSobj_nonCpG) # confirm
## get coverage
cov_nonCpg = getCoverage(BSobj_nonCpG)
mean_nonCpg_cov = rowMeans(cov_nonCpg)
quantile(mean_nonCpg_cov)

## change seqnames
r = rowRanges(BSobj_nonCpG)
seqlevels(r) = gsub("chr", "", seqlevels(r))
seqlevels(r)[which(seqlevels(r) == "M")] = "MT"
s_nonCpg = snpsByOverlaps(SNPlocs.Hsapiens.dbSNP144.GRCh37, r)
	
## what fraction of CpGs have snps?
length(s_nonCpg) / nrow(BSobj_nonCpG )
mean(s_nonCpg$alleles_as_ambig == "Y")

##check against our data
s_nonCpg$inSnpData = s_nonCpg$RefSNP_id %in% snpMap$name
table(s_nonCpg$inSnpData)
mean(s_nonCpg$inSnpData)
tt_non = table(s_nonCpg$alleles_as_ambig, s_nonCpg$inSnpData)
tt_non
prop.table(tt_non,1)

######################
## pheno-associated ##



sapply(snpOverlaps_non, function(x) mean(x$variable))

## make a little table
statsNon = data.frame(numCpGs =sapply(nonList, nrow))
statsNon$numSnps = lengths(snpOverlaps_non)
statsNon$numCTSnps = sapply(snpOverlaps_non, function(x) sum(x$alleles_as_ambig == "Y"))
statsNon$numWithMAF05 = sapply(snpOverlaps_non, function(x) sum(x$variable))
statsNon$prop = signif(statsNon$numWithMAF05 / statsNon$numCpGs,4)
