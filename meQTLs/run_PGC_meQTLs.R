###
library(bsseq)
library(MatrixEQTL)
library(readr)
library(stringr)
library(jaffelab)

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
snpMap$chrpos = paste0("chr", snpMap$CHR, ":", snpMap$POS)

### do PCA
oo = order(rowSds(meth),decreasing=TRUE)[1:1e6]
pca = prcomp(t(meth[oo,]))

pd$Sex = str_trim(pd$Sex)
mod = model.matrix(~mds$snpPC1 + mds$snpPC2 + mds$snpPC3 + pca$x[,1:5] + pd$Age + pd$Sex)
covs = SlicedData$new(t(mod[,-1]))

################
## subset SNPs to PGC
################

### load SZ GWAS clumps ###
pgc =read.delim("/users/ajaffe/Lieber/Projects/RNAseq/theControlEqtl/pgc/daner_PGC_SCZ52_0513a.gz.p4.clump.areator.sorted.1mhc.txt",
	as.is=TRUE, header=TRUE)
tmp = unlist(strsplit(pgc$LD.friends.0.1..p0.001,","))
tmp2 = strsplit(pgc$LD.friends.0.1..p0.001,",")

theSnps = data.frame(name = c(pgc$SNP, ss(tmp, "\\(")),
	R2 = as.numeric(c(rep(1,nrow(pgc)), ss(ss(tmp, "\\(",2),"/"))),
	Dist = as.numeric(c(rep(0,nrow(pgc)), 
		gsub(")", "", ss(ss(tmp, "\\(",2),"/",2),fixed=TRUE))))

# map back
theSnps$hitIndex=c(1:nrow(pgc), rep(seq(along=tmp2), times=sapply(tmp2,length)))
theSnps$Genes = pgc$genes.6.50kb.dist2index.[theSnps$hitIndex]
theSnps$pvalue = pgc$P[theSnps$hitIndex]

# drop those that don't map, and add coordinate info
library(GenomicRanges)
load("/users/ajaffe/Lieber/Projects/RNAseq/theControlEqtl/rdas/granges_pgc2.rda")
pgcFullGR$riskAllele = ifelse(pgcFullGR$or > 1, 
	pgcFullGR$A1, pgcFullGR$A2) # risk allele

## match up
mm = match(theSnps$name, names(pgcFullGR))
theSnps = theSnps[!is.na(mm),]
pgcStats = pgcFullGR[mm[!is.na(mm)]]
theSnps$chrpos = paste0(seqnames(pgcStats), ":", start(pgcStats))
theSnps$riskAllele = pgcStats$riskAllele

## drop those not with R^2 > 0.6
theSnps = theSnps[which(theSnps$R2 > 0.6),]

## just signif
pgcSig = read.delim("/users/ajaffe/Lieber/Projects/RNAseq/theControlEqtl/pgc/pgc2_128loci.txt")
pgcSig = pgcSig[order(pgcSig$Rank),]
mm = match(pgcSig$Index_SNP, names(pgcFullGR))
pgcStatsSig = pgcFullGR[mm]
pgcSig$chrpos = paste0(seqnames(pgcStatsSig), 
	":", start(pgcStatsSig))
theSnps$finalHitIndex = match(theSnps$chrpos, pgcSig$chrpos)

## check for being in LIBD
theSnps$inLibd = theSnps$chrpos %in% snpMap$chrpos  | 
	theSnps$name %in% snpMap$name
table(theSnps$inLibd, !is.na(theSnps$finalHitIndex))

snpMapSub = snpMap[snpMap$chrpos %in% theSnps$chrpos | 
	snpMap$name %in% theSnps$name),]
	
## only index
finalPgcSnps = theSnps[!is.na(theSnps$finalHitIndex),]

snpMatch = match(finalPgcSnps$chrpos, snpMap$chrpos)
snpMapSub = snpMap[snpMatch[!is.na(snpMatch)],]
snpSub = as.matrix(snp[snpMatch[!is.na(snpMatch)],])
snpMapSub$finalHitIndex = finalPgcSnps$finalHitIndex[
	match(snpMapSub$chrpos, finalPgcSnps$chrpos)]
#####################
## set up matrix eqtl run
theSnps = SlicedData$new(as.matrix(snpSub))

snpspos = snpMapSub[,c("SNP","CHR","POS")]
snpspos$CHR = paste0("chr",snpspos$CHR)
colnames(snpspos) = c("name","chr","pos")

### meth position
posMeth = as.data.frame(methMap)
posMeth$Name = paste0(posMeth$seqnames, ":", posMeth$start)
posMeth = posMeth[,c(6,1:3)]
names(posMeth)[2] = "chr"

rownames(meth) = posMeth$Name
theMeth = SlicedData$new(as.matrix(meth))
theMeth$ResliceCombined(sliceSize = 10000)

meMeth = Matrix_eQTL_main(snps=theSnps, gene = theMeth, 
	cvrt = covs, output_file_name.cis =  ".txt" ,
	pvOutputThreshold.cis = 1, pvOutputThreshold=0,
	snpspos = snpspos, genepos = posMeth, 
	useModel = modelLINEAR,	cisDist=200000,
	pvalue.hist = 100,min.pv.by.genesnp = TRUE)
save(meMeth, file="meQTLs_neurons_pgcSnps_cis.rda")

### reannotate
meqtls = meMeth$cis$eqtls
meqtls$snps = as.character(meqtls$snps)
meqtls$gene = as.character(meqtls$gene)
colnames(meqtls)[2] = "CpG"

## add annotation info
meqtls$finalHitIndex = snpMapSub$finalHitIndex[match(meqtls$snps, snpMapSub$SNP)]
meqtls$methChr = posMeth$chr[match(meqtls$CpG, posMeth$Name)]
meqtls$methPos = posMeth$start[match(meqtls$CpG, posMeth$Name)]
meqtls$snpPos = snpspos$pos[match(meqtls$snps, snpspos$name)]
meqtls$snpMinusMeth = meqtls$methPos - meqtls$snpPos

pgcList = split(meqtls, meqtls$finalHitIndex)

table(sapply(pgcList, function(x) sum(x$FDR < 0.05)) > 0)

## make some plots
pdf("meQTL_stat_plots_pgc.pdf",w=12, h= 7)
par(mar=c(5,6,3,2), cex.axis=1.6,cex.lab=1.6,cex.main = 1.6)
for(i in seq(along=pgcList)) {
	x = pgcList[[i]]
	x = x[order(x$methPos),]
	x$medStat = runmed(x$statistic, k=15)
	plot(x$methPos, x$statistic,xlab=x$methChr[1], pch=21, bg="grey",cex=0.6,
		ylab = "T-stats", main = paste("PGC Region", names(pgcList)[i]))
	lines(x$methPos, x$medStat,col="blue")
	abline(v=x$snpPos[1],col="red")
}
dev.off()


length(pgcList)