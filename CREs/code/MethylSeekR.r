library("bsseq")
library("GenomicFeatures")
library("MethylSeekR")
library("BSgenome")
library("BSgenome.Hsapiens.UCSC.hg19")


# Set random number generator seed to make results reproducible

set.seed(123)

## load phenotype data

load('/dcl01/lieber/ajaffe/lab/brain-epigenomics/bumphunting/BSobj_bsseqSmooth_Neuron_minCov_3.Rdata')

pd = pData(BSobj)
pd$Race[pd$Race== "CAUC "] <- 'CAUC'
pd$Sex[pd$Sex == " M"] <- 'M'
pd$RIN <- as.numeric(gsub(" ", "", pd$RIN))
pd$pg.DNA.nuclei.input <- as.numeric(pd$pg.DNA.nuclei.input)
pd$Reads <- as.numeric(pd$Reads)
pd$Percent.GreaterThan.Q30 <- as.numeric(pd$Percent.GreaterThan.Q30)


## load coordinates, methylation counts and total coverage for filtered CpGs

gr = granges(BSobj)
df = data.frame(gr)
cov = getCoverage(BSobj, type = "Cov")
m = getCoverage(BSobj, type = "M")
rm(BSobj,gr)


## Format info for methylseekR

CGlist = list()
for (i in 1:ncol(cov)) {
	CGlist[[i]] = data.frame(df[,c("seqnames","start","end")], T = cov[,i], M = m[,i])
}
names(CGlist) = colnames(cov)
CGlist = lapply(CGlist, function(x) makeGRangesFromDataFrame(x, keep.extra.columns=T))
	# each element in the list are the counts for one sample
rm(dr,cov,m)


## load methylation counts and total coverage for filtered CH

load('/dcl01/lieber/ajaffe/lab/brain-epigenomics/bsseq/bsobj_by_chr/allChrs_postNatal_cleaned_nonCG_noHomogenate_highCov.Rdata')

gr = granges(BSobj)
df = data.frame(gr)
cov = getCoverage(BSobj, type = "Cov")
m = getCoverage(BSobj, type = "M")
rm(BSobj,gr)


## Format info for methylseekR

CHlist = list()
for (i in 1:ncol(cov)) {
	CHlist[[i]] = data.frame(df[,c("seqnames","start","end")], T = cov[,i], M = m[,i])
}
names(CHlist) = colnames(cov)
CHlist = lapply(CHlist, function(x) makeGRangesFromDataFrame(x, keep.extra.columns=T))
	# each element in the list are the counts for one sample
rm(cov, df, m)

save(CGlist, CHlist, file="/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/methylSeekR_objects.rda")


## Plot the alpha values for one chromosome
	# α characterizes the distribution of methylation levels in sliding windows containing 100 consecutive CpGs along the genome

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/CREs/figures/AlphaDistribution_CG.pdf")
lapply(CGlist, function(x) plotAlphaDistributionOneChr(m= x, chr.sel="chr22", num.cores=1))
dev.off()

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/CREs/figures/AlphaDistribution_CH.pdf")
lapply(CHlist, function(x) plotAlphaDistributionOneChr(m= x, chr.sel="chr22", num.cores=1))
dev.off()


## Identify partially methylated domains (PMDs)

sLengths = seqlengths(Hsapiens)

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/CREs/figures/PMDsegments_methylSeekR_CG.pdf")
PMDsegments.CG = lapply(CGlist, function(x) segmentPMDs(m = x, chr.sel="chr22", seqLengths=sLengths, num.cores=2))
dev.off()

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/CREs/figures/PMDsegments_methylSeekR_CH.pdf")
PMDsegments.CH = lapply(CHlist, function(x) segmentPMDs(m = x, chr.sel="chr22", seqLengths=sLengths, num.cores=5))
dev.off()

save(PMDsegments.CG, PMDsegments.CH, file="/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/PMDs_methylSeekR.rda")


## Plot example PMDs

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/CREs/figures/PMDsegments_visualized_CG.pdf", width = 12, height = 12)
mapply(function(CG, PMD) plotPMDSegmentation(m = CG, segs = PMD), CGlist, PMDsegments.CG)
dev.off()

## Identify CpG Islands

session <- browserSession()
genome(session) <- "hg19" 
query <- ucscTableQuery(session, "cpgIslandExt") 
CpGislands.gr <- track(query) 
genome(CpGislands.gr) <- NA
CpGislands.gr <- suppressWarnings(resize(CpGislands.gr, 5000, fix="center"))


## Calculate FDR threshold for calling UMRs and LMRs

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/CREs/figures/calculateFDRthreshold_methylSeekR_CG.pdf")
stats.CG <- mapply(function(CG, PMD) calculateFDRs(m = CG, CGIs = CpGislands.gr, PMDs = PMD, num.cores=1), CGlist, PMDsegments.CG)
dev.off()
statcg = list()
for (i in 1:length(CGlist)) { statcg[[i]] =  stats.CG[,i] }
names(statcg) = names(CGlist)
stats.CG = statcg

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/CREs/figures/calculateFDRthreshold_methylSeekR_CH.pdf")
stats.CH <- lapply(CHlist, function(CH) calculateFDRs(m = CH, CGIs = CpGislands.gr, num.cores=5))
dev.off()
stats.CH

FDR.cutoff <- 10
m.sel <- 0.5

n.sel.CG = lapply(stats.CG, function(x) as.integer(names(x$FDRs[as.character(m.sel), ][x$FDRs[as.character(m.sel), ]<FDR.cutoff])[1]))
n.sel.CH = lapply(stats.CH, function(x) as.integer(names(x$FDRs[as.character(m.sel), ][x$FDRs[as.character(m.sel), ]<FDR.cutoff])[1]))

save(stats.CG, n.sel.CG, stats.CH, n.sel.CH, file="/media/Backup1_/amanda/NewFigures/setting_n_methylSeekR.rda")


# Recalculate with limited PMDs to >100kb
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/PMDs_methylSeekR_100kbLimit.rda")
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/setting_n_methylSeekR.rda")

names(total) = gsub("Postnatal-","", names(total))
pmds = total[names(total) %in% names(CGlist)]
identical(names(pmds), names(CGlist))

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/CREs/figures/calculateFDRthreshold_methylSeekR_CG.100kb.pdf")
stats.CG.100kb <- mapply(function(CG, PMD) calculateFDRs(m = CG, CGIs = CpGislands.gr, PMDs = PMD, num.cores=3), CGlist, pmds, SIMPLIFY = F)
dev.off()

n.sel.CG.100kb = lapply(stats.CG.100kb, function(x) as.integer(names(x$FDRs[as.character(m.sel), ][x$FDRs[as.character(m.sel), ]<FDR.cutoff])[1]))

save(stats.CG.100kb, n.sel.CG.100kb, stats.CG, n.sel.CG, stats.CH, n.sel.CH, file="/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/setting_n_methylSeekR.rda")


## Calculate UMRs and LMRs

UMRLMRsegments.CG <- mapply(function(CG,n,PMD) segmentUMRsLMRs(m = CG, meth.cutoff = m.sel, nCpG.cutoff= n, 
					 PMDs = PMD, num.cores=1, myGenomeSeq = Hsapiens, seqLengths = sLengths), CGlist, n.sel.CG, PMDsegments.CG)  

load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/UMRs_LMRs_methylSeekR.rda")

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/CREs/figures/UMRs_LMRs_methylseekr_100kbPMDs.pdf")
UMRLMRsegments.CG.100kb <- mapply(function(CG,n,PMD) segmentUMRsLMRs(m = CG, meth.cutoff = m.sel, nCpG.cutoff= n, 
                                                               PMDs = PMD, num.cores=1, myGenomeSeq = Hsapiens, seqLengths = sLengths), CGlist, n.sel.CG.100kb, pmds)  
dev.off()
save(UMRLMRsegments.CG,UMRLMRsegments.CG.100kb, file="/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/UMRs_LMRs_methylSeekR.rda")


## Call DMVs (parameters from Mo et al. 2015 supplement)

DMVs = lapply(UMRLMRsegments.CG, function(x) x[which(x$type=="UMR" & x$pmeth <=0.15)])
DMVs = lapply(DMVs, function(x) x[which(width(x)>=5000)])

load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/DMVs_hypoDMRs_methylSeekR.rda")

DMVs.100kb = lapply(UMRLMRsegments.CG.100kb, function(x) x[which(x$type=="UMR" & x$pmeth <=0.15)])
DMVs.100kb = lapply(DMVs.100kb, function(x) x[which(width(x)>=5000)])



## Call large hypo-DMRs

load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/DMR/DMR_objects.rda")

split = unlist(lapply(DMR, function(x) split(x, x$Dir)), recursive = F)
gr <- lapply(lapply(lapply(split, makeGRangesFromDataFrame), sortSeqlevels), sort)
gr2 <- lapply(lapply(UMRLMRsegments.CG, sortSeqlevels), sort)

for (i in 1:length(gr)) {
  write.table(gr[[i]], file=paste0("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/PWMEnrich/hypoDMRs_",names(split)[i],".tab"), 
              sep="\t", row.names = F,quote = F,col.names = F) }
for (i in 1:length(gr2)) {
  write.table(gr2[[i]], file=paste0("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/PWMEnrich/hypoDMRs_",names(UMRLMRsegments.CG)[i],".tab"), 
              sep="\t", row.names = F,quote = F,col.names = F) }
write.table(data.frame(c(names(split),names(UMRLMRsegments.CG))), 
            file="/dcl01/lieber/ajaffe/Amanda/sorted_nuclear_RNA/merged/hypoDMR_IDs.txt", 
            sep="\t", row.names = F,quote = F,col.names = F)

# run merge_hypoDMRs.sh

hypoDMRs = lapply(as.list(c(names(split),names(UMRLMRsegments.CG))), function(x) 
  read.table(paste0("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/PWMEnrich/",x,"_merged_hypoDMRs.tab"), 
             col.names = c("chromosome", "start","end")))
names(hypoDMRs) = c(names(split),names(UMRLMRsegments.CG))
elementNROWS(hypoDMRs)
hypoDMRs = lapply(hypoDMRs, makeGRangesFromDataFrame)
hypoDMRs = lapply(hypoDMRs, function(x) split(x, width(x)>=2000))

save(hypoDMRs, DMVs,DMVs.100kb, file="/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/DMVs_hypoDMRs_methylSeekR.rda")

              
              