library("bsseq")
library("GenomicFeatures")
library("MethylSeekR")
library("BSgenome")
library("BSgenome.Hsapiens.UCSC.hg19")


# Set random number generator seed to make results reproducible

set.seed(123)

## load phenotype data

load('/dcl01/lieber/ajaffe/lab/brain-epigenomics/bsseq/lister/BSobj_lister_minCov_3.Rdata')


## load coordinates, methylation counts and total coverage for filtered CpGs

gr = granges(BSobj)
df = data.frame(gr)
cov = getCoverage(BSobj, type = "Cov")
m = getCoverage(BSobj, type = "M")
rm(BSobj,gr)
pd = c("GSM1163695" = "fetal", "GSM1164631" = "16yr", "GSM1164632" = "25yr", "GSM1164630" = "12yr", "GSM1166274" = "5yr",
       "GSM1167004" = "35do", "GSM1167005" = "2yr", "GSM1173773" = "53yr_NeuN_pos", "GSM1173774" = "53yr_NeuN_neg",
       "GSM1173775" = "55yr_tissue", "GSM1173776" = "55yr_NeuN_pos", "GSM1173777" = "55yr_NeuN_neg", "GSM1173778" = "hues6")

# focus on sorted samples
cov = cov[,names(pd[grep("NeuN", pd)])]
m = m[,names(pd[grep("NeuN", pd)])]

## Format info for methylseekR

CGlist = list()
for (i in 1:ncol(cov)) {
  CGlist[[i]] = data.frame(df[,c("seqnames","start","end")], T = cov[,i], M = m[,i])
}
names(CGlist) = colnames(cov)
CGlist = lapply(CGlist, function(x) makeGRangesFromDataFrame(x, keep.extra.columns=T))
# each element in the list are the counts for one sample
rm(dr,cov,m)

save(CGlist, pd, file="/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/methylSeekR_objects_sortedLister.rda")


## Plot the alpha values for one chromosome
# α characterizes the distribution of methylation levels in sliding windows containing 100 consecutive CpGs along the genome

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/CREs/figures/AlphaDistribution_CG_sortedLister.pdf")
lapply(CGlist, function(x) plotAlphaDistributionOneChr(m= x, chr.sel="chr22", num.cores=1))
dev.off()


## Identify partially methylated domains (PMDs)

sLengths = seqlengths(Hsapiens)

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/CREs/figures/PMDsegments_methylSeekR_CG_sortedLister.pdf")
PMDsegments.CG.sortedLister = lapply(CGlist, function(x) segmentPMDs(m = x, chr.sel="chr22", seqLengths=sLengths, num.cores=2))
dev.off()

save(PMDsegments.CG.sortedLister, file="/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/PMDs_methylSeekR_sortedLister.rda")


## Identify CpG Islands

session <- browserSession()
genome(session) <- "hg19" 
query <- ucscTableQuery(session, "cpgIslandExt") 
CpGislands.gr <- track(query) 
genome(CpGislands.gr) <- NA
CpGislands.gr <- suppressWarnings(resize(CpGislands.gr, 5000, fix="center"))


## Calculate FDR threshold for calling UMRs and LMRs

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/CREs/figures/calculateFDRthreshold_methylSeekR_CG_sortedLister.pdf")
stats.CG <- mapply(function(CG, PMD) calculateFDRs(m = CG, CGIs = CpGislands.gr, PMDs = PMD, num.cores=1), CGlist, PMDsegments.CG.sortedLister)
dev.off()

statcg.sortedLister = list()
for (i in 1:length(CGlist)) { statcg.sortedLister[[i]] =  stats.CG[,i] }
names(statcg.sortedLister) = names(CGlist)
stats.CG.sortedLister = statcg.sortedLister

FDR.cutoff <- 10
m.sel <- 0.5

n.sel.CG.sortedLister = lapply(stats.CG.sortedLister, function(x) as.integer(names(x$FDRs[as.character(m.sel), ][x$FDRs[as.character(m.sel), ]<FDR.cutoff])[1]))

save(stats.CG.sortedLister, n.sel.CG.sortedLister, file="/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/setting_n_methylSeekR_sortedLister.rda")


## Calculate UMRs and LMRs

UMRLMRsegments.CG.sortedLister <- mapply(function(CG,n,PMD) segmentUMRsLMRs(m = CG, meth.cutoff = m.sel, nCpG.cutoff= n, 
                                                               PMDs = PMD, num.cores=1, myGenomeSeq = Hsapiens, seqLengths = sLengths), CGlist, n.sel.CG.sortedLister, PMDsegments.CG.sortedLister)  

## Call DMVs (parameters from Mo et al. 2015 supplement)

DMVs.sortedLister = lapply(UMRLMRsegments.CG.sortedLister, function(x) x[which(x$type=="UMR" & x$pmeth <=0.15)])
DMVs.sortedLister = lapply(DMVs.sortedLister, function(x) x[which(width(x)>=5000)])

save(UMRLMRsegments.CG.sortedLister,DMVs.sortedLister,
     file="/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/UMRs_LMRs_DMVs_methylSeekR_sortedLister.rda")
