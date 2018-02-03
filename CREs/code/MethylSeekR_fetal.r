library("bsseq")
library("GenomicFeatures")
library("MethylSeekR")
library("BSgenome")
library("BSgenome.Hsapiens.UCSC.hg19")


# Set random number generator seed to make results reproducible

set.seed(123)

## load phenotype data

load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/bumphunting/BSobj_bsseqSmooth_Neuron_minCov_3_prenatal.Rdata")

pd = pData(BSobj)

pd$avg.Cov
#[1]  9.73  9.87 11.40 10.48 11.45 10.31 11.07  9.54 10.48 11.64 11.16 12.97
#[13] 12.14 13.02 15.04 16.28 16.59 12.54 11.15 15.82
# seems high enough to me (10 minimum for this program)


## load coordinates, methylation counts and total coverage for filtered CpGs

gr = granges(BSobj)
df = data.frame(gr)
cov = getCoverage(BSobj, type = "Cov")
m = getCoverage(BSobj, type = "M")
rm(BSobj,gr)


## Format info for methylseekR

CGPrenlist = list()
for (i in 1:ncol(cov)) {
	CGPrenlist[[i]] = data.frame(df[,c("seqnames","start","end")], T = cov[,i], M = m[,i])
}
names(CGPrenlist) = colnames(cov)
CGPrenlist = lapply(CGPrenlist, function(x) makeGRangesFromDataFrame(x, keep.extra.columns=T))
	# each element in the list are the counts for one sample
rm(df,cov,m)

save(CGPrenlist, file="/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/methylSeekR_objects_prenatal.rda")


## Plot the alpha values for one chromosome
	# α characterizes the distribution of methylation levels in sliding windows containing 100 consecutive CpGs along the genome

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/CREs/figures/AlphaDistribution_CG_prenatal.pdf")
lapply(CGPrenlist, function(x) plotAlphaDistributionOneChr(m= x, chr.sel="chr22", num.cores=3))
dev.off()


## Identify partially methylated domains (PMDs)

sLengths = seqlengths(Hsapiens)

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/CREs/figures/PMDsegments_methylSeekR_CG_prenatal.pdf")
PMDsegments.CGpren = lapply(CGPrenlist, function(x) segmentPMDs(m = x, chr.sel="chr22", seqLengths=sLengths, num.cores=3))
dev.off()
save(PMDsegments.CGpren, file="/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/PMDs_methylSeekR_prenatal.rda")


## Identify CpG Islands

session <- browserSession()
genome(session) <- "hg19" 
query <- ucscTableQuery(session, "cpgIslandExt") 
CpGislands.gr <- track(query) 
genome(CpGislands.gr) <- NA
CpGislands.gr <- suppressWarnings(resize(CpGislands.gr, 5000, fix="center"))


## Calculate FDR threshold for calling UMRs and LMRs

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/CREs/figures/calculateFDRthreshold_methylSeekR_CGprenatal.pdf")
stats.CGpren <- mapply(function(CG, PMD) calculateFDRs(m = CG, CGIs = CpGislands.gr, PMDs = PMD, num.cores=1), CGPrenlist, PMDsegments.CGpren)
dev.off()
statcg = list()
for (i in 1:length(CGPrenlist)) { statcg[[i]] =  stats.CGpren[,i] }
names(statcg) = names(CGPrenlist)
stats.CGpren = statcg

FDR.cutoff <- 10
m.sel <- 0.5

n.sel.CGpren = lapply(stats.CGpren, function(x) as.integer(names(x$FDRs[as.character(m.sel), ][x$FDRs[as.character(m.sel), ]<FDR.cutoff])[1]))

save(stats.CGpren, n.sel.CGpren, file="/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/setting_n_methylSeekR_prenatal.rda")


# Recalculate with limited PMDs to >100kb
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/PMDs_methylSeekR_100kbLimit.rda")

names(total) = gsub("Prenatal-","", names(total))
pmds = total[names(total) %in% names(CGPrenlist)]
identical(names(pmds), names(CGPrenlist))

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/CREs/figures/calculateFDRthreshold_methylSeekR_CG.100kb_prenatal.pdf")
stats.CG.100kb.pren <- mapply(function(CG, PMD) calculateFDRs(m = CG, CGIs = CpGislands.gr, PMDs = PMD, num.cores=2), CGPrenlist, pmds)
dev.off()
statcg = list()
for (i in 1:length(CGPrenlist)) { statcg[[i]] =  stats.CG.100kb.pren[,i] }
names(statcg) = names(CGPrenlist)
stats.CG.100kb.pren = statcg

n.sel.CG.100kb.pren = lapply(stats.CG.100kb.pren, function(x) as.integer(names(x$FDRs[as.character(m.sel), ][x$FDRs[as.character(m.sel), ]<FDR.cutoff])[1]))

save(stats.CG.100kb.pren, n.sel.CG.100kb.pren, stats.CGpren, n.sel.CGpren, file="/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/setting_n_methylSeekR_prenatal.rda")


## Calculate UMRs and LMRs

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/CREs/figures/UMRs_LMRs_methylSeekR_prenatal.pdf")
UMRLMRsegments.CGpren <- mapply(function(CG,n,PMD) segmentUMRsLMRs(m = CG, meth.cutoff = m.sel, nCpG.cutoff= n, 
					 PMDs = PMD, num.cores=1, myGenomeSeq = Hsapiens, seqLengths = sLengths), CGPrenlist, n.sel.CGpren, PMDsegments.CGpren)  
dev.off()

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/CREs/figures/UMRs_LMRs_methylSeekR.100kb_prenatal.pdf")
UMRLMRsegments.CG.100kb.pren <- mapply(function(CG,n,PMD) segmentUMRsLMRs(m = CG, meth.cutoff = m.sel, nCpG.cutoff= n, 
                                                                   PMDs = PMD, num.cores=2, myGenomeSeq = Hsapiens, seqLengths = sLengths), CGPrenlist, n.sel.CG.100kb.pren, pmds)  
dev.off()

load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/UMRs_LMRs_methylSeekR_prenatal.rda")
save(UMRLMRsegments.CGpren,UMRLMRsegments.CG.100kb.pren, file="/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/UMRs_LMRs_methylSeekR_prenatal.rda")


## Call DMVs_pren (parameters from Mo et al. 2015 supplement)

DMVs_pren = lapply(UMRLMRsegments.CGpren, function(x) x[which(x$type=="UMR" & x$pmeth <=0.15)])
DMVs_pren = lapply(DMVs_pren, function(x) x[which(width(x)>=5000)])

DMVs_100kb.pren = lapply(UMRLMRsegments.CG.100kb.pren, function(x) x[which(x$type=="UMR" & x$pmeth <=0.15)])
DMVs_100kb.pren = lapply(DMVs_100kb.pren, function(x) x[which(width(x)>=5000)])


## Call large hypo-DMRs

gr2 <- lapply(lapply(UMRLMRsegments.CGpren, sortSeqlevels), sort)
for (i in 1:length(gr2)) {
  write.table(gr2[[i]], file=paste0("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/hypoDMRspren_",names(UMRLMRsegments.CGpren)[i],".tab"), 
              sep="\t", row.names = F,quote = F,col.names = F) }
write.table(data.frame(names(UMRLMRsegments.CGpren)), 
            file="/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/hypoDMR_prenatal_IDs.txt", 
            sep="\t", row.names = F,quote = F,col.names = F)


# run merge_hypoDMRspren.sh

hypoDMRspren = lapply(as.list(c(names(split),names(UMRLMRsegments.CGpren))), function(x) 
  read.table(paste0("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/",x,"_merged_hypoDMRspren.tab"), 
             col.names = c("chromosome", "start","end")))
names(hypoDMRspren) = c(names(split),names(UMRLMRsegments.CGpren))
elementNROWS(hypoDMRspren)
hypoDMRspren = lapply(hypoDMRspren, makeGRangesFromDataFrame)
hypoDMRspren = lapply(hypoDMRspren, function(x) split(x, width(x)>=2000))

load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/DMVs_hypoDMRs_methylSeekR_prenatal.rda")
save(hypoDMRspren, DMVs_pren, DMVs_100kb.pren, file="/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/DMVs_hypoDMRs_methylSeekR_prenatal.rda")