library(PWMEnrich)
library(GenomicRanges)
library(PWMEnrich.Hsapiens.background)
library(BSgenome.Hsapiens.UCSC.hg19)

load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/DMR/DMR_objects.rda")

setwd("/dcl01/lieber/ajaffe/lab/brain-epigenomics/PWMEnrich")


## Get sequences from different annotations in the correct format

for (i in 1:length(DMR)) { 
  rownames(DMR[[i]]) = DMR[[i]][,"regionID"]
  DMR[[i]][,"annotation"] = gsub("Other", "Intergenic", DMR[[i]][,"annotation"])
}

promoters = lapply(DMR, function(x) x[which(x$promoter=="Promoter" & x$sig=="FWER < 0.05" & x$width>7),]) # all DMRs that overlap a promoter
intergenic = lapply(DMR, function(x) x[which(x$annotation=="Intergenic" & x$sig=="FWER < 0.05" & x$width>=7),]) # all intergenic DMRs
introns = lapply(DMR, function(x) x[which(x$annotation=="Intron" & x$sig=="FWER < 0.05" & x$width>=7),]) #all DMRs that overlap exclusively introns

promotersgr = lapply(promoters, makeGRangesFromDataFrame)
intergenicgr = lapply(intergenic, makeGRangesFromDataFrame)
intronsgr = lapply(introns, makeGRangesFromDataFrame)
allgr = lapply(DMR, function(x) makeGRangesFromDataFrame(x[which(x$sig=="FWER < 0.05" & x$width>=7),]))

promoters_seq = lapply(promotersgr, function(x) getSeq(Hsapiens, x))
intergenic_seq = lapply(intergenicgr, function(x) getSeq(Hsapiens, x))
introns_seq = lapply(intronsgr, function(x) getSeq(Hsapiens, x))
all_seq = lapply(allgr, function(x) getSeq(Hsapiens, x))


### Run PWMEnrich

useBigMemoryPWMEnrich(TRUE)
registerCoresPWMEnrich(5)

# load the pre-compiled lognormal background computed using promoters
data(PWMLogn.hg19.MotifDb.Hsap)

promoters_res = intergenic_res = introns_res = list()
promoters_res$CellType = motifEnrichment(promoters_seq$CellType, PWMLogn.hg19.MotifDb.Hsap)
promoters_res$Age = motifEnrichment(promoters_seq$Age, PWMLogn.hg19.MotifDb.Hsap)
promoters_res$Interaction = motifEnrichment(promoters_seq$Interaction, PWMLogn.hg19.MotifDb.Hsap)

intergenic_res$CellType = motifEnrichment(intergenic_seq$CellType, PWMLogn.hg19.MotifDb.Hsap)
intergenic_res$Age = motifEnrichment(intergenic_seq$Age, PWMLogn.hg19.MotifDb.Hsap)
intergenic_res$Interaction = motifEnrichment(intergenic_seq$Interaction, PWMLogn.hg19.MotifDb.Hsap)

introns_res$CellType = motifEnrichment(introns_seq$CellType, PWMLogn.hg19.MotifDb.Hsap)
introns_res$Age = motifEnrichment(introns_seq$Age, PWMLogn.hg19.MotifDb.Hsap)
introns_res$Interaction = motifEnrichment(introns_seq$Interaction, PWMLogn.hg19.MotifDb.Hsap)

all_res = list()
all_res$CellType = motifEnrichment(all_seq$CellType, PWMLogn.hg19.MotifDb.Hsap)
all_res$Age = motifEnrichment(all_seq$Age, PWMLogn.hg19.MotifDb.Hsap)
all_res$Interaction = motifEnrichment(all_seq$Interaction, PWMLogn.hg19.MotifDb.Hsap)


save(all_res, promoters_res, intergenic_res, introns_res, file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/PWMEnrich/DMR_PWMEnrich_objects")


## Split by direction

promoters_split = intergenic_split = introns_split = list()
promoters = unlist(lapply(promoters, function(x) split(x, x$Dir)), recursive = F)
promotersgr = lapply(promoters, makeGRangesFromDataFrame)
promoters_seq = lapply(promotersgr, function(x) getSeq(Hsapiens, x))
promoters_split$CellType.pos = motifEnrichment(promoters_seq$CellType.pos, PWMLogn.hg19.MotifDb.Hsap)
promoters_split$CellType.neg = motifEnrichment(promoters_seq$CellType.neg, PWMLogn.hg19.MotifDb.Hsap)
promoters_split$Age.pos = motifEnrichment(promoters_seq$Age.pos, PWMLogn.hg19.MotifDb.Hsap)
promoters_split$Age.neg = motifEnrichment(promoters_seq$Age.neg, PWMLogn.hg19.MotifDb.Hsap)
promoters_split$Interaction.pos = motifEnrichment(promoters_seq$Interaction.pos, PWMLogn.hg19.MotifDb.Hsap)
promoters_split$Interaction.neg = motifEnrichment(promoters_seq$Interaction.neg, PWMLogn.hg19.MotifDb.Hsap)

intergenic = unlist(lapply(intergenic, function(x) split(x, x$Dir)), recursive = F)
intergenicgr = lapply(intergenic, makeGRangesFromDataFrame)
intergenic_seq = lapply(intergenicgr, function(x) getSeq(Hsapiens, x))
intergenic_split$CellType.pos = motifEnrichment(intergenic_seq$CellType.pos, PWMLogn.hg19.MotifDb.Hsap)
intergenic_split$CellType.neg = motifEnrichment(intergenic_seq$CellType.neg, PWMLogn.hg19.MotifDb.Hsap)
intergenic_split$Age.pos = motifEnrichment(intergenic_seq$Age.pos, PWMLogn.hg19.MotifDb.Hsap)
intergenic_split$Age.neg = motifEnrichment(intergenic_seq$Age.neg, PWMLogn.hg19.MotifDb.Hsap)
intergenic_split$Interaction.pos = motifEnrichment(intergenic_seq$Interaction.pos, PWMLogn.hg19.MotifDb.Hsap)
intergenic_split$Interaction.neg = motifEnrichment(intergenic_seq$Interaction.neg, PWMLogn.hg19.MotifDb.Hsap)

introns = unlist(lapply(introns, function(x) split(x, x$Dir)), recursive = F)
intronsgr = lapply(introns, makeGRangesFromDataFrame)
introns_seq = lapply(intronsgr, function(x) getSeq(Hsapiens, x))
introns_split$CellType.pos = motifEnrichment(introns_seq$CellType.pos, PWMLogn.hg19.MotifDb.Hsap)
introns_split$CellType.neg = motifEnrichment(introns_seq$CellType.neg, PWMLogn.hg19.MotifDb.Hsap)
introns_split$Age.pos = motifEnrichment(introns_seq$Age.pos, PWMLogn.hg19.MotifDb.Hsap)
introns_split$Age.neg = motifEnrichment(introns_seq$Age.neg, PWMLogn.hg19.MotifDb.Hsap)
introns_split$Interaction.pos = motifEnrichment(introns_seq$Interaction.pos, PWMLogn.hg19.MotifDb.Hsap)
introns_split$Interaction.neg = motifEnrichment(introns_seq$Interaction.neg, PWMLogn.hg19.MotifDb.Hsap)

all = lapply(DMR, function(x) x[which(x$sig=="FWER < 0.05" & x$width>=7),])
all = unlist(lapply(all, function(x) split(x, x$Dir)), recursive = F)
allgr = lapply(all, makeGRangesFromDataFrame)
all_seq = lapply(allgr, function(x) getSeq(Hsapiens, x))
all_split$CellType.pos = motifEnrichment(all_seq$CellType.pos, PWMLogn.hg19.MotifDb.Hsap)
all_split$CellType.neg = motifEnrichment(all_seq$CellType.neg, PWMLogn.hg19.MotifDb.Hsap)
all_split$Age.pos = motifEnrichment(all_seq$Age.pos, PWMLogn.hg19.MotifDb.Hsap)
all_split$Age.neg = motifEnrichment(all_seq$Age.neg, PWMLogn.hg19.MotifDb.Hsap)
all_split$Interaction.pos = motifEnrichment(all_seq$Interaction.pos, PWMLogn.hg19.MotifDb.Hsap)
all_split$Interaction.neg = motifEnrichment(all_seq$Interaction.neg, PWMLogn.hg19.MotifDb.Hsap)

load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/PWMEnrich/DMR_PWMEnrich_objects")
save(promoters_split, intergenic_split, introns_split, all_split, promoters_res, intergenic_res, introns_res,
     file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/PWMEnrich/DMR_PWMEnrich_objects")


### Extract results: group reports

promoters_groupReport = groupReport(promoters_res$Age)
promoters_sequenceReport = sapply(seq_along(promoters_seq$Age.pos), 
                                  function(i) sequenceReport(promoters_split$Age.pos, seq.id=i))

targ = promoters_sequenceReport[[1]]$target    
pvalMat = sapply(promoters_sequenceReport, function(x)    {
  x$p.value[match(targ, x$target)]
})
rownames(pvalMat) = targ
colnames(pvalMat) = names(promoters_seq$Age.pos)



## set of TFs as whole
grp = motifRankingForGroup(res)
s = cut(grp, c(0,1e-100,1e-20,1e-10, 1e-4, 0.05, 1),include.lowest=TRUE)
table(s)
length(unique(names(grp)[which(as.numeric(s)==1)]))
length(unique(names(grp)[which(as.numeric(s)%in% 1:2)]))

unique(names(grp)[which(as.numeric(s)==1)])
unique(names(grp)[which(as.numeric(s)==2)])

## plots
pdf("plots/top_motifs_all.pdf", h=9, w= 12)

### p-value for each sequencing by motif
rn = names(motifRankingForSequence(res,1))
resList = vector("list",nrow(sig))
for(i in seq(along=resList)) {
  if(i %% 1000 == 0) cat(".")
  resList[[i]] = motifRankingForSequence(res,i)[rn]
}

pMat = do.call("rbind",resList)
save(grp,pMat,file="rdas/motif_summaries.rda")

#### keep unique
load("rdas/motif_summaries.rda")
table(duplicated(t(pMat)))
pMat = pMat[,!duplicated(t(pMat))]

lMat = -log10(pMat)

# pca = prcomp(lMat)
# pcat = prcomp(t(lMat))

pdf("heatmap.pdf",h=16,w=8)
heatmap(zMat)
dev.off()

km = kmeans(t(zMat), 5)

##  cluster by TF
hc = hclust(dist(t(lMat)))
hc_cut = cutree(hc, k= 10)
palette(brewer.pal(12,"Paired"))
pdf("plots/dmrs_cluster_byTF.pdf", h = 5, w = 14)
myplclust(hc, lab.col=hc_cut,xlab="",hang=0.05,cex=1.1)
dev.off()

# cluster by DMR
hc2 = hclust(dist(lMat))
hc2_cut = cutree(hc2, k= 10)

pdf("plots/dmrs_cluster_byDMR.pdf", h = 5, w = 14)
palette(brewer.pal(12,"Paired"))
myplclust(hc2, lab.col=hc2_cut,xlab="",hang=0)
dev.off()

### go by cluster
sIndexes = splitit(hc2_cut)
nullgenes =  read.delim("/users/ajaffe/Lieber/Projects/450k/grant/ref_gene_hg19.txt", header=T,as.is=T)
genes = sig$annotation
goListByTF = lapply(sIndexes, function(x) dogo(genes[x],nullgenes[,2])[,-8])


pca= prcomp(t(zMat))








