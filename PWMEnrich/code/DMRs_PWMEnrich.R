library(PWMEnrich)
library(GenomicRanges)
library(PWMEnrich.Hsapiens.background)
library(BSgenome.Hsapiens.UCSC.hg19)
library(pheatmap)
library(RColorBrewer)

load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/DMR/DMR_objects.rda")


## Get sequences from different annotations in the correct format

for (i in 1:length(DMR)) { 
  rownames(DMR[[i]]) = DMR[[i]][,"regionID"]
  DMR[[i]][,"annotation"] = gsub("Other", "Intergenic", DMR[[i]][,"annotation"])
}

data = list(promoters = lapply(DMR, function(x) x[which(x$promoter=="Promoter" & x$sig=="FWER < 0.05" & x$width>29),]), # all DMRs that overlap a promoter
            intergenic = lapply(DMR, function(x) x[which(x$annotation=="Intergenic" & x$sig=="FWER < 0.05" & x$width>29),]), # all intergenic DMRs
            introns = lapply(DMR, function(x) x[which(x$annotation=="Intron" & x$sig=="FWER < 0.05" & x$width>29),]), #all DMRs that overlap exclusively introns
            all = lapply(DMR, function(x) x[which(x$sig=="FWER < 0.05" & x$width>29),]))
gr = lapply(data, function(x) lapply(x, makeGRangesFromDataFrame))
for (i in 1:length(gr)) { 
  for (j in 1:length(gr[[i]])) { 
    names(gr[[i]][[j]]) = data[[i]][[j]][,"regionID"] 
  } 
}
seq = lapply(gr, function(y) lapply(y, function(x) getSeq(Hsapiens, x)))


### Run PWMEnrich

#useBigMemoryPWMEnrich(TRUE)
registerCoresPWMEnrich(3)

# load the pre-compiled lognormal background computed using promoters
data(PWMLogn.hg19.MotifDb.Hsap)

promoters_int = motifEnrichment(promoters_seq$Interaction, PWMLogn.hg19.MotifDb.Hsap, verbose=F)
intergenic_int = motifEnrichment(intergenic_seq$Interaction, PWMLogn.hg19.MotifDb.Hsap, verbose=F)
introns_int = motifEnrichment(introns_seq$Interaction, PWMLogn.hg19.MotifDb.Hsap, verbose=F)
all_int = motifEnrichment(all_seq$Interaction, PWMLogn.hg19.MotifDb.Hsap, verbose=F)

save(promoters_int, intergenic_int, introns_int, all_int, file="/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/PWMEnrich/DMR_PWMEnrich_objects.rda")


## Split by direction

promoters_split = intergenic_split = introns_split = all_split = list()
split = lapply(data, function(y) unlist(lapply(y, function(x) split(x, x$Dir)), recursive = F))
gr = lapply(split, function(x) lapply(x, makeGRangesFromDataFrame))
for (i in 1:length(gr)) { 
  for (j in 1:length(gr[[i]])) { 
    names(gr[[i]][[j]]) = split[[i]][[j]][,"regionID"] 
  } 
}
seq = lapply(gr, function(y) lapply(y, function(x) getSeq(Hsapiens, x)))

all_split$CellType.pos = motifEnrichment(seq$all$CellType.pos, PWMLogn.hg19.MotifDb.Hsap, verbose=F)
all_split$CellType.neg = motifEnrichment(seq$all$CellType.neg, PWMLogn.hg19.MotifDb.Hsap, verbose=F)
all_split$Age.pos = motifEnrichment(seq$all$Age.pos, PWMLogn.hg19.MotifDb.Hsap, verbose=F)
all_split$Age.neg = motifEnrichment(seq$all$Age.neg, PWMLogn.hg19.MotifDb.Hsap, verbose=F)
all_split$Interaction.pos = motifEnrichment(seq$all$Interaction.pos, PWMLogn.hg19.MotifDb.Hsap, verbose=F)
all_split$Interaction.neg = motifEnrichment(seq$all$Interaction.neg, PWMLogn.hg19.MotifDb.Hsap, verbose=F)

promoters_split$CellType.pos = motifEnrichment(seq$promoters$CellType.pos, PWMLogn.hg19.MotifDb.Hsap, verbose=F)
promoters_split$CellType.neg = motifEnrichment(seq$promoters$CellType.neg, PWMLogn.hg19.MotifDb.Hsap, verbose=F)
promoters_split$Age.pos = motifEnrichment(seq$promoters$Age.pos, PWMLogn.hg19.MotifDb.Hsap, verbose=F)
promoters_split$Age.neg = motifEnrichment(seq$promoters$Age.neg, PWMLogn.hg19.MotifDb.Hsap, verbose=F)
promoters_split$Interaction.pos = motifEnrichment(seq$promoters$Interaction.pos, PWMLogn.hg19.MotifDb.Hsap, verbose=F)
promoters_split$Interaction.neg = motifEnrichment(seq$promoters$Interaction.neg, PWMLogn.hg19.MotifDb.Hsap, verbose=F)

intergenic_split$CellType.pos = motifEnrichment(seq$intergenic$CellType.pos, PWMLogn.hg19.MotifDb.Hsap, verbose=F)
intergenic_split$CellType.neg = motifEnrichment(seq$intergenic$CellType.neg, PWMLogn.hg19.MotifDb.Hsap, verbose=F)
intergenic_split$Age.pos = motifEnrichment(seq$intergenic$Age.pos, PWMLogn.hg19.MotifDb.Hsap, verbose=F)
intergenic_split$Age.neg = motifEnrichment(seq$intergenic$Age.neg, PWMLogn.hg19.MotifDb.Hsap, verbose=F)
intergenic_split$Interaction.pos = motifEnrichment(seq$intergenic$Interaction.pos, PWMLogn.hg19.MotifDb.Hsap, verbose=F)
intergenic_split$Interaction.neg = motifEnrichment(seq$intergenic$Interaction.neg, PWMLogn.hg19.MotifDb.Hsap, verbose=F)

introns_split$CellType.pos = motifEnrichment(seq$intergenic$CellType.pos, PWMLogn.hg19.MotifDb.Hsap, verbose=F)
introns_split$CellType.neg = motifEnrichment(seq$intergenic$CellType.neg, PWMLogn.hg19.MotifDb.Hsap, verbose=F)
introns_split$Age.pos = motifEnrichment(seq$intergenic$Age.pos, PWMLogn.hg19.MotifDb.Hsap, verbose=F)
introns_split$Age.neg = motifEnrichment(seq$intergenic$Age.neg, PWMLogn.hg19.MotifDb.Hsap, verbose=F)
introns_split$Interaction.pos = motifEnrichment(seq$intergenic$Interaction.pos, PWMLogn.hg19.MotifDb.Hsap, verbose=F)
introns_split$Interaction.neg = motifEnrichment(seq$intergenic$Interaction.neg, PWMLogn.hg19.MotifDb.Hsap, verbose=F)


save(promoters_split, intergenic_split, introns_split, all_split, promoters_int, intergenic_int, introns_int, all_int,
     file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/PWMEnrich/DMR_PWMEnrich_objects.rda")


## Cluster Interaction DMRs by TF enrichment to see how the groups materialize

load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/PWMEnrich/DMR_PWMEnrich_objects.rda")

all_int_seqRep = list()
for (i in 1:length(all_int$sequences)) {
	all_int_seqRep[[i]] = sequenceReport(all_int, seq.id=i)
}
targ = all_int_seqRep[[1]]$target    
pvalMat = do.call(cbind, lapply(all_int_seqRep, function(x) x$p.value[match(targ, x$target)]))
rownames(pvalMat) = targ
colnames(pvalMat) = names(all_int$sequences)

lMat = -log10(pvalMat)

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/PWMEnrich/figures/Interaction_DMR_TFenrichment_heatmap.pdf",h=16,w=8)
sampleDists <- dist(t(lMat))
sampleDistMatrix <- as.matrix(sampleDists)
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists,
         col=colors, main="Euclidean Distance Between\nInteraction DMRs by TF Enrichment")
sampleDists <- dist(lMat)
sampleDistMatrix <- as.matrix(sampleDists)
pheatmap(sampleDistMatrix,clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists,
         col=colors, main="Euclidean Distance Between\nTFs by Enrichment in Interaction DMRs")
dev.off()
# too large to be informative




pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/PWMEnrich/figures/Interaction_DMR_TFenrichment_cluster.pdf", h = 5, w = 14)
##  cluster by TF
hc = hclust(dist(t(lMat)))
hc_cut = cutree(hc, h= 40)
palette(brewer.pal(12,"Paired"))
myplclust(hc, lab.col=hc_cut,xlab="",hang=0.05,cex=1.1)
# cluster by DMR
hc2 = hclust(dist(lMat))
hc2_cut = cutree(hc2, k= 20)
myplclust(hc2, lab.col=hc2_cut,xlab="",hang=0)
dev.off()

### go by cluster
sIndexes = splitit(hc2_cut)
nullgenes =  read.delim("/users/ajaffe/Lieber/Projects/450k/grant/ref_gene_hg19.txt", header=T,as.is=T)
genes = sig$annotation
goListByTF = lapply(sIndexes, function(x) dogo(genes[x],nullgenes[,2])[,-8])



pca = prcomp(lMat)
pcat = prcomp(t(lMat))
ggplot(data = d, aes_string(x = "PC2", y = "PC3", color = "group")) + 
  geom_point(size = 3) + 
  ylab(paste0("PC3: ", round(percentVar[3] * 100), "% variance")) + 
  xlab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))

pca= prcomp(t(zMat))
km = kmeans(t(zMat), 5)


## Test for Differential TF binding

pos = c(all_split[grep("pos", names(all_split))], promoters_split[grep("pos", names(promoters_split))],
        introns_split[grep("pos", names(introns_split))], intergenic_split[grep("pos", names(intergenic_split))])
neg = c(all_split[grep("neg", names(all_split))], promoters_split[grep("neg", names(promoters_split))],
        introns_split[grep("neg", names(introns_split))], intergenic_split[grep("neg", names(intergenic_split))])
names(pos) = c("allAge","allInt","allCT","promCT","promAge","promInt","intronsCT","intronsAge","intronsInt","interCT","interAge","interInt")
pos = pos[!names(pos) %in% c("intronsAge","interAge")]
names(neg) = names(pos)
seqpos = c(allAge = seq$all$Age.pos, allInt = seq$all$Interaction.pos, allCT = seq$all$CellType.pos,
           promCT = seq$promoters$CellType.pos, promAge = seq$promoters$Age.pos, promInt = seq$promoters$Interaction.pos,
           intronsCT = seq$introns$CellType.pos, intronsAge = seq$introns$Age.pos, intronsInt = seq$introns$Interaction.pos,
           interCT = seq$intergenic$CellType.pos, interAge = seq$intergenic$Age.pos, interInt = seq$intergenic$Interaction.pos)
seqneg = c(allAge = seq$all$Age.neg, allInt = seq$all$Interaction.neg, allCT = seq$all$CellType.neg,
           promCT = seq$promoters$CellType.neg, promAge = seq$promoters$Age.neg, promInt = seq$promoters$Interaction.neg,
           intronsCT = seq$introns$CellType.neg, intronsInt = seq$introns$Interaction.neg,
           interCT = seq$intergenic$CellType.neg, interInt = seq$intergenic$Interaction.neg)

TFdiff = list()
for (i in 1:length(pos)) {
  TFdiff[[i]] = motifDiffEnrichment(sequences1 = seqpos[[i]], sequences2 = seqneg[[i]],
                                    res1 = pos[[i]], res2 = neg[[i]], PWMLogn.hg19.MotifDb.Hsap, verbose=FALSE)
}
names(TFdiff) = names(pos)

save(TFdiff, promoters_split, intergenic_split, introns_split, all_split, promoters_int, intergenic_int, introns_int, all_int,
     file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/PWMEnrich/DMR_PWMEnrich_objects.rda")




