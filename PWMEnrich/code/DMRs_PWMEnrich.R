library(PWMEnrich)
library(GenomicRanges)
library(PWMEnrich.Hsapiens.background)
library(BSgenome.Hsapiens.UCSC.hg19)
library(pheatmap)
library(RColorBrewer)

load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/DMR/DMR_objects.rda")


## Get sequences from different annotations in the correct format

for (i in 1:length(DMR)) { 
  DMR[[i]][,"annotation"] = gsub("Other", "Intergenic", DMR[[i]][,"annotation"])
}

data = list(promoters = lapply(DMR, function(x) x[which(x$promoter=="Promoter" & x$sig=="FWER < 0.05" & x$width>29),]), # all DMRs that overlap a promoter
            intergenic = lapply(DMR, function(x) x[which(x$annotation=="Intergenic" & x$sig=="FWER < 0.05" & x$width>29),]), # all intergenic DMRs
            introns = lapply(DMR, function(x) x[which(x$annotation=="Intron" & x$sig=="FWER < 0.05" & x$width>29),]), #all DMRs that overlap exclusively introns
            all = lapply(DMR, function(x) x[which(x$sig=="FWER < 0.05" & x$width>29),]))
gr = lapply(data, function(x) lapply(x, makeGRangesFromDataFrame))
gr = lapply(gr, function(x) lapply(x, reduce))
for (i in 1:length(gr)) { 
  for (j in 1:length(gr[[i]])) { 
    names(gr[[i]][[j]]) = paste0(as.character(seqnames(gr[[i]][[j]])), ":",as.character(start(gr[[i]][[j]])), "-",as.character(end(gr[[i]][[j]])))
  } 
}
seq = lapply(gr, function(y) lapply(y, function(x) getSeq(Hsapiens, x)))


### Run PWMEnrich

useBigMemoryPWMEnrich(TRUE)
registerCoresPWMEnrich(5)

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
gr = lapply(gr, function(x) lapply(x, reduce))
for (i in 1:length(gr)) { 
  for (j in 1:length(gr[[i]])) { 
    names(gr[[i]][[j]]) = paste0(as.character(seqnames(gr[[i]][[j]])), ":",as.character(start(gr[[i]][[j]])), "-",as.character(end(gr[[i]][[j]]))) 
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

introns_split$CellType.pos = motifEnrichment(seq$introns$CellType.pos, PWMLogn.hg19.MotifDb.Hsap, verbose=F)
introns_split$CellType.neg = motifEnrichment(seq$introns$CellType.neg, PWMLogn.hg19.MotifDb.Hsap, verbose=F)
introns_split$Age.pos = motifEnrichment(seq$introns$Age.pos, PWMLogn.hg19.MotifDb.Hsap, verbose=F)
introns_split$Age.neg = motifEnrichment(seq$introns$Age.neg, PWMLogn.hg19.MotifDb.Hsap, verbose=F)
introns_split$Interaction.pos = motifEnrichment(seq$introns$Interaction.pos, PWMLogn.hg19.MotifDb.Hsap, verbose=F)
introns_split$Interaction.neg = motifEnrichment(seq$introns$Interaction.neg, PWMLogn.hg19.MotifDb.Hsap, verbose=F)


save(promoters_split, intergenic_split, introns_split, all_split, promoters_int, intergenic_int, introns_int, all_int,
     file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/PWMEnrich/DMR_PWMEnrich_objects.rda")


## test each of the 6 kmeans clusters

load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/bumphunting/rda/limma_Neuron_CpGs_minCov_3_ageInfo_dmrs.Rdata")

names(dmrs) = dmrs$regionID
dmrs = dmrs[width(dmrs)>30]
int.groups = split(dmrs, dmrs$k6cluster_label)

intseq = lapply(int.groups, function(x) getSeq(Hsapiens, x))
int.kmeans = list()
int.kmeans$"1:G-N+" = motifEnrichment(intseq$"1:G-N+", PWMLogn.hg19.MotifDb.Hsap, verbose=F)
int.kmeans$"2:G0N+" = motifEnrichment(intseq$"2:G0N+", PWMLogn.hg19.MotifDb.Hsap, verbose=F)
int.kmeans$"3:G0N-" = motifEnrichment(intseq$"3:G0N-", PWMLogn.hg19.MotifDb.Hsap, verbose=F)
int.kmeans$"4:G+N0" = motifEnrichment(intseq$"4:G+N0", PWMLogn.hg19.MotifDb.Hsap, verbose=F)
int.kmeans$"5:G+N-" = motifEnrichment(intseq$"5:G+N-", PWMLogn.hg19.MotifDb.Hsap, verbose=F)
int.kmeans$"6:G-N0" = motifEnrichment(intseq$"6:G-N0", PWMLogn.hg19.MotifDb.Hsap, verbose=F)
    


## Test for Differential TF binding

pos = unlist(list(all = all_split[grep("pos", names(all_split))], prom = promoters_split[grep("pos", names(promoters_split))],
                  introns = introns_split[grep("pos", names(introns_split))], inter = intergenic_split[grep("pos", names(intergenic_split))]))
neg = unlist(list(all = all_split[grep("neg", names(all_split))], prom = promoters_split[grep("neg", names(promoters_split))],
                  introns = introns_split[grep("neg", names(introns_split))], inter = intergenic_split[grep("neg", names(intergenic_split))]))
pos = pos[order(names(pos))]
neg = neg[order(names(neg))]
names(pos) = gsub(".pos", "", names(pos))
names(neg) = gsub(".neg", "", names(neg))
pos = c(pos[names(pos)!="inter.Age"], pos[names(pos)=="inter.Age"])
match(names(pos), names(neg))

seq = unlist(seq, recursive = F)
seqpos = seq[grep("pos", names(seq))]
seqneg = seq[grep("neg", names(seq))]
seqpos = seqpos[order(names(seqpos))]
seqneg = seqneg[order(names(seqneg))]
names(seqpos) = gsub(".pos", "", names(seqpos))
names(seqneg) = gsub(".neg", "", names(seqneg))
seqpos = c(seqpos[names(seqpos)!="intergenic.Age"], seqpos[names(seqpos)=="intergenic.Age"])
match(names(seqpos), names(seqneg))
names(seqpos) = names(pos)
names(seqneg) = names(neg)

comps = Map(c, as.list(names(pos[names(pos)!="inter.Age"])), as.list(names(neg)))
names(comps) = c("allAge","allCT","allInt","interCT","interInt","intronsAge","intronsCT","intronsInt","promAge","promCT","promInt")
comps = c(comps,list(prom.intergenic.pos.CT = c("prom.CellType","inter.CellType"),
               prom.intron.pos.CT = c("prom.CellType","introns.CellType"),
               intron.intergenic.pos.CT = c("introns.CellType","inter.CellType"),
               prom.intergenic.neg.CT = c("prom.CellType", "inter.CellType"),
               prom.intron.neg.CT = c("prom.CellType", "introns.CellType"),
               intron.intergenic.neg.CT = c("introns.CellType", "inter.CellType"),
               prom.intergenic.pos.Age = c("prom.Age", "inter.Age"),
               prom.intron.pos.Age = c("prom.Age", "introns.Age"),
               intron.intergenic.pos.Age = c("introns.Age", "inter.Age"),
               prom.intron.neg.Age = c("prom.Age", "introns.Age"),
               prom.intergenic.pos.Int = c("prom.Interaction", "inter.Interaction"),
               prom.intron.pos.Int = c("prom.Interaction", "introns.Interaction"),
               intron.intergenic.pos.Int = c("introns.Interaction", "inter.Interaction"),
               prom.intergenic.neg.Int = c("prom.Interaction", "inter.Interaction"),
               prom.intron.neg.Int = c("prom.Interaction", "introns.Interaction"),
               intron.intergenic.neg.Int = c("introns.Interaction", "inter.Interaction")))

TFdiff = vector("list", length=length(comps))
for (i in 1:length(comps)) {
  if (i %in% c(1:11)) { 
    TFdiff[[i]] = motifDiffEnrichment(sequences1 = seqpos[[comps[[i]][1]]], sequences2 = seqneg[[comps[[i]][2]]],
                                      res1 = pos[[comps[[i]][1]]], res2 = neg[[comps[[i]][2]]], PWMLogn.hg19.MotifDb.Hsap, verbose=FALSE) } else
  if (i %in% c(12:14,18:20,22:24)) {
    TFdiff[[i]] = motifDiffEnrichment(sequences1 = seqpos[[comps[[i]][1]]], sequences2 = seqpos[[comps[[i]][2]]],
                                      res1 = pos[[comps[[i]][1]]], res2 = pos[[comps[[i]][2]]], PWMLogn.hg19.MotifDb.Hsap, verbose=FALSE) } else
  if (i %in% c(15:17,21,25:27)) {
    TFdiff[[i]] = motifDiffEnrichment(sequences1 = seqneg[[comps[[i]][1]]], sequences2 = seqneg[[comps[[i]][2]]],
                                      res1 = neg[[comps[[i]][1]]], res2 = neg[[comps[[i]][2]]], PWMLogn.hg19.MotifDb.Hsap, verbose=FALSE) }
}
names(TFdiff) = names(comps)

save(int.kmeans, TFdiff, promoters_split, intergenic_split, introns_split, all_split, promoters_int, intergenic_int, introns_int, all_int,
     file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/PWMEnrich/DMR_PWMEnrich_objects.rda")


