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

introns_split$CellType.pos = motifEnrichment(seq$introns$CellType.pos, PWMLogn.hg19.MotifDb.Hsap, verbose=F)
introns_split$CellType.neg = motifEnrichment(seq$introns$CellType.neg, PWMLogn.hg19.MotifDb.Hsap, verbose=F)
introns_split$Age.pos = motifEnrichment(seq$introns$Age.pos, PWMLogn.hg19.MotifDb.Hsap, verbose=F)
introns_split$Age.neg = motifEnrichment(seq$introns$Age.neg, PWMLogn.hg19.MotifDb.Hsap, verbose=F)
introns_split$Interaction.pos = motifEnrichment(seq$introns$Interaction.pos, PWMLogn.hg19.MotifDb.Hsap, verbose=F)
introns_split$Interaction.neg = motifEnrichment(seq$introns$Interaction.neg, PWMLogn.hg19.MotifDb.Hsap, verbose=F)


save(promoters_split, intergenic_split, introns_split, all_split, promoters_int, intergenic_int, introns_int, all_int,
     file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/PWMEnrich/DMR_PWMEnrich_objects.rda")


## Test for Differential TF binding

pos = c(all_split[grep("pos", names(all_split))], promoters_split[grep("pos", names(promoters_split))],
        introns_split[grep("pos", names(introns_split))], intergenic_split[grep("pos", names(intergenic_split))])
neg = c(all_split[grep("neg", names(all_split))], promoters_split[grep("neg", names(promoters_split))],
        introns_split[grep("neg", names(introns_split))], intergenic_split[grep("neg", names(intergenic_split))])
names(pos) = c("allAge","allInt","allCT","promCT","promAge","promInt","intronsCT","intronsAge","intronsInt","interCT","interAge","interInt")
names(neg) = c("allAge","allInt","allCT","promCT","promAge","promInt","intronsCT","intronsInt","intronsAge","interCT","interInt")
pos = pos[names(pos)!="interAge"]

seqpos = c(allAge = seq$all$Age.pos, 
           allInt = seq$all$Interaction.pos, 
           allCT = seq$all$CellType.pos,
           promCT = seq$promoters$CellType.pos, 
           promAge = seq$promoters$Age.pos, 
           promInt = seq$promoters$Interaction.pos,
           intronsCT = seq$introns$CellType.pos, 
           intronsAge = seq$introns$Age.pos, 
           intronsInt = seq$introns$Interaction.pos,
           interCT = seq$intergenic$CellType.pos, 
           interInt = seq$intergenic$Interaction.pos)
seqneg = c(allAge = seq$all$Age.neg, 
           allInt = seq$all$Interaction.neg, 
           allCT = seq$all$CellType.neg,
           promCT = seq$promoters$CellType.neg, 
           promAge = seq$promoters$Age.neg, 
           promInt = seq$promoters$Interaction.neg,
           intronsCT = seq$introns$CellType.neg, 
           intronsInt = seq$introns$Interaction.neg,
           intronsAge = seq$introns$Age.neg,
           interCT = seq$intergenic$CellType.neg, 
           interInt = seq$intergenic$Interaction.neg)

TFdiff = list()
for (i in 1:length(pos)) {
  TFdiff[[i]] = motifDiffEnrichment(sequences1 = seqpos[[i]], sequences2 = seqneg[[i]],
                                    res1 = pos[[i]], res2 = neg[[i]], PWMLogn.hg19.MotifDb.Hsap, verbose=FALSE)
}
names(TFdiff) = names(pos)
pos = c(all_split[grep("pos", names(all_split))], promoters_split[grep("pos", names(promoters_split))],
        introns_split[grep("pos", names(introns_split))], intergenic_split[grep("pos", names(intergenic_split))])
names(pos) = c("allAge","allInt","allCT","promCT","promAge","promInt","intronsCT","intronsAge","intronsInt","interCT","interAge","interInt")
seqpos = c(allAge = seq$all$Age.pos, allInt = seq$all$Interaction.pos, allCT = seq$all$CellType.pos, promCT = seq$promoters$CellType.pos, 
           promAge = seq$promoters$Age.pos, promInt = seq$promoters$Interaction.pos, intronsCT = seq$introns$CellType.pos, 
           intronsAge = seq$introns$Age.pos, intronsInt = seq$introns$Interaction.pos, interCT = seq$intergenic$CellType.pos,
           interAge = seq$intergenic$Age.pos, interInt = seq$intergenic$Interaction.pos)

TFdiff = c(TFdiff, 
           list(prom.intergenic.pos.CT = motifDiffEnrichment(sequences1 = seqpos$promCT, sequences2 = seqpos$interCT,res1 = pos$promCT, 
                                                             res2 = pos$interCT, PWMLogn.hg19.MotifDb.Hsap, verbose=FALSE),
                prom.intron.pos.CT = motifDiffEnrichment(sequences1 = seqpos$promCT, sequences2 = seqpos$intronsCT,res1 = pos$promCT, 
                                                         res2 = pos$intronsCT, PWMLogn.hg19.MotifDb.Hsap, verbose=FALSE),
                intron.intergenic.pos.CT = motifDiffEnrichment(sequences1 = seqpos$intronsCT, sequences2 = seqpos$interCT, res1 = pos$intronsCT, 
                                                               res2 = pos$interCT, PWMLogn.hg19.MotifDb.Hsap, verbose=FALSE),
                prom.intergenic.neg.CT = motifDiffEnrichment(sequences1 = seqneg$promCT, sequences2 = seqneg$interCT, res1 = neg$promCT, 
                                                             res2 = neg$interCT, PWMLogn.hg19.MotifDb.Hsap, verbose=FALSE),
                prom.intron.neg.CT = motifDiffEnrichment(sequences1 = seqneg$promCT, sequences2 = seqneg$intronsCT, res1 = neg$promCT, 
                                                         res2 = neg$intronsCT, PWMLogn.hg19.MotifDb.Hsap, verbose=FALSE),
                intron.intergenic.neg.CT = motifDiffEnrichment(sequences1 = seqneg$intronsCT, sequences2 = seqneg$interCT, res1 = neg$intronsCT, 
                                                               res2 = neg$interCT, PWMLogn.hg19.MotifDb.Hsap, verbose=FALSE),
                prom.intergenic.pos.Age = motifDiffEnrichment(sequences1 = seqpos$promAge, sequences2 = seqpos$interAge, res1 = pos$promAge, 
                                                              res2 = pos$interAge, PWMLogn.hg19.MotifDb.Hsap, verbose=FALSE),
                prom.intron.pos.Age = motifDiffEnrichment(sequences1 = seqpos$promAge, sequences2 = seqpos$intronsAge, res1 = pos$promAge, 
                                                          res2 = pos$intronsAge, PWMLogn.hg19.MotifDb.Hsap, verbose=FALSE),
                intron.intergenic.pos.Age = motifDiffEnrichment(sequences1 = seqpos$intronsAge, sequences2 = seqpos$interAge, res1 = pos$intronsAge, 
                                                                res2 = pos$interAge, PWMLogn.hg19.MotifDb.Hsap, verbose=FALSE),
                prom.intron.neg.Age = motifDiffEnrichment(sequences1 = seqneg$promAge, sequences2 = seqneg$intronsAge, res1 = neg$promAge, 
                                                          res2 = neg$intronsAge, PWMLogn.hg19.MotifDb.Hsap, verbose=FALSE),
                prom.intergenic.pos.Int = motifDiffEnrichment(sequences1 = seqpos$promInt, sequences2 = seqpos$interInt, res1 = pos$promInt, 
                                                              res2 = pos$interInt, PWMLogn.hg19.MotifDb.Hsap, verbose=FALSE),
                prom.intron.pos.Int = motifDiffEnrichment(sequences1 = seqpos$promInt, sequences2 = seqpos$intronsInt, res1 = pos$promInt, 
                                                          res2 = pos$intronsInt, PWMLogn.hg19.MotifDb.Hsap, verbose=FALSE),
                intron.intergenic.pos.Int = motifDiffEnrichment(sequences1 = seqpos$intronsInt, sequences2 = seqpos$interInt, res1 = pos$intronsInt, 
                                                                res2 = pos$interInt, PWMLogn.hg19.MotifDb.Hsap, verbose=FALSE),
                prom.intergenic.neg.Int = motifDiffEnrichment(sequences1 = seqneg$promInt, sequences2 = seqneg$interInt, res1 = neg$promInt, 
                                                              res2 = neg$interInt, PWMLogn.hg19.MotifDb.Hsap, verbose=FALSE),
                prom.intron.neg.Int = motifDiffEnrichment(sequences1 = seqneg$promInt, sequences2 = seqneg$intronsInt, res1 = neg$promInt, 
                                                          res2 = neg$intronsInt, PWMLogn.hg19.MotifDb.Hsap, verbose=FALSE),
                intron.intergenic.neg.Int = motifDiffEnrichment(sequences1 = seqneg$intronsInt, sequences2 = seqneg$interInt, res1 = neg$intronsInt, 
                                                                res2 = neg$interInt, PWMLogn.hg19.MotifDb.Hsap, verbose=FALSE)))

save(TFdiff, promoters_split, intergenic_split, introns_split, all_split, promoters_int, intergenic_int, introns_int, all_int,
     file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/PWMEnrich/DMR_PWMEnrich_objects.rda")


