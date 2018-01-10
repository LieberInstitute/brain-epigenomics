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

promoters = lapply(DMR, function(x) x[which(x$promoter=="Promoter" & x$sig=="FWER < 0.05" & x$width>=6),]) # all DMRs that overlap a promoter
intergenic = lapply(DMR, function(x) x[which(x$annotation=="Intergenic" & x$sig=="FWER < 0.05" & x$width>=6),]) # all intergenic DMRs
introns = lapply(DMR, function(x) x[which(x$annotation=="Intron" & x$sig=="FWER < 0.05" & x$width>=6),]) #all DMRs that overlap exclusively introns

promotersgr = lapply(promoters, makeGRangesFromDataFrame)
intergenicgr = lapply(intergenic, makeGRangesFromDataFrame)
intronsgr = lapply(introns, makeGRangesFromDataFrame)

promoters_seq = lapply(promotersgr, function(x) getSeq(Hsapiens, x))
intergenic_seq = lapply(intergenicgr, function(x) getSeq(Hsapiens, x))
introns_seq = lapply(intronsgr, function(x) getSeq(Hsapiens, x))


### Run PWMEnrich

useBigMemoryPWMEnrich(TRUE)
registerCoresPWMEnrich(5)

# load the pre-compiled lognormal background computed using promoters
data(PWMLogn.hg19.MotifDb.Hsap)

promoters_res = intergenic_res = introns_res = list()
promoters_res$CellType = motifEnrichment(promoters_seq$CellType, PWMLogn.hg19.MotifDb.Hsap)
promoters_res$Age = motifEnrichment(promoters_seq$Age, PWMLogn.hg19.MotifDb.Hsap)
promoters_res$Interaction = motifEnrichment(promoters_seq$Interaction, PWMLogn.hg19.MotifDb.Hsap)

save(promoters_res, file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/PWMEnrich/DMR_PWMEnrich_objects")

intergenic_res$CellType = motifEnrichment(intergenic_seq$CellType, PWMLogn.hg19.MotifDb.Hsap)
intergenic_res$Age = motifEnrichment(intergenic_seq$Age, PWMLogn.hg19.MotifDb.Hsap)
intergenic_res$Interaction = motifEnrichment(intergenic_seq$Interaction, PWMLogn.hg19.MotifDb.Hsap)

save(promoters_res, intergenic_res, file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/PWMEnrich/DMR_PWMEnrich_objects")

introns_res$CellType = motifEnrichment(introns_seq$CellType, PWMLogn.hg19.MotifDb.Hsap)
introns_res$Age = motifEnrichment(introns_seq$Age, PWMLogn.hg19.MotifDb.Hsap)
introns_res$Interaction = motifEnrichment(introns_seq$Interaction, PWMLogn.hg19.MotifDb.Hsap)

save(promoters_res, intergenic_res, introns_res, file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/PWMEnrich/DMR_PWMEnrich_objects")


lapply(promoters, function(x) min(x$width))
lapply(intergenic, function(x) min(x$width))
lapply(introns, function(x) min(x$width))




promoters_groupReport = groupReport(promoters_res$Age)
promoters_sequenceReport = sequenceReport(promoters_res$Age)









