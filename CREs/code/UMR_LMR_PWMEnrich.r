library(data.table)
library(ggplot2)
library(PWMEnrich)
library(GenomicRanges)
library(PWMEnrich.Hsapiens.background)
library(BSgenome.Hsapiens.UCSC.hg19)
library(pheatmap)
library(RColorBrewer)

load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/UMRs_LMRs_annotated.rda")
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/UMRs_LMRs_methylSeekR.rda")
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/UMRs_LMRs_methylSeekR_prenatal.rda")
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/DMR/DMR_objects.rda")


## Isolate UMRs and LMRs

segments = list(postnatal = UMRLMRsegments.CG.100kb, prenatal = UMRLMRsegments.CG.100kb.pren)
umrs = lapply(segments, function(x) lapply(x, function(y) y[y$type=="UMR"]))
lmrs = lapply(segments, function(x) lapply(x, function(y) y[y$type=="LMR"]))

# load regions

uall = list(AllUMR = Reduce(intersect, c(umrs$postnatal, umrs$prenatal)),
            PrenatalUMR = Reduce(intersect, umrs$prenatal),
            PostnatalUMR = Reduce(intersect, umrs$postnatal),
            NeuronsUMR = Reduce(intersect, umrs$postnatal[which(names(umrs$postnatal) %in% pd[pd$Cell.Type=="Neuron","Data.ID"])]),
            GliaUMR = Reduce(intersect, umrs$postnatal[which(names(umrs$postnatal) %in% pd[pd$Cell.Type=="Glia","Data.ID"])]))
lall = list(AllLMR = Reduce(intersect, c(lmrs$postnatal, lmrs$prenatal)),
            PrenatalLMR = Reduce(intersect, lmrs$prenatal),
            PostnatalLMR = Reduce(intersect, lmrs$postnatal),
            NeuronsLMR = Reduce(intersect, lmrs$postnatal[which(names(lmrs$postnatal) %in% pd[pd$Cell.Type=="Neuron","Data.ID"])]),
            GliaLMR = Reduce(intersect, lmrs$postnatal[which(names(lmrs$postnatal) %in% pd[pd$Cell.Type=="Glia","Data.ID"])]))
all = list(All = reduce(c(uall$AllUMR, lall$AllLMR)),
           Prenatal = reduce(c(uall$PrenatalUMR, lall$PrenatalLMR)),
           Postnatal = reduce(c(uall$PostnatalUMR, lall$PostnatalLMR)),
           Neurons = reduce(c(uall$NeuronsUMR, lall$NeuronsLMR)),
           Glia = reduce(c(uall$GliaUMR, lall$GliaLMR)))
for (i in 1:length(all)) {
names(all[[i]]) = paste0(seqnames(all[[i]]), ":", data.frame(ranges(all[[i]]))$start,"-", data.frame(ranges(all[[i]]))$end)
names(uall[[i]]) = paste0(seqnames(uall[[i]]), ":", data.frame(ranges(uall[[i]]))$start,"-", data.frame(ranges(uall[[i]]))$end)
names(lall[[i]]) = paste0(seqnames(lall[[i]]), ":", data.frame(ranges(lall[[i]]))$start,"-", data.frame(ranges(lall[[i]]))$end)
}


## Get sequences from different annotations in the correct format

seq = lapply(c(lapply(uall, function(x) x[which(width(x)>=30)]),
               lapply(lall, function(x) x[which(width(x)>=30)]),
               lapply(all, function(x) x[which(width(x)>=30)])), function(x) getSeq(Hsapiens, x))


### Run PWMEnrich

useBigMemoryPWMEnrich(TRUE)
registerCoresPWMEnrich(4)

# load the pre-compiled lognormal background computed using promoters

data(PWMLogn.hg19.MotifDb.Hsap)

for (i in 1:length(seq)) { 
  umrlmr_pwmenrich = motifEnrichment(seq[[i]], PWMLogn.hg19.MotifDb.Hsap, verbose=F)
  save(umrlmr_pwmenrich, 
       file = paste0("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/PWMEnrich/umrlmr/",names(seq)[i],"_PWMEnrich_object.rda"))
}
LMR_UMR_pwmenrich = list()
for (i in 1:length(seq)) { 
  load(paste0("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/PWMEnrich/umrlmr/",names(seq)[i],"_PWMEnrich_object.rda"))
  LMR_UMR_pwmenrich[[i]] = umrlmr_pwmenrich
  rm(umrlmr_pwmenrich)
}
names(LMR_UMR_pwmenrich) = names(seq)

save(LMR_UMR_pwmenrich, geneMap,pd, file="/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/PWMEnrich/UMR_LMR_PWMEnrich_object.rda")


## Test for Differential TF binding

ulTFdiff = mapply(function(u,l) motifDiffEnrichment(sequences1 = seq[[l]], sequences2 = seq[[u]], res1 = LMR_UMR_pwmenrich[[l]], 
                                                    res2 = LMR_UMR_pwmenrich[[u]], PWMLogn.hg19.MotifDb.Hsap, verbose=FALSE), c(1:5),c(6:10), SIMPLIFY = F)
names(ulTFdiff) = names(all)

ulTFdiff = c(ulTFdiff, 
           list(pre.post.UMR = motifDiffEnrichment(sequences1 = seq$PrenatalUMR, sequences2 = seq$PostnatalUMR,
                                                   res1 = LMR_UMR_pwmenrich$PrenatalUMR, res2 = LMR_UMR_pwmenrich$PostnatalUMR, PWMLogn.hg19.MotifDb.Hsap, verbose=FALSE),
                pre.neuron.UMR = motifDiffEnrichment(sequences1 = seq$PrenatalUMR, sequences2 = seq$NeuronsUMR,
                                                     res1 = LMR_UMR_pwmenrich$PrenatalUMR, res2 = LMR_UMR_pwmenrich$NeuronsUMR, PWMLogn.hg19.MotifDb.Hsap, verbose=FALSE),
                pre.glia.UMR = motifDiffEnrichment(sequences1 = seq$PrenatalUMR, sequences2 = seq$GliaUMR,
                                                   res1 = LMR_UMR_pwmenrich$PrenatalUMR, res2 = LMR_UMR_pwmenrich$GliaUMR, PWMLogn.hg19.MotifDb.Hsap, verbose=FALSE),
                neuron.glia.UMR = motifDiffEnrichment(sequences1 = seq$NeuronsUMR, sequences2 = seq$GliaUMR,
                                                      res1 = LMR_UMR_pwmenrich$NeuronsUMR, res2 = LMR_UMR_pwmenrich$GliaUMR, PWMLogn.hg19.MotifDb.Hsap, verbose=FALSE),
                pre.post.LMR = motifDiffEnrichment(sequences1 = seq$PrenatalLMR, sequences2 = seq$PostnatalLMR,
                                                   res1 = LMR_UMR_pwmenrich$PrenatalLMR, res2 = LMR_UMR_pwmenrich$PostnatalLMR, PWMLogn.hg19.MotifDb.Hsap, verbose=FALSE),
                pre.neuron.LMR = motifDiffEnrichment(sequences1 = seq$PrenatalLMR, sequences2 = seq$NeuronsLMR,
                                                     res1 = LMR_UMR_pwmenrich$PrenatalLMR, res2 = LMR_UMR_pwmenrich$NeuronsLMR, PWMLogn.hg19.MotifDb.Hsap, verbose=FALSE),
                pre.glia.LMR = motifDiffEnrichment(sequences1 = seq$PrenatalLMR, sequences2 = seq$GliaLMR,
                                                   res1 = LMR_UMR_pwmenrich$PrenatalLMR, res2 = LMR_UMR_pwmenrich$GliaLMR, PWMLogn.hg19.MotifDb.Hsap, verbose=FALSE),
                neuron.glia.LMR = motifDiffEnrichment(sequences1 = seq$NeuronsLMR, sequences2 = seq$GliaLMR,
                                                      res1 = LMR_UMR_pwmenrich$NeuronsLMR, res2 = LMR_UMR_pwmenrich$GliaLMR, PWMLogn.hg19.MotifDb.Hsap, verbose=FALSE)))

save(ulTFdiff, LMR_UMR_pwmenrich, geneMap,pd, file="/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/PWMEnrich/UMR_LMR_PWMEnrich_object.rda")


## Call binding motifs in 10kb region around TSSs for TF network analysis

umrs = unlist(umrs, recursive = F)
lmrs = unlist(lmrs, recursive = F)

TSSs = makeGRangesFromDataFrame(data.frame(seqnames = geneMap$Chr, start = geneMap$Start-5000, end = geneMap$Start+5000, strand = geneMap$Strand))
names(TSSs) = geneMap$gencodeID

uoo = lapply(umrs, function(x) findOverlaps(TSSs, x)) 
loo = lapply(lmrs, function(x) findOverlaps(TSSs, x)) 

umrs = mapply(function(u, oo) u[subjectHits(oo)], umrs, uoo, SIMPLIFY = F)
lmrs = mapply(function(l, oo) l[subjectHits(oo)], lmrs, loo, SIMPLIFY = F)

for (i in 1:length(umrs)) {
  names(umrs[[i]]) = paste0(seqnames(umrs[[i]]), ":", data.frame(ranges(umrs[[i]]))$start,"-", data.frame(ranges(umrs[[i]]))$end)
  names(lmrs[[i]]) = paste0(seqnames(lmrs[[i]]), ":", data.frame(ranges(lmrs[[i]]))$start,"-", data.frame(ranges(lmrs[[i]]))$end)
}

UMR.seq = lapply(lapply(umrs, function(x) x[which(width(x)>=30)]), function(x) getSeq(Hsapiens, x))
names(UMR.seq) = paste0("UMR.", names(UMR.seq))
LMR.seq = lapply(lapply(lmrs, function(x) x[which(width(x)>=30)]), function(x) getSeq(Hsapiens, x))
names(LMR.seq) = paste0("LMR.", names(LMR.seq))

rm(list=ls()[-which(ls() %in% c("LMR.seq","PWMLogn.hg19.MotifDb.Hsap", "geneMap", "pd"))])

for (i in 1:length(LMR.seq)) { 
  umrlmr_pwmenrich = motifEnrichment(LMR.seq[[i]], PWMLogn.hg19.MotifDb.Hsap, verbose=F)
  save(umrlmr_pwmenrich, 
       file = paste0("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/PWMEnrich/umrlmr/",names(LMR.seq)[i],"_TSS_TFs_bysample.rda"))
}

TSS_TFs_bysample = list()
for (i in 1:length(LMR.seq)) { 
  load(paste0("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/PWMEnrich/umrlmr/",names(LMR.seq)[i],"_TSS_TFs_bysample.rda"))
  TSS_TFs_bysample[[i]] = umrlmr_pwmenrich
  rm(umrlmr_pwmenrich)
}
names(TSS_TFs_bysample) = names(LMR.seq)

save(TSS_TFs_bysample, geneMap,pd, file="/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/PWMEnrich/TSS_TFs_bysample.rda")