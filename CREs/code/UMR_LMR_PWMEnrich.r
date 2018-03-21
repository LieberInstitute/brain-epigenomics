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
registerCoresPWMEnrich(6)

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
names(ulTFdiff) = c("All","Prenatal","Postnatal","Neurons","Glia")

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
                                                      res1 = LMR_UMR_pwmenrich$NeuronsLMR, res2 = LMR_UMR_pwmenrich$GliaLMR, PWMLogn.hg19.MotifDb.Hsap, verbose=FALSE),
                all.pre.UMR = motifDiffEnrichment(sequences1 = seq$AllUMR, sequences2 = seq$PrenatalUMR,
                                               res1 = LMR_UMR_pwmenrich$AllUMR, res2 = LMR_UMR_pwmenrich$PrenatalUMR, PWMLogn.hg19.MotifDb.Hsap, verbose=FALSE),
                all.post.UMR = motifDiffEnrichment(sequences1 = seq$AllUMR, sequences2 = seq$PostnatalUMR,
                                                   res1 = LMR_UMR_pwmenrich$AllUMR, res2 = LMR_UMR_pwmenrich$PostnatalUMR, PWMLogn.hg19.MotifDb.Hsap, verbose=FALSE),
                all.neuron.UMR = motifDiffEnrichment(sequences1 = seq$AllUMR, sequences2 = seq$NeuronsUMR,
                                                     res1 = LMR_UMR_pwmenrich$AllUMR, res2 = LMR_UMR_pwmenrich$NeuronsUMR, PWMLogn.hg19.MotifDb.Hsap, verbose=FALSE),
                all.glia.UMR = motifDiffEnrichment(sequences1 = seq$AllUMR, sequences2 = seq$GliaUMR,
                                                   res1 = LMR_UMR_pwmenrich$AllUMR, res2 = LMR_UMR_pwmenrich$GliaUMR, PWMLogn.hg19.MotifDb.Hsap, verbose=FALSE),
                all.pre.LMR = motifDiffEnrichment(sequences1 = seq$AllLMR, sequences2 = seq$PrenatalLMR,
                                                  res1 = LMR_UMR_pwmenrich$AllLMR, res2 = LMR_UMR_pwmenrich$PrenatalLMR, PWMLogn.hg19.MotifDb.Hsap, verbose=FALSE),
                all.post.LMR = motifDiffEnrichment(sequences1 = seq$AllLMR, sequences2 = seq$PostnatalLMR,
                                                   res1 = LMR_UMR_pwmenrich$AllLMR, res2 = LMR_UMR_pwmenrich$PostnatalLMR, PWMLogn.hg19.MotifDb.Hsap, verbose=FALSE),
                all.neuron.LMR = motifDiffEnrichment(sequences1 = seq$AllLMR, sequences2 = seq$NeuronsLMR,
                                                     res1 = LMR_UMR_pwmenrich$AllLMR, res2 = LMR_UMR_pwmenrich$NeuronsLMR, PWMLogn.hg19.MotifDb.Hsap, verbose=FALSE),
                all.glia.LMR = motifDiffEnrichment(sequences1 = seq$AllLMR, sequences2 = seq$GliaLMR,
                                                   res1 = LMR_UMR_pwmenrich$AllLMR, res2 = LMR_UMR_pwmenrich$GliaLMR, PWMLogn.hg19.MotifDb.Hsap, verbose=FALSE)))
                
            

save(ulTFdiff, LMR_UMR_pwmenrich, geneMap,pd, file="/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/PWMEnrich/UMR_LMR_PWMEnrich_object.rda")


## Call binding motifs in 10kb region around TSSs for TF network analysis

## Get PWMs

data(MotifDb.Hsap)
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/PWMEnrich/filtered_TF_list_expressionInfo.rda")
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/PWMEnrich/spliceRegions_PWMEnrich_objects.rda")
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/brainseq_pipeline/polyA_unstranded/rse_gene_polyA_dlpfc_n41.Rdata", verbose = T)

group = lapply(splice.exon, groupReport)
group = lapply(group, as.data.frame)
group = lapply(group, function(x) x[which(x$target %in% names(PostnataltargettogeneID)),])

pwms = MotifDb.Hsap[which(names(MotifDb.Hsap) %in% group$CpG.pos$id)] # all the motifs included in my other analyses

## Make new log-normal background on 2kb promoters

data(hg19.upstream2000)
lognPromoters = makePWMLognBackground(hg19.upstream2000, pwms, algorithm = "human", verbose=F)

## Filter regions to those around the TSS

geneMap = rowRanges(rse_gene)

# promoters, -5kb to +5kb
genePromoters = GRanges(seqnames(geneMap),
                        IRanges(start = ifelse(strand(geneMap) == "+",
                                               start(geneMap)-5000, end(geneMap)-5000),
                                end = ifelse(strand(geneMap) == "+",
                                             start(geneMap)+5000, end(geneMap)+5000)),
                        strand = strand(geneMap))
mcols(genePromoters) = mcols(geneMap)

lmrs = unlist(lmrs, recursive = F)
loo = lapply(lmrs, function(x) findOverlaps(genePromoters, x)) 
lmrs = mapply(function(l, oo) l[subjectHits(oo)], lmrs, loo, SIMPLIFY = F)
for (i in 1:length(lmrs)) {
  names(lmrs[[i]]) = paste0(seqnames(lmrs[[i]]), ":", data.frame(ranges(lmrs[[i]]))$start,"-", data.frame(ranges(lmrs[[i]]))$end)
}

LMR.seq = lapply(lapply(lmrs, function(x) x[which(width(x)>=30)]), function(x) getSeq(Hsapiens, x))

rm(list=ls()[-which(ls() %in% c("LMR.seq","lognPromoters", "geneMap", "pd"))])

save(LMR.seq,lognPromoters, geneMap, pd, file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/PWMEnrich/newPromoterBackground.rda")

for (i in 1:length(LMR.seq)) { 
  lmr_pwmenrich = motifEnrichment(LMR.seq[[i]], lognPromoters, verbose=F)
  save(lmr_pwmenrich, 
       file = paste0("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/PWMEnrich/umrlmr/",names(LMR.seq)[i],"_TSS_TFs_bysample.rda"))
}


TSS_TFs_bysample = list()
for (i in 1:length(LMR.seq)) { 
  load(paste0("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/PWMEnrich/umrlmr/",names(LMR.seq)[i],"_TSS_TFs_bysample.rda"))
  TSS_TFs_bysample[[i]] = lmr_pwmenrich
  rm(lmr_pwmenrich)
}
names(TSS_TFs_bysample) = names(LMR.seq)

save(TSS_TFs_bysample, geneMap,pd, file="/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/PWMEnrich/TSS_TFs_bysample.rda")