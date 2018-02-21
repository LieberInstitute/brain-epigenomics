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

