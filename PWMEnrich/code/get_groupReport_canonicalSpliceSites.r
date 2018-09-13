


load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/PWMEnrich/spliceRegions_PWMEnrich_objects_exonAdjustedBackgrounds.rda")

group = lapply(C.canonicalspliceBack, groupReport)
group = do.call(rbind, Map(cbind, lapply(group, as.data.frame), Group = as.list(names(group))))
group$FDR = p.adjust(group$p.value, method = "fdr")

write.csv(group, quote=F, file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/PWMEnrich/canonicalSpliceSite_enrichment.csv")