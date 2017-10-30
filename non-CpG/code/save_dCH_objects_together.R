library(GenomicRanges)
library(GenomicFeatures)

load("/dcl01/lieber/ajaffe/CellSorting/RNAseq_pipeline/rawCounts_CellSorting_July5_n12.rda")
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/non-CpG/limma_exploration_nonCG_highCov.Rdata")


## How many regions are dCH?

CHres = data.frame(as.data.frame(gr), lods.CellType = ebList$cell$lods[,"Cell.TypeNeuron"], Tstat.CellType = ebList$cell$t[,"Cell.TypeNeuron"], 
                   pval.CellType = ebList$cell$p.value[,"Cell.TypeNeuron"], padj.CellType = p.adjust(ebList$cell$p.value[,"Cell.TypeNeuron"], method = "fdr"),
                   lods.Age = ebList$age$lods[,"Age"], Tstat.Age = ebList$age$t[,"Age"], 
                   pval.Age = ebList$age$p.value[,"Age"], padj.Age = p.adjust(ebList$age$p.value[,"Age"], method = "fdr"),
                   lods.Interaction = ebList$interaction$lods[,"Age:Cell.TypeNeuron"], Tstat.Interaction = ebList$interaction$t[,"Age:Cell.TypeNeuron"], 
                   pval.Interaction = ebList$interaction$p.value[,"Age:Cell.TypeNeuron"], padj.Interaction = p.adjust(ebList$interaction$p.value[,"Age:Cell.TypeNeuron"], method = "fdr"))

# Annotate editing sites to features in the genome
txdb = loadDb("/dcl01/lieber/ajaffe/Amanda/annotation_objects/gencode.v25lift37.annotation.sqlite")
islands = read.table("/dcl01/lieber/ajaffe/Amanda/annotation_objects/cpgIslandExt.hg19.txt", sep="\t", header = T)
features = list(CDS = cdsBy(txdb, by="tx", use.names=T), Introns = intronsByTranscript(txdb, use.names=T), 
                UTR5 = fiveUTRsByTranscript(txdb, use.names=T), UTR3 = threeUTRsByTranscript(txdb, use.names=T))
features = lapply(features, function(x) unlist(x, recursive = TRUE, use.names = TRUE))
for (i in 1:length(features)){
  tmp = features[[i]]
  tmp$TxID = names(tmp)
  features[[i]] = tmp
}
features = c(features, islands = makeGRangesFromDataFrame(islands, keep.extra.columns = T, start.field = "chromStart", end.field = "chromEnd"),
             promoters = promoters(txdb, upstream=2000, downstream=200))
lapply(features, head)

grCH = makeGRangesFromDataFrame(CHres, keep.extra.columns = T)
annotation = lapply(features, function(y) findOverlaps(grCH, y))

grCH$rnum = 1:length(grCH)
grCH$cds = ifelse(grCH$rnum %in% queryHits(annotation[["CDS"]]), "CDS", NA)
grCH$intron = ifelse(grCH$rnum %in% queryHits(annotation[["Introns"]]), "Intron", NA)
grCH$UTR5 = ifelse(grCH$rnum %in% queryHits(annotation[["UTR5"]]), "UTR5", NA)
grCH$UTR3 = ifelse(grCH$rnum %in% queryHits(annotation[["UTR3"]]), "UTR3", NA)
grCH$islands = ifelse(grCH$rnum %in% queryHits(annotation[["islands"]]), "CpG_Island", "non-Island")
grCH$promoter = ifelse(grCH$rnum %in% queryHits(annotation[["promoters"]]), "Promoter", NA)
grCH$anno = paste0(grCH$cds,":",grCH$intron, ":", grCH$UTR5, ":", grCH$UTR3, ":", grCH$promoter)

CH = as.data.frame(grCH)
CH[which(CH$anno == "NA:NA:NA:NA:NA"),"annotation"] = "Other" 
CH[grep("CDS", CH$cds),"annotation"] = "CDS"
CH[which(is.na(CH$annotation) & CH$UTR5 == "UTR5"),"annotation"] = "UTR5"
CH[which(is.na(CH$annotation) & CH$UTR3 == "UTR3"),"annotation"] = "UTR3"
CH[which(is.na(CH$annotation) & CH$intron == "Intron"),"annotation"] = "Intron"
CH[which(is.na(CH$annotation) & CH$promoter == "Promoter"),"annotation"] = "Promoter"

# Mapping sites to nearest gene
geneMapGR = makeGRangesFromDataFrame(geneMap, start.field="Start",end.field="End",strand.field="Strand",keep=TRUE)
dA = distanceToNearest(grCH, geneMapGR)
CH$nearestSymbol = geneMapGR$Symbol[subjectHits(dA)]
CH$nearestID = names(geneMapGR)[subjectHits(dA)]
CH$distToGene = mcols(dA)$distance
CH$EntrezID = geneMapGR$EntrezID[subjectHits(dA)]
CH$regionID = paste0(CH$seqnames,":",CH$start,"-", CH$end)
CH$CT.sig = ifelse(CH$padj.CellType <=0.05, "FDR < 0.05", "FDR > 0.05")
CH$Age.sig = ifelse(CH$padj.Age <=0.05, "FDR < 0.05", "FDR > 0.05")
CH$Int.sig = ifelse(CH$padj.Interaction <=0.05, "FDR < 0.05", "FDR > 0.05")
CH$CT.dir = ifelse(CH$Tstat.CellType < 0, "neg", "pos")
CH$Age.dir = ifelse(CH$Tstat.Age < 0, "neg", "pos")
CH$Int.dir = ifelse(CH$Tstat.Interaction < 0, "neg", "pos")


save(CH, geneMap, pd, file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/non-CpG/CH_object.rda")

write.csv(CH[which(CH$padj.CellType <=0.001),], file="/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/non-CpG/nonCG_CellType_FDR_0.001.csv", quote=F)
write.csv(CH[which(CH$padj.Age <=0.001),], file="/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/non-CpG/nonCG_Age_FDR_0.001.csv", quote=F)
write.csv(CH[which(CH$padj.Interaction <=0.05),], file="/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/non-CpG/nonCG_Interaction_FDR_0.05.csv", quote=F)





