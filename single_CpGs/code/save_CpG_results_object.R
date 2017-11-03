library(bsseq)
library('limma')

load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/bumphunting/limma_Neuron_CpGs_minCov_3.Rdata")
load('/dcl01/lieber/ajaffe/lab/brain-epigenomics/bumphunting/BSobj_bsseqSmooth_Neuron_minCov_3.Rdata')
load("/dcl01/lieber/ajaffe/CellSorting/RNAseq_pipeline/rawCounts_CellSorting_July5_n12.rda")


## load phenotype data

pd = pData(BSobj)
pd$Race[pd$Race== "CAUC "] <- 'CAUC'
pd$Sex[pd$Sex == " M"] <- 'M'
pd$RIN <- as.numeric(gsub(" ", "", pd$RIN))
pd$pg.DNA.nuclei.input <- as.numeric(pd$pg.DNA.nuclei.input)
pd$Reads <- as.numeric(pd$Reads)
pd$Percent.GreaterThan.Q30 <- as.numeric(pd$Percent.GreaterThan.Q30)

## load methylation matrix for filtered CpGs
meth = getMeth(BSobj, type = 'raw')
gr = granges(BSobj)

summary(coef_interest)
summary(abs(coef_interest))

ebList = lapply(fits, ebayes)


## Make CpG results dataframe

CpG = data.frame(as.data.frame(gr), lods.CellType = ebList$cell$lods[,"Cell.TypeNeuron"], Tstat.CellType = ebList$cell$t[,"Cell.TypeNeuron"], 
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
features = c(features, islands = makeGRangesFromDataFrame(islands, keep.extra.columns = T, start.field = "CpGromStart", end.field = "CpGromEnd"),
             promoters = promoters(txdb, upstream=2000, downstream=200))
lapply(features, head)

grCpG = makeGRangesFromDataFrame(CpG, keep.extra.columns = T)
annotation = lapply(features, function(y) findOverlaps(grCpG, y))

grCpG$rnum = 1:length(grCpG)
grCpG$cds = ifelse(grCpG$rnum %in% queryHits(annotation[["CDS"]]), "CDS", NA)
grCpG$intron = ifelse(grCpG$rnum %in% queryHits(annotation[["Introns"]]), "Intron", NA)
grCpG$UTR5 = ifelse(grCpG$rnum %in% queryHits(annotation[["UTR5"]]), "5'UTR", NA)
grCpG$UTR3 = ifelse(grCpG$rnum %in% queryHits(annotation[["UTR3"]]), "3'UTR", NA)
grCpG$islands = ifelse(grCpG$rnum %in% queryHits(annotation[["islands"]]), "CpG Island", "non-Island")
grCpG$promoter = ifelse(grCpG$rnum %in% queryHits(annotation[["promoters"]]), "Promoter", NA)
grCpG$anno = paste0(grCpG$cds,":",grCpG$intron, ":", grCpG$UTR5, ":", grCpG$UTR3, ":", grCpG$promoter)

CpG = as.data.frame(grCpG)
CpG[which(CpG$anno == "NA:NA:NA:NA:NA"),"annotation"] = "Intergenic" 
CpG[grep("CDS", CpG$cds),"annotation"] = "CDS"
CpG[which(is.na(CpG$annotation) & CpG$UTR5 == "5'UTR"),"annotation"] = "5'UTR"
CpG[which(is.na(CpG$annotation) & CpG$UTR3 == "3'UTR"),"annotation"] = "3'UTR"
CpG[which(is.na(CpG$annotation) & CpG$intron == "Intron"),"annotation"] = "Intron"
CpG[which(is.na(CpG$annotation) & CpG$promoter == "Promoter"),"annotation"] = "Promoter"

# Mapping sites to nearest gene
geneMapGR = makeGRangesFromDataFrame(geneMap, start.field="Start",end.field="End",strand.field="Strand",keep=TRUE)
dA = distanceToNearest(grCpG, geneMapGR)
CpG$nearestSymbol = geneMapGR$Symbol[subjectHits(dA)]
CpG$nearestID = names(geneMapGR)[subjectHits(dA)]
CpG$distToGene = mcols(dA)$distance
CpG$EntrezID = geneMapGR$EntrezID[subjectHits(dA)]
CpG$regionID = paste0(CpG$seqnames,":",CpG$start,"-", CpG$end)
CpG$CT.sig = ifelse(CpG$padj.CellType <=0.05, "FDR < 0.05", "FDR > 0.05")
CpG$Age.sig = ifelse(CpG$padj.Age <=0.05, "FDR < 0.05", "FDR > 0.05")
CpG$Int.sig = ifelse(CpG$padj.Interaction <=0.05, "FDR < 0.05", "FDR > 0.05")
CpG$CT.dir = ifelse(CpG$Tstat.CellType < 0, "neg", "pos")
CpG$Age.dir = ifelse(CpG$Tstat.Age < 0, "neg", "pos")
CpG$Int.dir = ifelse(CpG$Tstat.Interaction < 0, "neg", "pos")


save(CpG, geneMap, pd, file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/single_CpGs/singleCpG_object.rda")