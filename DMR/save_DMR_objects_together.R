library(GenomicRanges)
library(GenomicFeatures)
library(data.table)
library(ggplot2)
library(clusterProfiler)
require(org.Hs.eg.db)

load("/dcl01/lieber/ajaffe/CellSorting/RNAseq_pipeline/rawCounts_CellSorting_July5_n12.rda")

### Cell Type ###

load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/bumphunting/bumps_bsseqSmooth_Neuron_cell_250_perm.Rdata")

## How many regions are DMRs? (model: ~ age + cell)

cell = bumps[[1]]
dim(cell) # 241924     14
dim(cell[which(cell$fwer<=0.05),]) # 11179    14

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

grcell = makeGRangesFromDataFrame(cell, keep.extra.columns = T)
annotation = lapply(features, function(y) findOverlaps(grcell, y))

grcell$rnum = 1:length(grcell)
grcell$cds = ifelse(grcell$rnum %in% queryHits(annotation[["CDS"]]), "CDS", NA)
grcell$intron = ifelse(grcell$rnum %in% queryHits(annotation[["Introns"]]), "Intron", NA)
grcell$UTR5 = ifelse(grcell$rnum %in% queryHits(annotation[["UTR5"]]), "UTR5", NA)
grcell$UTR3 = ifelse(grcell$rnum %in% queryHits(annotation[["UTR3"]]), "UTR3", NA)
grcell$islands = ifelse(grcell$rnum %in% queryHits(annotation[["islands"]]), "CpG_Island", "non-Island")
grcell$promoter = ifelse(grcell$rnum %in% queryHits(annotation[["promoters"]]), "Promoter", NA)
grcell$anno = paste0(grcell$cds,":",grcell$intron, ":", grcell$UTR5, ":", grcell$UTR3, ":", grcell$promoter)

cell = as.data.frame(grcell)
cell[which(cell$anno == "NA:NA:NA:NA:NA"),"annotation"] = "Other" 
cell[grep("CDS", cell$cds),"annotation"] = "CDS"
cell[which(is.na(cell$annotation) & cell$UTR5 == "UTR5"),"annotation"] = "UTR5"
cell[which(is.na(cell$annotation) & cell$UTR3 == "UTR3"),"annotation"] = "UTR3"
cell[which(is.na(cell$annotation) & cell$intron == "Intron"),"annotation"] = "Intron"
cell[which(is.na(cell$annotation) & cell$promoter == "Promoter"),"annotation"] = "Promoter"

# Mapping cell sites to nearest gene
geneMapGR = makeGRangesFromDataFrame(geneMap, start.field="Start",end.field="End",strand.field="Strand",keep=TRUE)
dA = distanceToNearest(grcell, geneMapGR)
cell$nearestSymbol = geneMapGR$Symbol[subjectHits(dA)]
cell$nearestID = names(geneMapGR)[subjectHits(dA)]
cell$distToGene = mcols(dA)$distance
cell$EntrezID = geneMapGR$EntrezID[subjectHits(dA)]
cell$regionID = paste0(cell$seqnames,":",cell$start,"-", cell$end)
cell$sig = ifelse(cell$fwer<=0.05, "FWER < 0.05", "FWER > 0.05")
cell$Dir = ifelse(cell$value<0, "neg", "pos")
dtcell = data.table(cell)

### Age ###

load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/bumphunting/bumps_bsseqSmooth_Neuron_age_250_perm.Rdata")

## How many regions are DMRs? (model: ~ cell + age)

age = bumps[[1]]
dim(age) # 52790     14
dim(age[which(age$fwer<=0.05),]) # 129    14

# Annotate editing sites to features in the genome
grage = makeGRangesFromDataFrame(age, keep.extra.columns = T)
annotation = lapply(features, function(y) findOverlaps(grage, y))

grage$rnum = 1:length(grage)
grage$cds = ifelse(grage$rnum %in% queryHits(annotation[["CDS"]]), "CDS", NA)
grage$intron = ifelse(grage$rnum %in% queryHits(annotation[["Introns"]]), "Intron", NA)
grage$UTR5 = ifelse(grage$rnum %in% queryHits(annotation[["UTR5"]]), "UTR5", NA)
grage$UTR3 = ifelse(grage$rnum %in% queryHits(annotation[["UTR3"]]), "UTR3", NA)
grage$islands = ifelse(grage$rnum %in% queryHits(annotation[["islands"]]), "CpG_Island", "non-Island")
grage$promoter = ifelse(grage$rnum %in% queryHits(annotation[["promoters"]]), "Promoter", NA)
grage$anno = paste0(grage$cds,":",grage$intron, ":", grage$UTR5, ":", grage$UTR3, ":", grage$promoter)

age = as.data.frame(grage)
age[which(age$anno == "NA:NA:NA:NA:NA"),"annotation"] = "Other" 
age[grep("CDS", age$cds),"annotation"] = "CDS"
age[which(is.na(age$annotation) & age$UTR5 == "UTR5"),"annotation"] = "UTR5"
age[which(is.na(age$annotation) & age$UTR3 == "UTR3"),"annotation"] = "UTR3"
age[which(is.na(age$annotation) & age$intron == "Intron"),"annotation"] = "Intron"
age[which(is.na(age$annotation) & age$promoter == "Promoter"),"annotation"] = "Promoter"

# Mapping age sites to nearest gene
dA = distanceToNearest(grage, geneMapGR)
age$nearestSymbol = geneMapGR$Symbol[subjectHits(dA)]
age$nearestID = names(geneMapGR)[subjectHits(dA)]
age$distToGene = mcols(dA)$distance
age$EntrezID = geneMapGR$EntrezID[subjectHits(dA)]
age$regionID = paste0(age$seqnames,":",age$start,"-", age$end)
age$sig = ifelse(age$fwer<=0.05, "FWER < 0.05", "FWER > 0.05")
age$Dir = ifelse(age$value<0, "neg", "pos")
dtage = data.table(age)

### Interaction ###

load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/bumphunting/bumps_bsseqSmooth_Neuron_interaction_250_perm.Rdata")


## How many regions are DMRs? (model: ~ age + cell + age:cell)

interaction = bumps[[1]]
dim(interaction) # 26282     14
dim(interaction[which(interaction$fwer<=0.05),]) # 2178    14

# Annotate editing sites to features in the genome
grinteraction = makeGRangesFromDataFrame(interaction, keep.extra.columns = T)
annotation = lapply(features, function(y) findOverlaps(grinteraction, y))

grinteraction$rnum = 1:length(grinteraction)
grinteraction$cds = ifelse(grinteraction$rnum %in% queryHits(annotation[["CDS"]]), "CDS", NA)
grinteraction$intron = ifelse(grinteraction$rnum %in% queryHits(annotation[["Introns"]]), "Intron", NA)
grinteraction$UTR5 = ifelse(grinteraction$rnum %in% queryHits(annotation[["UTR5"]]), "UTR5", NA)
grinteraction$UTR3 = ifelse(grinteraction$rnum %in% queryHits(annotation[["UTR3"]]), "UTR3", NA)
grinteraction$islands = ifelse(grinteraction$rnum %in% queryHits(annotation[["islands"]]), "CpG_Island", "non-Island")
grinteraction$promoter = ifelse(grinteraction$rnum %in% queryHits(annotation[["promoters"]]), "Promoter", NA)
grinteraction$anno = paste0(grinteraction$cds,":",grinteraction$intron, ":", grinteraction$UTR5, ":", grinteraction$UTR3, ":", grinteraction$promoter)

interaction = as.data.frame(grinteraction)
interaction[which(interaction$anno == "NA:NA:NA:NA:NA"),"annotation"] = "Other" 
interaction[grep("CDS", interaction$cds),"annotation"] = "CDS"
interaction[which(is.na(interaction$annotation) & interaction$UTR5 == "UTR5"),"annotation"] = "UTR5"
interaction[which(is.na(interaction$annotation) & interaction$UTR3 == "UTR3"),"annotation"] = "UTR3"
interaction[which(is.na(interaction$annotation) & interaction$intron == "Intron"),"annotation"] = "Intron"
interaction[which(is.na(interaction$annotation) & interaction$promoter == "Promoter"),"annotation"] = "Promoter"

# Mapping cell sites to nearest gene
dA = distanceToNearest(grinteraction, geneMapGR)
interaction$nearestSymbol = geneMapGR$Symbol[subjectHits(dA)]
interaction$nearestID = names(geneMapGR)[subjectHits(dA)]
interaction$distToGene = mcols(dA)$distance
interaction$EntrezID = geneMapGR$EntrezID[subjectHits(dA)]
interaction$regionID = paste0(interaction$seqnames,":",interaction$start,"-", interaction$end)
interaction$sig = ifelse(interaction$fwer<=0.05, "FWER < 0.05", "FWER > 0.05")
interaction$Dir = ifelse(interaction$value<0, "neg", "pos")
dtinteraction = data.table(interaction)

### All together ###

DMR = list(CellType = cell, Age = age, Interaction = interaction)
lapply(DMR, head)

## Fix pheno data
pd$Race[pd$Race== "CAUC "] <- 'CAUC'
pd$Sex[pd$Sex == " M"] <- 'M'
pd$RIN <- as.numeric(gsub(" ", "", pd$RIN))
pd$pg.DNA.nuclei.input <- as.numeric(pd$pg.DNA.nuclei.input)
pd$Reads <- as.numeric(pd$Reads)
pd$Percent.GreaterThan.Q30 <- as.numeric(pd$Percent.GreaterThan.Q30)


save(DMR, geneMap, pd, file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/DMR/DMR_objects.rda")

names = c("CellType", "Age", "CT_Age_Interaction")
for (i in (1:length(DMR))){
  tmp = DMR[[i]]
  write.csv(tmp[which(tmp$fwer<=0.05),], 
            file=paste0("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/DMR/", names[i], "/", names(DMR)[i],"_DMRs_fwer_0.05.csv"),
            quote=F)
}





