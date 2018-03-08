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
cell$regionID = paste0(cell$chr,":",cell$start,"-", cell$end)

# Mapping cell sites to nearest gene
grcell = makeGRangesFromDataFrame(cell, keep.extra.columns = T)
geneMapGR = makeGRangesFromDataFrame(geneMap, start.field="Start",end.field="End",strand.field="Strand",keep=TRUE)
dA = distanceToNearest(grcell, geneMapGR)
oo = findOverlaps(grcell, geneMapGR)
x = cbind(cell[queryHits(oo),], nearestSymbol = geneMapGR$Symbol[subjectHits(oo)], nearestID = names(geneMapGR)[subjectHits(oo)],
          EntrezID = geneMapGR$EntrezID[subjectHits(oo)], distToGene = 0)

cell$nearestSymbol = geneMapGR$Symbol[subjectHits(dA)]
cell$nearestID = names(geneMapGR)[subjectHits(dA)]
cell$EntrezID = geneMapGR$EntrezID[subjectHits(dA)]
cell$distToGene = mcols(dA)$distance
cell = rbind(x, cell[which(!cell$regionID %in% x$regionID),])


# Annotate DMRs to features in the genome
txdb = loadDb("/dcl01/lieber/ajaffe/Amanda/annotation_objects/gencode.v25lift37.annotation.sqlite")
islands = read.table("/dcl01/lieber/ajaffe/Amanda/annotation_objects/cpgIslandExt.hg19.txt", sep="\t", header = T)
rpmsk = read.table("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/RepeatMasker_genomewide_HG19.txt", header=T)
features = list(CDS = cdsBy(txdb, by="tx", use.names=T), Introns = intronsByTranscript(txdb, use.names=T), 
                UTR5 = fiveUTRsByTranscript(txdb, use.names=T), UTR3 = threeUTRsByTranscript(txdb, use.names=T))
features = lapply(features, function(x) unlist(x, recursive = TRUE, use.names = TRUE))
for (i in 1:length(features)){
  tmp = features[[i]]
  tmp$TxID = names(tmp)
  features[[i]] = tmp
}
features = c(features, 
             rpmskgr = makeGRangesFromDataFrame(rpmsk, seqnames.field="genoName",start.field="genoStart",end.field="genoEnd",keep.extra.columns=TRUE),
             islands = makeGRangesFromDataFrame(islands, keep.extra.columns = T, start.field = "chromStart", end.field = "chromEnd"),
             promoters = promoters(txdb, upstream=2000, downstream=200))
lapply(features, head)

grcell = makeGRangesFromDataFrame(cell, keep.extra.columns = T)
annotation = lapply(features, function(y) findOverlaps(grcell, y))

grcell$rnum = 1:length(grcell)
grcell$cds = ifelse(grcell$rnum %in% queryHits(annotation[["CDS"]]), "CDS", NA)
grcell$intron = ifelse(grcell$rnum %in% queryHits(annotation[["Introns"]]), "Intron", NA)
grcell$UTR5 = ifelse(grcell$rnum %in% queryHits(annotation[["UTR5"]]), "5'UTR", NA)
grcell$UTR3 = ifelse(grcell$rnum %in% queryHits(annotation[["UTR3"]]), "3'UTR", NA)
grcell$islands = ifelse(grcell$rnum %in% queryHits(annotation[["islands"]]), "CpG-Island", "non-Island")
grcell$promoter = ifelse(grcell$rnum %in% queryHits(annotation[["promoters"]]), "Promoter", NA)
grcell$anno = paste0(grcell$cds,":",grcell$intron, ":", grcell$UTR5, ":", grcell$UTR3, ":", grcell$promoter)
cell = as.data.frame(grcell)

x = data.frame(cell[queryHits(annotation$rpmskgr),], repeats = as.character(features$rpmskgr$repClass)[subjectHits(annotation$rpmskgr)])
cell = rbind(x, data.frame(cell[-unique(queryHits(annotation$rpmskgr)),], repeats = "No repeats"))
cell[which(cell$anno == "NA:NA:NA:NA:NA"),"annotation"] = "Intergenic" 
cell[grep("CDS", cell$cds),"annotation"] = "CDS"
cell[which(is.na(cell$annotation) & cell$UTR5 == "5'UTR"),"annotation"] = "5'UTR"
cell[which(is.na(cell$annotation) & cell$UTR3 == "3'UTR"),"annotation"] = "3'UTR"
cell[which(is.na(cell$annotation) & cell$intron == "Intron"),"annotation"] = "Intron"
cell[which(is.na(cell$annotation) & cell$promoter == "Promoter"),"annotation"] = "Promoter"

cell$sig = ifelse(cell$fwer<=0.05, "FWER < 0.05", "FWER > 0.05")
cell$Dir = ifelse(cell$value<0, "neg", "pos")
dtcell = data.table(cell)


### Age ###

load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/bumphunting/bumps_bsseqSmooth_Neuron_age_250_perm.Rdata")

## How many regions are DMRs? (model: ~ cell + age)

age = bumps[[1]]
dim(age) # 52790     14
dim(age[which(age$fwer<=0.05),]) # 129    14
age$regionID = paste0(age$chr,":",age$start,"-", age$end)


# Mapping age sites to nearest gene
grage = makeGRangesFromDataFrame(age, keep.extra.columns = T)
dA = distanceToNearest(grage, geneMapGR)
oo = findOverlaps(grage, geneMapGR)
x = cbind(age[queryHits(oo),], nearestSymbol = geneMapGR$Symbol[subjectHits(oo)], nearestID = names(geneMapGR)[subjectHits(oo)],
          EntrezID = geneMapGR$EntrezID[subjectHits(oo)], distToGene = 0)

age$nearestSymbol = geneMapGR$Symbol[subjectHits(dA)]
age$nearestID = names(geneMapGR)[subjectHits(dA)]
age$EntrezID = geneMapGR$EntrezID[subjectHits(dA)]
age$distToGene = mcols(dA)$distance
age = rbind(x, age[which(!age$regionID %in% x$regionID),])


# Annotate DMRs to features in the genome

grage = makeGRangesFromDataFrame(age, keep.extra.columns = T)
annotation = lapply(features, function(y) findOverlaps(grage, y))

grage$rnum = 1:length(grage)
grage$cds = ifelse(grage$rnum %in% queryHits(annotation[["CDS"]]), "CDS", NA)
grage$intron = ifelse(grage$rnum %in% queryHits(annotation[["Introns"]]), "Intron", NA)
grage$UTR5 = ifelse(grage$rnum %in% queryHits(annotation[["UTR5"]]), "5'UTR", NA)
grage$UTR3 = ifelse(grage$rnum %in% queryHits(annotation[["UTR3"]]), "3'UTR", NA)
grage$islands = ifelse(grage$rnum %in% queryHits(annotation[["islands"]]), "CpG-Island", "non-Island")
grage$promoter = ifelse(grage$rnum %in% queryHits(annotation[["promoters"]]), "Promoter", NA)
grage$anno = paste0(grage$cds,":",grage$intron, ":", grage$UTR5, ":", grage$UTR3, ":", grage$promoter)
age = as.data.frame(grage)

x = data.frame(age[queryHits(annotation$rpmskgr),], repeats = as.character(features$rpmskgr$repClass)[subjectHits(annotation$rpmskgr)])
age = rbind(x, data.frame(age[-unique(queryHits(annotation$rpmskgr)),], repeats = "No repeats"))
age[which(age$anno == "NA:NA:NA:NA:NA"),"annotation"] = "Intergenic" 
age[grep("CDS", age$cds),"annotation"] = "CDS"
age[which(is.na(age$annotation) & age$UTR5 == "5'UTR"),"annotation"] = "5'UTR"
age[which(is.na(age$annotation) & age$UTR3 == "3'UTR"),"annotation"] = "3'UTR"
age[which(is.na(age$annotation) & age$intron == "Intron"),"annotation"] = "Intron"
age[which(is.na(age$annotation) & age$promoter == "Promoter"),"annotation"] = "Promoter"

age$sig = ifelse(age$fwer<=0.05, "FWER < 0.05", "FWER > 0.05")
age$Dir = ifelse(age$value<0, "neg", "pos")
dtage = data.table(age)


### Interaction ###

load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/bumphunting/bumps_bsseqSmooth_Neuron_interaction_250_perm.Rdata")


## How many regions are DMRs? (model: ~ age + cell + age:cell)

interaction = bumps[[1]]
dim(interaction) # 26282     14
dim(interaction[which(interaction$fwer<=0.05),]) # 2178    14
interaction$regionID = paste0(interaction$chr,":",interaction$start,"-", interaction$end)

# Mapping interaction sites to nearest gene
grinteraction = makeGRangesFromDataFrame(interaction, keep.extra.columns = T)
dA = distanceToNearest(grinteraction, geneMapGR)
oo = findOverlaps(grinteraction, geneMapGR)
x = cbind(interaction[queryHits(oo),], nearestSymbol = geneMapGR$Symbol[subjectHits(oo)], nearestID = names(geneMapGR)[subjectHits(oo)],
          EntrezID = geneMapGR$EntrezID[subjectHits(oo)], distToGene = 0)

interaction$nearestSymbol = geneMapGR$Symbol[subjectHits(dA)]
interaction$nearestID = names(geneMapGR)[subjectHits(dA)]
interaction$EntrezID = geneMapGR$EntrezID[subjectHits(dA)]
interaction$distToGene = mcols(dA)$distance
interaction = rbind(x,interaction[which(!interaction$regionID %in% x$regionID),])


# Annotate DMRs to features in the genome

grinteraction = makeGRangesFromDataFrame(interaction, keep.extra.columns = T)
annotation = lapply(features, function(y) findOverlaps(grinteraction, y))

grinteraction$rnum = 1:length(grinteraction)
grinteraction$cds = ifelse(grinteraction$rnum %in% queryHits(annotation[["CDS"]]), "CDS", NA)
grinteraction$intron = ifelse(grinteraction$rnum %in% queryHits(annotation[["Introns"]]), "Intron", NA)
grinteraction$UTR5 = ifelse(grinteraction$rnum %in% queryHits(annotation[["UTR5"]]), "5'UTR", NA)
grinteraction$UTR3 = ifelse(grinteraction$rnum %in% queryHits(annotation[["UTR3"]]), "3'UTR", NA)
grinteraction$islands = ifelse(grinteraction$rnum %in% queryHits(annotation[["islands"]]), "CpG-Island", "non-Island")
grinteraction$promoter = ifelse(grinteraction$rnum %in% queryHits(annotation[["promoters"]]), "Promoter", NA)
grinteraction$anno = paste0(grinteraction$cds,":",grinteraction$intron, ":", grinteraction$UTR5, ":", grinteraction$UTR3, ":", grinteraction$promoter)
interaction = as.data.frame(grinteraction)

x = data.frame(interaction[queryHits(annotation$rpmskgr),], repeats = as.character(features$rpmskgr$repClass)[subjectHits(annotation$rpmskgr)])
interaction = rbind(x, data.frame(interaction[-unique(queryHits(annotation$rpmskgr)),], repeats = "No repeats"))
interaction[which(interaction$anno == "NA:NA:NA:NA:NA"),"annotation"] = "Intergenic" 
interaction[grep("CDS", interaction$cds),"annotation"] = "CDS"
interaction[which(is.na(interaction$annotation) & interaction$UTR5 == "5'UTR"),"annotation"] = "5'UTR"
interaction[which(is.na(interaction$annotation) & interaction$UTR3 == "3'UTR"),"annotation"] = "3'UTR"
interaction[which(is.na(interaction$annotation) & interaction$intron == "Intron"),"annotation"] = "Intron"
interaction[which(is.na(interaction$annotation) & interaction$promoter == "Promoter"),"annotation"] = "Promoter"

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





