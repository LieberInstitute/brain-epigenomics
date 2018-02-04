library(ggplot2)
library(GenomicRanges)
library(SummarizedExperiment)
library(PWMEnrich)


load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/PWMEnrich/DMR_PWMEnrich_objects.rda")


## Get homogenate gene level data 

load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/brainseq_pipeline/polyA_unstranded/rse_gene_polyA_dlpfc_n41.Rdata")

homRPKM = getRPKM(rse_gene)
hompd = colData(rse_gene)
geneMap = rowData(rse_gene)

## Extract TF genes

x = groupReport(all_split[[1]])
noW = unique(x$target)[-grep("UW", unique(x$target))]
length(noW) # 866
noW[-which(noW %in% geneMap$Symbol)]

ensIDs = c("HINFP1" = "ENSG00000172273", "MZF1_1-4" = "ENSG00000099326", "Myf" = "ENSG00000111046", "HIF1A::ARNT" = "ENSG00000100644", 
           "MZF1_5-13" = "ENSG00000099326", "MYC::MAX" = "ENSG00000136997", "RXRA::VDR" = "ENSG00000186350", "Trp53" = "ENSG00000141510", 
           "ZNF306" = "ENSG00000189298", "MIZF" = "ENSG00000172273", "EWSR1-FLI1" = "ENSG00000182944", "BHLHB3" = "ENSG00000123095",
           "BHLHB2" = "ENSG00000134107", "TLX1::NFIC" = "ENSG00000107807", "RORA_1" = "ENSG00000069667", "NR1H2::RXRA" = "ENSG00000131408", 
           "Trp73" = "ENSG00000078900", "JUND (var.2)" = "ENSG00000130522", "TAL1::TCF3" = "ENSG00000162367", "RXR::RAR_DR5" = "ENSG00000186350",
           "Pax6" = "ENSG00000007372", "DUX4" = "ENSG00000260596", "ZNF238" = "ENSG00000179456", "SMAD2::SMAD3::SMAD4" = "ENSG00000175387",
           "JUN (var.2)" = "ENSG00000177606", "NFE2::MAF" = "ENSG00000123405", "BATF::JUN" = "ENSG00000156127", "RAXL1" = "ENSG00000173976", 
           "CART1" = "ENSG00000180318", "ZNF435" = "ENSG00000196812", "TAL1::GATA1" = "ENSG00000162367",  "RORA_2" = "ENSG00000069667", 
           "Sox4" = "ENSG00000124766", "RBM8A" = "ENSG00000265241", "POU5F1P1" = "ENSG00000212993", "AP1" = "ENSG00000248394",
           "STAT2::STAT1" = "ENSG00000170581", "Oct-1" = "ENSG00000143190", "SMCR7L" = "ENSG00000100335", "HNF1B" = "ENSG00000275410",
           "GLTPD1" = "ENSG00000224051", "C9orf156" = "ENSG00000136932", "PRKRIR" = "ENSG00000137492", "C19orf40" = "ENSG00000131944", "MEF2BNB-MEF2B" = "ENSG00000254901")
# DUX4 not included in geneMap

geneIDs = c(geneMap[match(noW[which(noW %in% geneMap$Symbol)], geneMap$Symbol), "gencodeID"], 
            geneMap[match(ensIDs, geneMap$ensemblID), "gencodeID"])
length(unique(geneIDs)) # 839
TFhomRPKM = data.frame(homRPKM[which(rownames(homRPKM) %in% geneIDs),])
TFhomRPKM = as.matrix(TFhomRPKM)

## Test cutoffs of expression

quantile(homRPKM)
#0%          25%          50%          75%         100% 
#0.000000e+00 0.000000e+00 4.469787e-02 9.757915e-01 7.486556e+04 

cutoff = data.frame(geneID = rownames(TFhomRPKM), Threshold = NA)
for (i in 1:nrow(TFhomRPKM)) {
  cutoff[i,"Threshold"] = ifelse(max(TFhomRPKM[i,])>5, "TRUE","FALSE")
}
table(cutoff$Threshold=="TRUE")

no = geneMap[which(geneMap$gencodeID %in% cutoff$geneID[cutoff$Threshold=="FALSE"]),"Symbol"]

df = t(TFhomRPKM)

theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

geneIDs = c(geneMap$gencodeID[which(geneMap$Symbol %in% un)], geneMap$gencodeID[which(geneMap$ensemblID %in% ensIDs)]

load("/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/Consortium/rse_gene_BrainSeq_Phase1_hg19_TopHat2_EnsemblV75.rda")
load("/dcl01/lieber/ajaffe/CellSorting/RNAseq_pipeline/rawCounts_CellSorting_July5_n12.rda")



