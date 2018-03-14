library(ggplot2)
library(GenomicRanges)
library(SummarizedExperiment)
library(PWMEnrich)


load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/PWMEnrich/DMR_PWMEnrich_objects.rda")
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/sorted_nuclear_RNA/DE_limma_results_objects.rda")

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
prepostTFhomRPKM = data.frame(homRPKM[which(rownames(homRPKM) %in% geneIDs),])
prepostTFhomRPKM = as.matrix(prepostTFhomRPKM)

## Test cutoffs of expression

quantile(homRPKM)
#0%          25%          50%          75%         100% 
#0.000000e+00 0.000000e+00 4.469787e-02 9.757915e-01 7.486556e+04 

cutoff = data.frame(geneID = rownames(prepostTFhomRPKM), Threshold = NA)
for (i in 1:nrow(prepostTFhomRPKM)) {
  cutoff[i,"Threshold"] = ifelse(max(prepostTFhomRPKM[i,])>1, "TRUE","FALSE")
}
table(cutoff$Threshold=="TRUE")
#FALSE  TRUE 
#188   650 

# Excluded TFs
geneMap[which(geneMap$gencodeID %in% cutoff$geneID[cutoff$Threshold=="FALSE"]),"Symbol"]
#[1] "TP73"      "RUNX3"     "DMBX1"     "FOXD2"     "GLIS1"     "FOXD3"    
#[7] "BARHL2"    "ALX3"      "TBX15"     "HIST2H2AB" "LMX1A"     "TBX19"    
#[13] "ELF3"      "MYOG"      "IL24"      "IRF6"      "MIXL1"     "REL"      
#[19] "FIGLA"     "VAX2"      "NOTO"      "LBX2"      "EN1"       "EVX2"     
#[25] "HOXD13"    "HOXD12"    "HOXD11"    "HOXD8"     "HOXD3"     "BOLL"     
#[31] "FEV"       "PAX3"      "GBX2"      "HESX1"     "POU1F1"    "SOX14"    
#[37] "SHOX2"     "TP63"      "HMX1"      "NKX3-2"    "PHOX2B"    "NKX6-1"   
#[43] "FOSL1P1"   "POU4F2"    "IRX2"      "DDX4"      "FOXD1"     "PITX1"    
#[49] "POU4F3"    "CDX1"      "HAND1"     "FOXI1"     "PROP1"     "IRF4"     
#[55] "RREB1"     "TFAP2A"    "GCM2"      "HIST1H2BN" "TULP1"     "RUNX2"    
#[61] "TFAP2B"    "GCM1"      "DDX43"     "OLIG3"     "ESR1"      "PLG"      
#[67] "T"         "UNCX"      "MEOX2"     "HOXA1"     "HOXA2"     "HOXA10"   
#[73] "HOXA13"    "EVX1"      "TBX20"     "PGAM2"     "TFEC"      "PAX4"     
#[79] "KLF14"     "GBX1"      "EN2"       "MNX1"      "GATA4"     "SLC18A1"  
#[85] "NKX3-1"    "HNF4G"     "ESRP1"     "POU5F1B"   "MAFA"      "FOXH1"    
#[91] "PAX5"      "BARX1"     "LMX1B"     "GATA3"     "PTF1A"     "DRGX"     
#[97] "A1CF"      "NKX2-3"    "PAX2"      "TLX1"      "PITX3"     "HMX3"     
#[103] "HMX2"      "VENTX"     "MYOD1"     "ELF5"      "EHF"       "FGF19"    
#[109] "PHOX2A"    "POU2F3"    "BSX"       "PDE6H"     "VDR"       "HOXC13"   
#[115] "HOXC12"    "HOXC11"    "HOXC10"    "NFE2"      "MYF6"      "ALX1"     
#[121] "SPIC"      "TBX5"      "HNF1A"     "ZNF26"     "PDX1"      "CDX2"     
#[127] "POU4F1"    "CEBPE"     "NKX2-8"    "PAX9"      "FOXA1"     "ESR2"     
#[133] "ZNF410"    "VSX2"      "BATF"      "ESRRB"     "GSC"       "ONECUT1"  
#[139] "FOXB1"     "CELF6"     "ISL2"      "MCTP2"     "IRX5"      "ESRP2"    
#[145] "FOXL1"     "HNF1B"     "MEOX1"     "HOXB2"     "HOXB3"     "HOXB5"    
#[151] "HOXB9"     "HOXB13"    "DLX4"      "TBX4"      "RAX"       "NFATC1"   
#[157] "ONECUT3"   "RAX2"      "MEF2B"     "CEBPA"     "SPIB"      "DPRX"     
#[163] "DUXA"      "ZSCAN4"    "PAX1"      "VSX1"      "HNF4A"     "WISP2"    
#[169] "GATA5"     "U2AF1"     "GSC2"      "TBX1"      "ISX"       "SOX10"    
#[175] "SHOX"      "XG"        "YY2"       "DDX53"     "SSX3"      "FOXP3"    
#[181] "SSX2"      "AR"        "TGIF2LX"   "ESX1"      "ELF4"      "MAGEA8"   
#[187] "SRY"       "HSFY2"    

no = geneMap[which(geneMap$gencodeID %in% cutoff$geneID[cutoff$Threshold=="FALSE"]),"gencodeID"]


prepostTargettoGeneID = c(geneMap[match(noW[which(noW %in% geneMap$Symbol)], geneMap$Symbol), "gencodeID"], 
            geneMap[match(ensIDs, geneMap$ensemblID), "gencodeID"])
names(prepostTargettoGeneID) = c(geneMap[match(noW[which(noW %in% geneMap$Symbol)], geneMap$Symbol), "Symbol"], names(ensIDs))

prepostTargettoGeneID = prepostTargettoGeneID[!prepostTargettoGeneID %in% no]
prepostTFhomRPKM = prepostTFhomRPKM[-which(rownames(prepostTFhomRPKM) %in% no),]


## Get sorted nuclear RNA TF results

prepostTFnucres = nucRNAres[which(rownames(nucRNAres) %in% prepostTargettoGeneID),]



### Test cutoffs of expression in postnatal only

postRPKM = prepostTFhomRPKM[,which(colnames(prepostTFhomRPKM) %in% rownames(hompd[which(hompd$Age>0),]))]

cutoff = data.frame(geneID = rownames(postRPKM), Threshold = NA)
for (i in 1:nrow(postRPKM)) {
  cutoff[i,"Threshold"] = ifelse(max(postRPKM[i,])>1, "TRUE","FALSE")
}
table(cutoff$Threshold=="TRUE")
#FALSE  TRUE 
#   32   618 

# Additional excluded TFs
geneMap[which(geneMap$gencodeID %in% cutoff$geneID[cutoff$Threshold=="FALSE"]),"Symbol"]
# [1] "MORN1"   "E2F2"    "LHX4"    "LHX9"    "DTL"     "ZNF695"  "GLI2"   
# [8] "EOMES"   "GLYCTK"  "GSX2"    "TFAP2D"  "PTCD1"   "PLAG1"   "VAX1"   
# [15] "E2F8"    "ALX4"    "HMGA2"   "E2F7"    "GSX1"    "NRL"     "OTX2"   
# [22] "SCAND2P" "FOXC2"   "BRCA1"   "TBX21"   "ONECUT2" "FAAP24"  "ETV2"   
# [29] "MYBL2"   "UBE2V1"  "TFAP2C"  "BHLHE23"

no = geneMap[which(geneMap$gencodeID %in% cutoff$geneID[cutoff$Threshold=="FALSE"]),"gencodeID"]

PostnataltargettogeneID = prepostTargettoGeneID[-which(prepostTargettoGeneID %in% no)]

postRPKM = postRPKM[-which(rownames(postRPKM) %in% no),]


## Get sorted nuclear RNA TF results

TFnucresPostnatal = nucRNAres[which(rownames(nucRNAres) %in% PostnataltargettogeneID),]

save(PostnataltargettogeneID, postRPKM, TFnucresPostnatal,prepostTargettoGeneID, prepostTFhomRPKM, prepostTFnucres, hompd, geneMap,  
     file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/PWMEnrich/filtered_TF_list_expressionInfo.rda")



