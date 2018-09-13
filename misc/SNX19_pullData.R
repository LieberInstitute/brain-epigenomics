library(GenomicRanges)

## look for differences in DNAm in SNX19 locus

load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/DMR/DMR_objects.rda")
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/bumphunting/rda/limma_Neuron_CpGs_minCov_3_ageInfo_dmrs.Rdata")
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/UMRs_LMRs_annotated.rda")


DMRgr = lapply(DMR, function(x) makeGRangesFromDataFrame(x[which(x$sig=="FWER < 0.05"),], keep.extra.columns = T))

dmrs = split(dmrs, dmrs$k6cluster_label)
names(dmrs) = c("Gr1","Gr2","Gr3","Gr4","Gr5","Gr6")
oo = lapply(dmrs, function(x) findOverlaps(x, makeGRangesFromDataFrame(DMR$Interaction)))
dmrs = lapply(oo, function(x) DMR$Interaction[subjectHits(x),])

DMRgr = c(DMRgr, list("CellType.N+" = DMRgr$CellType[which(DMRgr$CellType$Dir=="pos" & DMRgr$CellType$sig=="FWER < 0.05"),],
                      "CellType.N-" = DMRgr$CellType[which(DMRgr$CellType$Dir=="neg" & DMRgr$CellType$sig=="FWER < 0.05"),]),
          lapply(dmrs, function(x) makeGRangesFromDataFrame(x[which(x$sig=="FWER < 0.05"),], keep.extra.columns = T)))
DMRgr = lapply(DMRgr, reduce)

oo = lapply(DMRgr, function(x) findOverlaps(x, GRanges("chr11:130700235-130800234"))) # No overlaps

oo = lapply(DMRgr, function(x) findOverlaps(x, GRanges("chr11:129718380-131718880")))
x = lapply(DMRgr, as.data.frame)
x = Map(cbind, x, regionID = lapply(x, function(y) paste0(y$seqnames,":",y$start,"-",y$end)))
x = mapply(function(data,subset) data[queryHits(subset),], x, oo, SIMPLIFY = F)
y = mapply(function(data, subset) data[which(data$regionID %in% subset$regionID),], c(DMR, list("CellType.N+" = DMR$CellType[which(DMR$CellType$Dir=="pos" & DMR$CellType$sig=="FWER < 0.05"),],
                                                                                                "CellType.N-" = DMR$CellType[which(DMR$CellType$Dir=="neg" & DMR$CellType$sig=="FWER < 0.05"),]),
                                                                                      lapply(dmrs, function(x) x[which(x$sig=="FWER < 0.05"),])), x, SIMPLIFY = F)
x = do.call(rbind, Map(cbind, y[elementNROWS(y)>0], Group = as.list(names(y[elementNROWS(y)>0]))))
rownames(x) = NULL
x$Group = gsub("Gr1", "G-:N+", x$Group)
x$Group = gsub("Gr2", "G0:N+", x$Group)
x$Group = gsub("Gr3", "G0:N-", x$Group)
x$Group = gsub("Gr4", "G+:N0", x$Group)
x$Group = gsub("Gr5", "G+:N-", x$Group)
x$Group = gsub("Gr6", "G-:N0", x$Group)
x = unique(x)
elementNROWS(split(x, x$Group))
# CellType CellType.N- CellType.N+       G+:N0 Interaction 
#       40           5          35           2           2
x = x[-which(x$Group %in% c("CellType", "Interaction")),]
write.csv(x, quote = F, file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/DMR/SNX19_DMRs_2MBwindow.csv")


### Check LMRs

oo = lapply(lDMR, function(x) findOverlaps(makeGRangesFromDataFrame(x), GRanges("chr11:130700235-130800234")))
x = mapply(function(data,subset) data[queryHits(subset),], lDMR, oo, SIMPLIFY = F)
x = do.call(rbind, Map(cbind, x[elementNROWS(x)>0], SampleID = as.list(names(x[elementNROWS(x)>0]))))
rownames(x) = NULL
x = unique(x)
x$CellType = pd[match(x$SampleID, rownames(pd)),"Cell.Type"]
x$age = pd[match(x$SampleID, rownames(pd)),"Age"]
x$Age.Bin = pd[match(x$SampleID, rownames(pd)),"Age.Bin"]
x = x[,!colnames(x) %in% c("tog","dmr")]
write.csv(x, quote = F, file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/DMR/SNX19_LMRs_100KBwindow.csv")

oo = lapply(lDMR, function(x) findOverlaps(makeGRangesFromDataFrame(x), GRanges("chr11:129718380-131718880")))
x = mapply(function(data,subset) data[queryHits(subset),], lDMR, oo, SIMPLIFY = F)
x = do.call(rbind, Map(cbind, x[elementNROWS(x)>0], SampleID = as.list(names(x[elementNROWS(x)>0]))))
rownames(x) = NULL
x = unique(x)
x$CellType = pd[match(x$SampleID, rownames(pd)),"Cell.Type"]
x$age = as.numeric(pd[match(x$SampleID, rownames(pd)),"Age"])
x$Age.Bin = pd[match(x$SampleID, rownames(pd)),"Age.Bin"]
x = x[,!colnames(x) %in% c("tog","dmr")]
write.csv(x, quote = F, file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/DMR/SNX19_LMRs_2MBwindow.csv")

## genes in the window

oo = findOverlaps(makeGRangesFromDataFrame(geneMap), GRanges("chr11:129718380-131718880"))
rownames(geneMap[queryHits(oo),])
x = read.csv("./Downloads/LIBD_methylation_explorer_selection_2018-05-07 21_30_02.csv")
head(x)
x = x[which(x$feature_id %in% c("ENSG00000151715.7_1","ENSG00000170322.14_2","ENSG00000170325.14_2","ENSG00000233220.3_1","ENSG00000255057.1_1","ENSG00000255262.3_1",
                                "ENSG00000255535.2_1","ENSG00000084234.16_2","ENSG00000244451.1_1","ENSG00000149418.10_1","ENSG00000196323.11_2","ENSG00000255358.1_1", 
                                "ENSG00000255220.1_1","ENSG00000175773.12_1","ENSG00000134917.9_1","ENSG00000166106.3_1","ENSG00000255155.1_1","ENSG00000236616.3_1", 
                                "ENSG00000175728.4_1","ENSG00000255352.1_1","ENSG00000254842.6_1","ENSG00000255455.2_1","ENSG00000120451.10_1","ENSG00000236129.1_1", 
                                "ENSG00000227125.1_1","ENSG00000237612.1_1","ENSG00000231698.2_1","ENSG00000182667.14_2","ENSG00000237654.5_1","ENSG00000236082.2_1", 
                                "ENSG00000224795.1_1","ENSG00000271624.1_1","ENSG00000267981.1","ENSG00000215565.2","ENSG00000272412.1","ENSG00000242673.2",   
                                "ENSG00000221516.1","ENSG00000252351.1")),]
write.csv(x, quote = F, file = "./Downloads/SNX19_splicing_2MBwindow.csv")




