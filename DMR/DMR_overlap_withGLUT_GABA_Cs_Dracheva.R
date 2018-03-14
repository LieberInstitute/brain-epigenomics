library(GenomicRanges)
library(bumphunter)
library(ggplot2)

load('/dcl01/lieber/ajaffe/lab/brain-epigenomics/bumphunting/BSobj_bsseqSmooth_Neuron_minCov_3.Rdata')
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/DMR/DMR_objects.rda")

# load dmCs by neuronal subtype
gabglu.450K = openxlsx::read.xlsx('/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/Kozlenkov_Dracheva_NucAcidsRes_2016_CpGmeth_450K.xlsx')
gabglu.rrbs = openxlsx::read.xlsx('/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/Kozlenkov_Dracheva_NucAcidsRes_2016_CpGmeth_ERRBS.xlsx')
gabglu.CHrrbs = openxlsx::read.xlsx('/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/Kozlenkov_Dracheva_NucAcidsRes_2016_CHmeth_ERRBS.xlsx')
colnames(gabglu.450K) = gabglu.450K[1,]
gabglu.450K = gabglu.450K[-1,]
gabglu.450K$FDR.GABA.GLU = as.numeric(gabglu.450K$FDR.GABA.GLU)
gabglu.450Kdelta.GABA.GLU = as.numeric(gabglu.450K$delta.GABA.GLU)
colnames(gabglu.rrbs) = gabglu.rrbs[1,]
gabglu.rrbs = gabglu.rrbs[-1,]
gabglu.rrbs$qvalue = as.numeric(gabglu.rrbs$qvalue)
colnames(gabglu.CHrrbs) = gabglu.CHrrbs[1,]
gabglu.CHrrbs = gabglu.CHrrbs[-1,]
gabglu.CHrrbs$qvalue = as.numeric(gabglu.CHrrbs$qvalue)

gabglu = list("450K" = makeGRangesFromDataFrame(gabglu.450K[which(gabglu.450K$FDR.GABA.GLU<=0.05),], seqnames.field="CHR",start.field="POSITION", end.field="POSITION", keep.extra.columns=T),
			  "rrbs" = makeGRangesFromDataFrame(gabglu.rrbs[which(gabglu.rrbs$qvalue<=0.05),], seqnames.field="CHR",start.field="POSITION", end.field="POSITION", keep.extra.columns=T),
			  "CHrrbs" = makeGRangesFromDataFrame(gabglu.CHrrbs[which(gabglu.CHrrbs$qvalue<=0.05),], seqnames.field="CHR",start.field="POSITION", end.field="POSITION", strand.field="STRAND", keep.extra.columns=T),
			  "450K.upGABA" = makeGRangesFromDataFrame(gabglu.450K[which(gabglu.450K$FDR.GABA.GLU<=0.05 & gabglu.450K$delta.GABA.GLU < 0),], seqnames.field="CHR",start.field="POSITION", end.field="POSITION", keep.extra.columns=T),
			  "rrbs.upGABA" = makeGRangesFromDataFrame(gabglu.rrbs[which(gabglu.rrbs$qvalue<=0.05 & gabglu.rrbs$DiffMeth=="DM.GABA"),], seqnames.field="CHR",start.field="POSITION", end.field="POSITION", keep.extra.columns=T),
			  "CHrrbs.upGABA" = makeGRangesFromDataFrame(gabglu.CHrrbs[which(gabglu.CHrrbs$qvalue<=0.05 & gabglu.CHrrbs$DiffMeth=="GABA.DM"),], 
			  											 seqnames.field="CHR",start.field="POSITION", end.field="POSITION", strand.field="STRAND", keep.extra.columns=T),
			  "450K.upGlut" = makeGRangesFromDataFrame(gabglu.450K[which(gabglu.450K$FDR.GABA.GLU<=0.05 & gabglu.450K$delta.GABA.GLU > 0),], seqnames.field="CHR",start.field="POSITION", end.field="POSITION", keep.extra.columns=T),
			  "rrbs.upGlut" = makeGRangesFromDataFrame(gabglu.rrbs[which(gabglu.rrbs$qvalue<=0.05 & gabglu.rrbs$DiffMeth=="DM.GLU"),], seqnames.field="CHR",start.field="POSITION", end.field="POSITION", keep.extra.columns=T),
			  "CHrrbs.upGlut" = makeGRangesFromDataFrame(gabglu.CHrrbs[which(gabglu.CHrrbs$qvalue<=0.05 & gabglu.CHrrbs$DiffMeth=="GLU.DM"),], seqnames.field="CHR",start.field="POSITION", end.field="POSITION", 				
			  											 strand.field="STRAND", keep.extra.columns=T))


# Identify all CpG clusters in the genome
gr <- granges(BSobj)
cl = clusterMaker( chr = as.character(seqnames(gr)), 
                   pos = start(gr),   maxGap = 1000)
gr.clusters = split(gr, cl)
gr.clusters = unlist(range(gr.clusters))
df.clusters = as.data.frame(gr.clusters)

# Find overlaps with DMRS in all three models
DMRgr = lapply(DMR, function(x) makeGRangesFromDataFrame(x[which(x$fwer<=0.05),], keep.extra.columns=T))

ooDMR = lapply(DMRgr, function(x) findOverlaps(gr.clusters, x))
ooGG = lapply(gabglu, function(x) findOverlaps(gr.clusters,x))

df.clusters$regionID = paste0(df.clusters$seqnames,":",df.clusters$start,"-",df.clusters$end)
df.clusters$rnum = 1:length(gr.clusters)
df.clusters$CellType = ifelse(df.clusters$rnum %in% queryHits(ooDMR $CellType), "CellType","no")
df.clusters$Age = ifelse(df.clusters$rnum %in% queryHits(ooDMR $Age), "Age","no")
df.clusters$Interaction = ifelse(df.clusters$rnum %in% queryHits(ooDMR $Interaction), "Interaction","no")
df.clusters$DMRs = paste(df.clusters$CellType, df.clusters$Age, df.clusters$Interaction, sep=":")
df.clusters$"450K" = ifelse(df.clusters$rnum %in% queryHits(ooGG$"450K"), "450K","no")
df.clusters$rrbs = ifelse(df.clusters$rnum %in% queryHits(ooGG$"rrbs"), "rrbs","no")
df.clusters$CHrrbs = ifelse(df.clusters$rnum %in% queryHits(ooGG$"CHrrbs"), "CHrrbs","no")
df.clusters$"450K.upGABA" = ifelse(df.clusters$rnum %in% queryHits(ooGG$"450K.upGABA"), "450K.upGABA","no")
df.clusters$rrbs.upGABA = ifelse(df.clusters$rnum %in% queryHits(ooGG$"rrbs.upGABA"), "rrbs.upGABA","no")
df.clusters$CHrrbs.upGABA = ifelse(df.clusters$rnum %in% queryHits(ooGG$"CHrrbs.upGABA"), "CHrrbs.upGABA","no")
df.clusters$"450K.upGlut" = ifelse(df.clusters$rnum %in% queryHits(ooGG$"450K.upGlut"), "450K.upGlut","no")
df.clusters$rrbs.upGlut = ifelse(df.clusters$rnum %in% queryHits(ooGG$"rrbs.upGlut"), "rrbs.upGlut","no")
df.clusters$CHrrbs.upGlut = ifelse(df.clusters$rnum %in% queryHits(ooGG$"CHrrbs.upGlut"), "CHrrbs.upGlut","no")


## make contingency tables

gg = c("450K","rrbs","CHrrbs","450K.upGABA","rrbs.upGABA","CHrrbs.upGABA","450K.upGlut","rrbs.upGlut","CHrrbs.upGlut")
tables = list()
for (i in 1:length(gg)){
	tables[[i]] = list(CellType = data.frame(c(nrow(df.clusters[df.clusters[,which(colnames(df.clusters)==gg[i])]!="no" & df.clusters$CellType=="CellType",]),
                                               nrow(df.clusters[df.clusters[,which(colnames(df.clusters)==gg[i])]!="no" & df.clusters$CellType=="no",])),
                                             c(nrow(df.clusters[df.clusters[,which(colnames(df.clusters)==gg[i])]=="no" & df.clusters$CellType=="CellType",]),
                                               nrow(df.clusters[df.clusters[,which(colnames(df.clusters)==gg[i])]=="no" & df.clusters$CellType=="no",])), row.names = c("YesCT","NoCT")),
              		    Age = data.frame(c(nrow(df.clusters[df.clusters[,which(colnames(df.clusters)==gg[i])]!="no" & df.clusters$Age=="Age",]),
                                           nrow(df.clusters[df.clusters[,which(colnames(df.clusters)==gg[i])]!="no" & df.clusters$Age=="no",])),
                          	             c(nrow(df.clusters[df.clusters[,which(colnames(df.clusters)==gg[i])]=="no" & df.clusters$Age=="Age",]),
                                           nrow(df.clusters[df.clusters[,which(colnames(df.clusters)==gg[i])]=="no" & df.clusters$Age=="no",])), row.names = c("YesAge","NoAge")),
                        Interaction = data.frame(c(nrow(df.clusters[df.clusters[,which(colnames(df.clusters)==gg[i])]!="no" & df.clusters$Interaction=="Interaction",]),
                                                   nrow(df.clusters[df.clusters[,which(colnames(df.clusters)==gg[i])]!="no" & df.clusters$Interaction=="no",])),
                                         	     c(nrow(df.clusters[df.clusters[,which(colnames(df.clusters)==gg[i])]=="no" & df.clusters$Interaction=="Interaction",]),
                                                   nrow(df.clusters[df.clusters[,which(colnames(df.clusters)==gg[i])]=="no" & df.clusters$Interaction=="no",])), row.names = c("YesInt","NoInt")))
	colnames(tables[[i]][["CellType"]]) = colnames(tables[[i]][["Age"]]) = colnames(tables[[i]][["Interaction"]]) = c(paste0("Yes",gg[i]),paste0("No",gg[i]))
}
names(tables) = gg

fisher = lapply(unlist(tables, recursive=F), fisher.test)
df = t(data.frame(pval = unlist(lapply(fisher, function(x) x$p.value)),
                  OR = unlist(lapply(fisher, function(x) x$estimate))))

#                                   pval        OR
#450K.CellType              0.000000e+00  30.45269
#450K.Age                  6.783447e-124 188.51915
#450K.Interaction           0.000000e+00  51.64631
#rrbs.CellType              0.000000e+00  30.88139
#rrbs.Age                  1.209423e-118 124.06565
#rrbs.Interaction           0.000000e+00  44.79600
#CHrrbs.CellType            0.000000e+00  22.13099
#CHrrbs.Age                1.475627e-100  74.02029
#CHrrbs.Interaction         0.000000e+00  30.28603
#450K.upGABA.CellType       0.000000e+00  27.02892
#450K.upGABA.Age           6.350689e-114 104.56721
#450K.upGABA.Interaction    0.000000e+00  40.42791
#rrbs.upGABA.CellType       0.000000e+00  30.93308
#rrbs.upGABA.Age           2.650823e-128 145.53876
#rrbs.upGABA.Interaction    0.000000e+00  47.52948
#CHrrbs.upGABA.CellType     0.000000e+00  24.70284
#CHrrbs.upGABA.Age         7.886822e-111  93.15982
#CHrrbs.upGABA.Interaction  0.000000e+00  35.08002
#450K.upGlut.CellType       0.000000e+00  31.09919
#450K.upGlut.Age           7.119570e-122 137.08810
#450K.upGlut.Interaction    0.000000e+00  51.45362
#rrbs.upGlut.CellType       0.000000e+00  33.37056
#rrbs.upGlut.Age           8.737987e-122 126.48668
#rrbs.upGlut.Interaction    0.000000e+00  46.91774
#CHrrbs.upGlut.CellType     0.000000e+00  24.26651
#CHrrbs.upGlut.Age          3.479409e-94  61.99506
#CHrrbs.upGlut.Interaction  0.000000e+00  32.79693

write.csv(t(df),file="/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/DMR/DMR_Kozlenkov_Dracheva_NucAcidsRes_2016_GABA-GLU_dmC_Overlap.csv", quote=F)



## Count the number of DMRs that contain sites differentially methylated by neuronal subtype

ooCT = lapply(gabglu, function(x) findOverlaps(reduce(makeGRangesFromDataFrame(DMR$CellType[which(DMR$CellType$fwer<=0.05),])),x))
ooAge = lapply(gabglu, function(x) findOverlaps(reduce(makeGRangesFromDataFrame(DMR$Age[which(DMR$Age$fwer<=0.05),])),x))
ooInt = lapply(gabglu, function(x) findOverlaps(reduce(makeGRangesFromDataFrame(DMR$Interaction[which(DMR$Interaction$fwer<=0.05),])),x))

df = data.frame(Group = names(ooCT), count = c(unlist(lapply(ooCT, function(x) length(unique(queryHits(x))))),
											           unlist(lapply(ooAge, function(x) length(unique(queryHits(x))))),
											           unlist(lapply(ooInt, function(x) length(unique(queryHits(x)))))), 
			   percent = c(round(unlist(lapply(ooCT, function(x) length(unique(queryHits(x)))))/11179*100,2),
			   			  round(unlist(lapply(ooAge, function(x) length(unique(queryHits(x)))))/129*100,2),
			   			  round(unlist(lapply(ooInt, function(x) length(unique(queryHits(x)))))/2178*100,2)),
			   model = c(rep.int("Cell Type",9),rep.int("Age",9),rep.int("Interaction",9)), 
			   list = rep.int(c("450K","rrbs","CHrrbs"),9),
			   dir = rep.int(c(rep.int("Both",3),rep.int("More In\nGABA",3),rep.int("More In\nGlut",3)),3))
			   
## Plot

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/DMR/figures/DMR_Kozlenkov_Dracheva_NucAcidsRes_2016_GABA-GLU_dmC_Overlap.pdf",width=10,height=12)
ggplot(df[df$dir=="Both",], aes(x = model, y = percent)) + geom_bar(stat = "identity") +
  geom_text(aes( label = count), vjust = -.5) +
  facet_grid(. ~ list) + theme_classic() +
  labs(fill="") +
  ylab("Percent") + ylim(0,100) +
  xlab("") +
  ggtitle("DMRs Overlapping Neuronal Subtype DMPs") +
  theme(axis.text.x=element_text(angle=45,hjust=1)) +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
ggplot(df[df$dir!="Both",], aes(x = dir, y = percent)) + geom_bar(stat = "identity") +
  geom_text(aes( label = count), vjust = -.5) +
  facet_grid(list ~ model, scales = "free") + theme_classic() +
  labs(fill="") +
  ylab("Percent") + ylim(0,100) +
  xlab("") +
  ggtitle("DMRs Overlapping Neuronal Subtype DMPs") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
dev.off()