library(GenomicRanges)
library(bumphunter)
library(ggplot2)

load('/dcl01/lieber/ajaffe/lab/brain-epigenomics/bumphunting/BSobj_bsseqSmooth_Neuron_minCov_3.Rdata')
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/DMR/DMR_objects.rda")
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/bumphunting/rda/limma_Neuron_CpGs_minCov_3_ageInfo_dmrs.Rdata")


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

gabglu = list("450K" = makeGRangesFromDataFrame(gabglu.450K[which(gabglu.450K$FDR.GABA.GLU<=0.05),], 
                                                seqnames.field="CHR",start.field="POSITION", end.field="POSITION", keep.extra.columns=T),
              "rrbs" = makeGRangesFromDataFrame(gabglu.rrbs[which(gabglu.rrbs$qvalue<=0.05),], 
                                                seqnames.field="CHR",start.field="POSITION", end.field="POSITION", keep.extra.columns=T),
              "CHrrbs" = makeGRangesFromDataFrame(gabglu.CHrrbs[which(gabglu.CHrrbs$qvalue<=0.05),], 
                                                  seqnames.field="CHR",start.field="POSITION", end.field="POSITION", strand.field="STRAND", keep.extra.columns=T),
              "450K.upGABA" = makeGRangesFromDataFrame(gabglu.450K[which(gabglu.450K$FDR.GABA.GLU<=0.05 & gabglu.450K$delta.GABA.GLU < 0),], 
                                                       seqnames.field="CHR",start.field="POSITION", end.field="POSITION", keep.extra.columns=T),
              "rrbs.upGABA" = makeGRangesFromDataFrame(gabglu.rrbs[which(gabglu.rrbs$qvalue<=0.05 & gabglu.rrbs$DiffMeth=="DM.GABA"),], 
                                                       seqnames.field="CHR",start.field="POSITION", end.field="POSITION", keep.extra.columns=T),
              "CHrrbs.upGABA" = makeGRangesFromDataFrame(gabglu.CHrrbs[which(gabglu.CHrrbs$qvalue<=0.05 & gabglu.CHrrbs$DiffMeth=="GABA.DM"),],
                                                         seqnames.field="CHR",start.field="POSITION", end.field="POSITION", strand.field="STRAND", keep.extra.columns=T),
              "450K.upGlut" = makeGRangesFromDataFrame(gabglu.450K[which(gabglu.450K$FDR.GABA.GLU<=0.05 & gabglu.450K$delta.GABA.GLU > 0),], 
                                                       seqnames.field="CHR",start.field="POSITION", end.field="POSITION", keep.extra.columns=T),
              "rrbs.upGlut" = makeGRangesFromDataFrame(gabglu.rrbs[which(gabglu.rrbs$qvalue<=0.05 & gabglu.rrbs$DiffMeth=="DM.GLU"),], 
                                                       seqnames.field="CHR",start.field="POSITION", end.field="POSITION", keep.extra.columns=T),
              "CHrrbs.upGlut" = makeGRangesFromDataFrame(gabglu.CHrrbs[which(gabglu.CHrrbs$qvalue<=0.05 & gabglu.CHrrbs$DiffMeth=="GLU.DM"),], 
                                                         seqnames.field="CHR",start.field="POSITION", end.field="POSITION", strand.field="STRAND", keep.extra.columns=T))


# Identify all CpG clusters in the genome
gr <- granges(BSobj)
cl = clusterMaker( chr = as.character(seqnames(gr)), 
                   pos = start(gr),   maxGap = 1000)
gr.clusters = split(gr, cl)
gr.clusters = unlist(range(gr.clusters))
df.clusters = as.data.frame(gr.clusters)

# Find overlaps with DMRS in all three models

dmrs = split(dmrs, dmrs$k6cluster_label)
names(dmrs) = c("Group 1 (G-N+)","Group 2 (G0N+)","Group 3 (G0N-)","Group 4 (G+N0)","Group 5 (G+N-)","Group 6 (G-N0)")
oo = lapply(dmrs, function(x) findOverlaps(x, makeGRangesFromDataFrame(DMR$Interaction)))
dmrs = lapply(oo, function(x) DMR$Interaction[subjectHits(x),])
CT = split(DMR$CellType, DMR$CellType$Dir)
names(CT) = c("Hypomethylated in Neurons", "Hypomethylated in Glia")
DMRgr = lapply(c(CT, DMR[names(DMR)!="CellType"], dmrs), function(x) makeGRangesFromDataFrame(x[which(x$fwer<=0.05),], keep.extra.columns = T))

ooDMR = lapply(DMRgr, function(x) findOverlaps(gr.clusters, x))
ooGG = lapply(gabglu, function(x) findOverlaps(gr.clusters,x))

df.clusters$regionID = paste0(df.clusters$seqnames,":",df.clusters$start,"-",df.clusters$end)
df.clusters$rnum = 1:length(gr.clusters)
df.clusters$Age = ifelse(df.clusters$rnum %in% queryHits(ooDMR$Age), "yes","no")
df.clusters$Interaction = ifelse(df.clusters$rnum %in% queryHits(ooDMR$Interaction), "yes","no")
df.clusters$"Hypomethylated in Neurons" = ifelse(df.clusters$rnum %in% queryHits(ooDMR$"Hypomethylated in Neurons"), "yes","no")
df.clusters$"Hypomethylated in Glia" = ifelse(df.clusters$rnum %in% queryHits(ooDMR$"Hypomethylated in Glia"), "yes","no")
df.clusters$"Group 1 (G-N+)" = ifelse(df.clusters$rnum %in% queryHits(ooDMR$"Group 1 (G-N+)"), "yes","no")
df.clusters$"Group 2 (G0N+)" = ifelse(df.clusters$rnum %in% queryHits(ooDMR$"Group 2 (G0N+)"), "yes","no")
df.clusters$"Group 3 (G0N-)" = ifelse(df.clusters$rnum %in% queryHits(ooDMR$"Group 3 (G0N-)"), "yes","no")
df.clusters$"Group 4 (G+N0)" = ifelse(df.clusters$rnum %in% queryHits(ooDMR$"Group 4 (G+N0)"), "yes","no")
df.clusters$"Group 5 (G+N-)" = ifelse(df.clusters$rnum %in% queryHits(ooDMR$"Group 5 (G+N-)"), "yes","no")
df.clusters$"Group 6 (G-N0)" = ifelse(df.clusters$rnum %in% queryHits(ooDMR$"Group 6 (G-N0)"), "yes","no")

df.clusters$"450K" = ifelse(df.clusters$rnum %in% queryHits(ooGG$"450K"), "yes","no")
df.clusters$rrbs = ifelse(df.clusters$rnum %in% queryHits(ooGG$"rrbs"), "yes","no")
df.clusters$CHrrbs = ifelse(df.clusters$rnum %in% queryHits(ooGG$"CHrrbs"), "yes","no")
df.clusters$"450K.upGABA" = ifelse(df.clusters$rnum %in% queryHits(ooGG$"450K.upGABA"), "yes","no")
df.clusters$rrbs.upGABA = ifelse(df.clusters$rnum %in% queryHits(ooGG$"rrbs.upGABA"), "yes","no")
df.clusters$CHrrbs.upGABA = ifelse(df.clusters$rnum %in% queryHits(ooGG$"CHrrbs.upGABA"), "yes","no")
df.clusters$"450K.upGlut" = ifelse(df.clusters$rnum %in% queryHits(ooGG$"450K.upGlut"), "yes","no")
df.clusters$rrbs.upGlut = ifelse(df.clusters$rnum %in% queryHits(ooGG$"rrbs.upGlut"), "yes","no")
df.clusters$CHrrbs.upGlut = ifelse(df.clusters$rnum %in% queryHits(ooGG$"CHrrbs.upGlut"), "yes","no")


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


load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/bumphunting/rda/limma_Neuron_CpGs_minCov_3_ageInfo_dmrs.Rdata")

dmrs = split(dmrs, dmrs$k6cluster_label)
oo = lapply(dmrs, function(x) findOverlaps(x, makeGRangesFromDataFrame(dtinteraction[sig=="FWER < 0.05",,])))
dtinteraction = dtinteraction[sig=="FWER < 0.05",,]
intclusters = lapply(oo, function(x) dtinteraction[subjectHits(x),,])

gabglu = gabglu[which(names(gabglu) %in% c("rrbs.upGABA","CHrrbs.upGABA","rrbs.upGlut","CHrrbs.upGlut"))]
ooInt = lapply(intclusters, function(i) lapply(gabglu, function(x) findOverlaps(makeGRangesFromDataFrame(i),x)))

df2 = mapply(function(oo, i) lapply(oo, function(x) rbind(cbind(i[unique(queryHits(x)),,], Hit = "Yes"), cbind(i[-unique(queryHits(x)),,], Hit = "No"))),
            ooInt, intclusters, SIMPLIFY = F)
df2 = lapply(df2, function(x) lapply(x, function(y) y[,length(unique(regionID)), by = "Hit"]))
df2 = lapply(df2, function(x) lapply(x, function(y) cbind(y, Perc = round(y$V1/sum(y$V1)*100,1))))
df2 = do.call(rbind, Map(cbind, lapply(df2, function(d) do.call(rbind, Map(cbind, d, CT = as.list(names(d))))), Cluster = as.list(names(df2))))
df2$dir[grep("GABA", df2$CT)] = "More In\nGABA"
df2$dir[grep("Glut", df2$CT)] = "More In\nGlut"
df2$list[-grep("CHrrbs", df2$CT)] = "CpG"
df2$list[grep("CHrrbs", df2$CT)] = "CpH"
df2$Cluster = factor(df2$Cluster, levels = c("1:G-N+", "2:G0N+", "3:G0N-", "4:G+N0", "5:G+N-", "6:G-N0"))


## Plot

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/DMR/figures/DMR_Kozlenkov_Dracheva_NucAcidsRes_2016_GABA-GLU_dmC_Overlap.pdf",width=12,height=6)
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
ggplot(df[df$dir!="Both" & df$list!="450K",], aes(x = dir, y = percent)) + geom_bar(stat = "identity") +
  geom_text(aes( label = count), vjust = -.5) +
  facet_grid(list ~ model, scales = "free") + theme_classic() +
  labs(fill="") +
  ylab("Percent") + ylim(0,100) +
  xlab("") +
  ggtitle("DMRs Overlapping Neuronal Subtype DMPs") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
ggplot(df2[df2$Hit=="Yes",], aes(x = dir, y = Perc, fill=Cluster)) + geom_bar(stat = "identity") +
  geom_text(aes( label = V1), vjust = -.5) + 
  scale_fill_brewer(8, palette="Dark2") +
  facet_grid(list ~ Cluster, scales = "free") + theme_classic() +
  labs(fill="") +
  ylab("Percent") + ylim(0,100) +
  xlab("") +
  ggtitle("DMRs Overlapping Neuronal Subtype DMPs") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20)) +
  theme(legend.position="none")
dev.off()

## find enrichment for interaction clusters

ooClus = lapply(dmrs, function(x) findOverlaps(gr.clusters, x))
df.clusters$Gr1 = ifelse(df.clusters$rnum %in% queryHits(ooClus[[1]]), "yes","no")
df.clusters$Gr2 = ifelse(df.clusters$rnum %in% queryHits(ooClus[[2]]), "yes","no")
df.clusters$Gr3 = ifelse(df.clusters$rnum %in% queryHits(ooClus[[3]]), "yes","no")
df.clusters$Gr4 = ifelse(df.clusters$rnum %in% queryHits(ooClus[[4]]), "yes","no")
df.clusters$Gr5 = ifelse(df.clusters$rnum %in% queryHits(ooClus[[5]]), "yes","no")
df.clusters$Gr6 = ifelse(df.clusters$rnum %in% queryHits(ooClus[[6]]), "yes","no")

tables = lapply(as.list(c("Gr1","Gr2","Gr3","Gr4","Gr5","Gr6")), function(x)
  lapply(as.list(c("rrbs.upGABA","CHrrbs.upGABA","rrbs.upGlut","CHrrbs.upGlut")), function(y)
    data.frame(Hit = c(nrow(df.clusters[which(df.clusters[,colnames(df.clusters)==x]=="yes" & df.clusters[,colnames(df.clusters)==y]!="no"),]),
                 nrow(df.clusters[which(df.clusters[,colnames(df.clusters)==x]=="no" & df.clusters[,colnames(df.clusters)==y]!="no"),])),
               noHit = c(nrow(df.clusters[which(df.clusters[,colnames(df.clusters)==x]=="yes" & df.clusters[,colnames(df.clusters)==y]=="no"),]),
                 nrow(df.clusters[which(df.clusters[,colnames(df.clusters)==x]=="no" & df.clusters[,colnames(df.clusters)==y]=="no"),])))))

fisher = lapply(tables, function(x) lapply(x, fisher.test))
or = lapply(fisher, function(x) lapply(x, function(y) data.frame(OR = y$estimate, P = y$p.value)))
or = do.call(rbind, Map(cbind, lapply(or, function(x) 
  do.call(rbind, Map(cbind, x, Group = as.list(c("rrbs.upGABA","CHrrbs.upGABA","rrbs.upGlut","CHrrbs.upGlut"))))),
  cluster = as.list(c("Gr1","Gr2","Gr3","Gr4","Gr5","Gr6"))))

write.csv(or,file="/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/DMR/interaction_kmeans_clusters_Dracheva_2016_GABA-GLU_dmC_Overlap.csv", quote=F)



  
  
  
  
  
  
