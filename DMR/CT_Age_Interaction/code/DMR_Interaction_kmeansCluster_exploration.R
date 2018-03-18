library(GenomicRanges)
library(GenomicFeatures)
library(data.table)
library(ggplot2)
library(clusterProfiler)
require(org.Hs.eg.db)
library(RColorBrewer)
library(bumphunter)


load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/bumphunting/rda/limma_Neuron_CpGs_minCov_3_ageInfo_dmrs.Rdata")
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/DMR/DMR_objects.rda")
load('/dcl01/lieber/ajaffe/lab/brain-epigenomics/bumphunting/BSobj_bsseqSmooth_Neuron_minCov_3.Rdata')

# Identify all CpG clusters in the genome
gr = granges(BSobj)
cl = clusterMaker( chr = as.character(seqnames(gr)), 
                   pos = start(gr),   maxGap = 1000)
gr.clusters = split(gr, cl)
gr.clusters = unlist(range(gr.clusters))
df.clusters = as.data.frame(gr.clusters)
df.clusters$regionID = paste0(df.clusters$seqnames,":",df.clusters$start,"-",df.clusters$end)
df.clusters$rnum = 1:length(gr.clusters)

dtinteraction = data.table(DMR$Interaction)

txdb = loadDb("/dcl01/lieber/ajaffe/Amanda/annotation_objects/gencode.v25lift37.annotation.sqlite")
islands = read.table("/dcl01/lieber/ajaffe/Amanda/annotation_objects/cpgIslandExt.hg19.txt", sep="\t", header = T)
rpmsk = read.table("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/RepeatMasker_genomewide_HG19.txt", header=T)
features = list(CDS = cdsBy(txdb, by="tx", use.names=T), Introns = intronsByTranscript(txdb, use.names=T), 
                UTR5 = fiveUTRsByTranscript(txdb, use.names=T), UTR3 = threeUTRsByTranscript(txdb, use.names=T))
features = lapply(features, function(x) unlist(x, recursive = TRUE, use.names = TRUE))
for (i in 1:length(features)){
  tmp = features[[i]]
  tmp$TxID = names(tmp)
  features[[i]] = tmp }
features = c(features, 
             rpmskgr = makeGRangesFromDataFrame(rpmsk, seqnames.field="genoName",start.field="genoStart",end.field="genoEnd",keep.extra.columns=TRUE),
             islands = makeGRangesFromDataFrame(islands, keep.extra.columns = T, start.field = "chromStart", end.field = "chromEnd"),
             promoters = promoters(txdb, upstream=2000, downstream=200))

df.clusters$Islands = ifelse(df.clusters$rnum %in% subjectHits(findOverlaps(features$islands, gr.clusters)), "island","no")
df.clusters$repeats = ifelse(df.clusters$rnum %in% subjectHits(findOverlaps(features$rpmskgr, gr.clusters)), "repeat","no")
df.clusters$promoters = ifelse(df.clusters$rnum %in% subjectHits(findOverlaps(features$promoters, gr.clusters)), "promoter","no")
df.clusters$genes = ifelse(df.clusters$rnum %in% subjectHits(findOverlaps(makeGRangesFromDataFrame(geneMap), gr.clusters)), "gene","no")


## Separate Interaction DMRs into kmeans clusters

dmrs = split(dmrs, dmrs$k6cluster_label)
oo = lapply(dmrs, function(x) findOverlaps(x, makeGRangesFromDataFrame(dtinteraction[sig=="FWER < 0.05",,])))
dtinteraction = dtinteraction[sig=="FWER < 0.05",,]
intclusters = lapply(oo, function(x) dtinteraction[subjectHits(x),,])
oo = lapply(intclusters, function(x) findOverlaps(gr.clusters, makeGRangesFromDataFrame(x)))
df.clusters = Map(cbind, list(df.clusters,df.clusters,df.clusters,df.clusters,df.clusters,df.clusters), 
                  Interaction = lapply(oo, function(x) ifelse(df.clusters$rnum %in% queryHits(x), "Interaction","no")))
names(df.clusters) = names(oo)


## How many of the 2178 DMRs fall within each cluster?
elementNROWS(lapply(intclusters, function(x) unique(x$regionID)))
# 1:G-N+ 2:G0N+ 3:G0N- 4:G+N0 5:G+N- 6:G-N0 
#    567    423    362     96    170    560 

(567+423+560)/2179 # 0.7113355


# Which groups are my two interaction example regions?
lapply(intclusters, function(x) x[which(x$regionID %in% c("chr18:74712838-74729551","chr20:10219386-10230921")),unique(nearestSymbol),])
#`6:G-N0` = MBP
#`3:G0N-` = SNAP25

#how about one from group 1?
intclusters$'1:G-N+'[which(distToGene==0 & nearestSymbol %in% c("CUX1", "TCF4", "HDAC4", "CACNA1C", "MEGF6","NOTCH3")),]

intclusters$AorB = ifelse()

## how many fall within CpG islands?

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/DMR/CT_Age_Interaction/figures/DMR_overlap_CpG_Islands_kmeans_Interaction_clusters.pdf",width=16, height = 6)
x = lapply(intclusters, function(x) x[,length(unique(regionID)), by = "islands"])
x = do.call(rbind, Map(cbind, x, perc = lapply(x, function(y) round(y$V1/sum(y$V1)*100,1)), Cluster = as.list(names(intclusters))))
x$islands = ifelse(x$islands=="CpG-Island", "CGI", "non-\nCGI")
ggplot(x, aes(x = islands, y = perc)) + geom_bar(stat = "identity") +
  geom_text(aes( label = V1), vjust = -.5) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(fill="") + theme_classic() +
  facet_grid(. ~ Cluster) +
  ylab("Percent") + ylim(0,100) + 
  xlab("") +
  ggtitle("DMRs Overlapping CpG Islands: FWER < 0.05") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
dev.off()

# Is there a relationship between being significantly DM and overlapping a CpG island?

df = lapply(df.clusters, function(x) data.frame(YesIsland = c(nrow(x[x$Islands=="island" & x$Interaction=="Interaction",]),
                                                              nrow(x[x$Islands=="island" & x$Interaction=="no",])),
                                                NoIsland = c(nrow(x[x$Islands=="no" & x$Interaction=="Interaction",]),
                                                             nrow(x[x$Islands=="no" & x$Interaction=="no",])), 
                                                row.names = c("YesInteraction","NoInteraction")))
rbind(OR = unlist(lapply(lapply(df, fisher.test), function(y) y$estimate)), p = unlist(lapply(lapply(df, fisher.test), function(y) y$p.value)))
#1:G-N+.odds ratio 2:G0N+.odds ratio 3:G0N-.odds ratio 4:G+N0.odds ratio
#OR      1.813241e+02       3.25006e+02      3.467903e+02      2.331777e+02
#p      3.542515e-262      9.71136e-126     5.772372e-193      9.680598e-50
#5:G+N-.odds ratio 6:G-N0.odds ratio
#OR      2.589743e+02      1.516762e+02
#p       9.724690e-45     6.289749e-285

df = data.frame(YesIsland = c(x[islands=="CGI" & Cluster %in% c("1:G-N+","2:G0N+", "6:G-N0"), sum(V1),],
                              x[islands=="CGI" & Cluster %in% c("3:G0N-","4:G+N0","5:G+N-"),sum(V1)]),
                NoIsland = c(x[islands!="CGI" & Cluster %in% c("1:G-N+","2:G0N+", "6:G-N0"), sum(V1),],
                             x[islands!="CGI" & Cluster %in% c("3:G0N-","4:G+N0","5:G+N-"),sum(V1)]), 
                             row.names = c("GroupA","GroupB"))
fisher.test(df)
#p-value < 2.2e-16
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.2273814 0.3888052
#sample estimates:
#  odds ratio 
#0.2974871 


## Repetitive elements

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/DMR/CT_Age_Interaction/figures/DMR_overlap_repetitiveElements_kmeans_Interaction_clusters.pdf",
    width = 18, height = 10)
x = lapply(intclusters, function(x) x[,length(unique(regionID)), by = "repeats"])
x = do.call(rbind, Map(cbind, x, perc = lapply(x, function(y) round(y$V1/sum(y$V1)*100,1)), Cluster = as.list(names(intclusters))))
x$repeats = factor(x$repeats, levels = c("No repeats","SINE","LINE","Simple_repeat","DNA","LTR","Low_complexity","Unknown","snRNA",
                                         "Other","DNA?","Satellite","RC","scRNA","tRNA","SINE?","srpRNA","RNA","rRNA"))
ggplot(x, aes(x = repeats, y = perc)) + geom_bar(stat = "identity") +
  geom_text(aes( label = V1), vjust = -.5) +
  theme(axis.text.x=element_text(angle=45, hjust=1)) +
  labs(fill="") + theme_classic() +
  facet_grid(Cluster ~ ., scales = "free_x") +
  ylab("Percent") + ylim(0,100) + 
  xlab("") +
  ggtitle("DMRs in repetitive elements: FWER < 0.05") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
dev.off()


# Is there a relationship between being significantly DM and overlapping a repeat?

df = lapply(df.clusters, function(x) data.frame(Yesrepeats = c(nrow(x[x$repeats=="repeat" & x$Interaction=="Interaction",]),
                                                              nrow(x[x$repeats=="repeat" & x$Interaction=="no",])),
                                                Norepeats = c(nrow(x[x$repeats=="no" & x$Interaction=="Interaction",]),
                                                             nrow(x[x$repeats=="no" & x$Interaction=="no",])), 
                                                row.names = c("YesInteraction","NoInteraction")))
rbind(OR = unlist(lapply(lapply(df, fisher.test), function(y) y$estimate)), p = unlist(lapply(lapply(df, fisher.test), function(y) y$p.value)))
#              1:G-N+            2:G0N+            3:G0N-            4:G+N0
#OR               Inf               Inf               Inf               Inf
#p       1.108879e-12      1.695593e-05      2.066534e-08         0.0120985
#              5:G+N-            6:G-N0
#OR               Inf               Inf
#p         0.02982786      1.810973e-14

fisher.test(data.frame(Yes = c(x[repeats=="No repeats" & Cluster %in% c("1:G-N+","2:G0N+", "6:G-N0"), sum(V1),],
                               x[repeats=="No repeats" & Cluster %in% c("3:G0N-","4:G+N0","5:G+N-"),sum(V1)]),
                       No = c(x[repeats!="No repeats" & Cluster %in% c("1:G-N+","2:G0N+", "6:G-N0"), sum(V1),],
                              x[repeats!="No repeats" & Cluster %in% c("3:G0N-","4:G+N0","5:G+N-"),sum(V1)]), 
                       row.names = c("GroupA","GroupB")))
#p-value = 0.07977
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.6990714 1.0253477
#sample estimates:
#  odds ratio 
#0.8457447 



# assign genomic features

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/DMR/CT_Age_Interaction/figures/DMR_annotation_kmeans_Interaction_clusters.pdf", width = 10)
x = lapply(intclusters, function(x) x[,length(unique(regionID)), by = "annotation"])
x = do.call(rbind, Map(cbind, x, perc = lapply(x, function(y) round(y$V1/sum(y$V1)*100,1)), Cluster = as.list(names(intclusters))))
x$annotation = factor(x$annotation, levels = c("Intergenic", "Promoter", "5'UTR", "CDS", "Intron", "3'UTR"))
ggplot(x, aes(x = annotation, y = perc)) + geom_bar(stat = "identity") +
  geom_text(aes( label = V1), vjust = -.5) +
  facet_grid(Cluster ~ ., scales = "free") +
  labs(fill="") + theme_classic() +
  ylab("Percent") + ylim(0,100) +
  xlab("") +
  ggtitle("DMR annotation: FWER < 0.05") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
dev.off()

# Is there a relationship between being significantly DM and overlapping a gene?

df = lapply(df.clusters, function(x) data.frame(Yesrepeats = c(nrow(x[x$genes=="gene" & x$Interaction=="Interaction",]),
                                                               nrow(x[x$genes=="gene" & x$Interaction=="no",])),
                                                Norepeats = c(nrow(x[x$genes=="no" & x$Interaction=="Interaction",]),
                                                              nrow(x[x$genes=="no" & x$Interaction=="no",])), 
                                                row.names = c("YesInteraction","NoInteraction")))
rbind(OR = unlist(lapply(lapply(df, fisher.test), function(y) y$estimate)), p = unlist(lapply(lapply(df, fisher.test), function(y) y$p.value)))
#              1:G-N+            2:G0N+            3:G0N-            4:G+N0
#OR               Inf               Inf               Inf               Inf
#p       1.081069e-60      5.139346e-26      2.658813e-39      3.782012e-11
#              5:G+N-            6:G-N0
#OR               Inf               Inf
#p        8.46862e-10      2.523357e-69


# Is there a relationship between being significantly DM and overlapping a promoter?

df = lapply(df.clusters, function(x) data.frame(Yesrepeats = c(nrow(x[x$promoters=="promoter" & x$Interaction=="Interaction",]),
                                                               nrow(x[x$promoters=="promoter" & x$Interaction=="no",])),
                                                Norepeats = c(nrow(x[x$promoters=="no" & x$Interaction=="Interaction",]),
                                                              nrow(x[x$promoters=="no" & x$Interaction=="no",])), 
                                                row.names = c("YesInteraction","NoInteraction")))
rbind(OR = unlist(lapply(lapply(df, fisher.test), function(y) y$estimate)), p = unlist(lapply(lapply(df, fisher.test), function(y) y$p.value)))
#              1:G-N+            2:G0N+            3:G0N-            4:G+N0
#OR      1.064465e+02      3.502212e+02      3.565863e+02      2.881842e+02
#p      7.366667e-180      5.629673e-86     1.076833e-130      8.651033e-36
#              5:G+N-            6:G-N0
#OR      2.518328e+02      8.575791e+01
#p       2.877309e-31     3.757413e-198

f = lapply(as.list(unique(x$annotation)), function(a) fisher.test(data.frame(Yes = c(x[annotation==as.character(a) & Cluster %in% c("1:G-N+","2:G0N+", "6:G-N0"), sum(V1),],
                                                                                     x[annotation==as.character(a) & Cluster %in% c("3:G0N-","4:G+N0","5:G+N-"),sum(V1)]),
                                                                             No = c(x[annotation!=as.character(a) & Cluster %in% c("1:G-N+","2:G0N+", "6:G-N0"), sum(V1),],
                                                                                    x[annotation!=as.character(a) & Cluster %in% c("3:G0N-","4:G+N0","5:G+N-"),sum(V1)]), 
                                                                             row.names = c("GroupA","GroupB"))))
names(f) = unique(x$annotation)
rbind(OR = unlist(lapply(f, function(y) y$estimate)), p = p.adjust(unlist(lapply(f, function(y) y$p.value)), method = "fdr"))
#   CDS.odds ratio Intron.odds ratio Promoter.odds ratio Intergenic.odds ratio
#OR   2.586352e-01        1.19864741        2.4905823308          4.943342e+00
#p    1.437108e-35        0.08747343        0.0001086105          1.573111e-23
#   5'UTR.odds ratio 3'UTR.odds ratio
#OR        1.3720357        0.6853964
#p         0.2971496        0.2971496

#compared to groups 3,4 and 5, groups 1,2 and 6 are: depleted in CDS, enriched in promoter and intergenic sequence


### Gene Ontology
entrez = lapply(intclusters, function(x) list(All = x[,list(na.omit(EntrezID)),], GenesPlusPromoters = x[annotation != "Intergenic",list(na.omit(EntrezID)),],
                                              Genes = x[distToGene==0,list(na.omit(EntrezID)),], Promoters = x[annotation == "Promoter",list(na.omit(EntrezID)),]))
entrez = lapply(entrez, function(x) lapply(x, function(y) as.character(unique(y$V1))))              

GeneUniverse = dtinteraction[,list(na.omit(EntrezID)),]
GeneUniverse = as.character(unique(GeneUniverse$V1))

# Find enriched Pathways via KEGG

keggList = lapply(entrez, function(y) lapply(y, function(x) enrichKEGG(x, organism="human", universe= GeneUniverse,
                                                                       minGSSize=5, pAdjustMethod="BH", qvalueCutoff=1)))
goList_MF = lapply(entrez, function(y) lapply(y, function(x) enrichGO(x, ont = "MF", OrgDb = org.Hs.eg.db, universe= GeneUniverse, 
                                                                      minGSSize=5, pAdjustMethod="BH", qvalueCutoff=1)))
goList_BP = lapply(entrez, function(y) lapply(y, function(x) enrichGO(x, ont = "BP", OrgDb = org.Hs.eg.db, universe= GeneUniverse, 
                                                                      minGSSize=5, pAdjustMethod="BH", qvalueCutoff=1)))
goList_CC = lapply(entrez, function(y) lapply(y, function(x) enrichGO(x, ont = "CC", OrgDb = org.Hs.eg.db, universe= GeneUniverse, 
                                                                      minGSSize=5, pAdjustMethod="BH", qvalueCutoff=1)))
goList_DO = lapply(entrez, function(y) lapply(y, function(x) enrichDO(x, ont = "DO", universe= GeneUniverse, minGSSize=5, pAdjustMethod="BH", 
                                                                      qvalueCutoff=1, readable=TRUE)))

# Compare the enriched terms between groups

entrez2 = mapply(function(a,b,c,d,e,f) list("1:G-N+" = a, "2:G0N+" = b, "3:G0N-" = c, "4:G+N0" = d, "5:G+N-" = e, "6:G-N0" = f), 
                 entrez[["1:G-N+"]], entrez[["2:G0N+"]], entrez[["3:G0N-"]], entrez[["4:G+N0"]], entrez[["5:G+N-"]], entrez[["6:G-N0"]], SIMPLIFY = F)


compareKegg = lapply(entrez2, function(x) compareCluster(x, fun="enrichKEGG", qvalueCutoff = 0.05, pvalueCutoff = 0.05))
compareBP = lapply(entrez2, function(x) compareCluster(x, fun="enrichGO", ont = "BP", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05))
compareMF = lapply(entrez2, function(x) compareCluster(x, fun="enrichGO",  ont = "MF", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05))
compareCC = lapply(entrez2, function(x) compareCluster(x, fun="enrichGO",  ont = "CC", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05))
compareDO = lapply(entrez2, function(x) compareCluster(x, fun="enrichDO",  ont = "DO", qvalueCutoff = 0.05, pvalueCutoff = 0.05))

compareBP.nolimit = lapply(entrez2, function(x) compareCluster(x, fun="enrichGO", ont = "BP", OrgDb = org.Hs.eg.db, pvalueCutoff = 1))


# save
save(compareKegg, compareBP, compareMF, compareCC, keggList, goList_BP, goList_MF, goList_CC, goList_DO,
     file="/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/DMR/CT_Age_Interaction/DMR_KEGG_GO_DO_objects_Interaction_kmeans_clusters.rda")

# plot compared results
pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/DMR/CT_Age_Interaction/figures/DMR_KEGG_GO_DO_plots_kmeans_Interaction_clusters.pdf", height = 20, width = 10)
mapply(function(comp, title) plot(comp, colorBy="p.adjust", showCategory = 450, title= paste0("KEGG Pathway: ", title)), compareKegg, as.list(names(compareKegg)), SIMPLIFY = F)
mapply(function(comp, title) plot(comp, colorBy="p.adjust", showCategory = 450, title= paste0("Biological Process GO: ", title)), compareBP, as.list(names(compareBP)), SIMPLIFY = F)
mapply(function(comp, title) plot(comp, colorBy="p.adjust", showCategory = 450, title= paste0("Molecular Function GO: ", title)), compareMF, as.list(names(compareMF)), SIMPLIFY = F)
mapply(function(comp, title) plot(comp, colorBy="p.adjust", showCategory = 450, title= paste0("Cellular Compartment GO: ", title)), compareCC, as.list(names(compareCC)), SIMPLIFY = F)
dev.off()

## Explore results

compareBP = lapply(compareBP, as.data.frame)
compareBP = lapply(compareBP, function(x) x[order(x$p.adjust),])
geneBP = split(compareBP$Genes, compareBP$Genes$Cluster)
lapply(geneBP,head)

selected = c("GO:0022604", "GO:0007272", "GO:0022010","GO:0007265", "GO:1901184",
             "GO:0048813","GO:0031346","GO:0007411","GO:0051493","GO:0098742",
             "GO:0050808","GO:0010975","GO:0098742","GO:0061351")
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/DMR/CT_Age_Interaction/DMR_KEGG_GO_DO_objects_Interaction_kmeans_clusters.rda")

## Plot select terms

plotExample = compareBP$Genes # clusterProfiler output
plotExample@compareClusterResult = plotExample@compareClusterResult[which(plotExample@compareClusterResult$ID %in% selected),]
plotExample@compareClusterResult$Description = c("Ras protein signal transduction","regulation of ERBB signaling pathway",
                                                 "dendrite morphogenesis","positive regulation of\ncell projection organization",
                                                 "regulation of cytoskeleton organization","axon guidance","regulation of neuron projection development",
                                                 "cell-cell adhesion via\nplasma-membrane adhesion molecules","regulation of neuron projection development",
                                                 "synapse organization","cell-cell adhesion via\nplasma-membrane adhesion molecules","neural precursor cell proliferation",
                                                 "regulation of cell morphogenesis","ensheathment of neurons","central nervous system myelination",
                                                 "positive regulation of\ncell projection organization","Ras protein signal transduction","dendrite morphogenesis",
                                                 "regulation of neuron projection development")

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/DMR/CT_Age_Interaction/figures/BP_selectedTerms_Interaction_kmeansClusters.pdf", height = 6, width=8.5)
plot(plotExample, colorBy="p.adjust", showCategory = 450, title= "Biological Process Enrichment")
dev.off()


## Width of DMRs

data.frame(mean = unlist(lapply(dmrs, function(x) mean(width(x)))), median = unlist(lapply(dmrs, function(x) median(width(x)))),
           sd = unlist(lapply(dmrs, function(x) sd(width(x)))), min = unlist(lapply(dmrs, function(x) min(width(x)))),
           max = unlist(lapply(dmrs, function(x) max(width(x)))))
#           mean median        sd min   max
#1:G-N+ 1428.115 1097.0 1214.0185   1 11120
#2:G0N+ 1013.125  846.0  766.3112   1  5321
#3:G0N- 1676.950 1248.5 1456.1044   1 12075
#4:G+N0 1163.823  883.0 1001.6641   1  5697
#5:G+N- 1277.535  925.0 1060.0442   1  5545
#6:G-N0 1467.284 1048.0 1591.8450   1 16714


