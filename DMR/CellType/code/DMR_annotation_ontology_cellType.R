library(GenomicRanges)
library(GenomicFeatures)
library(data.table)
library(ggplot2)
library(clusterProfiler)
require(org.Hs.eg.db)
library(RColorBrewer)
library(bumphunter)


load("/dcl01/lieber/ajaffe/CellSorting/RNAseq_pipeline/rawCounts_CellSorting_July5_n12.rda")
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

dtcell = data.table(DMR$CellType)
oo = findOverlaps(gr.clusters, makeGRangesFromDataFrame(dtcell[sig=="FWER < 0.05",,]))
length(unique(queryHits(oo)))
df.clusters$CellType = ifelse(df.clusters$rnum %in% queryHits(oo), "CellType","no")

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
lapply(features, head)

df.clusters$Islands = ifelse(df.clusters$rnum %in% subjectHits(findOverlaps(features$islands, gr.clusters)), "island","no")
df.clusters$repeats = ifelse(df.clusters$rnum %in% subjectHits(findOverlaps(features$rpmskgr, gr.clusters)), "repeat","no")
df.clusters$promoters = ifelse(df.clusters$rnum %in% subjectHits(findOverlaps(features$promoters, gr.clusters)), "promoter","no")
df.clusters$genes = ifelse(df.clusters$rnum %in% subjectHits(findOverlaps(makeGRangesFromDataFrame(geneMap), gr.clusters)), "gene","no")


## how many fall within CpG islands?

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/DMR/CellType/figures/DMR_overap_with_CpG-Islands_byCellType.pdf")
x = dtcell[sig=="FWER < 0.05",length(unique(regionID)), by = "islands"]
x$perc = round(x$V1/sum(x$V1)*100,2)
ggplot(x, aes(x = islands, y = V1)) + geom_bar(stat = "identity") +
  geom_text(aes( label = paste0(perc,"%")), vjust = -.5) +
  labs(fill="") + theme_classic() +
  ylab("Count") + 
  xlab("") +
  ggtitle("DMRs Overlapping CpG Islands: FWER < 0.05") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))

x = dtcell[,length(unique(regionID)), by = c("islands", "sig")]
x$perc = unlist(c(round(x[1:2,"V1"]/sum(x[1:2,"V1"])*100,2),
           round(x[3:4,"V1"]/sum(x[3:4,"V1"])*100,2)))
ggplot(x, aes(x = islands, y = V1)) + geom_bar(stat = "identity") +
  geom_text(aes( label = paste0(perc,"%")), vjust = -.5) +
  facet_grid(. ~ sig) + theme_classic() +
  labs(fill="") +
  ylab("Count") + 
  xlab("") +
  ggtitle("DMRs Overlapping CpG Islands") +
  theme(axis.text.x=element_text(angle=45,hjust=1)) +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
dev.off()

# Is there a relationship between being significantly DM and overlapping a CpG island?

fisher.test(data.frame(YesIsland = c(nrow(df.clusters[df.clusters$Islands=="island" & df.clusters$CellType=="CellType",]),
                                  nrow(df.clusters[df.clusters$Islands=="island" & df.clusters$CellType=="no",])),
                       NoIsland = c(nrow(df.clusters[df.clusters$Islands=="no" & df.clusters$CellType=="CellType",]),
                                  nrow(df.clusters[df.clusters$Islands=="no" & df.clusters$CellType=="no",])), row.names = c("YesCT","NoCT")))
# CpG islands are overrepresented in DMRs
#p-value < 2.2e-16
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  30.78423 33.99262
#sample estimates:
#  odds ratio 
#32.34714 

# How about Repetitive Elements?

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/DMR/CellType/figures/DMR_repetitiveElements_byCellType.pdf", width = 12)
x = dtcell[sig=="FWER < 0.05",length(unique(regionID)), by = "repeats"]
x$perc = round(x$V1/sum(x$V1)*100,2)
x$repeats = factor(x$repeats, levels = c("SINE","LINE","Simple_repeat","DNA","LTR","Low_complexity","No repeats","Unknown","snRNA",
                                         "Other","DNA?","Satellite","RC","scRNA","tRNA","SINE?","srpRNA","RNA","rRNA","LINE?","LTR?","Unknown?"))
ggplot(x, aes(x = repeats, y = V1)) + geom_bar(stat = "identity") +
  geom_text(aes( label = paste0(perc,"%")), vjust = -.5) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  labs(fill="") + theme_classic() +
  ylab("Count") + 
  xlab("") +
  ggtitle("DMRs in repetitive elements: FWER < 0.05") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
dev.off()


# Is there a relationship between being significantly DM and overlapping a repeat?

fisher.test(data.frame(Yesrepeats = c(nrow(df.clusters[df.clusters$repeats=="repeat" & df.clusters$CellType=="CellType",]),
                                     nrow(df.clusters[df.clusters$repeats=="repeat" & df.clusters$CellType=="no",])),
                       Norepeats = c(nrow(df.clusters[df.clusters$repeats=="no" & df.clusters$CellType=="CellType",]),
                                    nrow(df.clusters[df.clusters$repeats=="no" & df.clusters$CellType=="no",])), row.names = c("YesCT","NoCT")))
#p-value < 2.2e-16
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  1.903427 2.311114
#sample estimates:
#  odds ratio 
#2.095256 


# assign genomic features

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/DMR/CellType/figures/DMR_annotation_byCellType.pdf")
x = dtcell[sig=="FWER < 0.05",length(unique(regionID)), by = "annotation"]
x$perc = round(x$V1/sum(x$V1)*100,2)
ggplot(x, aes(x = annotation, y = V1)) + geom_bar(stat = "identity") +
  geom_text(aes( label = paste0(perc,"%")), vjust = -.5) +
  labs(fill="") + theme_classic() +
  ylab("Count") + 
  xlab("") +
  ggtitle("DMR annotation: FWER < 0.05") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))

x = dtcell[,length(unique(regionID)), by = c("annotation", "sig")]
x$perc = unlist(c(round(x[1:6,"V1"]/sum(x[1:6,"V1"])*100,2),
                  round(x[7:12,"V1"]/sum(x[7:12,"V1"])*100,2)))
ggplot(x, aes(x = annotation, y = V1)) + geom_bar(stat = "identity") +
  geom_text(aes( label = paste0(perc,"%")), vjust = -.5) +
  facet_grid(. ~ sig) + theme_classic() +
  labs(fill="") +
  ylab("Count") + 
  theme(axis.text.x=element_text(angle=45,hjust=1)) +
  xlab("") +
  ggtitle("DMR annotation") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
dev.off()

# Is there a relationship between being significantly DM and overlapping a gene?

fisher.test(data.frame(Yesgenes = c(nrow(df.clusters[df.clusters$genes=="gene" & df.clusters$CellType=="CellType",]),
                                      nrow(df.clusters[df.clusters$genes=="gene" & df.clusters$CellType=="no",])),
                       Nogenes = c(nrow(df.clusters[df.clusters$genes=="no" & df.clusters$CellType=="CellType",]),
                                     nrow(df.clusters[df.clusters$genes=="no" & df.clusters$CellType=="no",])), row.names = c("YesCT","NoCT")))
# Genes are overrepresented in DMRs
#p-value < 2.2e-16
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  4.900227 5.590851
#sample estimates:
#  odds ratio 
#5.23127 

# Is there a relationship between being significantly DM and overlapping a promoter?

fisher.test(data.frame(Yespromoters = c(nrow(df.clusters[df.clusters$promoters=="promoter" & df.clusters$CellType=="CellType",]),
                                    nrow(df.clusters[df.clusters$promoters=="promoter" & df.clusters$CellType=="no",])),
                       Nopromoters = c(nrow(df.clusters[df.clusters$promoters=="no" & df.clusters$CellType=="CellType",]),
                                   nrow(df.clusters[df.clusters$promoters=="no" & df.clusters$CellType=="no",])), row.names = c("YesCT","NoCT")))
# promoters by themselves aren't overrepresented in DMRs
#p-value < 2.2e-16
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  12.42343 13.64845
#sample estimates:
#  odds ratio 
#13.02401 


### Gene Ontology
entrez = list(All = dtcell[sig=="FWER < 0.05",list(na.omit(EntrezID)),], 
              GenesPlusPromoters = dtcell[annotation != "Intergenic" & sig=="FWER < 0.05",list(na.omit(EntrezID)),],
              Genes = dtcell[distToGene==0 & sig=="FWER < 0.05",list(na.omit(EntrezID)),],
              Promoters = dtcell[annotation == "Promoter" & sig=="FWER < 0.05",list(na.omit(EntrezID)),])
entrez.dir = list(All.pos = dtcell[value>0 & sig=="FWER < 0.05",list(na.omit(EntrezID)),], 
                  GenesPlusPromoters.pos = dtcell[value>0 & annotation != "Intergenic" & sig=="FWER < 0.05",list(na.omit(EntrezID)),],
                  Genes.pos = dtcell[value>0 & distToGene==0 & sig=="FWER < 0.05",list(na.omit(EntrezID)),],
                  Promoters.pos = dtcell[value>0 & annotation == "Promoter" & sig=="FWER < 0.05",list(na.omit(EntrezID)),],
                  All.neg = dtcell[value<0 & sig=="FWER < 0.05",list(na.omit(EntrezID)),], 
                  GenesPlusPromoters.neg = dtcell[value<0 & annotation != "Intergenic" & sig=="FWER < 0.05",list(na.omit(EntrezID)),],
                  Genes.neg = dtcell[value<0 & distToGene==0 & sig=="FWER < 0.05",list(na.omit(EntrezID)),],
                  Promoters.neg = dtcell[value<0 & annotation == "Promoter" & sig=="FWER < 0.05",list(na.omit(EntrezID)),])
entrez = lapply(entrez, function(x) as.character(unique(x$V1)))              
entrez.dir = lapply(entrez.dir, function(x) as.character(unique(x$V1)))       

GeneUniverse = dtcell[,list(na.omit(EntrezID)),]
GeneUniverse = as.character(unique(GeneUniverse$V1))

# Find enriched Pathways via KEGG
elementNROWS(entrez)
keggList = lapply(entrez, function(x) enrichKEGG(as.character(x), organism="human", universe= GeneUniverse, 
                                                   minGSSize=5, pAdjustMethod="BH", qvalueCutoff=1))
keggListdf = lapply(keggList, function(x) as.data.frame(x))
keggList.dir = lapply(entrez.dir, function(x) enrichKEGG(as.character(x), organism="human", universe= GeneUniverse, 
                                                   minGSSize=5, pAdjustMethod="BH", qvalueCutoff=1))
keggList.dir.df = lapply(keggList.dir, function(x) as.data.frame(x))

# Enriched Molecular Function GOs
goList_MF = lapply(entrez, function(x) enrichGO(as.character(x), ont = "MF", OrgDb = org.Hs.eg.db, 
                                                  universe= GeneUniverse, minGSSize=5, pAdjustMethod="BH",
                                                  qvalueCutoff=1))
goListdf_MF = lapply(goList_MF, function(x) as.data.frame(x))
goList_MF.dir = lapply(entrez.dir, function(x) enrichGO(as.character(x), ont = "MF", OrgDb = org.Hs.eg.db, 
                                                universe= GeneUniverse, minGSSize=5, pAdjustMethod="BH",
                                                qvalueCutoff=1))
goListdf_MF.dir = lapply(goList_MF.dir, function(x) as.data.frame(x))

# Biological Process GO enrichment
goList_BP = lapply(entrez, function(x) enrichGO(as.character(x), ont = "BP", OrgDb = org.Hs.eg.db, 
                                                universe= GeneUniverse, minGSSize=5, pAdjustMethod="BH",
                                                qvalueCutoff=1))
goListdf_BP = lapply(goList_BP, function(x) as.data.frame(x))
goList_BP.dir = lapply(entrez.dir, function(x) enrichGO(as.character(x), ont = "BP", OrgDb = org.Hs.eg.db, 
                                                        universe= GeneUniverse, minGSSize=5, pAdjustMethod="BH",
                                                        qvalueCutoff=1))
goListdf_BP.dir = lapply(goList_BP.dir, function(x) as.data.frame(x))

# Cellular Compartment GO enrichment
goList_CC = lapply(entrez, function(x) enrichGO(as.character(x), ont = "CC", OrgDb = org.Hs.eg.db, 
                                                universe= GeneUniverse, minGSSize=5, pAdjustMethod="BH",
                                                qvalueCutoff=1))
goListdf_CC = lapply(goList_CC, function(x) as.data.frame(x))
goList_CC.dir = lapply(entrez.dir, function(x) enrichGO(as.character(x), ont = "CC", OrgDb = org.Hs.eg.db, 
                                                        universe= GeneUniverse, minGSSize=5, pAdjustMethod="BH",
                                                        qvalueCutoff=1))
goListdf_CC.dir = lapply(goList_CC.dir, function(x) as.data.frame(x))

# Disease Ontology
goList_DO = lapply(entrez, function(x) enrichDO(as.character(x), ont = "DO", universe= GeneUniverse, 
                                                  minGSSize=5, pAdjustMethod="BH", qvalueCutoff=1, readable=TRUE))
goListdf_DO = lapply(goList_DO, function(x) as.data.frame(x))
goList_DO.dir = lapply(entrez.dir, function(x) enrichDO(as.character(x), ont = "DO", universe= GeneUniverse, 
                                                minGSSize=5, pAdjustMethod="BH", qvalueCutoff=1, readable=TRUE))
goListdf_DO.dir = lapply(goList_DO.dir, function(x) as.data.frame(x))

# Compare the enriched terms between 7 groups
# KEGG
compareKegg = compareCluster(entrez, fun="enrichKEGG", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
compareKegg.dir = compareCluster(entrez.dir, fun="enrichKEGG", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
# Biological Process
compareBP = compareCluster(entrez, fun="enrichGO", ont = "BP", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
compareBP.dir = compareCluster(entrez.dir, fun="enrichGO",  ont = "BP", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
# Molecular Function
compareMF = compareCluster(entrez, fun="enrichGO",  ont = "MF", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
compareMF.dir = compareCluster(entrez.dir, fun="enrichGO",  ont = "MF", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
# Cellular Component
compareCC = compareCluster(entrez, fun="enrichGO",  ont = "CC", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
compareCC.dir = compareCluster(entrez.dir, fun="enrichGO",  ont = "CC", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
# Disease Ontology
compareDO = compareCluster(entrez, fun="enrichDO",  ont = "DO", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
compareDO.dir = compareCluster(entrez.dir, fun="enrichDO",  ont = "DO", qvalueCutoff = 0.05, pvalueCutoff = 0.05)

# save, write to csv
save(compareKegg, compareKegg.dir, compareBP, compareBP.dir, compareMF, compareMF.dir, compareCC, compareCC.dir, compareDO, compareDO.dir,
     keggList, keggList.dir, goList_BP, goList_BP.dir, goList_MF, goList_MF.dir, goList_CC, goList_CC.dir, goList_DO, goList_DO.dir,
     keggListdf, keggList.dir.df, goListdf_BP, goListdf_BP.dir, goListdf_MF, goListdf_MF.dir, goListdf_CC, goListdf_CC.dir, goListdf_DO, goListdf_DO.dir,
     file="/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/DMR/CellType/DMR_KEGG_GO_DO_objects_byCellType.rda")

# plot compared results
pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/DMR/CellType/figures/DMR_KEGG_GO_DO_plots_byCellType.pdf", height = 120, width = 20)
plot(compareKegg, colorBy="p.adjust", showCategory = 800, title= "KEGG Pathway Enrichment")
plot(compareBP, colorBy="p.adjust", showCategory = 900, title= "Biological Process GO Enrichment")
plot(compareMF, colorBy="p.adjust", showCategory = 800, title= "Molecular Function GO Enrichment")
plot(compareCC, colorBy="p.adjust", showCategory = 800, title= "Cellular Compartment GO Enrichment")
plot(compareDO, colorBy="p.adjust", showCategory = 800, title= "Disease Ontology Enrichment")
plot(compareKegg.dir, colorBy="p.adjust", showCategory = 800, title= "KEGG Pathway Enrichment")
plot(compareBP.dir, colorBy="p.adjust", showCategory = 800, title= "Biological Process GO Enrichment")
plot(compareMF.dir, colorBy="p.adjust", showCategory = 800, title= "Molecular Function GO Enrichment")
plot(compareCC.dir, colorBy="p.adjust", showCategory = 800, title= "Cellular Compartment GO Enrichment")
plot(compareDO.dir, colorBy="p.adjust", showCategory = 800, title= "Disease Ontology Enrichment")
dev.off()


# Make more specific enrichment figure

x = as.data.frame(compareBP.dir)
x = split(x, x$Cluster)

pos = x$GenesPlusPromoters.pos[-which(x$GenesPlusPromoters.pos$Description %in% x$GenesPlusPromoters.neg$Description),]
neg = x$GenesPlusPromoters.neg[-which(x$GenesPlusPromoters.neg$Description %in% x$GenesPlusPromoters.pos$Description),]

toplot = rbind(pos[1:5,], neg[1:5,])
toplot$log10 = -log10(toplot$p.adjust)
toplot = toplot[,colnames(toplot) !="geneID"]
toplot$Description = factor(c(toplot$Description[1:8], "presynaptic process involved\nin chemical synaptic transmission", toplot$Description[10]),
                              levels = c("gliogenesis","glial cell differentiation","ensheathment of neurons","axon ensheathment","myelination",
                                         "neurotransmitter secretion","signal release from synapse",
                                         "neurotransmitter transport","presynaptic process involved\nin chemical synaptic transmission",
                                         "vesicle localization"))
toplot$Cluster = gsub("GenesPlusPromoters.neg", "Hypomethylated\nIn Neurons", toplot$Cluster)
toplot$Cluster = gsub("GenesPlusPromoters.pos", "Hypomethylated\nIn Glia", toplot$Cluster)
toplot$Cluster = factor(toplot$Cluster, levels = c("Hypomethylated\nIn Glia", "Hypomethylated\nIn Neurons"))


pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/DMR/CellType/figures/DMR_byCellType_BP_top12.pdf", width=8,height=8)
ggplot(toplot, aes(x=Description, y =log10)) + scale_fill_brewer(8, palette="Dark2") + 
  geom_bar(stat = 'identity', aes(fill = Cluster), position = 'dodge', col = 'transparent') +
  geom_text(aes(label=GeneRatio), hjust=1.1, color="black", position = position_dodge(1), size=5) +
  coord_flip() + theme_classic() +
  ylab("-log10(FDR)") + 
  xlab("") +
  ggtitle("Biological Process\nEnrichment: Cell Type DMRs") + 
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20), legend.title=element_blank()) + theme(legend.position="bottom")
dev.off()


## Width of DMRs

cellgr = reduce(makeGRangesFromDataFrame(DMR$CellType[which(DMR$CellType$sig=="FWER < 0.05"),], keep.extra.columns = T))
mean(width(cellgr)) # 2783.261
median(width(cellgr)) # 2049
sd(width(cellgr)) # 2885.05
min(width(cellgr)) # 1
max(width(cellgr)) # 72913

geneMapgr = makeGRangesFromDataFrame(geneMap)
geneMapgr$geneID = names(geneMapgr)
mean(width(geneMapgr)) # 29389.18
median(width(geneMapgr)) # 2602.5
min(width(geneMapgr)) # 8

tiles = tile(geneMapgr[width(geneMapgr)>=20], n = 20)
names(tiles) = names(geneMapgr[width(geneMapgr)>=20])

regions = makeGRangesFromDataFrame(DMR$CellType[which(DMR$CellType$sig=="FWER < 0.05" & DMR$CellType$distToGene==0),], keep.extra.columns = T)

tiles = tiles[names(tiles) %in% regions$nearestID]
tiles = as.list(tiles)
pos = lapply(tiles, function(x) findOverlaps(x, reduce(regions[which(regions$Dir=="pos")])))
neg = lapply(tiles, function(x) findOverlaps(x, reduce(regions[which(regions$Dir=="neg")])))

overlaps1 = overlaps2 = list()
for (i in 1:length(pos)) {
  if (length(queryHits(pos[[i]]))>0) {
    overlaps1[[i]] = data.frame(region = queryHits(pos[[i]]), strand = as.character(strand(tiles[[i]])[1])) } else {
      overlaps1[[i]] = data.frame(region = 0, strand = as.character(strand(tiles[[i]])[1])) }
  if (length(queryHits(neg[[i]]))>0) {
    overlaps2[[i]] = data.frame(region = queryHits(neg[[i]]), strand = as.character(strand(tiles[[i]])[1])) } else {
      overlaps2[[i]] = data.frame(region = 0, strand = as.character(strand(tiles[[i]])[1])) }
}
overlaps = rbind(do.call(rbind, Map(cbind, overlaps1, geneID = as.list(names(tiles)), Direction = "Hypermethylated\nin Neurons")),
                 do.call(rbind, Map(cbind, overlaps2, geneID = as.list(names(tiles)), Direction = "Hypomethylated\nin Neurons")))

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/DMR/CellType/figures/CellType_DMR_position_in_gene.pdf")
ggplot(overlaps, aes(region, fill = Direction)) + scale_fill_brewer(8, palette="Dark2") +
  geom_density(alpha = 1/2) + theme_classic() +
  labs(fill="") +
  ylab("Density") + xlab("") +
  ggtitle("Cell Type DMR Position in Gene") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20), legend.position="bottom", legend.title=element_blank())
dev.off()




