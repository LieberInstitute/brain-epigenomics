library(GenomicFeatures)
library(ggplot2)
library(GenomicRanges)
library(bumphunter)
library(data.table)
library(clusterProfiler)
require(org.Hs.eg.db)

load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/DMVs_hypoDMRs_methylSeekR.rda")
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/DMVs_hypoDMRs_methylSeekR_prenatal.rda")
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/DMR/DMR_objects.rda")

rm(DMVs,DMVs_pren)

dmvs = c(DMVs.100kb,DMVs_100kb.pren)

## How many are identified?

num = data.frame(num = unlist(elementNROWS(dmvs)), id = names(dmvs))
num[which(num$id %in% names(DMVs_100kb.pren)),"celltype"] = "Prenatal"
num[which(num$id %in% pd[pd$Cell.Type=="Neuron","Data.ID"]),"celltype"] = "Neuron"
num[which(num$id %in% pd[pd$Cell.Type=="Glia","Data.ID"]),"celltype"] = "Glia"
write.csv(num, quote=F, file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/DMV_number.csv")

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/CREs/figures/DMV_number_byAge_postnatal.pdf")

w = num[which(num$celltype %in% c("Neuron", "Glia")),]
w = cbind(w, Age = pd[match(w$id, pd$Data.ID), "Age"])
cor.test(x = w[which(w$celltype=="Neuron"),"Age"], y = w[which(w$celltype=="Neuron"),"num"])
#t = -0.82869, df = 22, p-value = 0.4162
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
#  -0.5395138  0.2467310
#sample estimates:
#  cor 
#-0.1739822 
cor.test(x = w[which(w$celltype=="Glia"),"Age"], y = w[which(w$celltype=="Glia"),"num"])
#t = -1.1274, df = 6, p-value = 0.3026
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
#  -0.8672563  0.4062756
#sample estimates:
#  cor 
#-0.4180941 

ggplot(w, aes(x = Age, y = num, colour = celltype)) + geom_point() +
  geom_smooth(method=lm, se=T, fullrange=TRUE) + scale_colour_brewer(8, palette="Dark2") + 
  labs(fill="") + theme_classic() +
  ylab("Number") + xlab("Age (Years)") +
  ggtitle("Number DMVs Identified by Age") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20), legend.title=element_blank()) + theme(legend.position="bottom")
ggplot(num, aes(x = celltype, y = num)) + geom_boxplot() +
  theme_classic() +
  labs(fill="") +
  ylab("Count") + xlab("") +
  ggtitle("Number of DMVs") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
dev.off()


## how long are these regions?

widths = do.call(rbind, Map(cbind, lapply(dmvs, function(y) data.frame(mean = mean(width(y)), median = median(width(y)), 
                                                                       sd = sd(width(y)), min = min(width(y)), max = max(width(y)))), id = as.list(names(dmvs))))
widths[which(widths$id %in% names(DMVs_100kb.pren)),"celltype"] = "Prenatal"
widths[which(widths$id %in% pd[pd$Cell.Type=="Neuron","Data.ID"]),"celltype"] = "Neuron"
widths[which(widths$id %in% pd[pd$Cell.Type=="Glia","Data.ID"]),"celltype"] = "Glia"
write.csv(widths, quote=F, file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/DMV_widths.csv")

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/CREs/figures/DMV_widths.pdf")

w = widths[which(widths$celltype %in% c("Neuron", "Glia")),]
w = cbind(w, Age = pd[match(w$id, pd$Data.ID), "Age"])
cor.test(x = w[which(w$celltype=="Neuron"),"Age"], y = w[which(w$celltype=="Neuron"),"median"])
#t = 0.25929, df = 22, p-value = 0.7978
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
#  -0.3561291  0.4486052
#sample estimates:
#  cor 
#0.05519736
cor.test(x = w[which(w$celltype=="Glia"),"Age"], y = w[which(w$celltype=="Glia"),"median"])
#t = -0.99175, df = 6, p-value = 0.3596
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
#  -0.8540912  0.4478129
#sample estimates:
#  cor 
#-0.3752859 

ggplot(w, aes(x = Age, y = median, colour = celltype)) + geom_point() +
  geom_smooth(method=lm, se=T, fullrange=TRUE) + scale_colour_brewer(8, palette="Dark2") + 
  labs(fill="") + theme_classic() +
  ylab("Bases") + xlab("Age (Years)") +
  ggtitle("Median DMV Widths by Age") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20), legend.title=element_blank()) + theme(legend.position="bottom")
ggplot(widths, aes(x = celltype, y = median)) + geom_boxplot() +
  theme_classic() +
  labs(fill="") +
  ylab("Bases") + xlab("") +
  ggtitle("Median DMV Widths") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
ggplot(widths, aes(x = celltype, y = mean)) + geom_boxplot() +
  theme_classic() +
  labs(fill="") +
  ylab("Bases") + xlab("") +
  ggtitle("Mean DMV Widths") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
dev.off()


## What proportion of the genome do they cover?

load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/PMDs_methylSeekR.rda")
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/PMDs_methylSeekR_prenatal.rda")
sLengths = c(lapply(PMDsegments.CG, reduce), lapply(PMDsegments.CGpren, reduce))
rm(PMDsegments.CGpren, PMDsegments.CG)
percGenome = data.frame(percent = unlist(mapply(function(p, t) round((sum(width(p))/sum(as.numeric(width(t)))*100),2), dmvs, sLengths)),
                        id = names(dmvs))
percGenome[which(percGenome$id %in% names(DMVs_100kb.pren)),"celltype"] = "Prenatal"
percGenome[which(percGenome$id %in% pd[pd$Cell.Type=="Neuron","Data.ID"]),"celltype"] = "Neuron"
percGenome[which(percGenome$id %in% pd[pd$Cell.Type=="Glia","Data.ID"]),"celltype"] = "Glia"
write.csv(percGenome, quote=F, file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/DMV_genomeCoverage.csv")

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/CREs/figures/DMV_percentGenome.pdf")

w = percGenome[which(percGenome$celltype %in% c("Neuron", "Glia")),]
w = cbind(w, Age = pd[match(w$id, pd$Data.ID), "Age"])
cor.test(x = w[which(w$celltype=="Neuron"),"Age"], y = w[which(w$celltype=="Neuron"),"percent"])
#t = -1.2851, df = 22, p-value = 0.2121
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
#  -0.6033294  0.1557545
#sample estimates:
#  cor 
#-0.264245 
cor.test(x = w[which(w$celltype=="Glia"),"Age"], y = w[which(w$celltype=="Glia"),"percent"])
#t = -1.321, df = 6, p-value = 0.2346
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
#  -0.8837414  0.3456068
#sample estimates:
#  cor 
#-0.4746662 

ggplot(w, aes(x = Age, y = percent, colour = celltype)) + geom_point() +
  geom_smooth(method=lm, se=T, fullrange=TRUE) + scale_colour_brewer(8, palette="Dark2") + 
  labs(fill="") +  theme_classic() +
  ylab("Percent") + xlab("Age (Years)") +
  ggtitle("Percent of Genome that is DMV by Age") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20), legend.title=element_blank()) + theme(legend.position="bottom")
ggplot(percGenome, aes(x = celltype, y = percent)) + geom_boxplot() +
  theme_classic() +
  labs(fill="") +
  ylab("Percent") + xlab("") +
  ggtitle("Percent Genome: DMV") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
dev.off()


### How much of the dmvs overlap?

# Master list of all bases in each setting
dmaster = reduce(do.call(getMethod(c, "GenomicRanges"), GRangesList(dmvs))) # 43995 total
c(mean=mean(width(dmaster)), median=median(width(dmaster)), sd=sd(width(dmaster)), min=min(width(dmaster)), max=max(width(dmaster)))/1000
#         mean    median        sd       min       max 
#kb:  13.33677   9.27300  15.81653   5.00000 451.64300 



## regions that are shared by samples

all = list(All = Reduce(intersect, dmvs),
           Prenatal = Reduce(intersect, DMVs_100kb.pren),
           Postnatal = Reduce(intersect, DMVs.100kb),
           Neurons = Reduce(intersect, DMVs.100kb[which(names(DMVs.100kb) %in% pd[pd$Cell.Type=="Neuron","Data.ID"])]),
           Glia = Reduce(intersect, DMVs.100kb[which(names(DMVs.100kb) %in% pd[pd$Cell.Type=="Glia","Data.ID"])]))
elementNROWS(all)
#      All  Prenatal Postnatal   Neurons      Glia 
#      318       559       399       588       554 

round((unlist(lapply(all, function(x) (sum(width(x))/sum(as.numeric(width(sLengths[[1]]))))*100)))/
        (sum(width(dmaster))/sum(as.numeric(width(sLengths[[1]])))*100),2)
#      All  Prenatal Postnatal   Neurons      Glia 
#     0.01      0.03      0.02      0.03      0.03 
# 1-3% of total bases in the DMV state are shared by all or most

shared = mapply(function(all,ind) lapply( ind, function(x) round(all / x *100,2)), lapply(all, function(x) (sum(width(x))/sum(as.numeric(width(sLengths[[1]])))*100)), 
                lapply( list(dmvs,DMVs_100kb.pren,DMVs.100kb,
                             DMVs.100kb[which(names(DMVs.100kb) %in% pd[pd$Cell.Type=="Neuron","Data.ID"])],
                             DMVs.100kb[which(names(DMVs.100kb) %in% pd[pd$Cell.Type=="Glia","Data.ID"])]),
                        function(y) lapply(y, function(x) sum(width(x))/sum(as.numeric(width(sLengths[[1]])))*100)), SIMPLIFY = F)
shared = do.call(rbind, Map(data.frame, percent = lapply(shared, function(x) unlist(x, recursive=F)), id = lapply(shared, names), group = as.list(names(shared))))

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/CREs/figures/percent_shared_dmvs_perSample_byBaseCoverage.pdf")
ggplot(shared, aes(x = group, y = percent)) + geom_boxplot() +
  labs(fill="") + theme_classic() +
  ylim(0,100) +
  ylab("Percent") + xlab("") +
  ggtitle("Percent Bases Shared Per Sample: DMV") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
dev.off()


## Find overlaps with DMRs, features, shared groups

txdb = loadDb("/dcl01/lieber/ajaffe/Amanda/annotation_objects/gencode.v25lift37.annotation.sqlite")
islands = read.table("/dcl01/lieber/ajaffe/Amanda/annotation_objects/cpgIslandExt.hg19.txt", sep="\t", header = T)
rpmsk = read.table("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/RepeatMasker_genomewide_HG19.txt", header=T)
features = c(genes = makeGRangesFromDataFrame(geneMap, keep.extra.columns = T), 
             rpmskgr = makeGRangesFromDataFrame(rpmsk, seqnames.field="genoName",start.field="genoStart",end.field="genoEnd",keep.extra.columns=TRUE),
             islands = makeGRangesFromDataFrame(islands, keep.extra.columns = T, start.field = "chromStart", end.field = "chromEnd"),
             promoters = promoters(txdb, upstream=2000, downstream=200))
lapply(features, head)

DMR = lapply(DMR, function(x) x[which(x$sig=="FWER < 0.05"),])
DMRgr = lapply(DMR, function(x) makeGRangesFromDataFrame(x, keep.extra.columns = T))

oo = lapply(dmvs, function(x) lapply(c(DMRgr, all, features), function(d) findOverlaps(x,d)))

dmvDMR = lapply(dmvs, as.data.frame)
dmvDMR = lapply(dmvDMR, function(x) data.frame(x, rnum = 1:nrow(x)))
dmvDMR = mapply(function(d,oo) data.frame(d, CT = ifelse(d$rnum %in% queryHits(oo$CellType), "CT", "no"),
                                          Age = ifelse(d$rnum %in% queryHits(oo$Age), "Age", "no"),
                                          Interaction = ifelse(d$rnum %in% queryHits(oo$Interaction), "Int", "no"),
                                          All = ifelse(d$rnum %in% queryHits(oo$All), "All", "no"),
                                          Prenatal = ifelse(d$rnum %in% queryHits(oo$Prenatal), "Prenatal", "no"),
                                          Postnatal = ifelse(d$rnum %in% queryHits(oo$Postnatal), "Postnatal", "no"),
                                          Neurons = ifelse(d$rnum %in% queryHits(oo$Neurons), "Neurons", "no"),
                                          Glia = ifelse(d$rnum %in% queryHits(oo$Glia), "Glia", "no"),
                                          genes = ifelse(d$rnum %in% queryHits(oo$genes), "Gene", NA),
                                          islands = ifelse(d$rnum %in% queryHits(oo$islands), "CpG-Island", "non-Island"),
                                          promoters = ifelse(d$rnum %in% queryHits(oo$promoters), "promoter", "no")), dmvDMR, oo, SIMPLIFY = F) 
dmvDMR = lapply(dmvDMR, function(d) data.frame(d, tog = paste(d$CT, d$Age, d$Interaction, d$All, d$Prenatal, 
                                                              d$Postnatal, d$Neurons, d$Glia, sep = ":"),
                                               dmr = paste(d$CT, d$Age, d$Interaction, sep = ":"),
                                               regionID = paste0(d$seqnames,":",d$start,"-", d$end)))
for (i in 1:length(dmvDMR)) {
  if (length(unique(queryHits(oo[[i]][["rpmskgr"]])))==nrow(dmvDMR[[i]])) {
    dmvDMR[[i]] = data.frame(dmvDMR[[i]][queryHits(oo[[i]][["rpmskgr"]]),], 
                             repeats = as.character(features$rpmskgr$repClass)[subjectHits(oo[[i]][["rpmskgr"]])]) } else {
                               dmvDMR[[i]] = rbind(data.frame(dmvDMR[[i]][queryHits(oo[[i]][["rpmskgr"]]),], 
                                                              repeats = as.character(features$rpmskgr$repClass)[subjectHits(oo[[i]][["rpmskgr"]])]), 
                                                   data.frame(dmvDMR[[i]][-unique(queryHits(oo[[i]][["rpmskgr"]])),], repeats = "No repeats"))                    
                             }
}

dmvDMRgr = lapply(dmvDMR, function(x) makeGRangesFromDataFrame(x, keep.extra.columns = T))
oo = lapply(dmvDMRgr, function(x) findOverlaps(x, features$genes))
dmvDMR = mapply(function(p, ov) rbind(data.frame(p[queryHits(ov),], nearestSymbol = features$genes$Symbol[subjectHits(ov)],
                                                 nearestID = names(features$genes)[subjectHits(ov)], EntrezID = features$genes$EntrezID[subjectHits(ov)]),
                                      data.frame(p[-unique(queryHits(ov)),], nearestSymbol = "NoGeneOverlap", nearestID = "NoGeneOverlap", EntrezID = "NoGeneOverlap")),
                dmvDMR, oo, SIMPLIFY = F)

postnatalpd = pd
save(dmvDMR, postnatalpd, geneMap, file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/dmvs_annotated.rda")


### Annotate to genomic features

load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/dmvs_annotated.rda")
load('/dcl01/lieber/ajaffe/lab/brain-epigenomics/bumphunting/BSobj_bsseqSmooth_Neuron_minCov_3.Rdata')

# Identify all CpG clusters in the genome
gr = granges(BSobj)
cl = clusterMaker( chr = as.character(seqnames(gr)), pos = start(gr),   maxGap = 1000)
gr.clusters = split(gr, cl)
gr.clusters = unlist(range(gr.clusters))
df.clusters = as.data.frame(gr.clusters)
df.clusters$regionID = paste0(df.clusters$seqnames,":",df.clusters$start,"-",df.clusters$end)
df.clusters$rnum = 1:length(gr.clusters)
df.clusters$Islands = ifelse(df.clusters$rnum %in% subjectHits(findOverlaps(features$islands, gr.clusters)), "island","no")
df.clusters$repeats = ifelse(df.clusters$rnum %in% subjectHits(findOverlaps(features$rpmskgr, gr.clusters)), "repeat","no")
df.clusters$genes = ifelse(df.clusters$rnum %in% subjectHits(findOverlaps(makeGRangesFromDataFrame(geneMap), gr.clusters)), "gene","no")
df.clusters$promoters = ifelse(df.clusters$rnum %in% subjectHits(findOverlaps(makeGRangesFromDataFrame(geneMap), gr.clusters)), "promoter","no")

oo = lapply(dmvDMR, function(y) findOverlaps(gr.clusters, makeGRangesFromDataFrame(y)))
df.clusters = lapply(oo, function(y) data.frame(df.clusters, overlap = ifelse(df.clusters$rnum %in% queryHits(y), "hit","no")))


## how many fall within CpG islands?

CpGIslands = lapply(dmvDMR, function(y) round(length(unique(y[which(y$islands=="CpG-Island"),,]$regionID))/length(unique(y$regionID))*100,2))
CpGIslands = data.frame(islands = unlist(CpGIslands), id = names(CpGIslands))
CpGIslands[!(CpGIslands$id %in% postnatalpd$Data.ID),"celltype"] = "Prenatal"
CpGIslands[which(CpGIslands$id %in% postnatalpd[which(postnatalpd$Cell.Type=="Neuron"),"Data.ID"]),"celltype"] = "Neuron"
CpGIslands[which(CpGIslands$id %in% postnatalpd[which(postnatalpd$Cell.Type=="Glia"),"Data.ID"]),"celltype"] = "Glia"
write.csv(CpGIslands, quote=F,file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/DMV_CpG_Island_Overlap.csv")

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/CREs/figures/DMV_overlap_with_CpG-Islands.pdf")
ggplot(CpGIslands, aes(x = celltype, y = islands)) + geom_boxplot() +
  theme_classic() +
  labs(fill="") +
  ylab("Percent") + ylim(0,100) + 
  xlab("") +
  ggtitle("Percent Overlapping CpG Islands") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20), legend.title=element_blank()) + theme(legend.position="bottom")
dev.off()


# How about Repetitive Elements?

dmvDMRdt = lapply(dmvDMR, data.table)
repeats = lapply(dmvDMRdt, function(y) y[,length(unique(regionID)), by = "repeats"])
repeats = lapply(repeats, function(y) data.frame(y, perc = round(y$V1/sum(y$V1)*100,2)))
repeats = do.call(rbind, Map(cbind, repeats, id = as.list(names(repeats))))
repeats$repeats = factor(repeats$repeats, levels = c("SINE","LINE","Simple_repeat","DNA","LTR","Low_complexity","No repeats","Unknown","snRNA",
                                                     "Other","DNA?","Satellite","RC","scRNA","tRNA","SINE?","srpRNA","RNA","rRNA","LINE?","LTR?","Unknown?"))
repeats[!(repeats$id %in% postnatalpd$Data.ID),"celltype"] = "Prenatal"
repeats[which(repeats$id %in% postnatalpd[which(postnatalpd$Cell.Type=="Neuron"),"Data.ID"]),"celltype"] = "Neuron"
repeats[which(repeats$id %in% postnatalpd[which(postnatalpd$Cell.Type=="Glia"),"Data.ID"]),"celltype"] = "Glia"
write.csv(repeats, quote=F,file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/DMV_repeats_Overlap.csv")

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/CREs/figures/DMV_overlap_with_repeats.pdf")
x = repeats[which(repeats$repeats %in% c("SINE","LINE","Simple_repeat","DNA","LTR","Low_complexity","No repeats")),]
x$repeats = droplevels(x$repeats)
ggplot(x, aes(x = celltype, y = perc, fill = repeats)) + 
  theme_classic() + geom_bar(position = "fill",stat = "identity") +
  labs(fill="") +
  ylab("Percent") + 
  xlab("") +
  ggtitle("Percent Overlapping Repeats: DMV") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20), legend.title=element_blank()) + theme(legend.position="bottom")
dev.off()


## Overlap with DMRs

dmr = lapply(dmvDMRdt, function(y) y[,length(unique(regionID)), by = "dmr"])
dmr = lapply(dmr, function(y) data.frame(y, perc = round(y$V1/sum(y$V1)*100,2)))
models = lapply(dmr, function(y) data.frame(CT = sum(data.frame(y)[grep("CT", y$dmr),"perc"]), 
                                            Age = sum(data.frame(y)[grep("Age", y$dmr),"perc"]),
                                            Int = sum(data.frame(y)[grep("Int", y$dmr),"perc"])))
dmr = do.call(rbind, Map(cbind, dmr, id = as.list(names(dmr))))
models = do.call(rbind, Map(cbind, models, id = as.list(names(models))))

dmr$dmr = factor(dmr$dmr, levels = c("no:no:no","CT:no:no","no:no:Int","CT:no:Int","no:Age:no","CT:Age:Int","no:Age:Int","CT:Age:no"))
dmr[!(dmr$id %in% postnatalpd$Data.ID),"celltype"] = "Prenatal"
dmr[which(dmr$id %in% postnatalpd[which(postnatalpd$Cell.Type=="Neuron"),"Data.ID"]),"celltype"] = "Neuron"
dmr[which(dmr$id %in% postnatalpd[which(postnatalpd$Cell.Type=="Glia"),"Data.ID"]),"celltype"] = "Glia"
models[!(models$id %in% postnatalpd$Data.ID),"celltype"] = "Prenatal"
models[which(models$id %in% postnatalpd[which(postnatalpd$Cell.Type=="Neuron"),"Data.ID"]),"celltype"] = "Neuron"
models[which(models$id %in% postnatalpd[which(postnatalpd$Cell.Type=="Glia"),"Data.ID"]),"celltype"] = "Glia"

write.csv(dmr, quote=F,file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/DMV_dmr_Overlap.csv")
write.csv(models, quote=F,file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/DMV_dmr_Overlap_byModel.csv")


pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/CREs/figures/DMV_overlap_with_dmr.pdf", width = 14)
ggplot(dmr, aes(x = celltype, y = perc, fill = dmr)) + 
  theme_classic() + geom_bar(position = "fill",stat = "identity") +
  labs(fill="") + 
  ylab("Percent") + 
  xlab("") +
  ggtitle("Percent Overlapping DMRs") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20), legend.title=element_blank()) + theme(legend.position="bottom")
ggplot(dmr[which(dmr$dmr=="no:no:no"),], aes(x = celltype, y = perc, fill = dmr)) + 
  theme_classic() + geom_boxplot() +
  labs(fill="") + 
  ylab("Percent") + 
  xlab("") +
  ggtitle("Percent Not Overlapping DMRs") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20), legend.title=element_blank()) + theme(legend.position="bottom")
ggplot(models, aes(x = celltype, y = CT)) + 
  theme_classic() + geom_boxplot() +
  labs(fill="") + 
  ylab("Percent") + 
  xlab("") +
  ggtitle("Percent Overlapping Cell Type DMRs") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20), legend.title=element_blank()) + theme(legend.position="bottom")
ggplot(models, aes(x = celltype, y = Age)) + 
  theme_classic() + geom_boxplot() +
  labs(fill="") + 
  ylab("Percent") + 
  xlab("") +
  ggtitle("Percent Overlapping Age DMRs") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20), legend.title=element_blank()) + theme(legend.position="bottom")
ggplot(models, aes(x = celltype, y = Int)) + 
  theme_classic() + geom_boxplot() +
  labs(fill="") + 
  ylab("Percent") + 
  xlab("") +
  ggtitle("Percent Overlapping Interaction DMRs") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20), legend.title=element_blank()) + theme(legend.position="bottom")
dev.off()


## association of existing DMVs with genes

geneoverlap = lapply(dmvDMR, function(y) round(length(unique(as.character(data.frame(y)[which(y$nearestID!="NoGeneOverlap"),"regionID"])))/
                                                 length(unique(as.character(y$regionID)))*100,2))
geneoverlap = data.frame(geneoverlap = unlist(geneoverlap), id = names(geneoverlap))
geneoverlap[!(geneoverlap$id %in% postnatalpd$Data.ID),"celltype"] = "Prenatal"
geneoverlap[which(geneoverlap$id %in% postnatalpd[which(postnatalpd$Cell.Type=="Neuron"),"Data.ID"]),"celltype"] = "Neuron"
geneoverlap[which(geneoverlap$id %in% postnatalpd[which(postnatalpd$Cell.Type=="Glia"),"Data.ID"]),"celltype"] = "Glia"
write.csv(geneoverlap, quote=F,file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/DMV_Gene_Overlap.csv")

numgenes = lapply(dmvDMRdt, function(y) y[,length(unique(nearestID)), by = "regionID"])
numgenes = do.call(rbind, Map(cbind, lapply(numgenes, function(y) 
           data.frame(mean = mean(y$V1), median = median(y$V1), sd = sd(y$V1), min = min(y$V1), max = max(y$V1))), id = as.list(names(numgenes))))
numgenes[!which(numgenes$id %in% postnatalpd$Data.ID),"celltype"] = "Prenatal"
numgenes[which(numgenes$id %in% postnatalpd[postnatalpd$Cell.Type=="Neuron","Data.ID"]),"celltype"] = "Neuron"
numgenes[which(numgenes$id %in% postnatalpd[postnatalpd$Cell.Type=="Glia","Data.ID"]),"celltype"] = "Glia"
write.csv(numgenes, quote=F, file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/DMV_numgenes.csv")

w = numgenes[which(numgenes$celltype %in% c("Neuron", "Glia")),]
w = cbind(w, Age = postnatalpd[match(w$id, postnatalpd$Data.ID), "Age"])
cor.test(x = w[which(w$celltype=="Neuron"),"Age"], y = w[which(w$celltype=="Neuron"),"mean"])
#t = -0.86357, df = 22, p-value = 0.3971
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
#  -0.5446817  0.2398456
#sample estimates:
#  cor 
#-0.1810702 
cor.test(x = w[which(w$celltype=="Glia"),"Age"], y = w[which(w$celltype=="Glia"),"mean"])
#t = -1.0309, df = 6, p-value = 0.3423
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
#  -0.8580426  0.4359092
#sample estimates:
#  cor 
#-0.3879235 

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/CREs/figures/DMV_overlap_with_Genes.pdf")
ggplot(geneoverlap, aes(x = celltype, y = geneoverlap)) + geom_boxplot() +
  theme_classic() +
  labs(fill="") +
  ylab("Percent") + ylim(0,100) + 
  xlab("") +
  ggtitle("Percent Overlapping Genes") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20), legend.title=element_blank()) + theme(legend.position="bottom")
ggplot(w, aes(x = Age, y = mean, colour = celltype)) + geom_point() +
  geom_smooth(method=lm, se=T, fullrange=TRUE) + scale_colour_brewer(8, palette="Dark2") + 
  labs(fill="") + theme_classic() +
  ylab("Bases") + xlab("Age (Years)") +
  ggtitle("Median DMV Gene Number by Age") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20), legend.title=element_blank()) + theme(legend.position="bottom")
ggplot(numgenes, aes(x = celltype, y = median)) + geom_boxplot() +
  theme_classic() +
  labs(fill="") +
  ylab("Bases") + xlab("") +
  ggtitle("Median DMV Gene Number") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
ggplot(numgenes, aes(x = celltype, y = mean)) + geom_boxplot() +
  theme_classic() +
  labs(fill="") +
  ylab("Bases") + xlab("") +
  ggtitle("Mean DMV Gene Number") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
dev.off()


promoverlap = lapply(dmvDMR, function(y) round(length(unique(as.character(data.frame(y)[which(y$promoters=="promoter"),"regionID"])))/
                                                 length(unique(as.character(y$regionID)))*100,2))
promoverlap = data.frame(promoverlap = unlist(promoverlap), id = names(promoverlap))
promoverlap[!(promoverlap$id %in% postnatalpd$Data.ID),"celltype"] = "Prenatal"
promoverlap[which(promoverlap$id %in% postnatalpd[which(postnatalpd$Cell.Type=="Neuron"),"Data.ID"]),"celltype"] = "Neuron"
promoverlap[which(promoverlap$id %in% postnatalpd[which(postnatalpd$Cell.Type=="Glia"),"Data.ID"]),"celltype"] = "Glia"
write.csv(promoverlap, quote=F,file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/DMV_prom_Overlap.csv")

numproms = lapply(dmvDMRdt, function(y) y[promoters=="promoter",length(unique(nearestID)), by = "regionID"])
numproms = do.call(rbind, Map(cbind, lapply(numproms, function(y) 
  data.frame(mean = mean(y$V1), median = median(y$V1), sd = sd(y$V1), min = min(y$V1), max = max(y$V1))), id = as.list(names(numproms))))
numproms[!which(numproms$id %in% postnatalpd$Data.ID),"celltype"] = "Prenatal"
numproms[which(numproms$id %in% postnatalpd[postnatalpd$Cell.Type=="Neuron","Data.ID"]),"celltype"] = "Neuron"
numproms[which(numproms$id %in% postnatalpd[postnatalpd$Cell.Type=="Glia","Data.ID"]),"celltype"] = "Glia"
write.csv(numproms, quote=F, file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/DMV_numproms.csv")

w = numproms[which(numproms$celltype %in% c("Neuron", "Glia")),]
w = cbind(w, Age = postnatalpd[match(w$id, postnatalpd$Data.ID), "Age"])
cor.test(x = w[which(w$celltype=="Neuron"),"Age"], y = w[which(w$celltype=="Neuron"),"mean"])
#t = -1.2137, df = 22, p-value = 0.2377
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
#  -0.5938903  0.1700746
#sample estimates:
#  cor 
#-0.2505089 
cor.test(x = w[which(w$celltype=="Glia"),"Age"], y = w[which(w$celltype=="Glia"),"mean"])
#t = -1.1737, df = 6, p-value = 0.285
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
#  -0.8714336  0.3918957
#sample estimates:
#  cor 
#-0.4321072 

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/CREs/figures/DMV_overlap_with_proms.pdf")
ggplot(promoverlap, aes(x = celltype, y = promoverlap)) + geom_boxplot() +
  theme_classic() +
  labs(fill="") +
  ylab("Percent") + ylim(0,100) + 
  xlab("") +
  ggtitle("Percent Overlapping Promoters") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20), legend.title=element_blank()) + theme(legend.position="bottom")
ggplot(w, aes(x = Age, y = mean, colour = celltype)) + geom_point() +
  geom_smooth(method=lm, se=T, fullrange=TRUE) + scale_colour_brewer(8, palette="Dark2") + 
  labs(fill="") + theme_classic() +
  ylab("Bases") + xlab("Age (Years)") +
  ggtitle("Median DMV Promoter Number by Age") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20), legend.title=element_blank()) + theme(legend.position="bottom")
ggplot(numproms, aes(x = celltype, y = median)) + geom_boxplot() +
  theme_classic() +
  labs(fill="") +
  ylab("Bases") + xlab("") +
  ggtitle("Median DMV Promoter Number") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
ggplot(numproms, aes(x = celltype, y = mean)) + geom_boxplot() +
  theme_classic() +
  labs(fill="") +
  ylab("Bases") + xlab("") +
  ggtitle("Mean DMV Promoter Number") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
dev.off()


## Identify the regions that are getting smaller or larger over age

dmvDMRdt = Map(cbind, dmvDMRdt, tog2 = lapply(dmvDMR, function(d) paste(d$All, d$Prenatal, d$Postnatal, d$Neurons, d$Glia, sep = ":")))
lapply(dmvDMRdt, function(x) x[,length(unique(regionID)),by="tog2"])

subset = do.call(rbind, Map(cbind, lapply(dmvDMRdt, function(x) data.frame(regionID = x$regionID, width = x$width, tog2 = x$tog2)), id = as.list(names(dmvDMRdt))))
subset[!which(subset$id %in% postnatalpd$Data.ID),"celltype"] = "Prenatal"
subset[which(subset$id %in% postnatalpd[postnatalpd$Cell.Type=="Neuron","Data.ID"]),"celltype"] = "Neuron"
subset[which(subset$id %in% postnatalpd[postnatalpd$Cell.Type=="Glia","Data.ID"]),"celltype"] = "Glia"
subset = unique(subset)
w = subset[which(subset$celltype %in% c("Neuron", "Glia")),]
w = cbind(w, Age = postnatalpd[match(w$id, postnatalpd$Data.ID), "Age"])
w = split(w, w$tog2)
elementNROWS(w)
#All:Prenatal:Postnatal:Neurons:Glia            9200
#no:no:no:Neurons:no                            3423
#no:no:no:no:Glia                               1746
#no:no:no:no:no                                15341
#no:no:Postnatal:Neurons:Glia                   2756
#no:Prenatal:no:Neurons:no                      1772
#no:Prenatal:no:no:Glia                         1630
#no:Prenatal:no:no:no                           1387
#no:Prenatal:Postnatal:Neurons:Glia                7
#no:Prenatal:no:Neurons:Glia                       1

w = lapply(w, data.table)
mean = lapply(w, function(x) data.frame(x[,mean(width), by = c("id","celltype","Age")]))
med = lapply(w, function(x) data.frame(x[,median(width), by = c("id","celltype","Age")]))

nwidth = lapply(mean[1:9], function(z) cor.test(z[which(z$celltype=="Neuron"),"Age"], z[which(z$celltype=="Neuron"),"V1"]))
gwidth = lapply(mean[1:8], function(z) cor.test(z[which(z$celltype=="Glia"),"Age"], z[which(z$celltype=="Glia"),"V1"]))
mnwidth = lapply(med[1:9], function(z) cor.test(z[which(z$celltype=="Neuron"),"Age"], z[which(z$celltype=="Neuron"),"V1"]))
mgwidth = lapply(med[1:8], function(z) cor.test(z[which(z$celltype=="Glia"),"Age"], z[which(z$celltype=="Glia"),"V1"]))

write.csv(rbind(data.frame(Category = names(nwidth), tstat = unlist(lapply(nwidth, function(x) x$statistic)),
                           parameter = unlist(lapply(nwidth, function(x) x$parameter)), p.value = unlist(lapply(nwidth, function(x) x$p.value)),
                           estimate = unlist(lapply(nwidth, function(x) x$estimate)), CellType = "Neurons", Comparison = "Mean"),
                data.frame(Category = names(gwidth), tstat = unlist(lapply(gwidth, function(x) x$statistic)),
                           parameter = unlist(lapply(gwidth, function(x) x$parameter)), p.value = unlist(lapply(gwidth, function(x) x$p.value)),
                           estimate = unlist(lapply(gwidth, function(x) x$estimate)), CellType = "Glia", Comparison = "Mean"),
                data.frame(Category = names(mnwidth), tstat = unlist(lapply(mnwidth, function(x) x$statistic)),
                           parameter = unlist(lapply(mnwidth, function(x) x$parameter)), p.value = unlist(lapply(mnwidth, function(x) x$p.value)),
                           estimate = unlist(lapply(mnwidth, function(x) x$estimate)), CellType = "Neurons", Comparison = "Median"),
                data.frame(Category = names(mgwidth), tstat = unlist(lapply(mgwidth, function(x) x$statistic)),
                           parameter = unlist(lapply(mgwidth, function(x) x$parameter)), p.value = unlist(lapply(mgwidth, function(x) x$p.value)),
                           estimate = unlist(lapply(mgwidth, function(x) x$estimate)), CellType = "Glia", Comparison = "Median")), quote = F,
          file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/shared_DMVs_correlation.csv")

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/CREs/figures/shared_DMV_width_correlation.pdf")
for (i in 1:length(mean)) {
g = ggplot(mean[[i]], aes(x = Age, y = V1, colour = celltype)) + geom_point() +
      geom_smooth(method=lm, se=T, fullrange=TRUE) + scale_colour_brewer(8, palette="Dark2") + 
      labs(fill="") + theme_classic() +
      ylab("Width") + xlab("Age (Years)") +
      ggtitle(paste0("Mean Width of Overlapping DMVs:\n",names(w)[i])) +
      theme(title = element_text(size = 20)) +
      theme(text = element_text(size = 20), legend.title=element_blank()) + theme(legend.position="bottom")
print(g)
g = ggplot(med[[i]], aes(x = Age, y = V1, colour = celltype)) + geom_point() +
  geom_smooth(method=lm, se=T, fullrange=TRUE) + scale_colour_brewer(8, palette="Dark2") + 
  labs(fill="") + theme_classic() +
  ylab("Width") + xlab("Age (Years)") +
  ggtitle(paste0("Median Width of Overlapping DMVs:\n",names(w)[i])) +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20), legend.title=element_blank()) + theme(legend.position="bottom")
print(g)
}
dev.off()


## Gene ontology of genes in different categories

subset = do.call(rbind, Map(cbind, lapply(dmvDMRdt, function(x) data.frame(regionID = x$regionID, EntrezID = x$EntrezID, tog2 = x$tog2)),
                            id = as.list(names(dmvDMRdt))))
subset[-which(subset$id %in% postnatalpd$Data.ID),"celltype"] = "Prenatal"
subset[which(subset$id %in% postnatalpd[postnatalpd$Cell.Type=="Neuron","Data.ID"]),"celltype"] = "Neuron"
subset[which(subset$id %in% postnatalpd[postnatalpd$Cell.Type=="Glia","Data.ID"]),"celltype"] = "Glia"
subset = unique(subset)
subset = split(subset, subset$tog2)
subset = lapply(subset, function(x) split(x$EntrezID, x$celltype))
subset = lapply(subset, function(x) lapply(x, function(y) na.omit(unique(as.character(y)))))
subset = lapply(subset, function(x) lapply(x, function(y) y[which(y!="NoGeneOverlap")]))

# Compare the enriched terms between 7 groups
# KEGG
compareKegg = lapply(subset, function(x) compareCluster(x, fun="enrichKEGG", qvalueCutoff = 0.05, pvalueCutoff = 0.05))
compareKegg.collapsed = compareCluster(unlist(subset, recursive = F), fun="enrichKEGG", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
# Biological Process
compareBP = lapply(subset, function(x) compareCluster(x, fun="enrichGO", ont = "BP", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05))
compareBP.collapsed = compareCluster(unlist(subset, recursive = F), fun="enrichGO", ont = "BP", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
# Molecular Function
compareMF = lapply(subset, function(x) compareCluster(x, fun="enrichGO",  ont = "MF", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05))
compareMF.collapsed = compareCluster(unlist(subset, recursive = F), fun="enrichGO",  ont = "MF", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
# Cellular Component
compareCC = lapply(subset, function(x) compareCluster(x, fun="enrichGO",  ont = "CC", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05))
compareCC.collapsed = compareCluster(unlist(subset, recursive = F), fun="enrichGO",  ont = "CC", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
# Disease Ontology
compareDO = lapply(subset, function(x) compareCluster(x, fun="enrichDO",  ont = "DO", qvalueCutoff = 0.05, pvalueCutoff = 0.05))
compareDO.collapsed = compareCluster(unlist(subset, recursive = F), fun="enrichDO",  ont = "DO", qvalueCutoff = 0.05, pvalueCutoff = 0.05)

# save object
save(compareKegg.collapsed, compareBP.collapsed, compareMF.collapsed, compareCC.collapsed, compareDO.collapsed,
     file="/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/DMV_KEGG_GO_DO_objects.rda")

# plot compared results
pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/CREs/figures/DMV_KEGG_GO_DO_plots.pdf", height = 80, width = 24)
plot(compareKegg.collapsed, colorBy="p.adjust", showCategory = 1000, title= "KEGG Pathway Enrichment")
plot(compareMF.collapsed, colorBy="p.adjust", showCategory = 1000, title= "Molecular Function GO Enrichment")
plot(compareCC.collapsed, colorBy="p.adjust", showCategory = 1000, title= "Cellular Compartment GO Enrichment")
plot(compareDO.collapsed, colorBy="p.adjust", showCategory = 1000, title= "Disease Ontology Enrichment")
dev.off()
pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/CREs/figures/DMV_BP_plot.pdf", height = 275, width = 24)
plot(compareBP.collapsed, colorBy="p.adjust", showCategory = 1500, title= "Biological Process GO Enrichment")
dev.off()


#### association with H3K27me3, H3K9me and other states

gr <- granges(BSobj)
cl = clusterMaker( chr = as.character(seqnames(gr)), 
                   pos = start(gr),   maxGap = 1000)
gr.clusters = split(gr, cl)
gr.clusters = unlist(range(gr.clusters))

# load roadmap states

load("/dcl01/lieber/ajaffe/PublicData/EpigenomeRoadmap/ChromHMM/chromHMM_15state_coverageDF.rda")
dat = read.delim("/dcl01/lieber/ajaffe/PublicData/EpigenomeRoadmap/ChromHMM/jul2013.roadmapData.qc - Consolidated_EpigenomeIDs_summary_Table.tsv",
                 as.is=TRUE)
colnames(dat)[2] = "EID"
dat = dat[dat$EID != "",]

# split out DMRs per background

dmvList = do.call(rbind, Map(cbind, lapply(dmvDMRdt, function(x) data.frame(x[,1:5], regionID = x$regionID, tog2 = x$tog2)), id = as.list(names(dmvDMRdt))))
dmvList[-which(dmvList$id %in% postnatalpd$Data.ID),"celltype"] = "Prenatal"
dmvList[which(dmvList$id %in% postnatalpd[postnatalpd$Cell.Type=="Neuron","Data.ID"]),"celltype"] = "Neuron"
dmvList[which(dmvList$id %in% postnatalpd[postnatalpd$Cell.Type=="Glia","Data.ID"]),"celltype"] = "Glia"
dmvList = unique(dmvList)
dmvList = unlist(lapply(split(dmvList, dmvList$tog2), function(x) split(x, x$celltype)), recursive = F)
dmvList = lapply(dmvList, makeGRangesFromDataFrame, keep=TRUE)
dmvList = GRangesList(lapply(dmvList, reduce))

bgWithDmv = endoapply(dmvList, function(x) {
  x = granges(x)
  m = c(x, gr.clusters)
  d = disjoin(m)
  d$inDMV = countOverlaps(d,x) > 0
  return(d)
})


#### split by chromosome 
DMVsByChr = lapply(bgWithDmv, function(x) split(x, seqnames(x)))

## autosomal
DMVs_bpOverlapByChr = lapply(DMVsByChr, function(x) mapply(function(regByChr, stateByChr) {
  cat(".")
  inDMV = regByChr[regByChr$inDMV]
  outDMV = regByChr[!regByChr$inDMV]
  
  enr = stateByChr[ranges(inDMV),]
  bg = stateByChr[ranges(outDMV),]
  
  list(enrTab = sapply(enr, table),
       bgTab  = sapply(bg, table))
}, x[paste0("chr",1:22)], stateCovDf, SIMPLIFY=FALSE))

## states
DMVs_bpOverlapEnr = lapply(DMVs_bpOverlapByChr, function(x) sapply(x, "[", "enrTab"))
bpOverlapArray = lapply(DMVs_bpOverlapEnr, function(x) array(unlist(x), dim = c(nrow(x[[1]]), ncol(x[[1]]), length(x))))
DMVs_enrTab = lapply(bpOverlapArray, function(x) apply(x, 1:2, sum))
for (i in 1:length(DMVs_enrTab)) {
dimnames(DMVs_enrTab[[i]]) = list(rownames(DMVs_bpOverlapEnr[[i]][[1]]), colnames(DMVs_bpOverlapEnr[[i]][[1]]))
}

# states in rest of genome	
DMVs_bpOverlapBg = lapply(DMVs_bpOverlapByChr, function(x) sapply(x, "[", "bgTab"))
bpOverlapArray = lapply(DMVs_bpOverlapBg, function(x) array(unlist(x), dim = c(nrow(x[[1]]), ncol(x[[1]]), length(x))))
DMVs_bgTab = lapply(bpOverlapArray, function(y) apply(y, 1:2, function(x) sum(as.numeric(x))))
for (i in 1:length(DMVs_bgTab)) {
dimnames(DMVs_bgTab[[i]]) = list(rownames(DMVs_bpOverlapBg[[i]][[1]]), colnames(DMVs_bpOverlapBg[[i]][[1]]))
}

# take ratio
DMVs_statTab = mapply(function(enr,bg) as.data.frame(t(prop.table(enr,2) / prop.table(bg,2))), DMVs_enrTab, DMVs_bgTab, SIMPLIFY = F)
DMVs_statTab = Map(cbind, DMVs_statTab, Sample = lapply(DMVs_statTab, function(x) dat$Standardized.Epigenome.name[match(rownames(x), dat$EID)]))
DMVs_statTab = lapply(DMVs_statTab, function(x) x[,c(16,1:15)])
names(DMVs_statTab) = c("A:Pr:Po:N:G-Glia","A:Pr:Po:N:G-Neuron","A:Pr:Po:N:G-Prenatal",":::N:-Glia",":::N:-Neuron",":::N:-Prenatal","::::G-Glia",
                        "::::G-Neuron","::::G-Prenatal","::::-Glia","::::-Neuron","::::-Prenatal","::Po:N:G-Glia","::Po:N:G-Neuron","::Po:N:G-Prenatal",
                        ":Pr::N:-Glia",":Pr::N:-Neuron",":Pr::N:-Prenatal",":Pr:::G-Glia",":Pr:::G-Neuron",":Pr:::G-Prenatal",":Pr:::-Glia",
                        ":Pr:::-Neuron",":Pr:::-Prenatal",":Pr:Po:N:G-Glia",":Pr:Po:N:G-Neuron",":Pr:Po:N:G-Prenatal",":Pr::N:G-Neuron")
save(DMVs_statTab, file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/DMV_chromHMM_15state_enrichment.rda")

dlpfcStats = do.call(rbind, lapply(DMVs_statTab, function(x) x["E073",-1]))
rownames(dlpfcStats) = names(DMVs_statTab)
dlpfcStats = t(signif(dlpfcStats,3))

tmpList = lapply(DMVs_statTab, function(x) vector("list", nrow(x)))
for (j in 1:length(tmpList)) {
for(i in 1:nrow(DMVs_statTab[[j]])) {
  tmp =  do.call(rbind, lapply(DMVs_statTab, function(x) x[i,-1]))
  rownames(tmp) = names(DMVs_statTab)
  tmpList[[j]][[i]] = t(signif(tmp,3)) / dlpfcStats
  }
  names(tmpList[[j]]) = rownames(DMVs_statTab[[j]])
}


### relative enrichments

VsDlpfc = vector("list",length(tmpList))
VsDlpfc = lapply(VsDlpfc, function(x) vector("list", ncol(tmpList[[1]][[1]])))
for (i in 1:length(tmpList)) {
  for (j in 1:ncol(tmpList[[1]][[i]])) {
    VsDlpfc[[i]][[j]] = sapply(tmpList[[i]], function(x) x[,j])
  }
  names(VsDlpfc[[i]]) = colnames(tmpList[[i]][[j]])
}
names(VsDlpfc) = names(tmpList)


## general
gen.enr = lapply(VsDlpfc, function(x) lapply(x, function(y) rowSums(y > 2)))
gen.enr = lapply(gen.enr, function(x) do.call(rbind, x))
gen.depl = lapply(VsDlpfc, function(x) lapply(x, function(y) rowSums(y < 0.5)))
gen.depl = lapply(gen.depl, function(x) do.call(rbind, x))

## brain
brain.enr = lapply(VsDlpfc, function(x) lapply(x, function(y) rowMeans(y[,grep("Brain", dat$Standardized.Epigenome.name)] > 2)))
brain.enr = lapply(brain.enr, function(x) do.call(rbind, x))
brain.depl = lapply(VsDlpfc, function(x) lapply(x, function(y) rowMeans(y[,grep("Brain", dat$Standardized.Epigenome.name)] < 0.5)))
brain.depl = lapply(brain.depl, function(x) do.call(rbind, x))

## not brain
notbrain.enr = lapply(VsDlpfc, function(x) lapply(x, function(y) rowMeans(y[,-grep("Brain", dat$Standardized.Epigenome.name)] > 2)))
notbrain.enr = lapply(notbrain.enr, function(x) do.call(rbind, x))
notbrain.depl = lapply(VsDlpfc, function(x) lapply(x, function(y) rowMeans(y[,-grep("Brain", dat$Standardized.Epigenome.name)] < 0.5)))
notbrain.depl = lapply(notbrain.depl, function(x) do.call(rbind, x))

atriumStats = do.call(rbind, lapply(DMVs_statTab, function(x) x["E104",-1]))
rownames(atriumStats) = names(DMVs_statTab)
atriumStats = t(signif(atriumStats,3))

## other tissues
otherList = list(cellDMRs_statTab, ageDMRs_statTab, intDMRs_statTab)
ageVsCell = ageDMRs_statTab[,-1] / cellDMRs_statTab[,-1]
cor(t(ageVsCell),t(ageVsCell["E073",]))