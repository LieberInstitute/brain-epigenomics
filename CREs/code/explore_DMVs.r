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
num$celltype = ifelse(num$id %in% pd$Data.ID, pd[match(num$id, pd$Data.ID),"Cell.Type"], "Prenatal")
write.csv(num, quote=F, file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/DMV_number.csv")

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/CREs/figures/DMV_number_byAge_postnatal.pdf")

w = num[which(num$celltype %in% c("Neuron", "Glia")),]
w = cbind(w, Age = pd[match(w$id, pd$Data.ID), "Age"])
cor.test(x = w[which(w$celltype=="Neuron"),"Age"], y = w[which(w$celltype=="Neuron"),"num"])
#t = -0.82869, df = 22, p-value = 0.4162, cor -0.1739822 

cor.test(x = w[which(w$celltype=="Glia"),"Age"], y = w[which(w$celltype=="Glia"),"num"])
#t = -1.1274, df = 6, p-value = 0.3026, cor -0.4180941 

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
  ylab("Count") + xlab("") + ylim(0,5500) +
  ggtitle("Number of DMVs") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
dev.off()

num = read.csv("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/DMV_number.csv")
range(num$num[which(num$celltype=="Prenatal")]) # 1098 5305
range(num$num[which(num$celltype=="Neuron")]) # 1113 1410
range(num$num[which(num$celltype=="Glia")]) # 915 1144

t.test(num[num$celltype=="Prenatal","num"],num[num$celltype!="Prenatal","num"])
#t = 6.0393, df = 19.216, p-value = 7.875e-06
#  mean of x mean of y 
# 2856.850  1164.469 



## how long are these regions?

widths = do.call(rbind, Map(cbind, lapply(dmvs, function(y) data.frame(mean = mean(width(y)), median = median(width(y)), 
                                                                       sd = sd(width(y)), min = min(width(y)), max = max(width(y)))), id = as.list(names(dmvs))))
widths$celltype = ifelse(widths$id %in% pd$Data.ID, pd[match(widths$id, pd$Data.ID),"Cell.Type"], "Prenatal")
write.csv(widths, quote=F, file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/DMV_widths.csv")

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/CREs/figures/DMV_widths.pdf")

w = widths[which(widths$celltype %in% c("Neuron", "Glia")),]
w = cbind(w, Age = pd[match(w$id, pd$Data.ID), "Age"])
cor.test(x = w[which(w$celltype=="Neuron"),"Age"], y = w[which(w$celltype=="Neuron"),"median"])
#t = 0.25929, df = 22, p-value = 0.7978, cor 0.05519736

cor.test(x = w[which(w$celltype=="Glia"),"Age"], y = w[which(w$celltype=="Glia"),"median"])
#t = -0.99175, df = 6, p-value = 0.3596, cor -0.3752859 

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
  ggtitle("Median DMV Widths") + ylim(0,8500) +
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

t.test(widths[widths$celltype=="Prenatal","median"],widths[widths$celltype!="Prenatal","median"])
#t = 6.7297, df = 20.132, p-value = 1.459e-06
#  mean of x mean of y 
# 7239.625  6665.750 
t.test(widths[widths$celltype=="Prenatal","median"],widths[widths$celltype=="Neuron","median"])
#t = 6.8666, df = 19.939, p-value = 1.153e-06
#  mean of x mean of y 
#7239.625  6655.500 
t.test(widths[widths$celltype=="Glia","median"],widths[widths$celltype=="Neuron","median"])
#t = 0.92004, df = 8.387, p-value = 0.3833
#  mean of x mean of y 
#6696.5    6655.5 
t.test(widths[widths$celltype=="Glia","median"],widths[widths$celltype=="Prenatal","median"])
#t = -5.7655, df = 25.452, p-value = 4.887e-06
#  mean of x mean of y 
#6696.500  7239.625 


## What proportion of the genome do they cover?

load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/PMDs_methylSeekR.rda")
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/PMDs_methylSeekR_prenatal.rda")
sLengths = c(lapply(PMDsegments.CG, reduce), lapply(PMDsegments.CGpren, reduce))
rm(PMDsegments.CGpren, PMDsegments.CG)
percGenome = data.frame(percent = unlist(mapply(function(p, t) round((sum(width(p))/sum(as.numeric(width(t)))*100),2), dmvs, sLengths)),
                        id = names(dmvs))
percGenome$celltype = ifelse(percGenome$id %in% pd$Data.ID, pd[match(percGenome$id, pd$Data.ID),"Cell.Type"], "Prenatal")

write.csv(percGenome, quote=F, file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/DMV_genomeCoverage.csv")

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/CREs/figures/DMV_percentGenome.pdf")

w = percGenome[which(percGenome$celltype %in% c("Neuron", "Glia")),]
w = cbind(w, Age = pd[match(w$id, pd$Data.ID), "Age"])
cor.test(x = w[which(w$celltype=="Neuron"),"Age"], y = w[which(w$celltype=="Neuron"),"percent"])
#t = -1.2851, df = 22, p-value = 0.2121, cor -0.264245 

cor.test(x = w[which(w$celltype=="Glia"),"Age"], y = w[which(w$celltype=="Glia"),"percent"])
#t = -1.321, df = 6, p-value = 0.2346, cor -0.4746662 

ggplot(w, aes(x = Age, y = percent, colour = celltype)) + geom_point() +
  geom_smooth(method=lm, se=T, fullrange=TRUE) + scale_colour_brewer(8, palette="Dark2") + 
  labs(fill="") +  theme_classic() +
  ylab("Percent") + xlab("Age (Years)") +
  ggtitle("Percent of Genome that is DMV by Age") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20), legend.title=element_blank()) + theme(legend.position="bottom")
ggplot(percGenome, aes(x = celltype, y = percent)) + geom_boxplot() +
  theme_classic() +
  labs(fill="") + ylim(0,2) +
  ylab("Percent") + xlab("") +
  ggtitle("Percent Genome: DMV") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
dev.off()

genome = data.table(read.csv("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/DMV_genomeCoverage.csv"))
genome[,mean(percent), by = c("celltype")]
#   celltype        V1
#1:     Glia 0.2925000
#2:   Neuron 0.3408333
#3: Prenatal 0.9515000

t.test(genome[genome$celltype=="Prenatal","percent"],genome[genome$celltype!="Prenatal","percent"])
#t = 5.9192, df = 19.101, p-value = 1.044e-05
#  mean of x mean of y 
#0.95150   0.32875 
(0.95150-0.32875)/0.32875 # 1.894297

t.test(genome[genome$celltype=="Prenatal","percent"],genome[genome$celltype=="Neuron","percent"])
#t = 5.8072, df = 19.064, p-value = 1.338e-05
#  mean of x mean of y 
#0.9515000 0.3408333 
t.test(genome[genome$celltype=="Glia","percent"],genome[genome$celltype=="Neuron","percent"])
#t = -4.676, df = 10.093, p-value = 0.0008513
#  mean of x mean of y 
#0.2925000 0.3408333 
t.test(genome[genome$celltype=="Prenatal","percent"],genome[genome$celltype=="Glia","percent"])
#t = 6.2471, df = 19.302, p-value = 4.966e-06
#  mean of x mean of y 
#0.9515    0.2925 



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

round((unlist(lapply(all, function(x) sum(width(x))))/sum(width(dmaster))*100),2)
#      All  Prenatal Postnatal   Neurons      Glia 
#     1.49      2.88      1.97      2.85      2.79 
# 1.5-3% of total bases in the DMV state are shared by all or most

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


data.table(shared)[,mean(percent), by = c("group")]
#       group       V1
#1:       All 20.18327
#2:  Prenatal 22.17000
#3: Postnatal 34.00156
#4:   Neurons 47.20042
#5:      Glia 54.08875


## How much of prenatal DMVs is represented in each postnatal sample?

pren = reduce(do.call(getMethod(c, "GenomicRanges"), GRangesList(dmvs[-which(names(dmvs) %in% pd$Data.ID)])))

sharedP = lapply(dmvs, function(x) Reduce(intersect, list(x, pren)))

prenperc = mapply(function(pren,ind) round(sum(width(reduce(pren))) / sum(width(reduce(ind))) *100,2), 
                   sharedP, dmvs, SIMPLIFY = F)

prenperc = data.frame(perc = unlist(prenperc))
prenperc$Person = rownames(prenperc)
prenperc$Age = pd[match(prenperc$Person, pd$Data.ID),"Age"]
prenperc$CellType = ifelse(prenperc$Person %in% pd$Data.ID, pd[match(prenperc$Person, pd$Data.ID),"Cell.Type"], "Prenatal")

write.table(prenperc, quote = F, file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/DMV_sharedWithPrenatal.csv")

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/CREs/figures/percent_shared_withPrenatal_DMV_perSample_byBaseCoverage.pdf", height = 4)
ggplot(prenperc[prenperc$CellType!="Prenatal",], aes(x = CellType, y = perc)) + geom_boxplot() +
  labs(fill="") + theme_classic() +
  ylim(0,100) +
  ylab("Percent") + xlab("") +
  ggtitle("Percent Bases Shared with Prenatal") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
ggplot(prenperc[prenperc$CellType!="Prenatal",], aes(x = Age, y = perc, colour = CellType)) + 
  geom_path() + geom_point() + ylim(0,100) +
  theme_classic() + scale_colour_brewer(8, palette="Dark2") +
  ylab("Percent") + xlab("Age") + 
  ggtitle("Percent Bases Shared with Prenatal") + 
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20), legend.title=element_blank()) + theme(legend.position="bottom")
dev.off()

cor.test(x = prenperc[prenperc$CellType=="Neuron",]$Age, y = prenperc[prenperc$CellType=="Neuron",]$perc)
#t = -3.2655, df = 22, p-value = 0.00354, cor -0.5713691 
cor.test(x = prenperc[prenperc$CellType=="Glia",,]$Age, y = prenperc[prenperc$CellType=="Glia",,]$perc)
#t = -2.1369, df = 6, p-value = 0.07648, cor -0.6573872 


data.table(prenperc)[,mean(perc),by="CellType"]
#   CellType        V1
#1:     Glia  97.29125
#2:   Neuron  95.51083
#3: Prenatal 100.00000


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

load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/bumphunting/rda/limma_Neuron_CpGs_minCov_3_ageInfo_dmrs.Rdata")
dmrs = split(dmrs, dmrs$k6cluster_label)
names(dmrs) = c("Gr1","Gr2","Gr3","Gr4","Gr5","Gr6")

oo = lapply(dmvs, function(x) lapply(c(DMRgr[names(DMRgr) %in% c("CellType","Age")], as.list(dmrs), all, features), function(d) findOverlaps(x,d)))

dmvDMR = lapply(dmvs, as.data.frame)
dmvDMR = lapply(dmvDMR, function(x) data.frame(x, rnum = 1:nrow(x)))
dmvDMR = mapply(function(d,oo) data.frame(d, CT = ifelse(d$rnum %in% queryHits(oo$CellType), "CT", "no"),
                                          Age = ifelse(d$rnum %in% queryHits(oo$Age), "Age", "no"),
                                          Gr1 = ifelse(d$rnum %in% queryHits(oo$Gr1), "Gr1", "no"),
                                          Gr2 = ifelse(d$rnum %in% queryHits(oo$Gr2), "Gr2", "no"),
                                          Gr3 = ifelse(d$rnum %in% queryHits(oo$Gr3), "Gr3", "no"),
                                          Gr4 = ifelse(d$rnum %in% queryHits(oo$Gr4), "Gr4", "no"),
                                          Gr5 = ifelse(d$rnum %in% queryHits(oo$Gr5), "Gr5", "no"),
                                          Gr6 = ifelse(d$rnum %in% queryHits(oo$Gr6), "Gr6", "no"),
                                          All = ifelse(d$rnum %in% queryHits(oo$All), "All", "no"),
                                          Prenatal = ifelse(d$rnum %in% queryHits(oo$Prenatal), "Prenatal", "no"),
                                          Postnatal = ifelse(d$rnum %in% queryHits(oo$Postnatal), "Postnatal", "no"),
                                          Neurons = ifelse(d$rnum %in% queryHits(oo$Neurons), "Neurons", "no"),
                                          Glia = ifelse(d$rnum %in% queryHits(oo$Glia), "Glia", "no"),
                                          genes = ifelse(d$rnum %in% queryHits(oo$genes), "Gene", NA),
                                          islands = ifelse(d$rnum %in% queryHits(oo$islands), "CpG-Island", "non-Island"),
                                          promoters = ifelse(d$rnum %in% queryHits(oo$promoters), "promoter", "no")), dmvDMR, oo, SIMPLIFY = F) 
dmvDMR = lapply(dmvDMR, function(d) data.frame(d, tog = paste(d$CT, d$Age, d$Gr1,d$Gr2,d$Gr3,d$Gr4,d$Gr5,d$Gr6, d$All, d$Prenatal, 
                                                              d$Postnatal, d$Neurons, d$Glia, sep = ":"),
                                               dmr = paste(d$CT, d$Age, d$Gr1,d$Gr2,d$Gr3,d$Gr4,d$Gr5,d$Gr6, sep = ":"),
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

## how many fall within CpG islands?

CpGIslands = lapply(dmvDMR, function(y) round(length(unique(y[which(y$islands=="CpG-Island"),,]$regionID))/length(unique(y$regionID))*100,2))
CpGIslands = data.frame(islands = unlist(CpGIslands), id = names(CpGIslands))
CpGIslands$celltype = ifelse(CpGIslands$id %in% pd$Data.ID, pd[match(CpGIslands$id, pd$Data.ID),"Cell.Type"], "Prenatal")
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

cgi = data.table(read.csv("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/DMV_CpG_Island_Overlap.csv"))
cgi[,mean(islands), by = "celltype"]
#   celltype       V1
#1:     Glia 93.37250
#2:   Neuron 90.88583
#3: Prenatal 86.89900
mean(cgi$islands) # 89.735


# How about Repetitive Elements?

dmvDMRdt = lapply(dmvDMR, data.table)
repeats = lapply(dmvDMRdt, function(y) y[,length(unique(regionID)), by = "repeats"])
repeats = lapply(repeats, function(y) data.frame(y, perc = round(y$V1/sum(y$V1)*100,2)))
repeats = do.call(rbind, Map(cbind, repeats, id = as.list(names(repeats))))
repeats$repeats = factor(repeats$repeats, levels = c("SINE","LINE","Simple_repeat","DNA","LTR","Low_complexity","No repeats","Unknown","snRNA",
                                                     "Other","DNA?","Satellite","RC","scRNA","tRNA","SINE?","srpRNA","RNA","rRNA","LINE?","LTR?","Unknown?"))
repeats$celltype = ifelse(repeats$id %in% pd$Data.ID, pd[match(repeats$id, pd$Data.ID),"Cell.Type"], "Prenatal")
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

repeats = data.table(read.csv("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/DMV_repeats_Overlap.csv"))
repeats[repeats %in% c("No repeats","SINE","Low_complexity","Simple_repeat","LINE"),mean(perc), by = c("repeats","celltype")]
#           repeats celltype       V1
# 1: Low_complexity     Glia 27.05000
# 2:  Simple_repeat     Glia 25.42625
# 3:           SINE     Glia 19.77125
# 4:           LINE     Glia 14.25000
# 5:     No repeats     Glia  0.13125
# 6:  Simple_repeat   Neuron 25.06875
# 7: Low_complexity   Neuron 26.29000
# 8:           SINE   Neuron 20.25458
# 9:           LINE   Neuron 14.61542
#10:     No repeats   Neuron  0.08250
#11:           LINE Prenatal 17.10400
#12: Low_complexity Prenatal 21.68050
#13:  Simple_repeat Prenatal 20.32200
#14:           SINE Prenatal 21.17250
#15:     No repeats Prenatal  0.04850


## Overlap with DMRs

dmr = lapply(dmvDMRdt, function(y) y[,length(unique(regionID)), by = "dmr"])
dmr = lapply(dmr, function(y) data.frame(y, perc = round(y$V1/sum(y$V1)*100,2)))
models = lapply(dmr, function(y) data.frame(CT = sum(data.frame(y)[grep("CT", y$dmr),"perc"]), 
                                            Age = sum(data.frame(y)[grep("Age", y$dmr),"perc"]),
                                            Gr1 = sum(data.frame(y)[grep("Gr1", y$dmr),"perc"]),
                                            Gr2 = sum(data.frame(y)[grep("Gr2", y$dmr),"perc"]),
                                            Gr3 = sum(data.frame(y)[grep("Gr3", y$dmr),"perc"]),
                                            Gr4 = sum(data.frame(y)[grep("Gr4", y$dmr),"perc"]),
                                            Gr5 = sum(data.frame(y)[grep("Gr5", y$dmr),"perc"]),
                                            Gr6 = sum(data.frame(y)[grep("Gr6", y$dmr),"perc"])))
dmr = do.call(rbind, Map(cbind, dmr, id = as.list(names(dmr))))
models = do.call(rbind, Map(cbind, models, id = as.list(names(models))))

dmr$celltype = ifelse(dmr$id %in% pd$Data.ID, pd[match(dmr$id, pd$Data.ID),"Cell.Type"], "Prenatal")
models$celltype = ifelse(models$id %in% pd$Data.ID, pd[match(models$id, pd$Data.ID),"Cell.Type"], "Prenatal")

dmr$NO = ifelse(dmr$dmr=="no:no:no:no:no:no:no:no", "No overlap", "Overlap")

write.csv(dmr, quote=F,file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/DMV_dmr_Overlap.csv")
write.csv(models, quote=F,file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/DMV_dmr_Overlap_byModel.csv")


pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/CREs/figures/DMV_overlap_with_dmr.pdf", width = 14, height = 6)
ggplot(dmr, aes(x = celltype, y = perc, fill = NO)) + 
  theme_classic() + geom_bar(position = "fill",stat = "identity") +
  labs(fill="") + 
  ylab("Proportion") + 
  xlab("") +
  ggtitle("Proportion Overlapping DMRs") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20), legend.title=element_blank()) + theme(legend.position="bottom")
ggplot(dmr[which(dmr$dmr=="no:no:no:no:no:no:no:no"),], aes(x = celltype, y = perc)) + 
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

x = melt(models[,!colnames(models) %in% c("CT","Age")])
x$variable = gsub("Gr1","1:G-N+", x$variable)
x$variable = gsub("Gr2","2:G0N+", x$variable)
x$variable = gsub("Gr3","3:G0N-", x$variable)
x$variable = gsub("Gr4","4:G+N0", x$variable)
x$variable = gsub("Gr5","5:G+N-", x$variable)
x$variable = gsub("Gr6","6:G-N0", x$variable)

ggplot(x, aes(x = celltype, y = value, fill = variable)) + 
  theme_classic() + geom_boxplot() +
  labs(fill="") + facet_grid(. ~ variable) +
  ylab("Percent") +  
  xlab("") + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_brewer(8, palette="Dark2") +
  ggtitle("Percent Overlapping DMRs") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20), legend.title=element_blank()) + theme(legend.position="bottom")
dev.off()

dmr = data.table(read.csv("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/DMV_dmr_Overlap.csv"))
models = data.table(read.csv("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/DMV_dmr_Overlap_byModel.csv"))
models[,list(Gr1=mean(Gr1),Gr2=mean(Gr2),Gr3=mean(Gr3),Gr4=mean(Gr4),Gr5=mean(Gr5),Gr6=mean(Gr6)),by="celltype"]
#   celltype     Gr1       Gr2      Gr3      Gr4     Gr5     Gr6
#1:     Glia 2.71625 1.6100000 0.297500 1.258750 0.31750 1.72125
#2:   Neuron 1.13625 0.7054167 1.562083 1.418333 0.83375 0.32125
#3: Prenatal 1.54000 1.5705000 0.215000 1.360000 0.42350 0.51550


# Are DMRs overlapping DMVs at the DMV periphery, or uniformly distributed?

dmvDMR = Map(cbind, lapply(dmvDMR, function(x) data.frame(x[,1:5], dmr = x$dmr, regionID = x$regionID)), id = as.list(names(dmvDMR)))
dmvDMR = lapply(dmvDMR, unique)
tiles = lapply(dmvDMR, function(x) tile(makeGRangesFromDataFrame(x[which(x$dmr!="no:no:no:no:no:no:no:no"),]), n = 100))
tiles = lapply(tiles, as.list)
ti = list(Prenatal = unlist(tiles[-which(names(tiles) %in% postnatalpd$Data.ID)], recursive = F),
          Neuron = unlist(tiles[which(names(tiles) %in% postnatalpd[which(postnatalpd$Cell.Type=="Neuron"),"Data.ID"])], recursive = F),
          Glia = unlist(tiles[which(names(tiles) %in% postnatalpd[which(postnatalpd$Cell.Type=="Glia"),"Data.ID"])], recursive = F))
DMR = lapply(DMR, function(x) x[which(x$sig=="FWER < 0.05"),])
DMR = lapply(DMR, function(x) reduce(makeGRangesFromDataFrame(x)))
CTov = lapply(ti, function(x) lapply(x, function(y) findOverlaps(y,DMR$CellType)))
Aov = lapply(ti, function(x) lapply(x, function(y) findOverlaps(y,DMR$Age)))
Iov = lapply(ti, function(x) lapply(x, function(y) findOverlaps(y,DMR$Interaction)))

CToverlaps = Ageoverlaps = Intoverlaps = list(list(),list(),list())
for (i in 1:length(ti)) {
  for (j in 1:length(ti[[i]])) {
    if (length(queryHits(CTov[[i]][[j]]))>0) {
      CToverlaps[[i]][[j]] = data.frame(region = queryHits(CTov[[i]][[j]])) } else {
        CToverlaps[[i]][[j]] = data.frame(region = 0) }
    if (length(queryHits(Aov[[i]][[j]]))>0) {
      Ageoverlaps[[i]][[j]] = data.frame(region = queryHits(Aov[[i]][[j]])) } else {
        Ageoverlaps[[i]][[j]] = data.frame(region = 0) }
    if (length(queryHits(Iov[[i]][[j]]))>0) {
      Intoverlaps[[i]][[j]] = data.frame(region = queryHits(Iov[[i]][[j]])) } else {
        Intoverlaps[[i]][[j]] = data.frame(region = 0) } 
  }
}
CToverlaps = do.call(rbind, Map(cbind, lapply(CToverlaps, function(x) do.call(rbind, x)), celltype = as.list(names(ti))))
Ageoverlaps = do.call(rbind, Map(cbind, lapply(Ageoverlaps, function(x) do.call(rbind, x)), celltype = as.list(names(ti))))
Intoverlaps = do.call(rbind, Map(cbind, lapply(Intoverlaps, function(x) do.call(rbind, x)), celltype = as.list(names(ti))))
overlaps = rbind(cbind(CToverlaps, model = "Cell Type"), cbind(Ageoverlaps, model = "Age"),
                 cbind(Intoverlaps, model = "Interaction"))

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/CREs/figures/DMR_position_in_DMV.pdf", width = 10)
ggplot(overlaps, aes(region, colour = model)) + geom_density() + theme_classic() +
  labs(fill="") + facet_grid(. ~ celltype) +
  ylab("Density") + xlab("% of DMV") +
  ggtitle("DMR Position in DMV") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20), legend.position="bottom", legend.title=element_blank())
dev.off()


## association of existing DMVs with genes

oo = lapply(dmvs, function(x) findOverlaps(makeGRangesFromDataFrame(geneMap), makeGRangesFromDataFrame(x))) 
dmvGene = mapply(function(p, ov) if (nrow(p)>length(unique(subjectHits(ov)))) {
  rbind(data.frame(p[subjectHits(ov),], nearestSymbol = geneMap[queryHits(ov),"Symbol"],
                   nearestID = geneMap[queryHits(ov),"gencodeID"], 
                   EntrezID = geneMap[queryHits(ov),"EntrezID"]),
        data.frame(p[-unique(subjectHits(ov)),], nearestSymbol = "NoGeneOverlap", 
                   nearestID = "NoGeneOverlap", EntrezID = "NoGeneOverlap")) } else {
                     data.frame(p[subjectHits(ov),], nearestSymbol = geneMap[queryHits(ov),"Symbol"],
                                nearestID = geneMap[queryHits(ov),"gencodeID"], EntrezID = geneMap[queryHits(ov),"EntrezID"])      
                   },lapply(dmvs, as.data.frame), oo, SIMPLIFY = F)
save(dmvGene, postnatalpd, geneMap, file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/DMVs_geneOverlap.rda")

dmvDMRdt = lapply(dmvGene, data.table)
dmvDMRdt = Map(cbind, dmvDMRdt, regionID = lapply(dmvDMRdt, function(x) paste0(x$seqnames,":",x$start,"-",x$end)))
geneoverlap = lapply(dmvDMRdt, function(y) round(length(unique(as.character(data.frame(y)[which(y$nearestID!="NoGeneOverlap"),"regionID"])))/
                                                 length(unique(as.character(y$regionID)))*100,2))
geneoverlap = data.frame(geneoverlap = unlist(geneoverlap), id = names(geneoverlap))
geneoverlap$celltype = ifelse(geneoverlap$id %in% pd$Data.ID, pd[match(geneoverlap$id, pd$Data.ID),"Cell.Type"], "Prenatal")
write.csv(geneoverlap, quote=F,file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/DMV_Gene_Overlap.csv")

numgenes = lapply(dmvDMRdt, function(y) y[,length(unique(nearestID)), by = "regionID"])
numgenes = do.call(rbind, Map(cbind, lapply(numgenes, function(y) 
           data.frame(mean = mean(y$V1), median = median(y$V1), sd = sd(y$V1), min = min(y$V1), max = max(y$V1))), id = as.list(names(numgenes))))
numgenes$celltype = ifelse(numgenes$id %in% pd$Data.ID, pd[match(numgenes$id, pd$Data.ID),"Cell.Type"], "Prenatal")
write.csv(numgenes, quote=F, file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/DMV_numgenes.csv")

geneoverlap = read.csv("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/DMV_Gene_Overlap.csv")
mean(geneoverlap$geneoverlap)

w = numgenes[which(numgenes$celltype %in% c("Neuron", "Glia")),]
w = cbind(w, Age = postnatalpd[match(w$id, postnatalpd$Data.ID), "Age"])
cor.test(x = w[which(w$celltype=="Neuron"),"Age"], y = w[which(w$celltype=="Neuron"),"mean"])
#t = -0.86357, df = 22, p-value = 0.3971, cor -0.1810702 
cor.test(x = w[which(w$celltype=="Glia"),"Age"], y = w[which(w$celltype=="Glia"),"mean"])
#t = -1.0309, df = 6, p-value = 0.3423, cor -0.3879235 

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
  ylab("Number") + xlab("Age (Years)") +
  ggtitle("Median DMV Gene Number by Age") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20), legend.title=element_blank()) + theme(legend.position="bottom")
ggplot(numgenes, aes(x = celltype, y = median)) + geom_boxplot() +
  theme_classic() +
  labs(fill="") +
  ylab("Number") + xlab("") +
  ggtitle("Median DMV Gene Number") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
ggplot(numgenes, aes(x = celltype, y = mean)) + geom_boxplot() +
  theme_classic() +
  labs(fill="") +
  ylab("Number") + xlab("") +
  ggtitle("Mean DMV Gene Number") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
dev.off()


library(GenomicFeatures)
txdb = loadDb("/dcl01/lieber/ajaffe/Amanda/annotation_objects/gencode.v25lift37.annotation.sqlite")
promoters = promoters(txdb, upstream=2000, downstream=200)
promoters = reduce(promoters)
oo = lapply(dmvGene, function(x) findOverlaps(promoters, makeGRangesFromDataFrame(x))) 
for (i in 1:length(dmvGene)) {
  dmvGene[[i]][unique(subjectHits(oo[[i]])),"promoters"] = "promoters"
  dmvGene[[i]][-unique(subjectHits(oo[[i]])),"promoters"] = "no"
}

dmvGenedt = lapply(dmvGene, data.table)
dmvGenedt = Map(cbind, dmvGenedt, regionID = lapply(dmvGenedt, function(x) paste0(x$seqnames,":",x$start,"-",x$end)))
promoverlap = lapply(dmvGenedt, function(y) round(length(unique(as.character(data.frame(y)[which(y$promoters=="promoters"),"regionID"])))/
                                                 length(unique(as.character(y$regionID)))*100,2))
promoverlap = data.frame(promoverlap = unlist(promoverlap), id = names(promoverlap))
promoverlap$celltype = ifelse(promoverlap$id %in% pd$Data.ID, pd[match(promoverlap$id, pd$Data.ID),"Cell.Type"], "Prenatal")
write.csv(promoverlap, quote=F,file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/DMV_prom_Overlap.csv")

numproms = lapply(dmvGenedt, function(y) y[promoters=="promoters",length(unique(nearestID)), by = "regionID"])
numproms = do.call(rbind, Map(cbind, lapply(numproms, function(y) 
  data.frame(mean = mean(y$V1), median = median(y$V1), sd = sd(y$V1), min = min(y$V1), max = max(y$V1))), id = as.list(names(numproms))))
numproms$celltype = ifelse(numproms$id %in% pd$Data.ID, pd[match(numproms$id, pd$Data.ID),"Cell.Type"], "Prenatal")
write.csv(numproms, quote=F, file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/DMV_numproms.csv")

mean(promoverlap$promoverlap) # 92.38346

w = numproms[which(numproms$celltype %in% c("Neuron", "Glia")),]
w = cbind(w, Age = postnatalpd[match(w$id, postnatalpd$Data.ID), "Age"])
cor.test(x = w[which(w$celltype=="Neuron"),"Age"], y = w[which(w$celltype=="Neuron"),"mean"])
#t = -1.2137, df = 22, p-value = 0.2377, cor -0.2505089 
cor.test(x = w[which(w$celltype=="Glia"),"Age"], y = w[which(w$celltype=="Glia"),"mean"])
#t = -1.1737, df = 6, p-value = 0.285, cor -0.4321072 

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
  ylab("Number") + xlab("Age (Years)") +
  ggtitle("Median DMV Promoter Number by Age") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20), legend.title=element_blank()) + theme(legend.position="bottom")
ggplot(numproms, aes(x = celltype, y = median)) + geom_boxplot() +
  theme_classic() +
  labs(fill="") +
  ylab("Number") + xlab("") +
  ggtitle("Median DMV Promoter Number") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
ggplot(numproms, aes(x = celltype, y = mean)) + geom_boxplot() +
  theme_classic() +
  labs(fill="") +
  ylab("Number") + xlab("") +
  ggtitle("Mean DMV Promoter Number") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
dev.off()


## transcription factor genes in our DMVs?

library(PWMEnrich.Hsapiens.background)
data(MotifDb.Hsap)
tfs = unlist(lapply(MotifDb.Hsap, function(x) x$name))
tfs = tfs[-grep("UW",tfs)]

not = tfs[!tfs %in% geneMap$Symbol]
not = na.omit(not[-grep("::", not, fixed=T)])

nomatch = c("SMCR7L","DUX4","MYF","MZF1_1-4","MZF1_5-13","RORA_1","RORA_2","MIZF","EWSR1-FLI1","ZNF238",
            "ZNF306","POU5F1P1","BHLHB2","BHLHB3","CART1","RAXL1","TRP53","TRP73","ZNF435","HNF1B")
nomatch[which(nomatch %in% geneMap$Symbol)]      
nomatch[which(nomatch %in% good)]      

tfs = unique(c(tfs[tfs %in% geneMap$Symbol], c("HINFP1","RBM8B"), unique(unlist(strsplit(not[grep("::", not, fixed=T)], "::", fixed=T)))))
tfs = geneMap[which(geneMap$Symbol %in% tfs),] # 830 genes

oo = lapply(dmvs, function(p) findOverlaps(p,makeGRangesFromDataFrame(tfs)))
hits = lapply(oo, function(x) data.frame(tfs = length(unique(subjectHits(x))), dmvs = length(unique(queryHits(x)))))
hits = do.call(rbind, Map(cbind, hits, id = as.list(names(dmvs)), num.dmvs = as.list(elementNROWS(dmvs))))
hits$celltype = ifelse(hits$id %in% pd$Data.ID, pd[match(hits$id, pd$Data.ID),"Cell.Type"], "Prenatal")
hits$DMV.perc = round(hits$dmvs / hits$num.dmvs * 100,1)
hits$TF.perc = round(hits$tfs / nrow(tfs) * 100,1)
write.csv(hits, quote=F,file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/DMV_tfsgenes_Overlap.csv")

hits  = read.csv("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/DMV_tfsgenes_Overlap.csv")
hits = data.table(hits)
hits[,median(tfs),by="celltype"]
#   celltype    V1
#1:     Glia 176.5
#2:   Neuron 183.0
#3: Prenatal 304.0


pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/CREs/figures/DMV_tfsgenes_Overlap.pdf")
ggplot(reshape2::melt(hits)[grep("perc", reshape2::melt(hits)$variable),], aes(x = celltype, y = value)) + geom_boxplot() +
  theme_classic() + facet_grid( ~ variable) +
  labs(fill="") +
  ylab("Percent") + 
  xlab("") + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Percent Genes and DMVs Overlapping") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20), legend.title=element_blank()) + theme(legend.position="bottom")
dev.off()


## Enrichment in transcription factor genes

geneuniverse = na.omit(unique(geneMap$gencodeID))
sets = na.omit(unique(tfs$gencodeID))
inDMV = lapply(dmvGene, function(x) na.omit(unique(x$nearestID)))
outDMV = lapply(inDMV, function(x) geneuniverse[!(geneuniverse %in% x)])

DMVenrich = mapply(function(yes, no) {
  DMV_OVERLAP = c( sum( yes %in% sets),sum(!(yes %in% sets)))
  NOT_DMV_OVERLAP= c(sum(no %in% sets), sum(!(no %in% sets)))
  enrich_table = cbind(DMV_OVERLAP, NOT_DMV_OVERLAP)
  res = fisher.test(enrich_table)
  dat=c(res$p.value, res$estimate)
  names(dat) <- c("P.Value","Odds.Ratio")
  return(dat)
}, inDMV, outDMV, SIMPLIFY =F)
DMVenrich = do.call(rbind, Map(cbind, lapply(DMVenrich, function(y) data.frame(P.Value = y["P.Value"], Odds.Ratio = y["Odds.Ratio"])),
                               id = as.list(names(DMVenrich))))
DMVenrich$celltype = ifelse(DMVenrich$id %in% pd$Data.ID, pd[match(DMVenrich$id, pd$Data.ID),"Cell.Type"], "Prenatal")
DMVenrich$FDR = p.adjust(DMVenrich$P.Value, method = "fdr")
DMVenrich$celltype = factor(DMVenrich$celltype, levels = c("Prenatal", "Neuron", "Glia"))
write.csv(DMVenrich,file="/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/DMV_tfsgenes_enrichment.csv",quote=F)

DMVenrich = read.csv("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/DMV_tfsgenes_enrichment.csv")

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/CREs/figures/DMV_tfsgenes_Enrichment.pdf", height = 4, width = 6)
ggplot(DMVenrich[DMVenrich$FDR<=0.05,], aes(x = celltype, y = Odds.Ratio)) + geom_boxplot() +
  theme_classic() +
  labs(fill="") +
  ylab("Odds Ratio") + ylim(0,15) +
  geom_hline(yintercept=1, linetype="dotted") +
  xlab("") + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Fisher Test for TF Genes in DMVs (FDR<0.05)") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20), legend.title=element_blank()) + theme(legend.position="bottom")
dev.off()



## Identify the regions that are getting smaller or larger over age

load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/dmvs_annotated.rda")
dmvDMRdt = Map(cbind, lapply(dmvDMR, data.table), tog2 = lapply(dmvDMR, function(d) paste(d$All, d$Prenatal, d$Postnatal, d$Neurons, d$Glia, sep = ":")))
lapply(dmvDMRdt, function(x) x[,length(unique(regionID)),by="tog2"])

subset = do.call(rbind, Map(cbind, lapply(dmvDMRdt, function(x) data.frame(regionID = x$regionID, width = x$width, tog2 = x$tog2)), id = as.list(names(dmvDMRdt))))
subset$celltype = ifelse(subset$id %in% pd$Data.ID, pd[match(subset$id, pd$Data.ID),"Cell.Type"], "Prenatal")
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


## Do same comparison but of shared or absent genes in DMVs

inDMV = lapply(inDMV, as.character)
inDMV.byCT = list(Prenatal = unique(unlist(inDMV[-which(names(inDMV) %in% pd$Data.ID)])),
                  Neuron = unique(unlist(inDMV[which(names(inDMV) %in% pd[pd$Cell.Type=="Neuron","Data.ID"])])),
                  Glia = unique(unlist(inDMV[which(names(inDMV) %in% pd[pd$Cell.Type=="Glia","Data.ID"])])))

elementNROWS(inDMV.byCT)
# Prenatal   Neuron     Glia 
#    16320     3556     2735 

CTcomps = list(PnotN = inDMV.byCT$Prenatal[-which(inDMV.byCT$Prenatal %in% inDMV.byCT$Neuron)],
               NnotP = inDMV.byCT$Neuron[-which(inDMV.byCT$Neuron %in% inDMV.byCT$Prenatal)],
               sharedPN = inDMV.byCT$Prenatal[which(inDMV.byCT$Prenatal %in% inDMV.byCT$Neuron)],
               PnotG = inDMV.byCT$Prenatal[-which(inDMV.byCT$Prenatal %in% inDMV.byCT$Glia)],
               GnotP = inDMV.byCT$Glia[-which(inDMV.byCT$Glia %in% inDMV.byCT$Prenatal)],
               sharedPG = inDMV.byCT$Prenatal[which(inDMV.byCT$Prenatal %in% inDMV.byCT$Glia)],
               GnotN = inDMV.byCT$Glia[-which(inDMV.byCT$Glia %in% inDMV.byCT$Neuron)],
               NnotG = inDMV.byCT$Neuron[-which(inDMV.byCT$Neuron %in% inDMV.byCT$Glia)],
               sharedGN = inDMV.byCT$Glia[which(inDMV.byCT$Glia %in% inDMV.byCT$Neuron)],
               allshared = Reduce(intersect, inDMV.byCT))
elementNROWS(CTcomps)
#PnotN     NnotP  sharedPN     PnotG     GnotP  sharedPG     GnotN     NnotG 
#12940       176      3380     13635        50      2685       618      1439 
#sharedGN allshared 
#    2117      2105 


inDMV.byAgeN = list(infant = unique(unlist(inDMV[which(names(inDMV) %in% pd[pd$Cell.Type=="Neuron" & pd$Age<1,"Data.ID"])])),
                   child = unique(unlist(inDMV[which(names(inDMV) %in% pd[pd$Cell.Type=="Neuron" & pd$Age>1 & pd$Age<=12,"Data.ID"])])), 
                   teen = unique(unlist(inDMV[which(names(inDMV) %in% pd[pd$Cell.Type=="Neuron" & pd$Age>12 & pd$Age<=17,"Data.ID"])])),
                   adult = unique(unlist(inDMV[which(names(inDMV) %in% pd[pd$Cell.Type=="Neuron" & pd$Age>17,"Data.ID"])])))
elementNROWS(inDMV.byAgeN)
#infant  child   teen  adult 
#  2576   2782   2825   2906 

AgecompsN = list(InotC = inDMV.byAgeN$infant[-which(inDMV.byAgeN$infant %in% inDMV.byAgeN$child)],
                CnotI = inDMV.byAgeN$child[-which(inDMV.byAgeN$child %in% inDMV.byAgeN$infant)],
                sharedIC = inDMV.byAgeN$infant[which(inDMV.byAgeN$infant %in% inDMV.byAgeN$child)],
                InotT = inDMV.byAgeN$infant[-which(inDMV.byAgeN$infant %in% inDMV.byAgeN$teen)],
                TnotI = inDMV.byAgeN$teen[-which(inDMV.byAgeN$teen %in% inDMV.byAgeN$infant)],
                sharedIT = inDMV.byAgeN$infant[which(inDMV.byAgeN$infant %in% inDMV.byAgeN$teen)],
                InotA = inDMV.byAgeN$infant[-which(inDMV.byAgeN$infant %in% inDMV.byAgeN$adult)],
                AnotI = inDMV.byAgeN$adult[-which(inDMV.byAgeN$adult %in% inDMV.byAgeN$infant)],
                sharedIA = inDMV.byAgeN$infant[which(inDMV.byAgeN$infant %in% inDMV.byAgeN$adult)],
                CnotT = inDMV.byAgeN$child[-which(inDMV.byAgeN$child %in% inDMV.byAgeN$teen)],
                TnotC = inDMV.byAgeN$teen[-which(inDMV.byAgeN$teen %in% inDMV.byAgeN$child)],
                sharedCT = inDMV.byAgeN$child[which(inDMV.byAgeN$child %in% inDMV.byAgeN$teen)],
                CnotA = inDMV.byAgeN$child[-which(inDMV.byAgeN$child %in% inDMV.byAgeN$adult)],
                AnotC = inDMV.byAgeN$adult[-which(inDMV.byAgeN$adult %in% inDMV.byAgeN$child)],
                sharedCA = inDMV.byAgeN$child[which(inDMV.byAgeN$child %in% inDMV.byAgeN$adult)],
                TnotA = inDMV.byAgeN$teen[-which(inDMV.byAgeN$teen %in% inDMV.byAgeN$adult)],
                AnotT = inDMV.byAgeN$adult[-which(inDMV.byAgeN$adult %in% inDMV.byAgeN$teen)],
                sharedTA = inDMV.byAgeN$teen[which(inDMV.byAgeN$teen %in% inDMV.byAgeN$adult)])
elementNROWS(AgecompsN)
#   InotC    CnotI sharedIC    InotT    TnotI sharedIT    InotA    AnotI 
#     373      579     2203      368      617     2208      356      686 
#sharedIA    CnotT    TnotC sharedCT    CnotA    AnotC sharedCA    TnotA 
#    2220      285      328     2497      284      408     2498      234 
#AnotT sharedTA 
#  315     2591 

inDMV.byAgeG = list(infant = unique(unlist(inDMV[which(names(inDMV) %in% pd[pd$Cell.Type=="Glia" & pd$Age<1,"Data.ID"])])),
                    child = unique(unlist(inDMV[which(names(inDMV) %in% pd[pd$Cell.Type=="Glia" & pd$Age>1 & pd$Age<=12,"Data.ID"])])), 
                    teen = unique(unlist(inDMV[which(names(inDMV) %in% pd[pd$Cell.Type=="Glia" & pd$Age>12 & pd$Age<=17,"Data.ID"])])),
                    adult = unique(unlist(inDMV[which(names(inDMV) %in% pd[pd$Cell.Type=="Glia" & pd$Age>17,"Data.ID"])])))
elementNROWS(inDMV.byAgeG)
#infant  child   teen  adult 
#  1936   2175   1617   2036 

AgecompsG = list(InotC = inDMV.byAgeG$infant[-which(inDMV.byAgeG$infant %in% inDMV.byAgeG$child)],
                 CnotI = inDMV.byAgeG$child[-which(inDMV.byAgeG$child %in% inDMV.byAgeG$infant)],
                 sharedIC = inDMV.byAgeG$infant[which(inDMV.byAgeG$infant %in% inDMV.byAgeG$child)],
                 InotT = inDMV.byAgeG$infant[-which(inDMV.byAgeG$infant %in% inDMV.byAgeG$teen)],
                 TnotI = inDMV.byAgeG$teen[-which(inDMV.byAgeG$teen %in% inDMV.byAgeG$infant)],
                 sharedIT = inDMV.byAgeG$infant[which(inDMV.byAgeG$infant %in% inDMV.byAgeG$teen)],
                 InotA = inDMV.byAgeG$infant[-which(inDMV.byAgeG$infant %in% inDMV.byAgeG$adult)],
                 AnotI = inDMV.byAgeG$adult[-which(inDMV.byAgeG$adult %in% inDMV.byAgeG$infant)],
                 sharedIA = inDMV.byAgeG$infant[which(inDMV.byAgeG$infant %in% inDMV.byAgeG$adult)],
                 CnotT = inDMV.byAgeG$child[-which(inDMV.byAgeG$child %in% inDMV.byAgeG$teen)],
                 TnotC = inDMV.byAgeG$teen[-which(inDMV.byAgeG$teen %in% inDMV.byAgeG$child)],
                 sharedCT = inDMV.byAgeG$child[which(inDMV.byAgeG$child %in% inDMV.byAgeG$teen)],
                 CnotA = inDMV.byAgeG$child[-which(inDMV.byAgeG$child %in% inDMV.byAgeG$adult)],
                 AnotC = inDMV.byAgeG$adult[-which(inDMV.byAgeG$adult %in% inDMV.byAgeG$child)],
                 sharedCA = inDMV.byAgeG$child[which(inDMV.byAgeG$child %in% inDMV.byAgeG$adult)],
                 TnotA = inDMV.byAgeG$teen[-which(inDMV.byAgeG$teen %in% inDMV.byAgeG$adult)],
                 AnotT = inDMV.byAgeG$adult[-which(inDMV.byAgeG$adult %in% inDMV.byAgeG$teen)],
                 sharedTA = inDMV.byAgeG$teen[which(inDMV.byAgeG$teen %in% inDMV.byAgeG$adult)])
elementNROWS(AgecompsG)
# InotC    CnotI sharedIC    InotT    TnotI sharedIT    InotA    AnotI 
#   300      539     1636      605      286     1331      438      538 
# sharedIA    CnotT    TnotC sharedCT    CnotA    AnotC sharedCA    TnotA 
#     1498      695      137     1480      408      269     1767      203 
# AnotT sharedTA 
#   622     1414 



## Gene ontology of genes in different categories

entrez = list(CellType = lapply(CTcomps, function(x) na.omit(unique(geneMap[which(geneMap$gencodeID %in% x),"EntrezID"]))),
              Age = lapply(AgecompsN, function(x) na.omit(unique(geneMap[which(geneMap$gencodeID %in% x),"EntrezID"]))),
              Age.Glia = lapply(AgecompsG, function(x) na.omit(unique(geneMap[which(geneMap$gencodeID %in% x),"EntrezID"]))))


# Compare the enriched terms between 7 groups

compareKegg = lapply(entrez, function(x) compareCluster(x, fun="enrichKEGG", qvalueCutoff = 0.05, pvalueCutoff = 0.05))
compareBP = lapply(entrez, function(x) compareCluster(x, fun="enrichGO", ont = "BP", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05))
compareMF = lapply(entrez, function(x) compareCluster(x, fun="enrichGO",  ont = "MF", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05))
compareCC = lapply(entrez, function(x) compareCluster(x, fun="enrichGO",  ont = "CC", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05))
compareDO = lapply(entrez, function(x) compareCluster(x, fun="enrichDO",  ont = "DO", qvalueCutoff = 0.05, pvalueCutoff = 0.05))

# save object
save(compareKegg, compareBP, compareMF, compareCC, compareDO, file="/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/DMV_KEGG_GO_DO_objects.rda")

# plot compared results
pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/CREs/figures/DMV_KEGG_GO_DO_plots.pdf", height = 160, width = 22)
for (i in 1:length(compareKegg)) {
  print(plot(compareKegg[[i]], colorBy="p.adjust", showCategory = 1500, title= paste0("KEGG Pathway Enrichment: ", names(compareKegg)[i])))
  print(plot(compareBP[[i]], colorBy="p.adjust", showCategory = 1500, title= paste0("Biological Process GO Enrichment: ", names(compareKegg)[i])))
  print(plot(compareMF[[i]], colorBy="p.adjust", showCategory = 1500, title= paste0("Molecular Function GO Enrichment: ", names(compareKegg)[i])))
  print(plot(compareCC[[i]], colorBy="p.adjust", showCategory = 1500, title= paste0("Cellular Compartment GO Enrichment: ", names(compareKegg)[i])))
  print(plot(compareDO[[i]], colorBy="p.adjust", showCategory = 1500, title= paste0("Disease Ontology Enrichment: ", names(compareKegg)[i])))
}
dev.off()


## Plot select terms

plotCellType = compareBP$CellType # clusterProfiler output
plotCellType@compareClusterResult = plotCellType@compareClusterResult[-grep("shared", as.character(plotCellType@compareClusterResult$Cluster)),]
plotAge = compareBP$Age # clusterProfiler output
plotAge@compareClusterResult = plotAge@compareClusterResult[-grep("shared", as.character(plotAge@compareClusterResult$Cluster)),]
plotAge.Glia = compareBP$Age.Glia # clusterProfiler output
plotAge.Glia@compareClusterResult = plotAge.Glia@compareClusterResult[-grep("shared", as.character(plotAge.Glia@compareClusterResult$Cluster)),]

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/CREs/figures/DMV_BP_noshared.pdf", height = 70, width=16)
plot(plotCellType, colorBy="p.adjust", showCategory = 450, title= "Biological Process Enrichment: Cell type")
plot(plotAge, colorBy="p.adjust", showCategory = 450, title= "Biological Process Enrichment: Development")
dev.off()

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/CREs/figures/DMV_BP_noshared_top10.pdf", height = 16, width=12)
plot(plotCellType, colorBy="p.adjust", showCategory = 10, title= "Biological Process Enrichment: Cell type")
plot(plotAge, colorBy="p.adjust", showCategory = 10, title= "Biological Process Enrichment: Neuronal Development")
plot(plotAge.Glia, colorBy="p.adjust", showCategory = 10, title= "Biological Process Enrichment: Glial Development")
dev.off()

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/CREs/figures/DMV_BP_noshared_top5.pdf", height = 8, width=12)
plot(plotCellType, colorBy="p.adjust", showCategory = 5, title= "Biological Process Enrichment: Cell type")
plot(plotAge, colorBy="p.adjust", showCategory = 5, title= "Biological Process Enrichment: Neuronal Development")
plot(plotAge.Glia, colorBy="p.adjust", showCategory = 5, title= "Biological Process Enrichment: Glial Development")
dev.off()


DMV.CTcomps = lapply(CTcomps, function(x) geneMap[match(x, geneMap$gencodeID),])
DMV.Agecomps = lapply(Agecomps, function(x) geneMap[match(x, geneMap$gencodeID),])
DMV.AgecompsG = lapply(AgecompsG, function(x) geneMap[match(x, geneMap$gencodeID),])


## Plot venn diagram of genes

library(VennDiagram)
venn.diagram(inDMV.byCT, "/dcl01/lieber/ajaffe/lab/brain-epigenomics/CREs/figures/venn_diagram_DMVgenes_byCellType.jpeg", 
             main="Genes within DMVs by Cell Type",
             col = "transparent",
             fill = c("lightpink2","cornflowerblue", "olivedrab2"),
             cat.col = c("palevioletred4", "darkblue", "olivedrab4"),
             alpha = 0.50, fontfamily = "Arial", fontface = "bold", cat.fontfamily = "Arial", margin=0.2)

venn.diagram(inDMV.byAgeN, "/dcl01/lieber/ajaffe/lab/brain-epigenomics/CREs/figures/venn_diagram_DMVgenes_byAge_inNeurons.jpeg", 
             main="Genes within DMVs by Age in Neurons",
             col = "transparent",
             fill = c("lightpink2","cornflowerblue", "olivedrab2", "khaki1"),
             cat.col = c("palevioletred4", "darkblue", "olivedrab4", "lightgoldenrod4"),
             alpha = 0.50, fontfamily = "Arial", fontface = "bold", cat.fontfamily = "Arial", margin=0.2)

venn.diagram(inDMV.byAgeG, "/dcl01/lieber/ajaffe/lab/brain-epigenomics/CREs/figures/venn_diagram_DMVgenes_byAge_inGlia.jpeg", 
             main="Genes within DMVs by Age in Glia",
             col = "transparent",
             fill = c("lightpink2","cornflowerblue", "olivedrab2", "khaki1"),
             cat.col = c("palevioletred4", "darkblue", "olivedrab4", "lightgoldenrod4"),
             alpha = 0.50, fontfamily = "Arial", fontface = "bold", cat.fontfamily = "Arial", margin=0.2)


## Associate with gene expression

load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/sorted_nuclear_RNA/DE_limma_results_objects.rda")
DMV.CTcomps = Map(cbind, DMV.CTcomps, lapply(DMV.CTcomps, function(x) nucRNAres[match(x$gencodeID, nucRNAres$gencodeID),]))
elementNROWS(DMV.CTcomps)
elementNROWS(lapply(DMV.CTcomps, function(x) unique(x$gencodeID)))
tstat.CT = list(PnotN = t.test(DMV.CTcomps$PnotN$Tstat.CellTypeNeuron, DMV.CTcomps$NnotP$Tstat.CellTypeNeuron),
                NnotP = t.test(DMV.CTcomps$PnotN$Tstat.CellTypeNeuron, DMV.CTcomps$NnotP$Tstat.CellTypeNeuron),
                PnotG = t.test(DMV.CTcomps$PnotG$Tstat.CellTypeNeuron, DMV.CTcomps$GnotP$Tstat.CellTypeNeuron),
                GnotP = t.test(DMV.CTcomps$PnotG$Tstat.CellTypeNeuron, DMV.CTcomps$GnotP$Tstat.CellTypeNeuron),
                GnotN = t.test(DMV.CTcomps$GnotN$Tstat.CellTypeNeuron, DMV.CTcomps$NnotG$Tstat.CellTypeNeuron),
                NnotG = t.test(DMV.CTcomps$GnotN$Tstat.CellTypeNeuron, DMV.CTcomps$NnotG$Tstat.CellTypeNeuron))
ct = data.frame(Tstat = unlist(lapply(tstat.CT, function(x) x$statistic)), mean1 = unlist(lapply(tstat.CT, function(x) x$estimate[1])),
                mean2 = unlist(lapply(tstat.CT, function(x) x$estimate[2])), pval = unlist(lapply(tstat.CT, function(x) x$p.value)), comps = names(tstat.CT), row.names = NULL)
# negative Tstat.CellTypeNeuron means higher expressed in glia
#       Tstat      mean1     mean2         pval       comps
#1  -6.889639  0.4861798  2.852575 1.274992e-10 PnotN.NnotP
#2   5.871658  0.8119004 -2.901443 4.781105e-07 PnotG.GnotP
#3 -22.430365 -2.0202821  2.841388 6.231489e-90 GnotN.NnotG
# Genes that are escaping the DMV state (likely accumulating DNAm) are higher expressed in the cell type in which the gene escapes

load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/homogenate_RNA/DE_limma_results_homogenateRNAseq.rda")
DMV.Agecomps = Map(cbind, DMV.Agecomps, lapply(DMV.Agecomps, function(x) postRNAres[match(x$gencodeID, postRNAres$gencodeID),]))
elementNROWS(DMV.Agecomps)
elementNROWS(lapply(DMV.Agecomps, function(x) unique(x$gencodeID)))
tstat = list(InotC = t.test(DMV.Agecomps$InotC$Tstat, DMV.Agecomps$CnotI$Tstat),
             CnotI = t.test(DMV.Agecomps$InotC$Tstat, DMV.Agecomps$CnotI$Tstat),
             InotT = t.test(DMV.Agecomps$InotT$Tstat, DMV.Agecomps$TnotI$Tstat),
             TnotI = t.test(DMV.Agecomps$InotT$Tstat, DMV.Agecomps$TnotI$Tstat),
             InotA = t.test(DMV.Agecomps$InotA$Tstat, DMV.Agecomps$AnotI$Tstat),
             AnotI = t.test(DMV.Agecomps$InotA$Tstat, DMV.Agecomps$AnotI$Tstat),
             CnotT = t.test(DMV.Agecomps$CnotT$Tstat, DMV.Agecomps$TnotC$Tstat),
             TnotC = t.test(DMV.Agecomps$CnotT$Tstat, DMV.Agecomps$TnotC$Tstat),
             CnotA = t.test(DMV.Agecomps$CnotA$Tstat, DMV.Agecomps$AnotC$Tstat),
             AnotC = t.test(DMV.Agecomps$CnotA$Tstat, DMV.Agecomps$AnotC$Tstat),
             TnotA = t.test(DMV.Agecomps$TnotA$Tstat, DMV.Agecomps$AnotT$Tstat),
             AnotT = t.test(DMV.Agecomps$TnotA$Tstat, DMV.Agecomps$AnotT$Tstat))
age = data.frame(Tstat = unlist(lapply(tstat, function(x) x$statistic)), mean1 = unlist(lapply(tstat, function(x) x$estimate[1])),
                 mean2 = unlist(lapply(tstat, function(x) x$estimate[2])), pval = unlist(lapply(tstat, function(x) x$p.value)), comps = names(tstat), row.names = NULL)

DMV.AgecompsG = Map(cbind, DMV.AgecompsG, lapply(DMV.AgecompsG, function(x) postRNAres[match(x$gencodeID, postRNAres$gencodeID),]))
elementNROWS(DMV.AgecompsG)
elementNROWS(lapply(DMV.AgecompsG, function(x) unique(x$gencodeID)))
tstat = list(InotC.CnotI = t.test(DMV.AgecompsG$InotC$Tstat, DMV.AgecompsG$CnotI$Tstat),
             InotT.TnotI = t.test(DMV.AgecompsG$InotT$Tstat, DMV.AgecompsG$TnotI$Tstat),
             InotA.AnotI = t.test(DMV.AgecompsG$InotA$Tstat, DMV.AgecompsG$AnotI$Tstat),
             CnotT.TnotC = t.test(DMV.AgecompsG$CnotT$Tstat, DMV.AgecompsG$TnotC$Tstat),
             CnotA.AnotC = t.test(DMV.AgecompsG$CnotA$Tstat, DMV.AgecompsG$AnotC$Tstat),
             TnotA.AnotT = t.test(DMV.AgecompsG$TnotA$Tstat, DMV.AgecompsG$AnotT$Tstat))
ageG = data.frame(Tstat = unlist(lapply(tstat, function(x) x$statistic)), mean1 = unlist(lapply(tstat, function(x) x$estimate[1])),
                  mean2 = unlist(lapply(tstat, function(x) x$estimate[2])), pval = unlist(lapply(tstat, function(x) x$p.value)), comps = names(tstat), row.names = NULL)

expr = rbind(cbind(ct, Model = "Cell Type"), cbind(age, Model = "Age: Neurons"), cbind(ageG, Model = "Age: Glia"))
expr$FDR = p.adjust(expr$pval, method = "fdr")

write.csv(expr, quote=F, file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/t.test_DMV_gene_expression.csv")


## How does this relate to DNAm at these genes?

## CpG

load('/dcl01/lieber/ajaffe/lab/brain-epigenomics/bumphunting/BSobj_bsseqSmooth_Neuron_minCov_3.Rdata')
postpd = pData(BSobj)
postpd$Race[postpd$Race== "CAUC "] <- 'CAUC'
postpd$Sex[postpd$Sex == " M"] <- 'M'
postpd$RIN <- as.numeric(gsub(" ", "", postpd$RIN))
postpd$pg.DNA.nuclei.input <- as.numeric(postpd$pg.DNA.nuclei.input)
postpd$Reads <- as.numeric(postpd$Reads)
postpd$Percent.GreaterThan.Q30 <- as.numeric(postpd$Percent.GreaterThan.Q30)
meth =getMeth(BSobj, type = 'raw')
methMap = granges(BSobj)

ids = lapply(c(DMV.CTcomps,DMV.Agecomps,DMV.AgecompsG), function(x) paste0(x$Chr,":",x$Start,"-",x$End,":",x$Strand))
ids = lapply(ids, function(x) GRanges(x[which(x != "NA:NA-NA:NA")]))
oo = lapply(ids, function(x) findOverlaps(x, methMap))

meanMeth = lapply(oo, function(x) do.call("rbind", lapply(split(subjectHits(x), factor(queryHits(x))), function(ii) colMeans(t(t(meth[ii,]))))))
meanMeth = lapply(meanMeth, reshape2::melt)

tstat.CT = list(PnotN.NnotP = t.test(meanMeth$GnotN[meanMeth$GnotN$Var2 %in% pd[pd$Cell.Type=="Neuron","Data.ID"],"value"],
                                     meanMeth$NnotG[meanMeth$NnotG$Var2 %in% pd[pd$Cell.Type=="Glia","Data.ID"],"value"]),
                PnotG.GnotP = t.test(meanMeth$GnotN[meanMeth$GnotN$Var2 %in% pd[pd$Cell.Type=="Glia","Data.ID"],"value"],
                                     meanMeth$NnotG[meanMeth$NnotG$Var2 %in% pd[pd$Cell.Type=="Glia","Data.ID"],"value"]),
                GnotN.NnotG = t.test(meanMeth$GnotN[meanMeth$GnotN$Var2 %in% pd[pd$Cell.Type=="Glia","Data.ID"],"value"],
                                     meanMeth$NnotG[meanMeth$NnotG$Var2 %in% pd[pd$Cell.Type=="Neuron","Data.ID"],"value"]))
ct = data.frame(Tstat = unlist(lapply(tstat.CT, function(x) x$statistic)), mean1 = unlist(lapply(tstat.CT, function(x) x$estimate[1])),
                mean2 = unlist(lapply(tstat.CT, function(x) x$estimate[2])), pval = unlist(lapply(tstat.CT, function(x) x$p.value)), comps = names(tstat.CT), row.names = NULL)

names(meanMeth)[29:46] = paste0(names(meanMeth)[29:46],".G")
tstat.N = list(InotC.CnotI = t.test(meanMeth$InotC[meanMeth$InotC$Var2 %in% pd[pd$Age<=1 & pd$Cell.Type=="Neuron","Data.ID"],"value"],
                                    meanMeth$CnotI[meanMeth$CnotI$Var2 %in% pd[pd$Age>1 & pd$Age<=12 & pd$Cell.Type=="Neuron","Data.ID"],"value"]),
               InotT.TnotI = t.test(meanMeth$InotT[meanMeth$InotT$Var2 %in% pd[pd$Age<=1 & pd$Cell.Type=="Neuron","Data.ID"],"value"],
                                    meanMeth$TnotI[meanMeth$TnotI$Var2 %in% pd[pd$Age>12 & pd$Age<=17 & pd$Cell.Type=="Neuron","Data.ID"],"value"]),
               InotA.AnotI = t.test(meanMeth$InotA[meanMeth$InotA$Var2 %in% pd[pd$Age<=1 & pd$Cell.Type=="Neuron","Data.ID"],"value"],
                                    meanMeth$AnotI[meanMeth$AnotI$Var2 %in% pd[pd$Age>17 & pd$Cell.Type=="Neuron","Data.ID"],"value"]),
               CnotT.TnotC = t.test(meanMeth$CnotT[meanMeth$CnotT$Var2 %in% pd[pd$Age>1 & pd$Age<=12 & pd$Cell.Type=="Neuron","Data.ID"],"value"],
                                    meanMeth$TnotC[meanMeth$TnotC$Var2 %in% pd[pd$Age>12 & pd$Age<=17 & pd$Cell.Type=="Neuron","Data.ID"],"value"]),
               CnotA.AnotC = t.test(meanMeth$CnotA[meanMeth$CnotA$Var2 %in% pd[pd$Age>1 & pd$Age<=12 & pd$Cell.Type=="Neuron","Data.ID"],"value"],
                                    meanMeth$AnotC[meanMeth$AnotC$Var2 %in% pd[pd$Age>17 & pd$Cell.Type=="Neuron","Data.ID"],"value"]),
               TnotA.AnotT = t.test(meanMeth$TnotA[meanMeth$TnotA$Var2 %in% pd[pd$Age>12 & pd$Age<=17 & pd$Cell.Type=="Neuron","Data.ID"],"value"],
                                    meanMeth$AnotT[meanMeth$AnotT$Var2 %in% pd[pd$Age>17 & pd$Cell.Type=="Neuron","Data.ID"],"value"]))
neuro = data.frame(Tstat = unlist(lapply(tstat.N, function(x) x$statistic)), mean1 = unlist(lapply(tstat.N, function(x) x$estimate[1])),
                   mean2 = unlist(lapply(tstat.N, function(x) x$estimate[2])), pval = unlist(lapply(tstat.N, function(x) x$p.value)), comps = names(tstat.N), row.names = NULL)

tstat.G = list(InotC.CnotI = t.test(meanMeth$InotC.G[meanMeth$InotC.G$Var2 %in% pd[pd$Age<=1 & pd$Cell.Type=="Glia","Data.ID"],"value"],
                                    meanMeth$CnotI.G[meanMeth$CnotI.G$Var2 %in% pd[pd$Age>1 & pd$Age<=12 & pd$Cell.Type=="Glia","Data.ID"],"value"]),
               InotT.TnotI = t.test(meanMeth$InotT.G[meanMeth$InotT.G$Var2 %in% pd[pd$Age<=1 & pd$Cell.Type=="Glia","Data.ID"],"value"],
                                    meanMeth$TnotI.G[meanMeth$TnotI.G$Var2 %in% pd[pd$Age>12 & pd$Age<=17 & pd$Cell.Type=="Glia","Data.ID"],"value"]),
               InotA.AnotI = t.test(meanMeth$InotA.G[meanMeth$InotA.G$Var2 %in% pd[pd$Age<=1 & pd$Cell.Type=="Glia","Data.ID"],"value"],
                                    meanMeth$AnotI.G[meanMeth$AnotI.G$Var2 %in% pd[pd$Age>17 & pd$Cell.Type=="Glia","Data.ID"],"value"]),
               CnotT.TnotC = t.test(meanMeth$CnotT.G[meanMeth$CnotT.G$Var2 %in% pd[pd$Age>1 & pd$Age<=12 & pd$Cell.Type=="Glia","Data.ID"],"value"],
                                    meanMeth$TnotC.G[meanMeth$TnotC.G$Var2 %in% pd[pd$Age>12 & pd$Age<=17 & pd$Cell.Type=="Glia","Data.ID"],"value"]),
               CnotA.AnotC = t.test(meanMeth$CnotA.G[meanMeth$CnotA.G$Var2 %in% pd[pd$Age>1 & pd$Age<=12 & pd$Cell.Type=="Glia","Data.ID"],"value"],
                                    meanMeth$AnotC.G[meanMeth$AnotC.G$Var2 %in% pd[pd$Age>17 & pd$Cell.Type=="Glia","Data.ID"],"value"]),
               TnotA.AnotT = t.test(meanMeth$TnotA.G[meanMeth$TnotA.G$Var2 %in% pd[pd$Age>12 & pd$Age<=17 & pd$Cell.Type=="Glia","Data.ID"],"value"],
                                    meanMeth$AnotT.G[meanMeth$AnotT.G$Var2 %in% pd[pd$Age>17 & pd$Cell.Type=="Glia","Data.ID"],"value"]))
gli = data.frame(Tstat = unlist(lapply(tstat.G, function(x) x$statistic)), mean1 = unlist(lapply(tstat.G, function(x) x$estimate[1])),
                 mean2 = unlist(lapply(tstat.G, function(x) x$estimate[2])), pval = unlist(lapply(tstat.G, function(x) x$p.value)), comps = names(tstat.G), row.names = NULL)

dfCG = rbind(cbind(neuro, Comp = "Neurons Age"), cbind(gli, Comp = "Glia Age"), cbind(ct, Comp = "Cell Type"))



## CpH

load('/dcl01/lieber/ajaffe/lab/brain-epigenomics/bsseq/bsobj_by_chr/allChrs_postNatal_cleaned_nonCG_noHomogenate_highCov.Rdata')
BSobj_ch = BSobj
meth_ch = getMeth(BSobj_ch, type = 'raw')
methMap_ch = granges(BSobj_ch)

ids = lapply(c(DMV.CTcomps,DMV.Agecomps,DMV.AgecompsG), function(x) paste0(x$Chr,":",x$Start,"-",x$End,":",x$Strand))
ids = lapply(ids, function(x) GRanges(x[which(x != "NA:NA-NA:NA")]))
oo_ch = lapply(ids, function(x) findOverlaps(x, methMap_ch))
meanMeth_ch = lapply(oo_ch, function(x) do.call("rbind", lapply(split(subjectHits(x), factor(queryHits(x), levels=1:length(unique(queryHits(x))))), 
                                                                function(ii) colMeans(t(t(meth_ch[ii,]))))))
meanMeth_ch = lapply(meanMeth_ch, reshape2::melt)


GnotN.NnotG = t.test(meanMeth_ch$GnotN[meanMeth_ch$GnotN$Var2 %in% pd[pd$Cell.Type=="Glia","Data.ID"],"value"],
                     meanMeth_ch$NnotG[meanMeth_ch$NnotG$Var2 %in% pd[pd$Cell.Type=="Neuron","Data.ID"],"value"])
ct = data.frame(Tstat = GnotN.NnotG$statistic, mean1 = GnotN.NnotG$estimate[1], mean2 = GnotN.NnotG$estimate[2], pval = GnotN.NnotG$p.value, comps = "GnotN.NnotG")

names(meanMeth_ch)[29:46] = paste0(names(meanMeth_ch)[29:46],".G")
tstat.N = list(InotC.CnotI = t.test(meanMeth_ch$InotC[meanMeth_ch$InotC$Var2 %in% pd[pd$Age<=1 & pd$Cell.Type=="Neuron","Data.ID"],"value"],
                                    meanMeth_ch$CnotI[meanMeth_ch$CnotI$Var2 %in% pd[pd$Age>1 & pd$Age<=12 & pd$Cell.Type=="Neuron","Data.ID"],"value"]),
               InotT.TnotI = t.test(meanMeth_ch$InotT[meanMeth_ch$InotT$Var2 %in% pd[pd$Age<=1 & pd$Cell.Type=="Neuron","Data.ID"],"value"],
                                    meanMeth_ch$TnotI[meanMeth_ch$TnotI$Var2 %in% pd[pd$Age>12 & pd$Age<=17 & pd$Cell.Type=="Neuron","Data.ID"],"value"]),
               InotA.AnotI = t.test(meanMeth_ch$InotA[meanMeth_ch$InotA$Var2 %in% pd[pd$Age<=1 & pd$Cell.Type=="Neuron","Data.ID"],"value"],
                                    meanMeth_ch$AnotI[meanMeth_ch$AnotI$Var2 %in% pd[pd$Age>17 & pd$Cell.Type=="Neuron","Data.ID"],"value"]),
               CnotT.TnotC = t.test(meanMeth_ch$CnotT[meanMeth_ch$CnotT$Var2 %in% pd[pd$Age>1 & pd$Age<=12 & pd$Cell.Type=="Neuron","Data.ID"],"value"],
                                    meanMeth_ch$TnotC[meanMeth_ch$TnotC$Var2 %in% pd[pd$Age>12 & pd$Age<=17 & pd$Cell.Type=="Neuron","Data.ID"],"value"]),
               CnotA.AnotC = t.test(meanMeth_ch$CnotA[meanMeth_ch$CnotA$Var2 %in% pd[pd$Age>1 & pd$Age<=12 & pd$Cell.Type=="Neuron","Data.ID"],"value"],
                                    meanMeth_ch$AnotC[meanMeth_ch$AnotC$Var2 %in% pd[pd$Age>17 & pd$Cell.Type=="Neuron","Data.ID"],"value"]),
               TnotA.AnotT = t.test(meanMeth_ch$TnotA[meanMeth_ch$TnotA$Var2 %in% pd[pd$Age>12 & pd$Age<=17 & pd$Cell.Type=="Neuron","Data.ID"],"value"],
                                    meanMeth_ch$AnotT[meanMeth_ch$AnotT$Var2 %in% pd[pd$Age>17 & pd$Cell.Type=="Neuron","Data.ID"],"value"]))
neuro = data.frame(Tstat = unlist(lapply(tstat.N, function(x) x$statistic)), mean1 = unlist(lapply(tstat.N, function(x) x$estimate[1])),
                   mean2 = unlist(lapply(tstat.N, function(x) x$estimate[2])), pval = unlist(lapply(tstat.N, function(x) x$p.value)), comps = names(tstat.N), row.names = NULL)

tstat.G = list(InotC.CnotI = t.test(meanMeth_ch$InotC.G[meanMeth_ch$InotC.G$Var2 %in% pd[pd$Age<=1 & pd$Cell.Type=="Glia","Data.ID"],"value"],
                                    meanMeth_ch$CnotI.G[meanMeth_ch$CnotI.G$Var2 %in% pd[pd$Age>1 & pd$Age<=12 & pd$Cell.Type=="Glia","Data.ID"],"value"]),
               InotT.TnotI = t.test(meanMeth_ch$InotT.G[meanMeth_ch$InotT.G$Var2 %in% pd[pd$Age<=1 & pd$Cell.Type=="Glia","Data.ID"],"value"],
                                    meanMeth_ch$TnotI.G[meanMeth_ch$TnotI.G$Var2 %in% pd[pd$Age>12 & pd$Age<=17 & pd$Cell.Type=="Glia","Data.ID"],"value"]),
               InotA.AnotI = t.test(meanMeth_ch$InotA.G[meanMeth_ch$InotA.G$Var2 %in% pd[pd$Age<=1 & pd$Cell.Type=="Glia","Data.ID"],"value"],
                                    meanMeth_ch$AnotI.G[meanMeth_ch$AnotI.G$Var2 %in% pd[pd$Age>17 & pd$Cell.Type=="Glia","Data.ID"],"value"]),
               CnotT.TnotC = t.test(meanMeth_ch$CnotT.G[meanMeth_ch$CnotT.G$Var2 %in% pd[pd$Age>1 & pd$Age<=12 & pd$Cell.Type=="Glia","Data.ID"],"value"],
                                    meanMeth_ch$TnotC.G[meanMeth_ch$TnotC.G$Var2 %in% pd[pd$Age>12 & pd$Age<=17 & pd$Cell.Type=="Glia","Data.ID"],"value"]),
               CnotA.AnotC = t.test(meanMeth_ch$CnotA.G[meanMeth_ch$CnotA.G$Var2 %in% pd[pd$Age>1 & pd$Age<=12 & pd$Cell.Type=="Glia","Data.ID"],"value"],
                                    meanMeth_ch$AnotC.G[meanMeth_ch$AnotC.G$Var2 %in% pd[pd$Age>17 & pd$Cell.Type=="Glia","Data.ID"],"value"]),
               TnotA.AnotT = t.test(meanMeth_ch$TnotA.G[meanMeth_ch$TnotA.G$Var2 %in% pd[pd$Age>12 & pd$Age<=17 & pd$Cell.Type=="Glia","Data.ID"],"value"],
                                    meanMeth_ch$AnotT.G[meanMeth_ch$AnotT.G$Var2 %in% pd[pd$Age>17 & pd$Cell.Type=="Glia","Data.ID"],"value"]))
gli = data.frame(Tstat = unlist(lapply(tstat.G, function(x) x$statistic)), mean1 = unlist(lapply(tstat.G, function(x) x$estimate[1])),
                 mean2 = unlist(lapply(tstat.G, function(x) x$estimate[2])), pval = unlist(lapply(tstat.G, function(x) x$p.value)), comps = names(tstat.G), row.names = NULL)


dfmeth = rbind(cbind(dfCG[,-1], Context = "CpG"), cbind(rbind(cbind(neuro, Comp = "Neurons Age"), cbind(gli, Comp = "Glia Age"), cbind(ct, Comp = "Cell Type")), Context = "CpH"))
dfmeth$FDR = p.adjust(dfmeth$pval, method = "fdr")

write.csv(dfmeth, quote=F, file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/t.test_DMV_gene_methylation.csv")


## Which ones are TFs, and what are they doing in terms of our TF analysis? 

DMV.CTcomps = Map(cbind, DMV.CTcomps, TF = lapply(DMV.CTcomps, function(x) ifelse(x$gencodeID %in% tfs$gencodeID, "TF","no")))
DMV.Agecomps = Map(cbind, DMV.Agecomps, TF = lapply(DMV.Agecomps, function(x) ifelse(x$gencodeID %in% tfs$gencodeID, "TF","no")))
DMV.AgecompsG = Map(cbind, DMV.AgecompsG, TF = lapply(DMV.AgecompsG, function(x) ifelse(x$gencodeID %in% tfs$gencodeID, "TF","no")))

save(DMV.CTcomps, DMV.Agecomps, DMV.AgecompsG, file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/DMV_gene_comps.rda")

elementNROWS(lapply(DMV.CTcomps, function(x) x[x$TF=="TF",]))
elementNROWS(lapply(DMV.Agecomps, function(x) x[x$TF=="TF",]))
elementNROWS(lapply(DMV.AgecompsG, function(x) x[x$TF=="TF",]))


## Write the lists for supplementary tables

x = lapply(DMV.CTcomps[-grep("shared", names(DMV.CTcomps))], function(y) y[,c(1:13,31,18,22,26,29)])
x = do.call(rbind, Map(cbind, x, DMV.Group = as.list(names(x))))
colnames(x) = c("Chr","Start","End","Strand","Length","gencodeID","ensemblID","gene_type","Symbol","EntrezID","Class","meanExprs","NumTx","TF","Coeff",
                "Tstat", "pval", "padj","DMV.Group") 

y = lapply(DMV.Agecomps[-grep("shared", names(DMV.Agecomps))], function(y) y[,c(1:13,15,17,19,21:22)])
y2 = lapply(DMV.AgecompsG[-grep("shared", names(DMV.AgecompsG))], function(y) y[,c(1:13,15,17,19,21:22)])
y = do.call(rbind, Map(cbind, y, DMV.Group = as.list(names(y))))
y2 = do.call(rbind, Map(cbind, y2, DMV.Group = as.list(names(y2))))
colnames(y) = colnames(y2) = colnames(x)

z = rbind(cbind(x, model = "Cell Type"), cbind(y, model = "Age: Neurons"), cbind(y2, model = "Age: Glia"))
z = z[,colnames(z)!="meanExprs"]
rownames(z) = NULL

write.csv(z[z$TF=="TF",], quote=F, file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/DMV_genes_annotated.csv")

length(unique(z$gencodeID)) # 15514
length(unique(z[z$TF=="TF","gencodeID"])) # 458



## What's the relationship between heterochromatin spreading and TFs in hypo-DMRs? 

## Roadmap state in different groups of DMV genes? Compare to parent DMV roadmap state in different groups





