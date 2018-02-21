library(GenomicFeatures)
library(ggplot2)
library(GenomicRanges)
library(clusterProfiler)
require(org.Hs.eg.db)


load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/PMDs_methylSeekR.rda")
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/PMDs_methylSeekR_prenatal.rda")
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/PMDs_methylSeekR_sortedLister.rda")
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/DMR/DMR_objects.rda")


# rename Lister samples

listernames = c("GSM1173773" = "ListerNeuron1", "GSM1173774" = "ListerGlia1", "GSM1173776" = "ListerNeuron2", "GSM1173777" = "ListerGlia2")
names(PMDsegments.CG.sortedLister) = listernames[match(names(PMDsegments.CG.sortedLister), names(listernames))]
names(PMDsegments.CG) = paste0("Postnatal-", names(PMDsegments.CG))
names(PMDsegments.CGpren) = paste0("Prenatal-", names(PMDsegments.CGpren))

# What proportion of the genome is covered by each sample?

sLengths = c(lapply(PMDsegments.CG.sortedLister, reduce), lapply(PMDsegments.CG, reduce), 
             lapply(PMDsegments.CGpren, reduce))
pmds_total = c(lapply(PMDsegments.CG.sortedLister, function(x) reduce(x[x$type=="PMD"])),
               lapply(PMDsegments.CG, function(x) reduce(x[x$type=="PMD"])),
               lapply(PMDsegments.CGpren, function(x) reduce(x[x$type=="PMD"])))

propPMDs = mapply(function(s, l) round(sum(as.numeric(width(s)))/sum(as.numeric(width(l))), 2), pmds_total, sLengths)
propPMDs[order(propPMDs)] # 


## What proportion of called PMDs are >100kb?

pmds = c(lapply(PMDsegments.CG.sortedLister, function(x) x[x$type=="PMD"]),lapply(PMDsegments.CG, function(x) x[x$type=="PMD"]),
         lapply(PMDsegments.CGpren, function(x) x[x$type=="PMD"]))

propPMDs = round(elementNROWS(lapply(pmds, function(x) x[which(width(x)>100000)]))/elementNROWS(pmds), 2)
propPMDs[order(propPMDs)] # 3-10% of the reported PMDs are >100kb
# So PMDs aren't really a thing in brain cells...


## is there a difference in prenatal and postnatal?

num = elementNROWS(lapply(pmds, function(x) x[which(width(x)>100000)]))
t.test(num[grep("Pre",names(num))], num[grep("Post",names(num))]) # t = -0.17579, df = 45.961, p-value = 0.8612


## Create new objects reflecting size limit

pmds = c(lapply(PMDsegments.CG.sortedLister, function(x) x[x$type=="PMD"]),lapply(PMDsegments.CG, function(x) x[x$type=="PMD"]),
         lapply(PMDsegments.CGpren, function(x) x[x$type=="PMD"]))
nonpmds = c(lapply(PMDsegments.CG.sortedLister, function(x) x[x$type=="notPMD"]),lapply(PMDsegments.CG, function(x) x[x$type=="notPMD"]),
            lapply(PMDsegments.CGpren, function(x) x[x$type=="notPMD"]))
newnon = lapply(pmds, function(x) x[which(width(x)<100000)])

merge = mapply(function(gr1, gr2) reduce(c(gr1,gr2)), nonpmds, newnon)
oo1 = mapply(function(m, gr1) findOverlaps(m, gr1), merge, nonpmds)
oo2 = mapply(function(m, gr2) findOverlaps(m, gr2), merge, newnon)
grm1 = mapply(function(gr1, o1) split(gr1[subjectHits(o1)], queryHits(o1)), nonpmds, oo1)
grm2 = mapply(function(gr2, o2) split(gr2[subjectHits(o2)], queryHits(o2)), newnon, oo2)

merge2 = rep( list(list()), length(merge)) 
for (j in 1:length(merge)) {
  for (i in 1:length(merge[[j]])) {
    if (as.character(i) %in% names(grm2[[j]])) {
    merge2[[j]][[i]] = c(grm1[[j]][[as.character(i)]], grm2[[j]][[as.character(i)]]) } else {
    merge2[[j]][[i]] = grm1[[j]][[as.character(i)]]
  } } }
nCG = lapply(merge2, function(y) lapply(y, function(x) sum(x$nCG)))
merge = Map(cbind, lapply(merge, as.data.frame), type = "notPMD", nCG = lapply(nCG, unlist))
merge = lapply(merge, function(x) makeGRangesFromDataFrame(x, keep.extra.columns = T))

pmds = lapply(pmds, function(x) x[which(width(x)>=100000)])

total = mapply(function(non, pmd) c(non, pmd), merge, pmds)
total = lapply(total, function(x) x[order(x)])

table(unlist(lapply(c(PMDsegments.CG.sortedLister, PMDsegments.CG, PMDsegments.CGpren), function(x) sum(x$nCG)))==18664892)
table(unlist(lapply(total, function(x) sum(x$nCG)))==18664892)
names(unlist(lapply(total, function(x) sum(x$nCG))))[unlist(lapply(total, function(x) sum(x$nCG)))!=18664892 & 
                                                       unlist(lapply(total, function(x) sum(x$nCG)))!=18308447]
# ??? "Postnatal-WGC052317L" "Postnatal-WGC052318L" "Postnatal-WGC052311L" "Postnatal-WGC059603L"

sLengths = c(lapply(PMDsegments.CG.sortedLister, reduce), lapply(PMDsegments.CG, reduce), 
             lapply(PMDsegments.CGpren, reduce))
red = lapply(total, reduce)
table(elementNROWS(red)==25) # all true
table(unlist(mapply(function(x,y) identical(x,y), red, sLengths))=="TRUE") # all true

save(total, file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/PMDs_methylSeekR_100kbLimit.rda")

################################################################
#### Analysis of regions that still meet PMD definition criteria
################################################################

pmds = c(lapply(PMDsegments.CG, function(x) x[x$type=="PMD"]), lapply(PMDsegments.CGpren, function(x) x[x$type=="PMD"]))
pmds = lapply(pmds, function(x) x[which(width(x)>=100000)])

## How many are identified?

num = data.frame(num = unlist(elementNROWS(pmds)), id = names(pmds))
num[which(num$id %in% names(PMDsegments.CGpren)),"celltype"] = "Prenatal"
num[which(num$id %in% pd[pd$Cell.Type=="Neuron","Data.ID"]),"celltype"] = "Neuron"
num[which(num$id %in% pd[pd$Cell.Type=="Glia","Data.ID"]),"celltype"] = "Glia"
write.csv(num, quote=F, file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/PMD_number.csv")

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/CREs/figures/PMD_number_byAge_postnatal.pdf")

w = num[which(num$celltype %in% c("Neuron", "Glia")),]
w = cbind(w, Age = pd[match(w$id, pd$Data.ID), "Age"])
cor.test(x = w[which(w$celltype=="Neuron"),"Age"], y = w[which(w$celltype=="Neuron"),"num"])
#t = 1.506, df = 22, p-value = 0.1463
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
#  -0.1114205  0.6312632
#sample estimates:
#  cor 
#0.3057172 
cor.test(x = w[which(w$celltype=="Glia"),"Age"], y = w[which(w$celltype=="Glia"),"num"])
#t = -0.39397, df = 6, p-value = 0.7072
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
#  -0.7765714  0.6146563
#sample estimates:
#  cor 
#-0.1587964

ggplot(w, aes(x = Age, y = num, colour = celltype)) + geom_point() +
  geom_smooth(method=lm, se=T, fullrange=TRUE) + scale_colour_brewer(8, palette="Dark2") + 
  labs(fill="") + theme_classic() +
  ylab("Number") + xlab("Age (Years)") +
  ggtitle("Number PMDs Identified by Age") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20), legend.title=element_blank()) + theme(legend.position="bottom")
ggplot(num, aes(x = celltype, y = num)) + geom_boxplot() +
  theme_classic() +
  labs(fill="") +
  ylab("Count") + xlab("") +
  ggtitle("Number of PMDs") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
dev.off()



## how long are these regions?

widths = do.call(rbind, Map(cbind, lapply(pmds, function(y) data.frame(mean = mean(width(y)), median = median(width(y)), 
                                             sd = sd(width(y)), min = min(width(y)), max = max(width(y)))), id = as.list(names(pmds))))
widths[which(widths$id %in% names(PMDsegments.CGpren)),"celltype"] = "Prenatal"
widths[which(widths$id %in% pd[pd$Cell.Type=="Neuron","Data.ID"]),"celltype"] = "Neuron"
widths[which(widths$id %in% pd[pd$Cell.Type=="Glia","Data.ID"]),"celltype"] = "Glia"
write.csv(widths, quote=F, file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/PMD_widths.csv")

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/CREs/figures/PMD_widths.pdf")

w = widths[which(widths$celltype %in% c("Neuron", "Glia")),]
w = cbind(w, Age = pd[match(w$id, pd$Data.ID), "Age"])
cor.test(x = w[which(w$celltype=="Neuron"),"Age"], y = w[which(w$celltype=="Neuron"),"median"])
#t = 3.4531, df = 22, p-value = 0.002265
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
#  0.2490229 0.8039808
#sample estimates:
#  cor 
#0.5928635
cor.test(x = w[which(w$celltype=="Glia"),"Age"], y = w[which(w$celltype=="Glia"),"median"])
#t = 0.4732, df = 6, p-value = 0.6528
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
#  -0.5944511  0.7889042
#sample estimates:
#  cor 
#0.1896757

ggplot(w, aes(x = Age, y = median, colour = celltype)) + geom_point() +
  geom_smooth(method=lm, se=T, fullrange=TRUE) + scale_colour_brewer(8, palette="Dark2") + 
  labs(fill="") + theme_classic() +
  ylab("Bases") + xlab("Age (Years)") +
  ggtitle("Median PMD Widths by Age") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20), legend.title=element_blank()) + theme(legend.position="bottom")
ggplot(widths, aes(x = celltype, y = median)) + geom_boxplot() +
  theme_classic() +
  labs(fill="") +
  ylab("Bases") + xlab("") +
  ggtitle("Median PMD Widths") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
ggplot(widths, aes(x = celltype, y = mean)) + geom_boxplot() +
  theme_classic() +
  labs(fill="") +
  ylab("Bases") + xlab("") +
  ggtitle("Mean PMD Widths") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
dev.off()


## What proportion of the genome do they cover?

sLengths = c(lapply(PMDsegments.CG, reduce), lapply(PMDsegments.CGpren, reduce))
percGenome = data.frame(percent = unlist(mapply(function(p, t) round((sum(width(p))/sum(as.numeric(width(t)))*100),2), pmds, sLengths)),
                        id = names(pmds))
percGenome[which(percGenome$id %in% names(PMDsegments.CGpren)),"celltype"] = "Prenatal"
percGenome[which(percGenome$id %in% pd[pd$Cell.Type=="Neuron","Data.ID"]),"celltype"] = "Neuron"
percGenome[which(percGenome$id %in% pd[pd$Cell.Type=="Glia","Data.ID"]),"celltype"] = "Glia"
write.csv(percGenome, quote=F, file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/PMD_genomeCoverage.csv")

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/CREs/figures/PMD_percentGenome.pdf")

w = percGenome[which(percGenome$celltype %in% c("Neuron", "Glia")),]
w = cbind(w, Age = pd[match(w$id, pd$Data.ID), "Age"])
cor.test(x = w[which(w$celltype=="Neuron"),"Age"], y = w[which(w$celltype=="Neuron"),"percent"])
#t = -1.0978, df = 22, p-value = 0.2841
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
#  -0.5781465  0.1932627
#sample estimates:
#  cor 
#-0.2279016 
cor.test(x = w[which(w$celltype=="Glia"),"Age"], y = w[which(w$celltype=="Glia"),"percent"])
#t = -0.13414, df = 6, p-value = 0.8977
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
#  -0.7311795  0.6760418
#sample estimates:
#  cor 
#-0.05468004 

ggplot(w, aes(x = Age, y = percent, colour = celltype)) + geom_point() +
  geom_smooth(method=lm, se=T, fullrange=TRUE) + scale_colour_brewer(8, palette="Dark2") + 
  labs(fill="") +  theme_classic() +
  ylab("Percent") + xlab("Age (Years)") +
  ggtitle("Percent of Genome that is PMD by Age") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20), legend.title=element_blank()) + theme(legend.position="bottom")
ggplot(percGenome, aes(x = celltype, y = percent)) + geom_boxplot() +
  theme_classic() +
  labs(fill="") +
  ylab("Count") + xlab("") +
  ggtitle("Percent Genome: PMD") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
dev.off()


### How much of the PMDs overlap?

# Master list of all bases in each setting
pmaster = reduce(do.call(getMethod(c, "GenomicRanges"), GRangesList(pmds))) # 43995 total
c(mean=mean(width(pmaster)), median=median(width(pmaster)), sd=sd(width(pmaster)), min=min(width(pmaster)), max=max(width(pmaster)))/1000
#         mean     median         sd        min        max 
#kb:  377.0875   160.1570  2170.9122   100.0660 59373.5660 


## regions that are shared by samples

all = list(All = Reduce(intersect, pmds),
           Prenatal = Reduce(intersect, lapply(PMDsegments.CGpren, function(x) x[which(x$type=="PMD" & width(x)>=100000)])),
           Postnatal = Reduce(intersect, lapply(PMDsegments.CG, function(x) x[which(x$type=="PMD" & width(x)>=100000)])),
           Neurons = Reduce(intersect, lapply(PMDsegments.CG, function(x) x[which(x$type=="PMD" & width(x)>=100000)])
                            [which(names(PMDsegments.CG) %in% pd[pd$Cell.Type=="Neuron","Data.ID"])]),
           Glia = Reduce(intersect, lapply(PMDsegments.CG, function(x) x[which(x$type=="PMD" & width(x)>=100000)])
                                           [which(names(PMDsegments.CG) %in% pd[pd$Cell.Type=="Glia","Data.ID"])]))
elementNROWS(all)
#      All  Prenatal Postnatal   Neurons      Glia 
#      105       137        98       189       277 

round((unlist(lapply(all, function(x) (sum(width(x))/sum(as.numeric(width(sLengths[[1]]))))*100)))/
        (sum(width(pmaster))/sum(as.numeric(width(sLengths[[1]])))*100),2)
#      All  Prenatal Postnatal   Neurons      Glia 
#     0.25      0.27      0.25      0.26      0.30 
# 25-30% of total bases in the PMD state are shared by all or most

shared = mapply(function(all,ind) lapply( ind, function(x) round(all / x *100,2)), lapply(all, function(x) (sum(width(x))/sum(as.numeric(width(sLengths[[1]])))*100)), 
                 lapply( list(pmds, 
                              lapply(PMDsegments.CGpren, function(x) x[which(x$type=="PMD" & width(x)>=100000)]), 
                              lapply(PMDsegments.CG, function(x) x[which(x$type=="PMD" & width(x)>=100000)]),
                              lapply(PMDsegments.CG, function(x) x[which(x$type=="PMD" & width(x)>=100000)])[which(names(PMDsegments.CG) %in% pd[pd$Cell.Type=="Neuron","Data.ID"])],
                              lapply(PMDsegments.CG, function(x) x[which(x$type=="PMD" & width(x)>=100000)])[which(names(PMDsegments.CG) %in% pd[pd$Cell.Type=="Glia","Data.ID"])]),
                         function(y) lapply(y, function(x) sum(width(x))/sum(as.numeric(width(sLengths[[1]])))*100)), SIMPLIFY = F)
shared = do.call(rbind, Map(data.frame, percent = lapply(shared, function(x) unlist(x, recursive=F)), id = lapply(shared, names), group = as.list(names(shared))))

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/CREs/figures/percent_shared_PMDs_perSample_byBaseCoverage.pdf")
ggplot(shared, aes(x = group, y = percent)) + geom_boxplot() +
  labs(fill="") + theme_classic() +
  ylim(0,100) +
  ylab("Percent") + xlab("") +
  ggtitle("Percent Bases Shared Per Sample: PMD") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
dev.off()


## Find overlaps with DMRs, features, shared groups

islands = read.table("/dcl01/lieber/ajaffe/Amanda/annotation_objects/cpgIslandExt.hg19.txt", sep="\t", header = T)
rpmsk = read.table("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/RepeatMasker_genomewide_HG19.txt", header=T)
features = c(genes = makeGRangesFromDataFrame(geneMap, keep.extra.columns = T), 
             rpmskgr = makeGRangesFromDataFrame(rpmsk, seqnames.field="genoName",start.field="genoStart",end.field="genoEnd",keep.extra.columns=TRUE),
             islands = makeGRangesFromDataFrame(islands, keep.extra.columns = T, start.field = "chromStart", end.field = "chromEnd"))
lapply(features, head)

DMR = lapply(DMR, function(x) x[which(x$sig=="FWER < 0.05"),])
DMRgr = lapply(DMR, function(x) makeGRangesFromDataFrame(x, keep.extra.columns = T))

oo = lapply(pmds, function(x) lapply(c(DMRgr, all, features), function(d) findOverlaps(x,d)))

pmdDMR = lapply(pmds, as.data.frame)
pmdDMR = lapply(pmdDMR, function(x) data.frame(x, rnum = 1:nrow(x)))
pmdDMR = mapply(function(d,oo) data.frame(d, CT = ifelse(d$rnum %in% queryHits(oo$CellType), "CT", "no"),
                                        Age = ifelse(d$rnum %in% queryHits(oo$Age), "Age", "no"),
                                        Interaction = ifelse(d$rnum %in% queryHits(oo$Interaction), "Int", "no"),
                                        All = ifelse(d$rnum %in% queryHits(oo$All), "All", "no"),
                                        Prenatal = ifelse(d$rnum %in% queryHits(oo$Prenatal), "Prenatal", "no"),
                                        Postnatal = ifelse(d$rnum %in% queryHits(oo$Postnatal), "Postnatal", "no"),
                                        Neurons = ifelse(d$rnum %in% queryHits(oo$Neurons), "Neurons", "no"),
                                        Glia = ifelse(d$rnum %in% queryHits(oo$Glia), "Glia", "no"),
                                        genes = ifelse(d$rnum %in% queryHits(oo$genes), "Gene", NA),
                                        islands = ifelse(d$rnum %in% queryHits(oo$islands), "CpG-Island", "non-Island")), pmdDMR, oo, SIMPLIFY = F) 
pmdDMR = lapply(pmdDMR, function(d) data.frame(d, tog = paste(d$CT, d$Age, d$Interaction, d$All, d$Prenatal, 
                                                          d$Postnatal, d$Neurons, d$Glia, sep = ":"),
                                           dmr = paste(d$CT, d$Age, d$Interaction, sep = ":"),
                                           regionID = paste0(d$seqnames,":",d$start,"-", d$end),
                                           tog2 = paste(d$All, d$Prenatal, d$Postnatal, d$Neurons, d$Glia, sep = ":")))
for (i in 1:length(pmdDMR)) {
  if (length(unique(queryHits(oo[[i]][["rpmskgr"]])))==nrow(pmdDMR[[i]])) {
  pmdDMR[[i]] = data.frame(pmdDMR[[i]][queryHits(oo[[i]][["rpmskgr"]]),], 
                               repeats = as.character(features$rpmskgr$repClass)[subjectHits(oo[[i]][["rpmskgr"]])]) } else {
  pmdDMR[[i]] = rbind(data.frame(pmdDMR[[i]][queryHits(oo[[i]][["rpmskgr"]]),], 
                                 repeats = as.character(features$rpmskgr$repClass)[subjectHits(oo[[i]][["rpmskgr"]])]), 
                      data.frame(pmdDMR[[i]][-unique(queryHits(oo[[i]][["rpmskgr"]])),], repeats = "No repeats"))                    
                    }
}

pmdDMRgr = lapply(pmdDMR, function(x) makeGRangesFromDataFrame(x, keep.extra.columns = T))
oo = lapply(pmdDMRgr, function(x) findOverlaps(x, makeGRangesFromDataFrame(geneMap)))

postnatalpd = pd
save(pmdDMR, postnatalpd, geneMap,oo, file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/PMDs_annotated.rda")


### Annotate to genomic features

load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/PMDs_annotated.rda")
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

oo = lapply(pmdDMR, function(y) findOverlaps(gr.clusters, makeGRangesFromDataFrame(y)))
df.clusters = lapply(oo, function(y) data.frame(df.clusters, overlap = ifelse(df.clusters$rnum %in% queryHits(y), "hit","no")))


## how many fall within CpG islands?

CpGIslands = lapply(pmdDMR, function(y) round(length(unique(y[which(y$islands=="CpG-Island"),,]$regionID))/length(unique(y$regionID))*100,2))
CpGIslands = data.frame(islands = unlist(CpGIslands), id = names(CpGIslands))
CpGIslands[!(CpGIslands$id %in% postnatalpd$Data.ID),"celltype"] = "Prenatal"
CpGIslands[which(CpGIslands$id %in% postnatalpd[which(postnatalpd$Cell.Type=="Neuron"),"Data.ID"]),"celltype"] = "Neuron"
CpGIslands[which(CpGIslands$id %in% postnatalpd[which(postnatalpd$Cell.Type=="Glia"),"Data.ID"]),"celltype"] = "Glia"
write.csv(CpGIslands, quote=F,file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/PMD_CpG_Island_Overlap.csv")

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/CREs/figures/PMD_overlap_with_CpG-Islands.pdf")
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

pmdDMRdt = lapply(pmdDMR, data.table)
repeats = lapply(pmdDMRdt, function(y) y[,length(unique(regionID)), by = "repeats"])
repeats = lapply(repeats, function(y) data.frame(y, perc = round(y$V1/sum(y$V1)*100,2)))
repeats = do.call(rbind, Map(cbind, repeats, id = as.list(names(repeats))))
repeats$repeats = factor(repeats$repeats, levels = c("SINE","LINE","Simple_repeat","DNA","LTR","Low_complexity","No repeats","Unknown","snRNA",
                                                     "Other","DNA?","Satellite","RC","scRNA","tRNA","SINE?","srpRNA","RNA","rRNA","LINE?","LTR?","Unknown?"))
repeats[!(repeats$id %in% postnatalpd$Data.ID),"celltype"] = "Prenatal"
repeats[which(repeats$id %in% postnatalpd[which(postnatalpd$Cell.Type=="Neuron"),"Data.ID"]),"celltype"] = "Neuron"
repeats[which(repeats$id %in% postnatalpd[which(postnatalpd$Cell.Type=="Glia"),"Data.ID"]),"celltype"] = "Glia"
write.csv(repeats, quote=F,file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/PMD_repeats_Overlap.csv")

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/CREs/figures/PMD_overlap_with_repeats.pdf")
x = repeats[which(repeats$repeats %in% c("SINE","LINE","Simple_repeat","DNA","LTR","Low_complexity","No repeats")),]
x$repeats = droplevels(x$repeats)
ggplot(x, aes(x = celltype, y = perc, fill = repeats)) + 
  theme_classic() + geom_bar(position = "fill",stat = "identity") +
  labs(fill="") +
  ylab("Percent") + 
  xlab("") +
  ggtitle("Percent Overlapping Repeats: PMD") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20), legend.title=element_blank()) + theme(legend.position="bottom")
dev.off()


## Overlap with DMRs

dmr = lapply(pmdDMRdt, function(y) y[,length(unique(regionID)), by = "dmr"])
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

write.csv(dmr, quote=F,file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/PMD_dmr_Overlap.csv")
write.csv(models, quote=F,file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/PMD_dmr_Overlap_byModel.csv")


pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/CREs/figures/PMD_overlap_with_dmr.pdf", width = 14)
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

# Are DMRs overlapping PMDs at the PMD periphery, or uniformly distributed?

pmdDMR = Map(cbind, lapply(pmdDMR, function(x) data.frame(x[,1:5], dmr = x$dmr, regionID = x$regionID)), id = as.list(names(pmdDMR)))
pmdDMR = lapply(pmdDMR, unique)
tiles = lapply(pmdDMR, function(x) tile(makeGRangesFromDataFrame(x[which(x$dmr!="no:no:no"),]), n = 100))
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
        CToverlaps[[i]][[j]] = data.frame(region = 0)
      } } }
for (i in 1:length(ti)) {
  for (j in 1:length(ti[[i]])) {
    if (length(queryHits(Aov[[i]][[j]]))>0) {
      Ageoverlaps[[i]][[j]] = data.frame(region = queryHits(Aov[[i]][[j]])) } else {
        Ageoverlaps[[i]][[j]] = data.frame(region = 0)
      } } }
for (i in 1:length(ti)) {
  for (j in 1:length(ti[[i]])) {
    if (length(queryHits(Iov[[i]][[j]]))>0) {
      Intoverlaps[[i]][[j]] = data.frame(region = queryHits(Iov[[i]][[j]])) } else {
        Intoverlaps[[i]][[j]] = data.frame(region = 0)
      } } }
CToverlaps = do.call(rbind, Map(cbind, lapply(CToverlaps, function(x) do.call(rbind, x)), celltype = as.list(names(ti))))
Ageoverlaps = do.call(rbind, Map(cbind, lapply(Ageoverlaps, function(x) do.call(rbind, x)), celltype = as.list(names(ti))))
Intoverlaps = do.call(rbind, Map(cbind, lapply(Intoverlaps, function(x) do.call(rbind, x)), celltype = as.list(names(ti))))
overlaps = rbind(cbind(CToverlaps, model = "Cell Type"), cbind(Ageoverlaps, model = "Age"),
                 cbind(Intoverlaps, model = "Interaction"))

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/CREs/figures/DMR_position_in_PMD.pdf", width = 10)
ggplot(overlaps, aes(region, colour = model)) + geom_density() + theme_classic() +
  labs(fill="") + facet_grid(. ~ celltype) +
  ylab("Density") + xlab("% in PMD") +
  ggtitle("DMR Position in Gene") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20), legend.position="bottom", legend.title=element_blank())
dev.off()


## Overlap with Lister samples

load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/PMDs_methylSeekR_sortedLister.rda")
listerpd = c("GSM1173773" = "53yr_NeuN_pos", "GSM1173774" = "53yr_NeuN_neg",
             "GSM1173776" = "55yr_NeuN_pos", "GSM1173777" = "55yr_NeuN_neg")
listerpmds = lapply(PMDsegments.CG.sortedLister, function(x) x[which(x$type=="PMD" & width(x)>100000)])
names(listerpmds) = listerpd[match(names(listerpmds),names(listerpd))]

# Total bases
celltype = list(Neuron = pmdDMR[which(names(pmdDMR) %in% postnatalpd[postnatalpd$Cell.Type=="Neuron","Data.ID"])],
                Glia = pmdDMR[which(names(pmdDMR) %in% postnatalpd[postnatalpd$Cell.Type=="Glia","Data.ID"])])
celltype = lapply(celltype, function(y) lapply(y, makeGRangesFromDataFrame, keep=TRUE))
celltype = lapply(celltype, function(x) reduce(do.call(getMethod(c, "GenomicRanges"), GRangesList(x))))
lister = list(Neuron = reduce(do.call(getMethod(c, "GenomicRanges"), GRangesList(listerpmds[grep("pos", names(listerpmds))]))),
              Glia = reduce(do.call(getMethod(c, "GenomicRanges"), GRangesList(listerpmds[grep("neg", names(listerpmds))]))))
mapply(function(c,l) c( OurOverlap = round(sum(width(Reduce(intersect, list(c,l))))/sum(width(c))*100,2),
                        ListerOverlap = round(sum(width(Reduce(intersect, list(c,l))))/sum(width(l))*100,2)), celltype, lister, SIMPLIFY = F)
# Neuron
#OurOverlap ListerOverlap 
#91.96         96.14 
# Glia
#OurOverlap ListerOverlap 
#70.66         97.06 


## association of existing PMDs with genes

oo = lapply(pmdDMR, function(x) findOverlaps(makeGRangesFromDataFrame(geneMap), makeGRangesFromDataFrame(x))) 
pmdDMR = mapply(function(p, ov) if (nrow(p)>length(unique(subjectHits(ov)))) {
                  rbind(data.frame(p[subjectHits(ov),], nearestSymbol = geneMap[queryHits(ov),"Symbol"],
                                                 nearestID = geneMap[queryHits(ov),"gencodeID"], 
                                                 EntrezID = geneMap[queryHits(ov),"EntrezID"]),
                                      data.frame(p[-unique(subjectHits(ov)),], nearestSymbol = "NoGeneOverlap", 
                                                 nearestID = "NoGeneOverlap", EntrezID = "NoGeneOverlap")) } else {
                        data.frame(p[subjectHits(ov),], nearestSymbol = geneMap[queryHits(ov),"Symbol"],
                                   nearestID = geneMap[queryHits(ov),"gencodeID"], EntrezID = geneMap[queryHits(ov),"EntrezID"])      
                                                 },pmdDMR, oo, SIMPLIFY = F)

pmdDMRdt = lapply(pmdDMR, data.table)
geneoverlap = lapply(pmdDMR, function(y) round(length(unique(as.character(data.frame(y)[which(y$nearestID!="NoGeneOverlap"),"regionID"])))/
                                                 length(unique(as.character(y$regionID)))*100,2))
geneoverlap = data.frame(geneoverlap = unlist(geneoverlap), id = names(geneoverlap))
geneoverlap[!(geneoverlap$id %in% postnatalpd$Data.ID),"celltype"] = "Prenatal"
geneoverlap[which(geneoverlap$id %in% postnatalpd[which(postnatalpd$Cell.Type=="Neuron"),"Data.ID"]),"celltype"] = "Neuron"
geneoverlap[which(geneoverlap$id %in% postnatalpd[which(postnatalpd$Cell.Type=="Glia"),"Data.ID"]),"celltype"] = "Glia"
write.csv(geneoverlap, quote=F,file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/PMD_Gene_Overlap.csv")

numgenes = lapply(pmdDMRdt, function(y) y[,length(unique(nearestID)), by = "regionID"])
numgenes = do.call(rbind, Map(cbind, lapply(numgenes, function(y) 
  data.frame(mean = mean(y$V1), median = median(y$V1), sd = sd(y$V1), min = min(y$V1), max = max(y$V1))), id = as.list(names(numgenes))))
numgenes[!which(numgenes$id %in% postnatalpd$Data.ID),"celltype"] = "Prenatal"
numgenes[which(numgenes$id %in% postnatalpd[postnatalpd$Cell.Type=="Neuron","Data.ID"]),"celltype"] = "Neuron"
numgenes[which(numgenes$id %in% postnatalpd[postnatalpd$Cell.Type=="Glia","Data.ID"]),"celltype"] = "Glia"
write.csv(numgenes, quote=F, file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/PMD_numgenes.csv")

w = numgenes[which(numgenes$celltype %in% c("Neuron", "Glia")),]
w = cbind(w, Age = postnatalpd[match(w$id, postnatalpd$Data.ID), "Age"])
cor.test(x = w[which(w$celltype=="Neuron"),"Age"], y = w[which(w$celltype=="Neuron"),"mean"])
#t = -2.7944, df = 22, p-value = 0.01057
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
#  -0.7585919 -0.1366317
#sample estimates:
#  cor 
#-0.5118186 
cor.test(x = w[which(w$celltype=="Glia"),"Age"], y = w[which(w$celltype=="Glia"),"mean"])
#t = 0.50224, df = 6, p-value = 0.6334
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
#  -0.5868793  0.7932557
#sample estimates:
#  cor 
#0.2008616 

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/CREs/figures/PMD_overlap_with_Genes.pdf")
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
  ggtitle("Median PMD Gene Number by Age") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20), legend.title=element_blank()) + theme(legend.position="bottom")
ggplot(numgenes, aes(x = celltype, y = median)) + geom_boxplot() +
  theme_classic() +
  labs(fill="") +
  ylab("Bases") + xlab("") +
  ggtitle("Median PMD Gene Number") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
ggplot(numgenes, aes(x = celltype, y = mean)) + geom_boxplot() +
  theme_classic() +
  labs(fill="") +
  ylab("Bases") + xlab("") +
  ggtitle("Mean PMD Gene Number") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
dev.off()

library(GenomicFeatures)
txdb = loadDb("/dcl01/lieber/ajaffe/Amanda/annotation_objects/gencode.v25lift37.annotation.sqlite")
promoters = promoters(txdb, upstream=2000, downstream=200)
oo = lapply(pmdDMR, function(x) findOverlaps(promoters, makeGRangesFromDataFrame(x))) 
for (i in 1:length(pmdDMR)) {
pmdDMR[[i]][subjectHits(oo[[i]]),"promoters"] = "promoters"
pmdDMR[[i]][-subjectHits(oo[[i]]),"promoters"] = "no"
}
lapply(pmdDMR, function(x) unique(x$promoters))

pmdDMRdt = lapply(pmdDMR, data.table)
promoverlap = lapply(pmdDMR, function(y) round(length(unique(as.character(data.frame(y)[which(y$promoters=="promoters"),"regionID"])))/
                                                 length(unique(as.character(y$regionID)))*100,2))
promoverlap = data.frame(promoverlap = unlist(promoverlap), id = names(promoverlap))
promoverlap[!(promoverlap$id %in% postnatalpd$Data.ID),"celltype"] = "Prenatal"
promoverlap[which(promoverlap$id %in% postnatalpd[which(postnatalpd$Cell.Type=="Neuron"),"Data.ID"]),"celltype"] = "Neuron"
promoverlap[which(promoverlap$id %in% postnatalpd[which(postnatalpd$Cell.Type=="Glia"),"Data.ID"]),"celltype"] = "Glia"
write.csv(promoverlap, quote=F,file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/PMD_prom_Overlap.csv")

numproms = lapply(pmdDMRdt, function(y) y[promoters=="promoters",length(unique(nearestID)), by = "regionID"])
numproms = do.call(rbind, Map(cbind, lapply(numproms, function(y) 
  data.frame(mean = mean(y$V1), median = median(y$V1), sd = sd(y$V1), min = min(y$V1), max = max(y$V1))), id = as.list(names(numproms))))
numproms[!which(numproms$id %in% postnatalpd$Data.ID),"celltype"] = "Prenatal"
numproms[which(numproms$id %in% postnatalpd[postnatalpd$Cell.Type=="Neuron","Data.ID"]),"celltype"] = "Neuron"
numproms[which(numproms$id %in% postnatalpd[postnatalpd$Cell.Type=="Glia","Data.ID"]),"celltype"] = "Glia"
write.csv(numproms, quote=F, file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/PMD_numproms.csv")

w = numproms[which(numproms$celltype %in% c("Neuron", "Glia")),]
w = cbind(w, Age = postnatalpd[match(w$id, postnatalpd$Data.ID), "Age"])
cor.test(x = w[which(w$celltype=="Neuron"),"Age"], y = w[which(w$celltype=="Neuron"),"mean"])
#t = -2.953, df = 22, p-value = 0.007351
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
#  -0.7705711 -0.1648164
#sample estimates:
#  cor 
#-0.5327899 
cor.test(x = w[which(w$celltype=="Glia"),"Age"], y = w[which(w$celltype=="Glia"),"mean"])
#t = 0.31921, df = 6, p-value = 0.7604
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
#  -0.6330993  0.7642994
#sample estimates:
#  cor 
#0.1292242 

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/CREs/figures/PMD_overlap_with_proms.pdf")
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
  ggtitle("Median PMD Promoter Number by Age") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20), legend.title=element_blank()) + theme(legend.position="bottom")
ggplot(numproms, aes(x = celltype, y = median)) + geom_boxplot() +
  theme_classic() +
  labs(fill="") +
  ylab("Bases") + xlab("") +
  ggtitle("Median PMD Promoter Number") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
ggplot(numproms, aes(x = celltype, y = mean)) + geom_boxplot() +
  theme_classic() +
  labs(fill="") +
  ylab("Bases") + xlab("") +
  ggtitle("Mean PMD Promoter Number") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
dev.off()


## Identify the regions that are getting smaller or larger over age

load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/PMDs_annotated.rda")
pmdDMRdt = lapply(pmdDMR, data.table)
lapply(pmdDMRdt, function(x) x[,length(unique(regionID)),by="tog2"])

subset = do.call(rbind, Map(cbind, lapply(pmdDMRdt, function(x) data.frame(regionID = x$regionID, width = x$width, tog2 = x$tog2)), id = as.list(names(pmdDMRdt))))
subset[!which(subset$id %in% postnatalpd$Data.ID),"celltype"] = "Prenatal"
subset[which(subset$id %in% postnatalpd[postnatalpd$Cell.Type=="Neuron","Data.ID"]),"celltype"] = "Neuron"
subset[which(subset$id %in% postnatalpd[postnatalpd$Cell.Type=="Glia","Data.ID"]),"celltype"] = "Glia"
subset = unique(subset)
w = subset[which(subset$celltype %in% c("Neuron", "Glia")),]
w = cbind(w, Age = postnatalpd[match(w$id, postnatalpd$Data.ID), "Age"])
w = split(w, w$tog2)
data.frame(elementNROWS(w)[order(elementNROWS(w), decreasing=T)])
#no:no:no:no:no                                 9606
#All:Prenatal:Postnatal:Neurons:Glia            2017
#no:no:no:Neurons:no                            1945
#no:no:no:no:Glia                               1420
#no:Prenatal:no:no:Glia                          595
#no:Prenatal:no:no:no                            133
#no:no:Postnatal:Neurons:Glia                    129
#no:Prenatal:no:Neurons:Glia                       1

w = lapply(w, data.table)
mean = lapply(w, function(x) data.frame(x[,mean(width), by = c("id","celltype","Age")]))
med = lapply(w, function(x) data.frame(x[,median(width), by = c("id","celltype","Age")]))

nwidth = lapply(mean[1:7], function(z) cor.test(z[which(z$celltype=="Neuron"),"Age"], z[which(z$celltype=="Neuron"),"V1"]))
gwidth = lapply(mean[1:7], function(z) cor.test(z[which(z$celltype=="Glia"),"Age"], z[which(z$celltype=="Glia"),"V1"]))
mnwidth = lapply(med[1:7], function(z) cor.test(z[which(z$celltype=="Neuron"),"Age"], z[which(z$celltype=="Neuron"),"V1"]))
mgwidth = lapply(med[1:7], function(z) cor.test(z[which(z$celltype=="Glia"),"Age"], z[which(z$celltype=="Glia"),"V1"]))

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
          file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/shared_PMDs_correlation.csv")

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/CREs/figures/shared_PMD_width_correlation.pdf")
for (i in 1:length(mean)) {
  g = ggplot(mean[[i]], aes(x = Age, y = V1, colour = celltype)) + geom_point() +
    geom_smooth(method=lm, se=T, fullrange=TRUE) + scale_colour_brewer(8, palette="Dark2") + 
    labs(fill="") + theme_classic() +
    ylab("Width") + xlab("Age (Years)") +
    ggtitle(paste0("Mean Width of Overlapping PMDs:\n",names(w)[i])) +
    theme(title = element_text(size = 20)) +
    theme(text = element_text(size = 20), legend.title=element_blank()) + theme(legend.position="bottom")
  print(g)
  g = ggplot(med[[i]], aes(x = Age, y = V1, colour = celltype)) + geom_point() +
    geom_smooth(method=lm, se=T, fullrange=TRUE) + scale_colour_brewer(8, palette="Dark2") + 
    labs(fill="") + theme_classic() +
    ylab("Width") + xlab("Age (Years)") +
    ggtitle(paste0("Median Width of Overlapping PMDs:\n",names(w)[i])) +
    theme(title = element_text(size = 20)) +
    theme(text = element_text(size = 20), legend.title=element_blank()) + theme(legend.position="bottom")
  print(g)
}
dev.off()


## Gene ontology of genes in different categories

p = lapply(pmdDMR, function(x) data.frame(x[,1:5], tog2 = x$tog2, dmr = x$dmr, regionID = x$regionID))
p = lapply(p, function(x) unique(x))
oo = lapply(p, function(x) findOverlaps(makeGRangesFromDataFrame(geneMap), makeGRangesFromDataFrame(x))) 
p = mapply(function(q, ov) if (nrow(q)>length(unique(subjectHits(ov)))) {
  rbind(data.frame(q[subjectHits(ov),], nearestSymbol = geneMap[queryHits(ov),"Symbol"],
                   nearestID = geneMap[queryHits(ov),"gencodeID"], 
                   EntrezID = geneMap[queryHits(ov),"EntrezID"]),
        data.frame(q[-unique(subjectHits(ov)),], nearestSymbol = "NoGeneOverlap", 
                   nearestID = "NoGeneOverlap", EntrezID = "NoGeneOverlap")) } else {
                     data.frame(q[subjectHits(ov),], nearestSymbol = geneMap[queryHits(ov),"Symbol"],
                                nearestID = geneMap[queryHits(ov),"gencodeID"], EntrezID = geneMap[queryHits(ov),"EntrezID"])      
                   },p, oo, SIMPLIFY = F)

pmdDMRdt = lapply(p, data.table)
subset = do.call(rbind, Map(cbind, lapply(pmdDMRdt, function(x) data.frame(regionID = x$regionID, EntrezID = x$EntrezID, tog2 = x$tog2)),
                            id = as.list(names(pmdDMRdt))))
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
save(compareKegg.collapsed, compareBP.collapsed, compareMF, compareMF.collapsed, compareCC.collapsed, compareDO.collapsed,
     file="/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/PMD_KEGG_GO_DO_objects.rda")

# plot compared results
pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/CREs/figures/PMD_KEGG_GO_DO_plots.pdf", height = 60, width = 24)
plot(compareKegg.collapsed, colorBy="p.adjust", showCategory = 1000, title= "KEGG Pathway Enrichment")
plot(compareBP.collapsed, colorBy="p.adjust", showCategory = 1500, title= "Biological Process GO Enrichment")
lapply(compareMF, function(x) plot(x, colorBy="p.adjust", showCategory = 1000, title= "Molecular Function GO Enrichment"))
plot(compareMF.collapsed, colorBy="p.adjust", showCategory = 1000, title= "Molecular Function GO Enrichment")
plot(compareCC.collapsed, colorBy="p.adjust", showCategory = 1000, title= "Cellular Compartment GO Enrichment")
plot(compareDO.collapsed, colorBy="p.adjust", showCategory = 1000, title= "Disease Ontology Enrichment")
dev.off()

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/CREs/figures/PMD_BP_plot.pdf", height = 275, width = 24)
plot(compareBP.collapsed, colorBy="p.adjust", showCategory = 1500, title= "Biological Process GO Enrichment")
dev.off()


#### association with H3K27me3, H3K9me and other states

library(bumphunter)
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

pmdList = do.call(rbind, Map(cbind, lapply(pmdDMRdt, function(x) data.frame(x[,1:5], regionID = x$regionID, tog2 = x$tog2)), id = as.list(names(pmdDMRdt))))
pmdList[-which(pmdList$id %in% postnatalpd$Data.ID),"celltype"] = "Prenatal"
pmdList[which(pmdList$id %in% postnatalpd[postnatalpd$Cell.Type=="Neuron","Data.ID"]),"celltype"] = "Neuron"
pmdList[which(pmdList$id %in% postnatalpd[postnatalpd$Cell.Type=="Glia","Data.ID"]),"celltype"] = "Glia"
pmdList = unique(pmdList)
pmdList = unlist(lapply(split(pmdList, pmdList$tog2), function(x) split(x, x$celltype)), recursive = F)
pmdList = lapply(pmdList, makeGRangesFromDataFrame, keep=TRUE)
pmdList = GRangesList(lapply(pmdList, reduce))

bgWithPMD = endoapply(pmdList, function(x) {
  x = granges(x)
  m = c(x, gr.clusters)
  d = disjoin(m)
  d$inPMD = countOverlaps(d,x) > 0
  return(d)
})


#### split by chromosome 
PMDsByChr = lapply(bgWithPMD, function(x) split(x, seqnames(x)))

## autosomal
PMDs_bpOverlapByChr = lapply(PMDsByChr, function(x) mapply(function(regByChr, stateByChr) {
  cat(".")
  inPMD = regByChr[regByChr$inPMD]
  outPMD = regByChr[!regByChr$inPMD]
  
  enr = stateByChr[ranges(inPMD),]
  bg = stateByChr[ranges(outPMD),]
  
  list(enrTab = sapply(enr, table),
       bgTab  = sapply(bg, table))
}, x[paste0("chr",1:22)], stateCovDf, SIMPLIFY=FALSE))

## states
PMDs_bpOverlapEnr = lapply(PMDs_bpOverlapByChr, function(x) sapply(x, "[", "enrTab"))
bpOverlapArray = lapply(PMDs_bpOverlapEnr, function(x) array(unlist(x), dim = c(nrow(x[[1]]), ncol(x[[1]]), length(x))))
PMDs_enrTab = lapply(bpOverlapArray, function(x) apply(x, 1:2, sum))
for (i in 1:length(PMDs_enrTab)) {
  dimnames(PMDs_enrTab[[i]]) = list(rownames(PMDs_bpOverlapEnr[[i]][[1]]), colnames(PMDs_bpOverlapEnr[[i]][[1]]))
}

# states in rest of genome	
PMDs_bpOverlapBg = lapply(PMDs_bpOverlapByChr, function(x) sapply(x, "[", "bgTab"))
bpOverlapArray = lapply(PMDs_bpOverlapBg, function(x) array(unlist(x), dim = c(nrow(x[[1]]), ncol(x[[1]]), length(x))))
PMDs_bgTab = lapply(bpOverlapArray, function(y) apply(y, 1:2, function(x) sum(as.numeric(x))))
for (i in 1:length(PMDs_bgTab)) {
  dimnames(PMDs_bgTab[[i]]) = list(rownames(PMDs_bpOverlapBg[[i]][[1]]), colnames(PMDs_bpOverlapBg[[i]][[1]]))
}

# take ratio
PMDs_statTab = mapply(function(enr,bg) as.data.frame(t(prop.table(enr,2) / prop.table(bg,2))), PMDs_enrTab, PMDs_bgTab, SIMPLIFY = F)
PMDs_statTab = Map(cbind, PMDs_statTab, Sample = lapply(PMDs_statTab, function(x) dat$Standardized.Epigenome.name[match(rownames(x), dat$EID)]))
PMDs_statTab = lapply(PMDs_statTab, function(x) x[,c(16,1:15)])
names(PMDs_statTab) = c("A:Pr:Po:N:G-Glia","A:Pr:Po:N:G-Neuron","A:Pr:Po:N:G-Prenatal",":::N:-Glia",":::N:-Neuron",":::N:-Prenatal","::::G-Glia",
                        "::::G-Neuron","::::G-Prenatal","::::-Glia","::::-Neuron","::::-Prenatal","::Po:N:G-Glia","::Po:N:G-Neuron","::Po:N:G-Prenatal",
                        ":Pr:::G-Glia",":Pr:::G-Neuron",":Pr:::G-Prenatal",":Pr:::-Glia",
                        ":Pr:::-Neuron",":Pr:::-Prenatal",":Pr::N:G-Neuron")
save(PMDs_statTab, file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/PMD_chromHMM_15state_enrichment.rda")

dlpfcStats = do.call(rbind, lapply(PMDs_statTab, function(x) x["E073",-1]))
rownames(dlpfcStats) = names(PMDs_statTab)
dlpfcStats = t(signif(dlpfcStats,3))

tmpList = lapply(PMDs_statTab, function(x) vector("list", nrow(x)))
for (j in 1:length(tmpList)) {
  for(i in 1:nrow(PMDs_statTab[[j]])) {
    tmp =  do.call(rbind, lapply(PMDs_statTab, function(x) x[i,-1]))
    rownames(tmp) = names(PMDs_statTab)
    tmpList[[j]][[i]] = t(signif(tmp,3)) / dlpfcStats
  }
  names(tmpList[[j]]) = rownames(PMDs_statTab[[j]])
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

atriumStats = do.call(rbind, lapply(PMDs_statTab, function(x) x["E104",-1]))
rownames(atriumStats) = names(PMDs_statTab)
atriumStats = t(signif(atriumStats,3))

## other tissues
otherList = list(cellDMRs_statTab, ageDMRs_statTab, intDMRs_statTab)
ageVsCell = ageDMRs_statTab[,-1] / cellDMRs_statTab[,-1]
cor(t(ageVsCell),t(ageVsCell["E073",]))


## Plot example PMDs

library("MethylSeekR")
library("BSgenome")
library("BSgenome.Hsapiens.UCSC.hg19")

load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/PMDs_methylSeekR.rda")
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/methylSeekR_objects_prenatal.rda")
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/methylSeekR_objects.rda")

names(total) = gsub("Postnatal-","", names(total))
names(total) = gsub("Prenatal-","", names(total))
pmds = total[-grep("Lister",names(total))]

table(names(c(CGlist,CGPrenlist))==names(pmds))

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/CREs/figures/PMDsegments_visualized_CG_100kb.pdf", width = 12, height = 12)
mapply(function(CG, PMD) plotPMDSegmentation(m = CG, segs = PMD), c(CGlist,CGPrenlist), pmds)
dev.off()





