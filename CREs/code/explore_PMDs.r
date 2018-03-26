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

grm1 = lapply(grm1, as.list)
grm2 = lapply(grm2, as.list)

sumCGgrm1 = lapply(grm1, function(x) lapply(x, function(y) sum(y$nCG)))
sumCGgrm2 = lapply(grm2, function(x) lapply(x, function(y) sum(y$nCG)))

grm1 = lapply(grm1, function(g) lapply(g, function(x) x[1]))
for (i in 1:length(grm1)) { for (j in 1:length(grm1[[i]])) { grm1[[i]][[j]]$nCG = sumCGgrm1[[i]][[j]] } }
grm2 = lapply(grm2, function(g) lapply(g, function(x) x[1]))
for (i in 1:length(grm2)) { for (j in 1:length(grm2[[i]])) { grm2[[i]][[j]]$nCG = sumCGgrm2[[i]][[j]] } }

grm1 = lapply(grm1, function(x) do.call(getMethod(c, "GenomicRanges"), GRangesList(x)))
grm2 = lapply(grm2, function(x) do.call(getMethod(c, "GenomicRanges"), GRangesList(x)))

merge = mapply(function(gr1, gr2) reduce(c(gr1,gr2)), nonpmds, newnon)
oo1 = mapply(function(m, gr1) findOverlaps(m, gr1), merge, grm1)
for (i in 1:length(merge)) { 
  merge[[i]]$type = "notPMD"
  if (length(subjectHits(oo1[[i]]))==length(merge[[i]])) {
    tmp = merge[[i]][queryHits(oo1[[i]])]
    tmp$nCG1 = grm1[[i]][subjectHits(oo1[[i]])]$nCG
    merge[[i]] = tmp } else {
      tmp = merge[[i]][queryHits(oo1[[i]])]
      tmp1 = merge[[i]][-unique(queryHits(oo1[[i]]))]
      tmp$nCG1 = grm1[[i]][subjectHits(oo1[[i]])]$nCG
      tmp1$nCG1 = 0
      merge[[i]] = c(tmp, tmp1)
    }
}
oo2 = mapply(function(m, gr2) findOverlaps(m, gr2), merge, grm2)
for (i in 1:length(merge)) { 
  tmp = merge[[i]][queryHits(oo2[[i]])]
  tmp2 = merge[[i]][-unique(queryHits(oo2[[i]]))]
  tmp$nCG2 = grm2[[i]][subjectHits(oo2[[i]])]$nCG
  tmp2$nCG2 = 0
  merge[[i]] = c(tmp, tmp2)
}

for (i in 1:length(merge)) { 
  merge[[i]]$nCG = rowSums(as.data.frame(mcols(merge[[i]])[,colnames(mcols(merge[[i]])) %in% c("nCG1","nCG2")]))
  mcols(merge[[i]]) = mcols(merge[[i]])[,colnames(mcols(merge[[i]])) %in% c("type","nCG")]
}

pmds = lapply(pmds, function(x) x[which(width(x)>=100000)])

total = mapply(function(non, pmd) c(non, pmd), merge, pmds)
total = lapply(total, function(x) x[order(x)])

table(unlist(lapply(c(PMDsegments.CG.sortedLister, PMDsegments.CG, PMDsegments.CGpren), function(x) sum(x$nCG)))==18664892)
table(unlist(lapply(total, function(x) sum(x$nCG)))==18664892)

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


## Is this independent of centromeres, telomeres?

centtelo = read.table("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/centromeres_telomeres_hg19.txt", header=T)
centtelo = makeGRangesFromDataFrame(centtelo, keep.extra.columns = T)

coo = lapply(pmds, function(x) findOverlaps(x, centtelo[which(centtelo$type=="centromere")]))
too = lapply(pmds, function(x) findOverlaps(x, centtelo[which(centtelo$type=="telomere")]))
soo = lapply(pmds, function(x) findOverlaps(x, centtelo[which(centtelo$type=="short_arm")]))
aoo = lapply(pmds, function(x) findOverlaps(x, centtelo[which(centtelo$type %in% c("clone","contig"))]))
pmd.centtelo = data.frame(sampleID = names(pmds), num.PMDS = elementNROWS(pmds), 
                          num.centromere = unlist(lapply(coo, function(x) length(unique(queryHits(x))))),
                          num.telomere = unlist(lapply(too, function(x) length(unique(queryHits(x))))),
                          num.shortarm = unlist(lapply(soo, function(x) length(unique(queryHits(x))))),
                          num.clone.contig = unlist(lapply(aoo, function(x) length(unique(queryHits(x))))))
pmd.centtelo$perc.centromere = round(pmd.centtelo$num.centromere / pmd.centtelo$num.PMDS * 100, 1)
pmd.centtelo$perc.telomere = round(pmd.centtelo$num.telomere / pmd.centtelo$num.PMDS * 100, 1)
pmd.centtelo$perc.shortarm = round(pmd.centtelo$num.shortarm / pmd.centtelo$num.PMDS * 100, 1)
pmd.centtelo$perc.clone.contig = round(pmd.centtelo$num.clone.contig / pmd.centtelo$num.PMDS * 100, 1)
pmd.centtelo$sum.perc = rowSums(pmd.centtelo[,which(colnames(pmd.centtelo) %in% c("perc.centromere","perc.telomere","perc.shortarm","perc.clone.contig"))])
range(pmd.centtelo$sum.perc) # 10.8 - 39.2% of PMDs identified overlap with a gap from the genome table browser

## To be safe, remove these going forward

oo = lapply(pmds, function(x) findOverlaps(x, centtelo[which(centtelo$type != "heterochromatin")]))
pmds = mapply(function(p,o) p[-queryHits(o)], pmds, oo, SIMPLIFY = F)


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
#t = 1.9391, df = 22, p-value = 0.06543, cor 0.3820491

cor.test(x = w[which(w$celltype=="Glia"),"Age"], y = w[which(w$celltype=="Glia"),"num"])
#t = -0.367, df = 6, p-value = 0.7262, cor -0.1481735 

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

num = read.csv("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/PMD_number.csv")
range(num$num[which(num$celltype=="Prenatal")]) # 237 548
range(num$num[which(num$celltype=="Neuron")]) # 245 404
range(num$num[which(num$celltype=="Glia")]) # 471 1013

t.test(num[num$celltype=="Prenatal","num"],num[num$celltype!="Prenatal","num"])
#t = -0.82741, df = 47.026, p-value = 0.4122
#mean of x mean of y 
#403.0500  433.7812 

t.test(num[num$celltype=="Glia","num"],num[num$celltype!="Glia","num"])
#t = 4.9442, df = 7.3556, p-value = 0.001444
#  mean of x mean of y 
#698.2500  371.7273 


## how long are these regions?

widths = do.call(rbind, Map(cbind, lapply(pmds, function(y) data.frame(mean = mean(width(y)), median = median(width(y)), 
                                             sd = sd(width(y)), min = min(width(y)), max = max(width(y)))), id = as.list(names(pmds))))
widths[which(widths$id %in% names(PMDsegments.CGpren)),"celltype"] = "Prenatal"
widths[which(widths$id %in% pd[pd$Cell.Type=="Neuron","Data.ID"]),"celltype"] = "Neuron"
widths[which(widths$id %in% pd[pd$Cell.Type=="Glia","Data.ID"]),"celltype"] = "Glia"
write.csv(widths, quote=F, file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/PMD_widths.csv")

widths = read.csv("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/PMD_widths.csv")
range(widths$median[which(widths$celltype=="Prenatal")]) # 130562.5 138176.0
range(widths$median[which(widths$celltype=="Neuron")]) # 136230 158870
range(widths$median[which(widths$celltype=="Glia")]) # 147797 164566


pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/CREs/figures/PMD_widths.pdf")

w = widths[which(widths$celltype %in% c("Neuron", "Glia")),]
w = cbind(w, Age = pd[match(w$id, pd$Data.ID), "Age"])
cor.test(x = w[which(w$celltype=="Neuron"),"Age"], y = w[which(w$celltype=="Neuron"),"median"])
#t = 3.3111, df = 22, p-value = 0.003177, cor 0.576706 

cor.test(x = w[which(w$celltype=="Glia"),"Age"], y = w[which(w$celltype=="Glia"),"median"])
#t = 0.52835, df = 6, p-value = 0.6162, cor 0.2108469 

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

t.test(widths[widths$celltype=="Prenatal","median"],widths[widths$celltype!="Prenatal","median"])
#t = -13.579, df = 40.241, p-value < 2.2e-16
#  mean of x mean of y 
#134570.0  151129.1 
t.test(widths[widths$celltype=="Prenatal","median"],widths[widths$celltype=="Neuron","median"])
#t = -12.42, df = 30.638, p-value = 1.709e-13
#  mean of x mean of y 
#134570.0  149289.4 
t.test(widths[widths$celltype=="Glia","median"],widths[widths$celltype=="Neuron","median"])
#t = 2.9443, df = 10.524, p-value = 0.01394
#  mean of x mean of y 
#156648.0  149289.4 
t.test(widths[widths$celltype=="Glia","median"],widths[widths$celltype=="Prenatal","median"])
#t = 9.6256, df = 7.5883, p-value = 1.615e-05
#  mean of x mean of y 
#156648    134570 


## What proportion of the genome do they cover?

sLengths = c(lapply(PMDsegments.CG, reduce), lapply(PMDsegments.CGpren, reduce))
percGenome = data.frame(percent = unlist(mapply(function(p, t) round((sum(width(p))/sum(as.numeric(width(t)))*100),2), pmds, sLengths)),
                        id = names(pmds))
percGenome[which(percGenome$id %in% names(PMDsegments.CGpren)),"celltype"] = "Prenatal"
percGenome[which(percGenome$id %in% pd[pd$Cell.Type=="Neuron","Data.ID"]),"celltype"] = "Neuron"
percGenome[which(percGenome$id %in% pd[pd$Cell.Type=="Glia","Data.ID"]),"celltype"] = "Glia"
write.csv(percGenome, quote=F, file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/PMD_genomeCoverage.csv")

range(percGenome[which(percGenome$celltype=="Neuron"),"percent"]) # 1.47 3.30

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/CREs/figures/PMD_percentGenome.pdf")

w = percGenome[which(percGenome$celltype %in% c("Neuron", "Glia")),]
w = cbind(w, Age = pd[match(w$id, pd$Data.ID), "Age"])
cor.test(x = w[which(w$celltype=="Neuron"),"Age"], y = w[which(w$celltype=="Neuron"),"percent"])
#t = 0.90213, df = 22, p-value = 0.3768, cor 0.1888729 

cor.test(x = w[which(w$celltype=="Glia"),"Age"], y = w[which(w$celltype=="Glia"),"percent"])
#t = -0.1845, df = 6, p-value = 0.8597, cor -0.07510928 

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

genome = data.table(read.csv("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/PMD_genomeCoverage.csv"))
genome[,mean(percent), by = c("celltype")]
#celltype     V1
#1:     Glia 4.8425
#2:   Neuron 2.2350
#3: Prenatal 2.4855

t.test(genome[genome$celltype=="Prenatal","percent"],genome[genome$celltype!="Prenatal","percent"])
#t = -1.5283, df = 38.828, p-value = 0.1345
#  mean of x mean of y 
#2.485500  2.886875 
t.test(genome[genome$celltype=="Prenatal","percent"],genome[genome$celltype=="Neuron","percent"])
#t = 1.9502, df = 41.617, p-value = 0.05791
#  mean of x mean of y 
#2.4855    2.2350 
t.test(genome[genome$celltype=="Glia","percent"],genome[genome$celltype=="Neuron","percent"])
#t = 4.9806, df = 7.4389, p-value = 0.001338
#  mean of x mean of y 
#4.8425    2.2350 
t.test(genome[genome$celltype=="Prenatal","percent"],genome[genome$celltype=="Glia","percent"])
#t = -4.5021, df = 7.4387, p-value = 0.002397
#  mean of x mean of y 
#2.4855    4.8425 


### How much of the PMDs overlap?

# Master list of all bases in each setting
pmaster = reduce(do.call(getMethod(c, "GenomicRanges"), GRangesList(pmds))) # 2874 total
c(mean=mean(width(pmaster)), median=median(width(pmaster)), sd=sd(width(pmaster)), min=min(width(pmaster)), max=max(width(pmaster)))/1000
#       mean     median         sd        min        max 
#kb:239.7188   160.3775   343.9905   100.0460 11078.1050 


## regions that are shared by samples

all = list(All = Reduce(intersect, pmds),
           Prenatal = Reduce(intersect, lapply(PMDsegments.CGpren, function(x) x[which(x$type=="PMD" & width(x)>=100000)])),
           Postnatal = Reduce(intersect, lapply(PMDsegments.CG, function(x) x[which(x$type=="PMD" & width(x)>=100000)])),
           Neurons = Reduce(intersect, lapply(PMDsegments.CG, function(x) x[which(x$type=="PMD" & width(x)>=100000)])
                            [which(names(PMDsegments.CG) %in% pd[pd$Cell.Type=="Neuron","Data.ID"])]),
           Glia = Reduce(intersect, lapply(PMDsegments.CG, function(x) x[which(x$type=="PMD" & width(x)>=100000)])
                                           [which(names(PMDsegments.CG) %in% pd[pd$Cell.Type=="Glia","Data.ID"])]))
elementNROWS(all)
#All  Prenatal Postnatal   Neurons      Glia 
# 17       137        98       189       277 

round((unlist(lapply(all, function(x) sum(width(x))))/sum(width(pmaster))*100),2)
#      All  Prenatal Postnatal   Neurons      Glia 
#     0.38     40.30     38.12     39.95     44.95 
# 38-45% of total bases in the PMD state are shared by all in a cell type but not across all

shared = mapply(function(all,ind) lapply( ind, function(x) round(all / x *100,2)), lapply(all, function(x) sum(width(x))), 
                 lapply( list(pmds, 
                              lapply(PMDsegments.CGpren, function(x) x[which(x$type=="PMD" & width(x)>=100000)]), 
                              lapply(PMDsegments.CG, function(x) x[which(x$type=="PMD" & width(x)>=100000)]),
                              lapply(PMDsegments.CG, function(x) x[which(x$type=="PMD" & width(x)>=100000)])[which(names(PMDsegments.CG) %in% pd[pd$Cell.Type=="Neuron","Data.ID"])],
                              lapply(PMDsegments.CG, function(x) x[which(x$type=="PMD" & width(x)>=100000)])[which(names(PMDsegments.CG) %in% pd[pd$Cell.Type=="Glia","Data.ID"])]),
                         function(y) lapply(y, function(x) sum(width(x)))), SIMPLIFY = F)
shared = do.call(rbind, Map(data.frame, percent = lapply(shared, function(x) unlist(x, recursive=F)), id = lapply(shared, names), 
                            group = as.list(names(shared))))

write.table(shared, quote = F, file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/PMD_shared_byGroup.csv")


pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/CREs/figures/percent_shared_PMDs_perSample_byBaseCoverage.pdf")
ggplot(shared, aes(x = group, y = percent)) + geom_boxplot() +
  labs(fill="") + theme_classic() +
  ylim(0,100) +
  ylab("Percent") + xlab("") +
  ggtitle("Percent Bases Shared Per Sample: PMD") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
dev.off()

data.table(shared)[,mean(percent), by = c("group")]
#group        V1
#1:       All  3.464231
#2:  Prenatal 65.415000
#3: Postnatal 68.937187
#4:   Neurons 76.955417
#5:      Glia 65.411250


## How much of prenatal PMDs is represented in each postnatal sample?
  
pren = reduce(do.call(getMethod(c, "GenomicRanges"), GRangesList(pmds[-which(names(pmds) %in% pd$Data.ID)])))

sharedP = lapply(pmds, function(x) Reduce(intersect, list(x, pren)))

prenpercP = mapply(function(pren,ind) round(sum(width(reduce(pren))) / sum(width(reduce(ind))) *100,2), 
                   sharedP, pmds, SIMPLIFY = F)

prenpercP = data.frame(perc = unlist(prenpercP))
prenpercP$Person = rownames(prenpercP)
prenpercP$Age = pd[match(prenpercP$Person, pd$Data.ID),"Age"]
prenpercP$CellType = ifelse(prenpercP$Person %in% pd$Data.ID, pd[match(prenpercP$Person, pd$Data.ID),"Cell.Type"], "Prenatal")

write.table(prenpercP, quote = F, file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/PMD_sharedWithPrenatal.csv")

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/CREs/figures/percent_shared_withPrenatal_PMD_perSample_byBaseCoverage.pdf", width = 4,height = 4)
ggplot(prenpercP[prenpercP$CellType!="Prenatal",], aes(x = CellType, y = perc)) + geom_boxplot() +
  labs(fill="") + theme_classic() +
  ylim(0,100) +
  ylab("Percent") + xlab("") +
  ggtitle("Percent Bases Shared\nwith Prenatal") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
ggplot(prenpercP[prenpercP$CellType!="Prenatal",], aes(x = Age, y = perc, colour = CellType)) + 
  geom_path() + geom_point() + ylim(0,100) +
  theme_classic() + scale_colour_brewer(8, palette="Dark2") +
  ylab("Percent") + xlab("Age") + 
  ggtitle("Percent Bases Shared\nwith Prenatal") + 
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20), legend.title=element_blank()) + theme(legend.position="bottom")
dev.off()

cor.test(x = prenpercP[prenpercP$CellType=="Neuron",]$Age, y = prenpercP[prenpercP$CellType=="Neuron",]$perc)
#t = -6.2697, df = 22, p-value = 2.61e-06, cor -0.800725 
cor.test(x = prenpercP[prenpercP$CellType=="Glia",,]$Age, y = prenpercP[prenpercP$CellType=="Glia",,]$perc)
#t = -0.23369, df = 6, p-value = 0.823 cor -0.09497095 

data.table(prenpercP)[,mean(perc),by="CellType"]
#   CellType       V1
#1:     Glia  57.4575
#2:   Neuron  42.6600
#3: Prenatal 100.0000
(57.4575-42.6600)/57.4575 # 25.8%


## Find overlaps with DMRs, features, shared groups

islands = read.table("/dcl01/lieber/ajaffe/Amanda/annotation_objects/cpgIslandExt.hg19.txt", sep="\t", header = T)
rpmsk = read.table("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/RepeatMasker_genomewide_HG19.txt", header=T)
features = c(genes = makeGRangesFromDataFrame(geneMap, keep.extra.columns = T), 
             rpmskgr = makeGRangesFromDataFrame(rpmsk, seqnames.field="genoName",start.field="genoStart",end.field="genoEnd",keep.extra.columns=TRUE),
             islands = makeGRangesFromDataFrame(islands, keep.extra.columns = T, start.field = "chromStart", end.field = "chromEnd"))
lapply(features, head)

DMR = lapply(DMR, function(x) x[which(x$sig=="FWER < 0.05"),])
DMRgr = lapply(DMR, function(x) makeGRangesFromDataFrame(x, keep.extra.columns = T))

load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/bumphunting/rda/limma_Neuron_CpGs_minCov_3_ageInfo_dmrs.Rdata")
dmrs = split(dmrs, dmrs$k6cluster_label)
names(dmrs) = c("Gr1","Gr2","Gr3","Gr4","Gr5","Gr6")

oo = lapply(pmds, function(x) lapply(c(DMRgr[names(DMRgr) %in% c("CellType","Age")], as.list(dmrs), all, features), function(d) findOverlaps(x,d)))

pmdDMR = lapply(pmds, as.data.frame)
pmdDMR = lapply(pmdDMR, function(x) data.frame(x, rnum = 1:nrow(x)))
pmdDMR = mapply(function(d,oo) data.frame(d, CT = ifelse(d$rnum %in% queryHits(oo$CellType), "CT", "no"),
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
                                        islands = ifelse(d$rnum %in% queryHits(oo$islands), "CpG-Island", "non-Island")), pmdDMR, oo, SIMPLIFY = F) 
pmdDMR = lapply(pmdDMR, function(d) data.frame(d, tog = paste(d$CT, d$Age, d$Gr1,d$Gr2,d$Gr3,d$Gr4,d$Gr5,d$Gr6, d$All, d$Prenatal, 
                                                          d$Postnatal, d$Neurons, d$Glia, sep = ":"),
                                           dmr = paste(d$CT, d$Age, d$Gr1,d$Gr2,d$Gr3,d$Gr4,d$Gr5,d$Gr6, sep = ":"),
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
save(pmdDMR, postnatalpd, geneMap, file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/PMDs_annotated.rda")


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

cgi = data.table(read.csv("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/PMD_CpG_Island_Overlap.csv"))
cgi[,mean(islands), by = "celltype"]
#   celltype       V1
#1:     Glia 18.34250
#2:   Neuron 40.67375
#3: Prenatal 39.84250
mean(cgi$islands) # 36.91846


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

repeats = data.table(read.csv("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/PMD_repeats_Overlap.csv"))
repeats[repeats %in% c("No repeats","SINE","Low_complexity","Simple_repeat","LINE"),mean(perc), by = c("repeats","celltype")]
#          repeats celltype       V1
#1:           SINE     Glia 14.06125
#2:  Simple_repeat     Glia 14.04250
#3: Low_complexity     Glia 14.04000
#4:           LINE     Glia 14.06250
#5:           SINE   Neuron 13.86250
#6:           LINE   Neuron 13.86250
#7: Low_complexity   Neuron 13.85167
#8:  Simple_repeat   Neuron 13.86250
#9:           LINE Prenatal 14.14800
#10: Low_complexity Prenatal 14.12450
#11:           SINE Prenatal 14.14450
#12:  Simple_repeat Prenatal 14.13850


## Overlap with DMRs

dmr = lapply(pmdDMRdt, function(y) y[,length(unique(regionID)), by = "dmr"])
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

dmr[!(dmr$id %in% postnatalpd$Data.ID),"celltype"] = "Prenatal"
dmr[which(dmr$id %in% postnatalpd[which(postnatalpd$Cell.Type=="Neuron"),"Data.ID"]),"celltype"] = "Neuron"
dmr[which(dmr$id %in% postnatalpd[which(postnatalpd$Cell.Type=="Glia"),"Data.ID"]),"celltype"] = "Glia"
models[!(models$id %in% postnatalpd$Data.ID),"celltype"] = "Prenatal"
models[which(models$id %in% postnatalpd[which(postnatalpd$Cell.Type=="Neuron"),"Data.ID"]),"celltype"] = "Neuron"
models[which(models$id %in% postnatalpd[which(postnatalpd$Cell.Type=="Glia"),"Data.ID"]),"celltype"] = "Glia"

dmr$NO = ifelse(dmr$dmr=="no:no:no:no:no:no:no:no", "No overlap", "Overlap")

write.csv(dmr, quote=F,file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/PMD_dmr_Overlap.csv")
write.csv(models, quote=F,file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/PMD_dmr_Overlap_byModel.csv")


pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/CREs/figures/PMD_overlap_with_dmr.pdf", width = 14, height = 6)
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

dmr = data.table(read.csv("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/PMD_dmr_Overlap.csv"))
models = data.table(read.csv("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/PMD_dmr_Overlap_byModel.csv"))
models[,list(Gr1=mean(Gr1),Gr2=mean(Gr2),Gr3=mean(Gr3),Gr4=mean(Gr4),Gr5=mean(Gr5),Gr6=mean(Gr6)),by="celltype"]
#   celltype     Gr1      Gr2      Gr3       Gr4      Gr5      Gr6
#1:     Glia 2.03875 1.555000 0.500000 0.0800000 0.246250 1.210000
#2:   Neuron 4.03000 1.874583 5.615417 0.7433333 2.807917 2.345833
#3: Prenatal 2.80850 3.211500 1.211500 0.0540000 0.766000 1.937500


# Are DMRs overlapping PMDs at the PMD periphery, or uniformly distributed?

pmdDMR = Map(cbind, lapply(pmdDMR, function(x) data.frame(x[,1:5], dmr = x$dmr, regionID = x$regionID)), id = as.list(names(pmdDMR)))
pmdDMR = lapply(pmdDMR, unique)
tiles = lapply(pmdDMR, function(x) tile(makeGRangesFromDataFrame(x[which(x$dmr!="no:no:no:no:no:no:no:no"),]), n = 100))
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

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/CREs/figures/DMR_position_in_PMD.pdf", width = 10)
ggplot(overlaps, aes(region, colour = model)) + geom_density() + theme_classic() +
  labs(fill="") + facet_grid(. ~ celltype) +
  ylab("Density") + xlab("% of PMD") +
  ggtitle("DMR Position in PMD") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20), legend.position="bottom", legend.title=element_blank())
dev.off()


## Overlap with Lister samples

load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/PMDs_methylSeekR_sortedLister.rda")
listerpd = c("GSM1173773" = "53yr_NeuN_pos", "GSM1173774" = "53yr_NeuN_neg",
             "GSM1173776" = "55yr_NeuN_pos", "GSM1173777" = "55yr_NeuN_neg")
listerpmds = lapply(PMDsegments.CG.sortedLister, function(x) x[which(x$type=="PMD" & width(x)>100000)])
names(listerpmds) = listerpd[match(names(listerpmds),names(listerpd))]
oo = lapply(listerpmds, function(x) findOverlaps(x, centtelo[which(centtelo$type != "heterochromatin")]))
listerpmds = mapply(function(p,o) p[-queryHits(o)], listerpmds, oo, SIMPLIFY = F)


# Total bases
celltype = list(LIBDneuron = pmdDMR[which(names(pmdDMR) %in% postnatalpd[postnatalpd$Cell.Type=="Neuron","Data.ID"])],
                LIBDglia = pmdDMR[which(names(pmdDMR) %in% postnatalpd[postnatalpd$Cell.Type=="Glia","Data.ID"])],
                LIBDprenatal = pmdDMR[-which(names(pmdDMR) %in% postnatalpd[postnatalpd$Cell.Type %in% c("Neuron","Glia"),"Data.ID"])])
celltype = lapply(celltype, function(y) lapply(y, makeGRangesFromDataFrame, keep=TRUE))
celltype = lapply(celltype, function(x) reduce(do.call(getMethod(c, "GenomicRanges"), GRangesList(x))))
lister = list(ListerNeuron = reduce(do.call(getMethod(c, "GenomicRanges"), GRangesList(listerpmds[grep("pos", names(listerpmds))]))),
              ListerGlia = reduce(do.call(getMethod(c, "GenomicRanges"), GRangesList(listerpmds[grep("neg", names(listerpmds))]))))

overlap = lapply(celltype, function(cell) lapply(lister, function(list) c( OurOverlap = round(sum(width(Reduce(intersect, list(cell,list))))/sum(width(cell))*100,1),
       ListerOverlap = round(sum(width(Reduce(intersect, list(cell,list))))/sum(width(list))*100,1))))
overlap = lapply(overlap, function(x) do.call(rbind, Map(cbind, lapply(x, function(y) data.frame(OurOverlap = y["OurOverlap"], ListerOverlap = y["ListerOverlap"])),
                                                         Lister = as.list(names(x)))))
overlap = do.call(rbind, Map(cbind, overlap, LIBD = as.list(names(celltype))))
overlap = reshape2::melt(overlap)
overlap$LIBD = gsub("LIBD", "", overlap$LIBD)
overlap$Lister = gsub("ListerN", "n", overlap$Lister)
overlap$Lister = gsub("ListerG", "g", overlap$Lister)
overlap$Lister = gsub("ListerP", "p", overlap$Lister)
overlap$variable = gsub("Our", "LIBD\n", overlap$variable)
overlap$variable = gsub("Lister", "Lister\n", overlap$variable)
colnames(overlap) = c("Lister", "LIBD", "Comparison", "Percent")
overlap$Lister = factor(overlap$Lister, levels = c("glia","neuron"))

write.csv(overlap, quote=F,file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/PMD_Lister_Overlap.csv")

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/CREs/figures/PMD_Lister_Overlap.pdf", width = 10)
ggplot(overlap, aes(x = Comparison, y = Percent, fill = Lister)) + 
  theme_classic() + geom_bar(stat = "identity") +
  labs(fill="") + facet_grid(Lister ~ LIBD) +
  ylab("Percent") + ylim(0,100) + xlab("") + 
  scale_fill_brewer(8, palette="Dark2") +
  ggtitle("Percent PMD Bases Replicated") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20), legend.title=element_blank()) + theme(legend.position="bottom")
dev.off()


## Do PMDs correspond to CpH deserts?

deserts = read.table("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CH_deserts_ListerEtal.txt", header=T)
deserts$Chr = paste0("chr", deserts$Chr)
deserts = makeGRangesFromDataFrame(deserts)

# percent PMDs overlapping deserts per sample
oo = lapply(pmds, function(x) findOverlaps(x, deserts))
desertstats = data.frame(id = names(pmds), perc.overlapping.PMD = unlist(mapply(function(p,ov) round(length(unique(queryHits(ov))) / length(p) * 100, 1), pmds, oo, SIMPLIFY = F)),
                         perc.bases.PMD = unlist(lapply(pmds, function(x) round(sum(width(Reduce(intersect, list(x, deserts))))/sum(width(x))*100,1))),
                         perc.overlapping.Desert = unlist(lapply(oo, function(ov) round(length(unique(subjectHits(ov))) / length(deserts) * 100, 1))),
                         perc.bases.Desert = unlist(lapply(pmds, function(x) round(sum(width(Reduce(intersect, list(x, deserts))))/sum(width(deserts))*100,1))))
desertstats[!which(desertstats$id %in% postnatalpd$Data.ID),"celltype"] = "Prenatal"
desertstats[which(desertstats$id %in% postnatalpd[postnatalpd$Cell.Type=="Neuron","Data.ID"]),"celltype"] = "Neuron"
desertstats[which(desertstats$id %in% postnatalpd[postnatalpd$Cell.Type=="Glia","Data.ID"]),"celltype"] = "Glia"
desertstats$age = postnatalpd[match(desertstats$id, postnatalpd$Data.ID),"Age"]

write.csv(desertstats, quote=F,file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/PMD_CpH_desert_Overlap.csv")


## association of existing PMDs with genes

oo = lapply(pmds, function(x) findOverlaps(makeGRangesFromDataFrame(geneMap), makeGRangesFromDataFrame(x))) 
pmdGene = mapply(function(p, ov) if (nrow(p)>length(unique(subjectHits(ov)))) {
                  rbind(data.frame(p[subjectHits(ov),], nearestSymbol = geneMap[queryHits(ov),"Symbol"],
                                                 nearestID = geneMap[queryHits(ov),"gencodeID"], 
                                                 EntrezID = geneMap[queryHits(ov),"EntrezID"]),
                                      data.frame(p[-unique(subjectHits(ov)),], nearestSymbol = "NoGeneOverlap", 
                                                 nearestID = "NoGeneOverlap", EntrezID = "NoGeneOverlap")) } else {
                        data.frame(p[subjectHits(ov),], nearestSymbol = geneMap[queryHits(ov),"Symbol"],
                                   nearestID = geneMap[queryHits(ov),"gencodeID"], EntrezID = geneMap[queryHits(ov),"EntrezID"])      
                                                 },lapply(pmds, as.data.frame), oo, SIMPLIFY = F)
save(pmdGene, postnatalpd, geneMap, file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/PMDs_geneOverlap.rda")

pmdDMRdt = lapply(pmdGene, data.table)
geneoverlap = lapply(pmdGene, function(y) round(length(unique(as.character(data.frame(y)[which(y$nearestID!="NoGeneOverlap"),"regionID"])))/
                                                 length(unique(as.character(y$regionID)))*100,2))
geneoverlap = data.frame(geneoverlap = unlist(geneoverlap), id = names(geneoverlap))
geneoverlap[!(geneoverlap$id %in% postnatalpd$Data.ID),"celltype"] = "Prenatal"
geneoverlap[which(geneoverlap$id %in% postnatalpd[which(postnatalpd$Cell.Type=="Neuron"),"Data.ID"]),"celltype"] = "Neuron"
geneoverlap[which(geneoverlap$id %in% postnatalpd[which(postnatalpd$Cell.Type=="Glia"),"Data.ID"]),"celltype"] = "Glia"
write.csv(geneoverlap, quote=F,file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/PMD_Gene_Overlap.csv")

numgenes = lapply(pmdDMRdt, function(y) y[,length(unique(nearestID)), by = "regionID"])
numgenes = do.call(rbind, Map(cbind, lapply(numgenes, function(y) 
  data.frame(mean = mean(y$V1), median = median(y$V1), sd = sd(y$V1), min = min(y$V1), max = max(y$V1))), id = as.list(names(numgenes))))
numgenes[-which(numgenes$id %in% postnatalpd$Data.ID),"celltype"] = "Prenatal"
numgenes$celltype = ifelse(numgenes$id %in% postnatalpd$Data.ID, postnatalpd[match(numgenes$id, postnatalpd$Data.ID),"Cell.Type"],"Prenatal")
numgenes[which(numgenes$id %in% postnatalpd[postnatalpd$Cell.Type=="Glia","Data.ID"]),"celltype"] = "Glia"
write.csv(numgenes, quote=F, file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/PMD_numgenes.csv")


w = numgenes[which(numgenes$celltype %in% c("Neuron", "Glia")),]
w = cbind(w, Age = postnatalpd[match(w$id, postnatalpd$Data.ID), "Age"])
cor.test(x = w[which(w$celltype=="Neuron"),"Age"], y = w[which(w$celltype=="Neuron"),"mean"])
#t = -2.0171, df = 22, p-value = 0.05605, cor -0.3950679 

cor.test(x = w[which(w$celltype=="Glia"),"Age"], y = w[which(w$celltype=="Glia"),"mean"])
#t = 0.69919, df = 6, p-value = 0.5106, cor 0.2744798 

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
  ylab("Number") + xlab("Age (Years)") +
  ggtitle("Median PMD Gene Number by Age") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20), legend.title=element_blank()) + theme(legend.position="bottom")
ggplot(numgenes, aes(x = celltype, y = median)) + geom_boxplot() +
  theme_classic() +
  labs(fill="") +
  ylab("Number") + xlab("") +
  ggtitle("Median PMD Gene Number") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
ggplot(numgenes, aes(x = celltype, y = mean)) + geom_boxplot() +
  theme_classic() +
  labs(fill="") +
  ylab("Number") + xlab("") +
  ggtitle("Mean PMD Gene Number") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
dev.off()

library(GenomicFeatures)
txdb = loadDb("/dcl01/lieber/ajaffe/Amanda/annotation_objects/gencode.v25lift37.annotation.sqlite")
promoters = promoters(txdb, upstream=2000, downstream=200)
promoters = reduce(promoters)
oo = lapply(pmdGene, function(x) findOverlaps(promoters, makeGRangesFromDataFrame(x))) 
for (i in 1:length(pmdGene)) {
  pmdGene[[i]][unique(subjectHits(oo[[i]])),"promoters"] = "promoters"
  pmdGene[[i]][-unique(subjectHits(oo[[i]])),"promoters"] = "no"
}

pmdGenedt = lapply(pmdGene, data.table)
pmdGenedt = Map(cbind, pmdGenedt, regionID = lapply(pmdGenedt, function(x) paste0(x$seqnames,":",x$start,"-",x$end)))
promoverlap = lapply(pmdGenedt, function(y) round(length(unique(as.character(data.frame(y)[which(y$promoters=="promoters"),"regionID"])))/
                                                 length(unique(as.character(y$regionID)))*100,2))
promoverlap = data.frame(promoverlap = unlist(promoverlap), id = numgenes$id)
promoverlap[!(promoverlap$id %in% postnatalpd$Data.ID),"celltype"] = "Prenatal"
promoverlap[which(promoverlap$id %in% postnatalpd[which(postnatalpd$Cell.Type=="Neuron"),"Data.ID"]),"celltype"] = "Neuron"
promoverlap[which(promoverlap$id %in% postnatalpd[which(postnatalpd$Cell.Type=="Glia"),"Data.ID"]),"celltype"] = "Glia"
write.csv(promoverlap, quote=F,file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/PMD_prom_Overlap.csv")

numproms = lapply(pmdGenedt, function(y) y[promoters=="promoters",length(unique(nearestID)), by = "regionID"])
numproms = do.call(rbind, Map(cbind, lapply(numproms, function(y) 
  data.frame(mean = mean(y$V1), median = median(y$V1), sd = sd(y$V1), min = min(y$V1), max = max(y$V1))), id = as.list(numgenes$id)))
numproms[-which(numproms$id %in% postnatalpd$Data.ID),"celltype"] = "Prenatal"
numproms[which(numproms$id %in% postnatalpd[postnatalpd$Cell.Type=="Neuron","Data.ID"]),"celltype"] = "Neuron"
numproms[which(numproms$id %in% postnatalpd[postnatalpd$Cell.Type=="Glia","Data.ID"]),"celltype"] = "Glia"
write.csv(numproms, quote=F, file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/PMD_numproms.csv")


## Overlap with promoters

w = numproms[which(numproms$celltype %in% c("Neuron", "Glia")),]
w = cbind(w, Age = postnatalpd[match(w$id, postnatalpd$Data.ID), "Age"])
cor.test(x = w[which(w$celltype=="Neuron"),"Age"], y = w[which(w$celltype=="Neuron"),"mean"])
#t = -2.4206, df = 22, p-value = 0.0242, cor -0.4586086 

cor.test(x = w[which(w$celltype=="Glia"),"Age"], y = w[which(w$celltype=="Glia"),"mean"])
#t = 0.46925, df = 6, p-value = 0.6555, cor o.1881491 


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
  ylab("Number") + xlab("Age (Years)") +
  ggtitle("Median PMD Promoter Number by Age") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20), legend.title=element_blank()) + theme(legend.position="bottom")
ggplot(numproms, aes(x = celltype, y = median)) + geom_boxplot() +
  theme_classic() +
  labs(fill="") +
  ylab("Number") + xlab("") +
  ggtitle("Median PMD Promoter Number") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
ggplot(numproms, aes(x = celltype, y = mean)) + geom_boxplot() +
  theme_classic() +
  labs(fill="") +
  ylab("Number") + xlab("") +
  ggtitle("Mean PMD Promoter Number") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
dev.off()


## SHSY5Y PMD genes in ours?

shsy5y = readxl::read_excel("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/PMD_genes_Schroeder2010_table3.xls")
shsy5y = as.data.frame(shsy5y)
shsy5y = list(Nhmd = na.omit(shsy5y$'N-HMD'), Lhmd = na.omit(shsy5y$"L-HMD"), Bpmd = na.omit(shsy5y$"B-PMD"))
shsy5y = lapply(shsy5y, function(x) geneMap[which(geneMap$Symbol %in% x),])

oo = lapply(pmds, function(p) lapply(shsy5y, function(sh) findOverlaps(p,makeGRangesFromDataFrame(sh))))
hits = lapply(oo, function(x) list(SHSY5Y = unlist(lapply(x, function(y) length(unique(subjectHits(y))))),
                                   PMDs = unlist(lapply(x, function(y) length(unique(queryHits(y)))))))
hits = do.call(rbind, Map(cbind, lapply(hits, function(x) do.call(rbind, Map(cbind, x, denominator = as.list(names(x))))), 
                          id = as.list(names(pmds)), num.pmds = as.list(elementNROWS(pmds))))
hits = data.frame(hits, group = rownames(hits))
hits$V1 = as.numeric(as.character(hits$V1))
hits$num.pmds = as.numeric(as.character(hits$num.pmds))
hits$percent = NA
hits[which(hits$denominator == "PMDs"),"percent"] = round(hits[which(hits$denominator == "PMDs"),"V1"]/hits[which(hits$denominator == "PMDs"),"num.pmds"] * 100,1)
hits[which(hits$denominator == "SHSY5Y" & hits$group == "Nhmd"),"percent"] = round(hits[which(hits$denominator == "SHSY5Y" & hits$group == "Nhmd"),"V1"]/2230 * 100,1)
hits[which(hits$denominator == "SHSY5Y" & hits$group == "Lhmd"),"percent"] = round(hits[which(hits$denominator == "SHSY5Y" & hits$group == "Lhmd"),"V1"]/613 * 100,1)
hits[which(hits$denominator == "SHSY5Y" & hits$group == "Bpmd"),"percent"] = round(hits[which(hits$denominator == "SHSY5Y" & hits$group == "Bpmd"),"V1"]/1706 * 100,1)
hits$celltype = ifelse(hits$id %in% pd$Data.ID, pd[match(hits$id, pd$Data.ID),"Cell.Type"], "Prenatal")
colnames(hits)[1] = "count"
write.csv(hits, quote=F,file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/PMD_SHSY5Ygenes_Overlap.csv")

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/CREs/figures/PMD_SHSY5Ygenes_Overlap.pdf")
ggplot(hits, aes(x = celltype, y = percent)) + geom_boxplot() +
  theme_classic() + facet_grid(denominator ~ group) +
  labs(fill="") +
  ylab("Percent") + ylim(0,50) + 
  xlab("") + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Percent Genes and PMDs Overlapping") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20), legend.title=element_blank()) + theme(legend.position="bottom")
dev.off()


## Enrichment in SHSY5Y PMD genes

geneuniverse = na.omit(unique(geneMap$gencodeID))
sets = lapply(shsy5y, function(x) na.omit(unique(x$gencodeID)))
inPMD = lapply(pmdGene, function(x) na.omit(unique(x$nearestID)))
outPMD = lapply(inPMD, function(x) geneuniverse[!(geneuniverse %in% x)])

PMDenrich = mapply(function(yes, no) lapply(sets, function(x) {
  PMD_OVERLAP = c( sum( yes %in% x),sum(!(yes %in% x)))
  NOT_PMD_OVERLAP= c(sum(no %in% x), sum(!(no %in% x)))
  enrich_table = cbind(PMD_OVERLAP, NOT_PMD_OVERLAP)
  res = fisher.test(enrich_table)
  dat=c(res$p.value, res$estimate)
  names(dat) <- c("P.Value","Odds.Ratio")
  return(dat)
}), inPMD, outPMD, SIMPLIFY =F)
PMDenrich = do.call(rbind, Map(cbind, lapply(PMDenrich, function(x) do.call(rbind, Map(cbind, lapply(x, function(y) data.frame(P.Value = y["P.Value"], Odds.Ratio = y["Odds.Ratio"])), 
                                                                                       SHSY5Y = as.list(names(x))))), id = as.list(names(PMDenrich))))
PMDenrich$celltype = ifelse(PMDenrich$id %in% pd$Data.ID, pd[match(PMDenrich$id, pd$Data.ID),"Cell.Type"], "Prenatal")
PMDenrich$FDR = p.adjust(PMDenrich$P.Value, method = "fdr")
PMDenrich$celltype = factor(PMDenrich$celltype, levels = c("Prenatal", "Neuron", "Glia"))
PMDenrich$SHSY5Y = gsub("Bpmd", "Both-PMD", PMDenrich$SHSY5Y)
PMDenrich$SHSY5Y = gsub("Nhmd","Neuron-HMD", PMDenrich$SHSY5Y)
PMDenrich$SHSY5Y = gsub("Lhmd", "Lung-HMD", PMDenrich$SHSY5Y)
write.csv(PMDenrich,file="/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/PMD_SHSY5Ygenes_enrichment.csv",quote=F)

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/CREs/figures/PMD_SHSY5Ygenes_Enrichment.pdf", height = 4, width = 8.5)
ggplot(PMDenrich[PMDenrich$FDR<=0.05,], aes(x = SHSY5Y, y = Odds.Ratio)) + geom_boxplot() +
  theme_classic() + facet_grid(. ~ celltype) +
  labs(fill="") +
  ylab("Odds Ratio") + ylim(0,10) +
  geom_hline(yintercept=1, linetype="dotted") +
  xlab("") + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Fisher Test for PMD Gene Sets (FDR<0.05)") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20), legend.title=element_blank()) + theme(legend.position="bottom")
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
#no:no:no:no:no                                           9243
#no:no:no:Neurons:no                                      1906
#no:no:no:no:Glia                                         1363
#All:Prenatal:Postnatal:Neurons:Glia                       494
#no:Prenatal:no:no:Glia                                    391
#no:Prenatal:Postnatal:Neurons:Glia                        267
#no:no:Postnatal:Neurons:Glia                              128
#no:Prenatal:no:no:no                                       89

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


## Do same comparison but of shared or absent genes in PMDs

inPMD = lapply(inPMD, as.character)
inPMD.byCT = list(Prenatal = unique(unlist(inPMD[-which(names(inPMD) %in% pd$Data.ID)])),
                  Neuron = unique(unlist(inPMD[which(names(inPMD) %in% pd[pd$Cell.Type=="Neuron","Data.ID"])])),
                  Glia = unique(unlist(inPMD[which(names(inPMD) %in% pd[pd$Cell.Type=="Glia","Data.ID"])])))

elementNROWS(inPMD.byCT)
# Prenatal   Neuron     Glia 
#     8373     3744     6200 

CTcomps = list(PnotN = inPMD.byCT$Prenatal[-which(inPMD.byCT$Prenatal %in% inPMD.byCT$Neuron)],
               NnotP = inPMD.byCT$Neuron[-which(inPMD.byCT$Neuron %in% inPMD.byCT$Prenatal)],
               sharedPN = inPMD.byCT$Prenatal[which(inPMD.byCT$Prenatal %in% inPMD.byCT$Neuron)],
               PnotG = inPMD.byCT$Prenatal[-which(inPMD.byCT$Prenatal %in% inPMD.byCT$Glia)],
               GnotP = inPMD.byCT$Glia[-which(inPMD.byCT$Glia %in% inPMD.byCT$Prenatal)],
               sharedPG = inPMD.byCT$Prenatal[which(inPMD.byCT$Prenatal %in% inPMD.byCT$Glia)],
               GnotN = inPMD.byCT$Glia[-which(inPMD.byCT$Glia %in% inPMD.byCT$Neuron)],
               NnotG = inPMD.byCT$Neuron[-which(inPMD.byCT$Neuron %in% inPMD.byCT$Glia)],
               sharedGN = inPMD.byCT$Glia[which(inPMD.byCT$Glia %in% inPMD.byCT$Neuron)],
               allshared = Reduce(intersect, inPMD.byCT))
# PnotN     NnotP  sharedPN     PnotG     GnotP  sharedPG     GnotN     NnotG 
#  6248      1619      2125      4358      2185      4015      4442      1986 
# sharedGN allshared 
#     1758      1576 

inPMD.byAge = list(infant = unique(unlist(inPMD[which(names(inPMD) %in% pd[pd$Cell.Type=="Neuron" & pd$Age<1,"Data.ID"])])),
                child = unique(unlist(inPMD[which(names(inPMD) %in% pd[pd$Cell.Type=="Neuron" & pd$Age>1 & pd$Age<=12,"Data.ID"])])), 
                teen = unique(unlist(inPMD[which(names(inPMD) %in% pd[pd$Cell.Type=="Neuron" & pd$Age>12 & pd$Age<=17,"Data.ID"])])),
                adult = unique(unlist(inPMD[which(names(inPMD) %in% pd[pd$Cell.Type=="Neuron" & pd$Age>17,"Data.ID"])])))
elementNROWS(inPMD.byAge)
# infant  child   teen  adult 
#   1896   2465   2561   2894 

Agecomps = list(InotC = inPMD.byAge$infant[-which(inPMD.byAge$infant %in% inPMD.byAge$child)],
               CnotI = inPMD.byAge$child[-which(inPMD.byAge$child %in% inPMD.byAge$infant)],
               sharedIC = inPMD.byAge$infant[which(inPMD.byAge$infant %in% inPMD.byAge$child)],
               InotT = inPMD.byAge$infant[-which(inPMD.byAge$infant %in% inPMD.byAge$teen)],
               TnotI = inPMD.byAge$teen[-which(inPMD.byAge$teen %in% inPMD.byAge$infant)],
               sharedIT = inPMD.byAge$infant[which(inPMD.byAge$infant %in% inPMD.byAge$teen)],
               InotA = inPMD.byAge$infant[-which(inPMD.byAge$infant %in% inPMD.byAge$adult)],
               AnotI = inPMD.byAge$adult[-which(inPMD.byAge$adult %in% inPMD.byAge$infant)],
               sharedIA = inPMD.byAge$infant[which(inPMD.byAge$infant %in% inPMD.byAge$adult)],
               CnotT = inPMD.byAge$child[-which(inPMD.byAge$child %in% inPMD.byAge$teen)],
               TnotC = inPMD.byAge$teen[-which(inPMD.byAge$teen %in% inPMD.byAge$child)],
               sharedCT = inPMD.byAge$child[which(inPMD.byAge$child %in% inPMD.byAge$teen)],
               CnotA = inPMD.byAge$child[-which(inPMD.byAge$child %in% inPMD.byAge$adult)],
               AnotC = inPMD.byAge$adult[-which(inPMD.byAge$adult %in% inPMD.byAge$child)],
               sharedCA = inPMD.byAge$child[which(inPMD.byAge$child %in% inPMD.byAge$adult)],
               TnotA = inPMD.byAge$teen[-which(inPMD.byAge$teen %in% inPMD.byAge$adult)],
               AnotT = inPMD.byAge$adult[-which(inPMD.byAge$adult %in% inPMD.byAge$teen)],
               sharedTA = inPMD.byAge$teen[which(inPMD.byAge$teen %in% inPMD.byAge$adult)])
elementNROWS(Agecomps)
# InotC    CnotI sharedIC    InotT    TnotI sharedIT    InotA    AnotI 
#   595     1164     1301      566     1231     1330      585     1583 
# sharedIA    CnotT    TnotC sharedCT    CnotA    AnotC sharedCA    TnotA 
#     1311      357      453     2108      264      693     2201      241 
# AnotT sharedTA 
#   574     2320 


## Gene ontology of genes in different categories

entrez = list(CellType = lapply(CTcomps, function(x) na.omit(unique(geneMap[which(geneMap$gencodeID %in% x),"EntrezID"]))),
              Age = lapply(Agecomps, function(x) na.omit(unique(geneMap[which(geneMap$gencodeID %in% x),"EntrezID"]))))


# Compare the enriched terms between 7 groups

compareKegg = lapply(entrez, function(x) compareCluster(x, fun="enrichKEGG", qvalueCutoff = 0.05, pvalueCutoff = 0.05))
compareBP = lapply(entrez, function(x) compareCluster(x, fun="enrichGO", ont = "BP", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05))
compareMF = lapply(entrez, function(x) compareCluster(x, fun="enrichGO",  ont = "MF", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05))
compareCC = lapply(entrez, function(x) compareCluster(x, fun="enrichGO",  ont = "CC", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05))
compareDO = lapply(entrez, function(x) compareCluster(x, fun="enrichDO",  ont = "DO", qvalueCutoff = 0.05, pvalueCutoff = 0.05))

# save object
save(compareKegg, compareBP, compareMF, compareCC, compareDO, file="/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/PMD_KEGG_GO_DO_objects.rda")

# plot compared results
pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/CREs/figures/PMD_KEGG_GO_DO_plots.pdf", height = 40, width = 24)
for (i in 1:length(compareKegg)) {
  print(plot(compareKegg[[i]], colorBy="p.adjust", showCategory = 1000, title= paste0("KEGG Pathway Enrichment: ", names(entrez)[i])))
  print(plot(compareBP[[i]], colorBy="p.adjust", showCategory = 1500, title= paste0("Biological Process GO Enrichment: ", names(entrez)[i])))
  print(plot(compareMF[[i]], colorBy="p.adjust", showCategory = 1000, title= paste0("Molecular Function GO Enrichment: ", names(entrez)[i])))
  print(plot(compareCC[[i]], colorBy="p.adjust", showCategory = 1000, title= paste0("Cellular Compartment GO Enrichment: ", names(entrez)[i])))
  print(plot(compareDO[[i]], colorBy="p.adjust", showCategory = 1000, title= paste0("Disease Ontology Enrichment ", names(entrez)[i])))
}
dev.off()


PMD.CTcomps = lapply(CTcomps, function(x) geneMap[match(x, geneMap$gencodeID),])
PMD.Agecomps = lapply(Agecomps, function(x) geneMap[match(x, geneMap$gencodeID),])

save(PMD.CTcomps, PMD.Agecomps, file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/PMD_gene_comps.rda")



## Plot example PMDs

library("MethylSeekR")
library("BSgenome")
library("BSgenome.Hsapiens.UCSC.hg19")

load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/PMDs_methylSeekR.rda")
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/PMDs_methylSeekR_sortedLister.rda")
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/PMDs_methylSeekR_prenatal.rda")
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/PMDs_methylSeekR_100kbLimit.rda", verbose = T)
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/methylSeekR_objects_prenatal.rda", verbose = T)
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/methylSeekR_objects.rda", verbose = T)

names(total) = gsub("Postnatal-","", names(total))
names(total) = gsub("Prenatal-","", names(total))
pmds = total[-grep("GSM",names(total))]
table(names(c(CGlist,CGPrenlist))==names(pmds))

## Create new objects reflecting Gap overlapping PMD removal

pmds = lapply(total, function(x) x[x$type=="PMD"])
nonpmds = lapply(total, function(x) x[x$type=="notPMD"])
oo = lapply(pmds, function(x) findOverlaps(x, centtelo[which(centtelo$type != "heterochromatin")]))
newnon = mapply(function(p,o) p[unique(queryHits(o))], pmds, oo, SIMPLIFY = F)

merge = mapply(function(gr1, gr2) reduce(c(gr1,gr2)), nonpmds, newnon, SIMPLIFY = F)
oo1 = mapply(function(m, gr1) findOverlaps(m, gr1), merge, nonpmds)
oo2 = mapply(function(m, gr2) findOverlaps(m, gr2), merge, newnon)
grm1 = mapply(function(gr1, o1) split(gr1[subjectHits(o1)], queryHits(o1)), nonpmds, oo1)
grm2 = mapply(function(gr2, o2) split(gr2[subjectHits(o2)], queryHits(o2)), newnon, oo2)

grm1 = lapply(grm1, as.list)
grm2 = lapply(grm2, as.list)
sumCGgrm1 = lapply(grm1, function(x) lapply(x, function(y) sum(y$nCG)))
sumCGgrm2 = lapply(grm2, function(x) lapply(x, function(y) sum(y$nCG)))
grm1 = lapply(grm1, function(g) lapply(g, function(x) x[1]))
for (i in 1:length(grm1)) { for (j in 1:length(grm1[[i]])) { grm1[[i]][[j]]$nCG = sumCGgrm1[[i]][[j]] } }
grm2 = lapply(grm2, function(g) lapply(g, function(x) x[1]))
for (i in 1:length(grm2)) { for (j in 1:length(grm2[[i]])) { grm2[[i]][[j]]$nCG = sumCGgrm2[[i]][[j]] } }
grm1 = lapply(grm1, function(x) do.call(getMethod(c, "GenomicRanges"), GRangesList(x)))
grm2 = lapply(grm2, function(x) do.call(getMethod(c, "GenomicRanges"), GRangesList(x)))

merge = mapply(function(gr1, gr2) reduce(c(gr1,gr2)), nonpmds, newnon, SIMPLIFY = F)
oo1 = mapply(function(m, gr1) findOverlaps(m, gr1), merge, grm1)
for (i in 1:length(merge)) { 
  merge[[i]]$type = "notPMD"
  if (length(subjectHits(oo1[[i]]))==length(merge[[i]])) {
    tmp = merge[[i]][queryHits(oo1[[i]])]
    tmp$nCG1 = grm1[[i]][subjectHits(oo1[[i]])]$nCG
    merge[[i]] = tmp } else {
      tmp = merge[[i]][queryHits(oo1[[i]])]
      tmp1 = merge[[i]][-unique(queryHits(oo1[[i]]))]
      tmp$nCG1 = grm1[[i]][subjectHits(oo1[[i]])]$nCG
      tmp1$nCG1 = 0
      merge[[i]] = c(tmp, tmp1)
    }
}
oo2 = mapply(function(m, gr2) findOverlaps(m, gr2), merge, grm2)
for (i in 1:length(merge)) { 
  tmp = merge[[i]][queryHits(oo2[[i]])]
  tmp2 = merge[[i]][-unique(queryHits(oo2[[i]]))]
  tmp$nCG2 = grm2[[i]][subjectHits(oo2[[i]])]$nCG
  tmp2$nCG2 = 0
  merge[[i]] = c(tmp, tmp2)
}

for (i in 1:length(merge)) { 
  merge[[i]]$nCG = rowSums(as.data.frame(mcols(merge[[i]])[,colnames(mcols(merge[[i]])) %in% c("nCG1","nCG2")]))
  mcols(merge[[i]]) = mcols(merge[[i]])[,colnames(mcols(merge[[i]])) %in% c("type","nCG")]
}

oo = lapply(pmds, function(x) findOverlaps(x, centtelo[which(centtelo$type != "heterochromatin")]))
pmds = mapply(function(p,o) p[-unique(queryHits(o))], pmds, oo, SIMPLIFY = F)
table(elementNROWS(pmds)==elementNROWS(lapply(pmds, function(x) x[which(width(x)>100000)])))

total.nogaps = mapply(function(non, pmd) c(non, pmd), merge, pmds)
total.nogaps = lapply(total.nogaps, function(x) x[order(x)])

table(unlist(lapply(c(PMDsegments.CG.sortedLister, PMDsegments.CG, PMDsegments.CGpren), function(x) sum(x$nCG)))==18664892)
table(unlist(lapply(total.nogaps, function(x) sum(x$nCG)))==18664892)

sLengths = c(lapply(PMDsegments.CG.sortedLister, reduce), lapply(PMDsegments.CG, reduce), 
             lapply(PMDsegments.CGpren, reduce))
red = lapply(total.nogaps, reduce)
table(elementNROWS(red)==25) # all true
table(unlist(mapply(function(x,y) identical(x,y), red, sLengths))=="TRUE") # all true

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/CREs/figures/PMDsegments_visualized_CG_100kb.pdf", width = 12, height = 12)
mapply(function(CG, PMD) plotPMDSegmentation(m = CG, segs = PMD), c(CGlist,CGPrenlist), total.nogaps)
dev.off()

save(total.nogaps, file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/PMDs_100kb_noGaps.rda")