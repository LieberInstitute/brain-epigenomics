library("BSgenome")
library("BSgenome.Hsapiens.UCSC.hg19")
library(jaffelab)
library(GenomicRanges)
library(GenomicFeatures)
library(data.table)
library(ggplot2)
library(clusterProfiler)
require(org.Hs.eg.db)


load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/UMRs_LMRs_methylSeekR.rda")
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/UMRs_LMRs_methylSeekR_prenatal.rda")
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/UMRs_LMRs_DMVs_methylSeekR_sortedLister.rda")
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/DMR/DMR_objects.rda")
load('/dcl01/lieber/ajaffe/lab/brain-epigenomics/bumphunting/BSobj_bsseqSmooth_Neuron_minCov_3.Rdata')


## Isolate UMRs and LMRs

segments = list(postnatal = UMRLMRsegments.CG, postnatal.100kb = UMRLMRsegments.CG.100kb, prenatal = UMRLMRsegments.CGpren,
                prenatal.100kb = UMRLMRsegments.CG.100kb.pren, Lister = UMRLMRsegments.CG.sortedLister, 
                Lister.100kb = UMRLMRsegments.CG.100kb.SL)
elementNROWS(segments)
lapply(segments, function(x) mean(elementNROWS(x))) # limiting PMDs to 100 kb increases the # of UMRs and LMRs found

umrs = lapply(segments, function(x) lapply(x, function(y) y[y$type=="UMR"]))
lapply(umrs, function(x) mean(elementNROWS(x))) # more comparable

lmrs = lapply(segments, function(x) lapply(x, function(y) y[y$type=="LMR"]))
unlist(lapply(lmrs, function(x) mean(elementNROWS(x)))) # much bigger difference in number identified
# postnatal  prenatal    Lister 
#  59580.97  27371.60  54445.25 

# I'm going to stick with the 100kb defined PMD removal UMRs and LMRs

umrs = umrs[grep("100", names(umrs))]
lmrs = lmrs[grep("100", names(lmrs))]
names(umrs) = names(lmrs) = gsub(".100kb", "", names(umrs))
unlist(lapply(lmrs, function(x) mean(elementNROWS(x))))
# postnatal  prenatal    Lister 
#  59580.97  27371.60  54445.25 
unlist(lapply(umrs, function(x) mean(elementNROWS(x))))
# postnatal  prenatal    Lister 
#  20082.50  21344.00  19722.75 


## how many regions are identified?

unum = do.call(rbind, lapply(umrs, function(x) data.frame(num = unlist(elementNROWS(x)))))
unum$id = gsub("^.*\\.","",rownames(unum))
unum[grep("prenatal", rownames(unum)),"celltype"] = "Prenatal"
unum[grep("Lister", rownames(unum)),"celltype"] = c("Neuron","Glia","Neuron","Glia")
unum[which(unum$id %in% pd[pd$Cell.Type=="Neuron","Data.ID"]),"celltype"] = "Neuron"
unum[which(unum$id %in% pd[pd$Cell.Type=="Glia","Data.ID"]),"celltype"] = "Glia"
unum$rep = c(rep.int("Discovery",52), rep.int("Replication",4))

lnum = do.call(rbind, lapply(lmrs, function(x) data.frame(num = unlist(elementNROWS(x)))))
lnum$id = gsub("^.*\\.","",rownames(lnum))
lnum[grep("prenatal", rownames(lnum)),"celltype"] = "Prenatal"
lnum[grep("Lister", rownames(lnum)),"celltype"] = c("Neuron","Glia","Neuron","Glia")
lnum[which(lnum$id %in% pd[pd$Cell.Type=="Neuron","Data.ID"]),"celltype"] = "Neuron"
lnum[which(lnum$id %in% pd[pd$Cell.Type=="Glia","Data.ID"]),"celltype"] = "Glia"
lnum$rep = c(rep.int("Discovery",52), rep.int("Replication",4))
write.csv(rbind(data.frame(unum, category = "UMR"),data.frame(lnum, category = "LMR")), quote=F,
          file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/UMR_LMR_number.csv")


pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/CREs/figures/UMR_LMR_number_byAge_postnatal.pdf")
w = unum[which(unum$celltype %in% c("Neuron", "Glia")),]
w = cbind(w, Age = pd[match(w$id, pd$Data.ID), "Age"])
cor.test(x = w[which(w$celltype=="Neuron"),"Age"], y = w[which(w$celltype=="Neuron"),"num"])
# t = -1.6681, df = 22, p-value = 0.1095
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
#  -0.65054434  0.07899163
#sample estimates:
#  cor 
#-0.3350825 

cor.test(x = w[which(w$celltype=="Glia"),"Age"], y = w[which(w$celltype=="Glia"),"num"])
# t = -0.64527, df = 6, p-value = 0.5426
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
#  -0.8134002  0.5483710
#sample estimates:
#  cor 
#-0.2547386 

ggplot(w, aes(x = Age, y = num, colour = celltype)) + geom_point() +
  geom_smooth(method=lm, se=T, fullrange=TRUE) + scale_colour_brewer(8, palette="Dark2") + 
  labs(fill="") + theme_classic() +
  ylab("Number") + xlab("Age (Years)") +
  ggtitle("Number UMRs Identified by Age") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20), legend.title=element_blank()) + theme(legend.position="bottom")

w = lnum[which(lnum$celltype %in% c("Neuron", "Glia")),]
w = cbind(w, Age = pd[match(w$id, pd$Data.ID), "Age"])
cor.test(x = w[which(w$celltype=="Neuron"),"Age"], y = w[which(w$celltype=="Neuron"),"num"])
# t = -7.6017, df = 22, p-value = 1.364e-07
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
#  -0.9338403 -0.6816540
#sample estimates:
#  cor 
#-0.8510355 

cor.test(x = w[which(w$celltype=="Glia"),"Age"], y = w[which(w$celltype=="Glia"),"num"])
# t = 0.25125, df = 6, p-value = 0.81
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
#  -0.6493243  0.7525961
#sample estimates:
#  cor 
#0.1020367 

ggplot(w, aes(x = Age, y = num, colour = celltype)) + geom_point() +
  geom_smooth(method=lm, se=T, fullrange=TRUE) + scale_colour_brewer(8, palette="Dark2") + 
  labs(fill="") +  theme_classic() +
  ylab("Number") + xlab("Age (Years)") +
  ggtitle("Number of LMRs Identified by Age") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20), legend.title=element_blank()) + theme(legend.position="bottom")

dev.off()


## how long are these regions?

uwidths = do.call(rbind, lapply(umrs, function(x) do.call(rbind, lapply(x, function(y) 
  data.frame(mean = mean(width(y)), median = median(width(y)), sd = sd(width(y)), min = min(width(y)), max = max(width(y)))))))
uwidths$id = gsub("^.*\\.","",rownames(uwidths))
uwidths[grep("prenatal", rownames(uwidths)),"celltype"] = "Prenatal"
uwidths[grep("Lister", rownames(uwidths)),"celltype"] = c("Neuron","Glia","Neuron","Glia")
uwidths[which(uwidths$id %in% pd[pd$Cell.Type=="Neuron","Data.ID"]),"celltype"] = "Neuron"
uwidths[which(uwidths$id %in% pd[pd$Cell.Type=="Glia","Data.ID"]),"celltype"] = "Glia"
uwidths$rep = c(rep.int("Discovery",52), rep.int("Replication",4))

lwidths = do.call(rbind, lapply(lmrs, function(x) do.call(rbind, lapply(x, function(y) 
  data.frame(mean = mean(width(y)), median = median(width(y)), sd = sd(width(y)), min = min(width(y)), max = max(width(y)))))))
lwidths$id = gsub("^.*\\.","",rownames(lwidths))
lwidths[grep("prenatal", rownames(lwidths)),"celltype"] = "Prenatal"
lwidths[grep("Lister", rownames(lwidths)),"celltype"] = c("Neuron","Glia","Neuron","Glia")
lwidths[which(lwidths$id %in% pd[pd$Cell.Type=="Neuron","Data.ID"]),"celltype"] = "Neuron"
lwidths[which(lwidths$id %in% pd[pd$Cell.Type=="Glia","Data.ID"]),"celltype"] = "Glia"
lwidths$rep = c(rep.int("Discovery",52), rep.int("Replication",4))
write.csv(rbind(data.frame(uwidths, category = "UMR"),data.frame(lwidths, category = "LMR")), quote=F,
          file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/UMR_LMR_widths.csv")


pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/CREs/figures/UMR_LMR_widths_byAge_postnatal.pdf")
w = uwidths[which(uwidths$celltype %in% c("Neuron", "Glia")),]
w = cbind(w, Age = pd[match(w$id, pd$Data.ID), "Age"])
cor.test(x = w[which(w$celltype=="Neuron"),"Age"], y = w[which(w$celltype=="Neuron"),"median"])
#t = -3.2393, df = 22, p-value = 0.003766
#95 percent confidence interval:
#  -0.7904663 -0.2139193
#sample estimates:
#  cor 
#-0.568277 

cor.test(x = w[which(w$celltype=="Glia"),"Age"], y = w[which(w$celltype=="Glia"),"median"])
#t = -1.5712, df = 6, p-value = 0.1672
#95 percent confidence interval:
#  -0.901571  0.265945
#sample estimates:
#  cor 
#-0.5399092 

ggplot(w, aes(x = Age, y = median, colour = celltype)) + geom_point() +
  geom_smooth(method=lm, se=T, fullrange=TRUE) + scale_colour_brewer(8, palette="Dark2") + 
  labs(fill="") + theme_classic() +
  ylab("Bases") + xlab("Age (Years)") +
  ggtitle("Median UMR Widths by Age") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20), legend.title=element_blank()) + theme(legend.position="bottom")

w = lwidths[which(lwidths$celltype %in% c("Neuron", "Glia")),]
w = cbind(w, Age = pd[match(w$id, pd$Data.ID), "Age"])
cor.test(x = w[which(w$celltype=="Neuron"),"Age"], y = w[which(w$celltype=="Neuron"),"median"])
#t = -3.2923, df = 22, p-value = 0.003322
#95 percent confidence interval:
#  -0.7939182 -0.2227460
#sample estimates:
#  cor 
#-0.5745193
cor.test(x = w[which(w$celltype=="Glia"),"Age"], y = w[which(w$celltype=="Glia"),"median"])
#t = -0.045923, df = 6, p-value = 0.9649
#95 percent confidence interval:
#  -0.7139867  0.6951099
#sample estimates:
#  cor 
#-0.01874469 

ggplot(w, aes(x = Age, y = median, colour = celltype)) + geom_point() +
  geom_smooth(method=lm, se=T, fullrange=TRUE) + scale_colour_brewer(8, palette="Dark2") + 
  labs(fill="") +  theme_classic() +
  ylab("Bases") + xlab("Age (Years)") +
  ggtitle("Median LMR Widths by Age") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20), legend.title=element_blank()) + theme(legend.position="bottom")

dev.off()


pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/CREs/figures/UMR_LMR_widths.pdf")
ggplot(uwidths, aes(x = celltype, y = median)) + geom_boxplot() +
  facet_grid(. ~ rep, scales = "free_x") +  theme_classic() +
  labs(fill="") +
  ylab("Count") + xlab("") +
  ggtitle("Median UMR Widths") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
ggplot(lwidths, aes(x = celltype, y = median)) + geom_boxplot() +
  facet_grid(. ~ rep, scales = "free_x") +  theme_classic() +
  labs(fill="") + ylab("Count") + xlab("") +
  ggtitle("Median LMR Widths") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
ggplot(uwidths, aes(x = celltype, y = mean)) + geom_boxplot() +
  facet_grid(. ~ rep, scales = "free_x") +  theme_classic() +
  labs(fill="") +
  ylab("Count") + xlab("") +
  ggtitle("Mean UMR Widths") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
ggplot(lwidths, aes(x = celltype, y = mean)) + geom_boxplot() +
  facet_grid(. ~ rep, scales = "free_x") +  theme_classic() +
  labs(fill="") + ylab("Count") + xlab("") +
  ggtitle("Mean LMR Widths") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
dev.off()


## What proportion of the genome do they cover?

sLengths = seqlengths(Hsapiens)
sLengths = sLengths[-c(grep("gl", names(sLengths)), grep("hap", names(sLengths)))]

ugenome = lapply(umrs, function(x) unlist(lapply(x, function(y) sum(width(reduce(y)))/sum(as.numeric(sLengths))*100)))
ugenome = do.call(rbind, Map(cbind, lapply(ugenome, function(x) data.frame(id = names(x), percent = x)), 
              rep = as.list(c(rep.int("Discovery",2),"Replication"))))
ugenome[grep("prenatal", rownames(ugenome)),"celltype"] = "Prenatal"
ugenome[grep("Lister", rownames(ugenome)),"celltype"] = c("Neuron","Glia","Neuron","Glia")
ugenome[which(ugenome$id %in% pd[pd$Cell.Type=="Neuron","Data.ID"]),"celltype"] = "Neuron"
ugenome[which(ugenome$id %in% pd[pd$Cell.Type=="Glia","Data.ID"]),"celltype"] = "Glia"

lgenome = lapply(lmrs, function(x) unlist(lapply(x, function(y) sum(width(reduce(y)))/sum(as.numeric(sLengths))*100)))
lgenome = do.call(rbind, Map(cbind, lapply(lgenome, function(x) data.frame(id = names(x), percent = x)), 
                             rep = as.list(c(rep.int("Discovery",2),"Replication"))))
lgenome[grep("prenatal", rownames(lgenome)),"celltype"] = "Prenatal"
lgenome[grep("Lister", rownames(lgenome)),"celltype"] = c("Neuron","Glia","Neuron","Glia")
lgenome[which(lgenome$id %in% pd[pd$Cell.Type=="Neuron","Data.ID"]),"celltype"] = "Neuron"
lgenome[which(lgenome$id %in% pd[pd$Cell.Type=="Glia","Data.ID"]),"celltype"] = "Glia"
write.csv(rbind(data.frame(ugenome, category = "UMR"),data.frame(lgenome, category = "LMR")), quote=F,
          file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/UMR_LMR_genomeCoverage.csv")


pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/CREs/figures/UMR_LMR_percentGenome_byAge_postnatal.pdf")
w = ugenome[which(ugenome$celltype %in% c("Neuron", "Glia")),]
w = cbind(w, Age = pd[match(w$id, pd$Data.ID), "Age"])
cor.test(x = w[which(w$celltype=="Neuron"),"Age"], y = w[which(w$celltype=="Neuron"),"percent"])
#t = -2.2398, df = 22, p-value = 0.03554
#95 percent confidence interval:
#  -0.71075736 -0.03330204
#sample estimates:
#  cor 
#-0.4309099 

cor.test(x = w[which(w$celltype=="Glia"),"Age"], y = w[which(w$celltype=="Glia"),"percent"])
#t = -1.0744, df = 6, p-value = 0.3239
#95 percent confidence interval:
#  -0.8622795  0.4226219
#sample estimates:
#  cor 
#-0.4016735 

ggplot(w, aes(x = Age, y = percent, colour = celltype)) + geom_point() +
  geom_smooth(method=lm, se=T, fullrange=TRUE) + scale_colour_brewer(8, palette="Dark2") + 
  labs(fill="") +  theme_classic() +
  ylab("Percent") + xlab("Age (Years)") +
  ggtitle("Percent of Genome that is UMR by Age") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20), legend.title=element_blank()) + theme(legend.position="bottom")

w = lgenome[which(lgenome$celltype %in% c("Neuron", "Glia")),]
w = cbind(w, Age = pd[match(w$id, pd$Data.ID), "Age"])
cor.test(x = w[which(w$celltype=="Neuron"),"Age"], y = w[which(w$celltype=="Neuron"),"percent"])
#data:  w[which(w$celltype == "Neuron"), "Age"] and w[which(w$celltype == "Neuron"), "percent"]
#t = -7.1489, df = 22, p-value = 3.624e-07
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
#  -0.9268808 -0.6529278
#sample estimates:
#  cor 
#-0.836104 

cor.test(x = w[which(w$celltype=="Glia"),"Age"], y = w[which(w$celltype=="Glia"),"percent"])
#t = 0.25282, df = 6, p-value = 0.8088
#95 percent confidence interval:
#  -0.6489549  0.7528727
#sample estimates:
#  cor 
#0.1026683 

ggplot(w, aes(x = Age, y = percent, colour = celltype)) + geom_point() +
  geom_smooth(method=lm, se=T, fullrange=TRUE) + scale_colour_brewer(8, palette="Dark2") + 
  labs(fill="") +  theme_classic() +
  ylab("Percent") + xlab("Age (Years)") +
  ggtitle("Percent of Genome that is LMR by Age") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20), legend.title=element_blank()) + theme(legend.position="bottom")

dev.off()


pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/CREs/figures/UMR_LMR_percent_genome.pdf")
ggplot(ugenome, aes(x = celltype, y = percent)) + geom_boxplot() +
  facet_grid(. ~ rep, scales = "free_x") +  theme_classic() +
  labs(fill="") +
  ylab("Count") + xlab("") +
  ggtitle("Percent Genome: UMR") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
ggplot(lgenome, aes(x = celltype, y = percent)) + geom_boxplot() +
  facet_grid(. ~ rep, scales = "free_x") +  theme_classic() +
  labs(fill="") + ylab("Count") + xlab("") +
  ggtitle("Percent Genome: LMR") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
dev.off()


### How much of the UMRs and LMRs overlap?

# Master list of all bases in each setting
umaster = reduce(do.call(getMethod(c, "GenomicRanges"), c(GRangesList(umrs$postnatal), GRangesList(umrs$prenatal)))) # 43995 total
lmaster = reduce(do.call(getMethod(c, "GenomicRanges"), c(GRangesList(lmrs$postnatal), GRangesList(lmrs$prenatal)))) # 238576 total
do.call(rbind, lapply(list(umaster,lmaster), function(x) { c(mean=mean(width(x)), median=median(width(x)), sd=sd(width(x)), 
                                                             min=min(width(x)), max=max(width(x))) }))
#     mean median        sd min    max
# 7584.632   4822 12218.260 221 461743
# 1073.481    762  1064.274   8  52517


## regions that are shared by samples

uall = list("All" = Reduce(intersect, c(umrs$postnatal, umrs$prenatal)),
            "Prenatal" = Reduce(intersect, umrs$prenatal),
            "Postnatal" = Reduce(intersect, umrs$postnatal),
            "Neurons" = Reduce(intersect, umrs$postnatal[which(names(umrs$postnatal) %in% pd[pd$Cell.Type=="Neuron","Data.ID"])]),
            "Glia" = Reduce(intersect, umrs$postnatal[which(names(umrs$postnatal) %in% pd[pd$Cell.Type=="Glia","Data.ID"])]))
elementNROWS(uall)
#         all  allPrenatal allPostnatal   allNeurons      allGlia 
#       12900        13185        17098        17525        18427

round((unlist(lapply(uall, function(x) (sum(width(x))/sum(as.numeric(sLengths)))*100)))/
        (sum(width(umaster))/sum(as.numeric(sLengths))*100),2)
# all  allPrenatal allPostnatal   allNeurons      allGlia 
# 0.07         0.08         0.09         0.10         0.11 
# 7-11% of total bases in UMR state shared by all or most

ushared = mapply(function(all,ind) lapply( ind, function(x) round(all / x *100,2)), unlist(lapply(uall, function(x) (sum(width(reduce(x)))/sum(as.numeric(sLengths))*100))), 
                  lapply(list(c(umrs$postnatal, umrs$prenatal), umrs$prenatal, umrs$postnatal, umrs$postnatal[which(names(umrs$postnatal) %in% pd[pd$Cell.Type=="Neuron","Data.ID"])],
                   umrs$postnatal[which(names(umrs$postnatal) %in% pd[pd$Cell.Type=="Glia","Data.ID"])]), function(y) lapply(y, function(x) sum(width(reduce(x)))/sum(as.numeric(sLengths))*100)))
ushared = data.frame(perc = unlist(ushared))
ushared$Group = gsub("\\..*","", rownames(ushared))

lall = list("All" = Reduce(intersect, c(lmrs$postnatal, lmrs$prenatal)),
            "Prenatal" = Reduce(intersect, lmrs$prenatal),
            "Postnatal" = Reduce(intersect, lmrs$postnatal),
            "Neurons" = Reduce(intersect, lmrs$postnatal[which(names(lmrs$postnatal) %in% pd[pd$Cell.Type=="Neuron","Data.ID"])]),
            "Glia" = Reduce(intersect, lmrs$postnatal[which(names(lmrs$postnatal) %in% pd[pd$Cell.Type=="Glia","Data.ID"])]))
elementNROWS(lall)
#      All  Prenatal Postnatal   Neurons      Glia 
#       14        41      6162     11499     26041 

round((unlist(lapply(lall, function(x) (sum(width(x))/sum(as.numeric(sLengths)))*100)))/
        (sum(width(lmaster))/sum(as.numeric(sLengths))*100),2)
#  all  allPrenatal allPostnatal   allNeurons      allGlia 
# 0.00         0.00         0.01         0.02         0.05  
# Very little of total bases in LMR state shared by all or most

lshared = mapply(function(all,ind) lapply( ind, function(x) round(all / x *100,2)), unlist(lapply(lall, function(x) (sum(width(reduce(x)))/sum(as.numeric(sLengths))*100))), 
                 lapply(list(c(lmrs$postnatal, lmrs$prenatal), lmrs$prenatal, lmrs$postnatal, lmrs$postnatal[which(names(lmrs$postnatal) %in% pd[pd$Cell.Type=="Neuron","Data.ID"])],
                             lmrs$postnatal[which(names(lmrs$postnatal) %in% pd[pd$Cell.Type=="Glia","Data.ID"])]), function(y) lapply(y, function(x) sum(width(reduce(x)))/sum(as.numeric(sLengths))*100)))
lshared = data.frame(perc = unlist(lshared))
lshared$Group = gsub("\\..*","", rownames(lshared))

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/CREs/figures/percent_shared_UMRs_LMRs_perSample_byBaseCoverage.pdf")
ggplot(ushared, aes(x = Group, y = perc)) + geom_boxplot() +
  labs(fill="") + theme_classic() +
  ylim(0,100) +
  ylab("Percent") + xlab("") +
  ggtitle("Percent Bases Shared Per Sample: UMR") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
ggplot(lshared, aes(x = Group, y = perc)) + geom_boxplot() +
  labs(fill="") + ylab("Percent") + xlab("") +
  ylim(0,100) + theme_classic() +
  ggtitle("Percent Bases Shared Per Sample: LMR") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
dev.off()


## Find overlaps with DMRs, features, shared groups

txdb = loadDb("/dcl01/lieber/ajaffe/Amanda/annotation_objects/gencode.v25lift37.annotation.sqlite")
islands = read.table("/dcl01/lieber/ajaffe/Amanda/annotation_objects/cpgIslandExt.hg19.txt", sep="\t", header = T)
rpmsk = read.table("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/RepeatMasker_genomewide_HG19.txt", header=T)
features = list(CDS = cdsBy(txdb, by="tx", use.names=T), Introns = intronsByTranscript(txdb, use.names=T), 
                UTR5 = fiveUTRsByTranscript(txdb, use.names=T), UTR3 = threeUTRsByTranscript(txdb, use.names=T))
features = lapply(features, function(x) unlist(x, recursive = TRUE, use.names = TRUE))
for (i in 1:length(features)){
  tmp = features[[i]]
  tmp$TxID = names(tmp)
  features[[i]] = tmp
}
features = c(features, 
             rpmskgr = makeGRangesFromDataFrame(rpmsk, seqnames.field="genoName",start.field="genoStart",end.field="genoEnd",keep.extra.columns=TRUE),
             islands = makeGRangesFromDataFrame(islands, keep.extra.columns = T, start.field = "chromStart", end.field = "chromEnd"),
             promoters = promoters(txdb, upstream=2000, downstream=200))
lapply(features, head)

DMR = lapply(DMR, function(x) x[which(x$sig=="FWER < 0.05"),])
DMRgr = lapply(DMR, function(x) makeGRangesFromDataFrame(x, keep.extra.columns = T))

loo = lapply(c(lmrs$postnatal, lmrs$prenatal), function(x) lapply(c(DMRgr, lall, features), function(d) findOverlaps(x,d)))
uoo = lapply(c(umrs$postnatal, umrs$prenatal), function(x) lapply(c(DMRgr, uall, features), function(d) findOverlaps(x,d)))

uDMR = lapply(c(umrs$postnatal, umrs$prenatal), as.data.frame)
uDMR = lapply(uDMR, function(x) data.frame(x, rnum = 1:nrow(x)))
uDMR = mapply(function(d,oo) data.frame(d, CT = ifelse(d$rnum %in% queryHits(oo$CellType), "CT", "no"),
                                        Age = ifelse(d$rnum %in% queryHits(oo$Age), "Age", "no"),
                                        Interaction = ifelse(d$rnum %in% queryHits(oo$Interaction), "Int", "no"),
                                        All = ifelse(d$rnum %in% queryHits(oo$All), "All", "no"),
                                        Prenatal = ifelse(d$rnum %in% queryHits(oo$Prenatal), "Prenatal", "no"),
                                        Postnatal = ifelse(d$rnum %in% queryHits(oo$Postnatal), "Postnatal", "no"),
                                        Neurons = ifelse(d$rnum %in% queryHits(oo$Neurons), "Neurons", "no"),
                                        Glia = ifelse(d$rnum %in% queryHits(oo$Glia), "Glia", "no"),
                                        cds = ifelse(d$rnum %in% queryHits(oo$CDS), "CDS", NA),
                                        intron = ifelse(d$rnum %in% queryHits(oo$Introns), "Intron", NA),
                                        UTR5 = ifelse(d$rnum %in% queryHits(oo$UTR5), "5'UTR", NA),
                                        UTR3 = ifelse(d$rnum %in% queryHits(oo$UTR3), "3'UTR", NA),
                                        islands = ifelse(d$rnum %in% queryHits(oo$islands), "CpG-Island", "non-Island"),
                                        promoter = ifelse(d$rnum %in% queryHits(oo$promoters), "Promoter", NA)), uDMR, uoo, SIMPLIFY = F) 
uDMR = lapply(uDMR, function(d) data.frame(d, tog = paste(d$CT, d$Age, d$Interaction, d$All, d$Prenatal, 
                                                          d$Postnatal, d$Neurons, d$Glia, sep = ":"),
                                           dmr = paste(d$CT, d$Age, d$Interaction, sep = ":"),
                                           regionID = paste0(d$seqnames,":",d$start,"-", d$end),
                                           anno = paste0(d$cds,":",d$intron, ":", d$UTR5, ":", d$UTR3, ":", d$promoter)))
for (i in 1:length(uDMR)) {
  uDMR[[i]] = rbind(data.frame(uDMR[[i]][queryHits(uoo[[i]][["rpmskgr"]]),], 
                               repeats = as.character(features$rpmskgr$repClass)[subjectHits(uoo[[i]][["rpmskgr"]])]), 
                    data.frame(uDMR[[i]][-unique(queryHits(uoo[[i]][["rpmskgr"]])),], repeats = "No repeats"))
  uDMR[[i]][which(uDMR[[i]][,"anno"] == "NA:NA:NA:NA:NA"),"annotation"] = "Intergenic" 
  uDMR[[i]][grep("CDS", uDMR[[i]][,"cds"]),"annotation"] = "CDS"
  uDMR[[i]][which(is.na(uDMR[[i]][,"annotation"]) & uDMR[[i]][,"UTR5"] == "5'UTR"),"annotation"] = "5'UTR"
  uDMR[[i]][which(is.na(uDMR[[i]][,"annotation"]) & uDMR[[i]][,"UTR3"] == "3'UTR"),"annotation"] = "3'UTR"
  uDMR[[i]][which(is.na(uDMR[[i]][,"annotation"]) & uDMR[[i]][,"intron"] == "Intron"),"annotation"] = "Intron"
  uDMR[[i]][which(is.na(uDMR[[i]][,"annotation"]) & uDMR[[i]][,"promoter"] == "Promoter"),"annotation"] = "Promoter"
}

uDMRgr = lapply(uDMR, function(x) makeGRangesFromDataFrame(x, keep.extra.columns = T))
geneMapGR = makeGRangesFromDataFrame(geneMap, start.field="Start",end.field="End",strand.field="Strand",keep=TRUE)
uoo = lapply(uDMRgr, function(x) findOverlaps(x, geneMapGR))

dA = mapply(function(d,oo) distanceToNearest(d[-unique(queryHits(oo)),], geneMapGR), uDMRgr, uoo, SIMPLIFY = F)

for (i in 1:length(uDMR)) {
  uDMR[[i]] = rbind(data.frame(uDMR[[i]][queryHits(uoo[[i]]),], nearestSymbol = geneMapGR$Symbol[subjectHits(uoo[[i]])], 
                               nearestID = names(geneMapGR)[subjectHits(uoo[[i]])],
                               EntrezID = geneMapGR$EntrezID[subjectHits(uoo[[i]])], distToGene = 0), 
                    data.frame(uDMR[[i]][-unique(queryHits(uoo[[i]])),], nearestSymbol = geneMapGR$Symbol[subjectHits(dA[[i]])], 
                               nearestID = names(geneMapGR)[subjectHits(dA[[i]])],
                               EntrezID = geneMapGR$EntrezID[subjectHits(dA[[i]])], distToGene = mcols(dA[[i]])$distance))
}
  

lDMR = lapply(c(lmrs$postnatal, lmrs$prenatal), as.data.frame)
lDMR = lapply(lDMR, function(x) data.frame(x, rnum = 1:nrow(x)))
lDMR = mapply(function(d,oo) data.frame(d, CT = ifelse(d$rnum %in% queryHits(oo$CellType), "CT", "no"),
                                        Age = ifelse(d$rnum %in% queryHits(oo$Age), "Age", "no"),
                                        Interaction = ifelse(d$rnum %in% queryHits(oo$Interaction), "Int", "no"),
                                        All = ifelse(d$rnum %in% queryHits(oo$All), "All", "no"),
                                        Prenatal = ifelse(d$rnum %in% queryHits(oo$Prenatal), "Prenatal", "no"),
                                        Postnatal = ifelse(d$rnum %in% queryHits(oo$Postnatal), "Postnatal", "no"),
                                        Neurons = ifelse(d$rnum %in% queryHits(oo$Neurons), "Neurons", "no"),
                                        Glia = ifelse(d$rnum %in% queryHits(oo$Glia), "Glia", "no"),
                                        cds = ifelse(d$rnum %in% queryHits(oo$CDS), "CDS", NA),
                                        intron = ifelse(d$rnum %in% queryHits(oo$Introns), "Intron", NA),
                                        UTR5 = ifelse(d$rnum %in% queryHits(oo$UTR5), "5'UTR", NA),
                                        UTR3 = ifelse(d$rnum %in% queryHits(oo$UTR3), "3'UTR", NA),
                                        islands = ifelse(d$rnum %in% queryHits(oo$islands), "CpG-Island", "non-Island"),
                                        promoter = ifelse(d$rnum %in% queryHits(oo$promoters), "Promoter", NA)), lDMR, loo, SIMPLIFY = F) 
lDMR = lapply(lDMR, function(d) data.frame(d, tog = paste(d$CT, d$Age, d$Interaction, d$All, d$Prenatal, 
                                                          d$Postnatal, d$Neurons, d$Glia, sep = ":"),
                                           dmr = paste(d$CT, d$Age, d$Interaction, sep = ":"),
                                           regionID = paste0(d$seqnames,":",d$start,"-", d$end),
                                           anno = paste0(d$cds,":",d$intron, ":", d$UTR5, ":", d$UTR3, ":", d$promoter)))
for (i in 1:length(lDMR)) {
  lDMR[[i]] = rbind(data.frame(lDMR[[i]][queryHits(loo[[i]][["rpmskgr"]]),], 
                               repeats = as.character(features$rpmskgr$repClass)[subjectHits(loo[[i]][["rpmskgr"]])]), 
                    data.frame(lDMR[[i]][-unique(queryHits(loo[[i]][["rpmskgr"]])),], repeats = "No repeats"))
  lDMR[[i]][which(lDMR[[i]][,"anno"] == "NA:NA:NA:NA:NA"),"annotation"] = "Intergenic" 
  lDMR[[i]][grep("CDS", lDMR[[i]][,"cds"]),"annotation"] = "CDS"
  lDMR[[i]][which(is.na(lDMR[[i]][,"annotation"]) & lDMR[[i]][,"UTR5"] == "5'UTR"),"annotation"] = "5'UTR"
  lDMR[[i]][which(is.na(lDMR[[i]][,"annotation"]) & lDMR[[i]][,"UTR3"] == "3'UTR"),"annotation"] = "3'UTR"
  lDMR[[i]][which(is.na(lDMR[[i]][,"annotation"]) & lDMR[[i]][,"intron"] == "Intron"),"annotation"] = "Intron"
  lDMR[[i]][which(is.na(lDMR[[i]][,"annotation"]) & lDMR[[i]][,"promoter"] == "Promoter"),"annotation"] = "Promoter"
}

lDMRgr = lapply(lDMR, function(x) makeGRangesFromDataFrame(x, keep.extra.columns = T))
geneMapGR = makeGRangesFromDataFrame(geneMap, start.field="Start",end.field="End",strand.field="Strand",keep=TRUE)
loo = lapply(lDMRgr, function(x) findOverlaps(x, geneMapGR))

dA = mapply(function(d,oo) distanceToNearest(d[-unique(queryHits(oo)),], geneMapGR), lDMRgr, loo, SIMPLIFY = F)

for (i in 1:length(lDMR)) {
  lDMR[[i]] = rbind(data.frame(lDMR[[i]][queryHits(loo[[i]]),], nearestSymbol = geneMapGR$Symbol[subjectHits(loo[[i]])], 
                               nearestID = names(geneMapGR)[subjectHits(loo[[i]])],
                               EntrezID = geneMapGR$EntrezID[subjectHits(loo[[i]])], distToGene = 0), 
                    data.frame(lDMR[[i]][-unique(queryHits(loo[[i]])),], nearestSymbol = geneMapGR$Symbol[subjectHits(dA[[i]])], 
                               nearestID = names(geneMapGR)[subjectHits(dA[[i]])],
                               EntrezID = geneMapGR$EntrezID[subjectHits(dA[[i]])], distToGene = mcols(dA[[i]])$distance))
}

postnatalpd = pd
save(uDMR, lDMR, postnatalpd, geneMap, file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/UMRs_LMRs_annotated.rda")


## Annotate to genomic features

load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/UMRs_LMRs_annotated.rda")

# Identify all CpG clusters in the genome
gr = granges(BSobj)
cl = clusterMaker( chr = as.character(seqnames(gr)), 
                   pos = start(gr),   maxGap = 1000)
gr.clusters = split(gr, cl)
gr.clusters = unlist(range(gr.clusters))
df.clusters = as.data.frame(gr.clusters)
df.clusters$regionID = paste0(df.clusters$seqnames,":",df.clusters$start,"-",df.clusters$end)
df.clusters$rnum = 1:length(gr.clusters)
df.clusters$Islands = ifelse(df.clusters$rnum %in% subjectHits(findOverlaps(features$islands, gr.clusters)), "island","no")
df.clusters$repeats = ifelse(df.clusters$rnum %in% subjectHits(findOverlaps(features$rpmskgr, gr.clusters)), "repeat","no")
df.clusters$promoters = ifelse(df.clusters$rnum %in% subjectHits(findOverlaps(features$promoters, gr.clusters)), "promoter","no")
df.clusters$genes = ifelse(df.clusters$rnum %in% subjectHits(findOverlaps(makeGRangesFromDataFrame(geneMap), gr.clusters)), "gene","no")

dt = list(UMR = lapply(uDMR, data.table), LMR = lapply(lDMR, data.table))

oo = lapply(dt, function(x) lapply(x, function(y) findOverlaps(gr.clusters, makeGRangesFromDataFrame(y))))
df.clusters = lapply(oo, function(x) lapply(x, function(y) 
  data.frame(df.clusters, overlap = ifelse(df.clusters$rnum %in% queryHits(y), "hit","no"))))


## how many fall within CpG islands?

CpGIslands = lapply(dt, function(x) lapply(x, function(y) round(length(unique(y[which(y$islands=="CpG-Island"),,]$regionID))/length(unique(y$regionID))*100,2)))
CpGIslands = lapply(CpGIslands, function(x) data.frame(islands = unlist(x), id = names(x)))
CpGIslands = do.call(rbind, Map(cbind, CpGIslands, mr = as.list(names(CpGIslands))))
CpGIslands[!(CpGIslands$id %in% postnatalpd$Data.ID),"celltype"] = "Prenatal"
CpGIslands[which(CpGIslands$id %in% postnatalpd[which(postnatalpd$Cell.Type=="Neuron"),"Data.ID"]),"celltype"] = "Neuron"
CpGIslands[which(CpGIslands$id %in% postnatalpd[which(postnatalpd$Cell.Type=="Glia"),"Data.ID"]),"celltype"] = "Glia"
write.csv(CpGIslands, quote=F,file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/UMR_LMR_CpG_Island_Overlap.csv")

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/CREs/figures/UMR_LMR_overlap_with_CpG-Islands.pdf")
ggplot(CpGIslands, aes(x = celltype, y = islands)) + geom_boxplot() +
  theme_classic() + facet_grid(. ~ mr) +
  labs(fill="") +
  ylab("Percent") + ylim(0,100) + 
  xlab("") +
  ggtitle("Percent Overlapping CpG Islands") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20), legend.title=element_blank()) + theme(legend.position="bottom")
dev.off()


# How about Repetitive Elements?

repeats = lapply(dt, function(x) lapply(x, function(y) y[,length(unique(regionID)), by = "repeats"]))
repeats = lapply(repeats, function(x) lapply(x, function(y) data.frame(y, perc = round(y$V1/sum(y$V1)*100,2))))
repeats = lapply(repeats, function(x) do.call(rbind, Map(cbind, x, id = as.list(names(x)))))
repeats = do.call(rbind, Map(cbind, repeats, mr = as.list(names(repeats))))
repeats$repeats = factor(repeats$repeats, levels = c("SINE","LINE","Simple_repeat","DNA","LTR","Low_complexity","No repeats","Unknown","snRNA",
                                         "Other","DNA?","Satellite","RC","scRNA","tRNA","SINE?","srpRNA","RNA","rRNA","LINE?","LTR?","Unknown?"))
repeats[!(repeats$id %in% postnatalpd$Data.ID),"celltype"] = "Prenatal"
repeats[which(repeats$id %in% postnatalpd[which(postnatalpd$Cell.Type=="Neuron"),"Data.ID"]),"celltype"] = "Neuron"
repeats[which(repeats$id %in% postnatalpd[which(postnatalpd$Cell.Type=="Glia"),"Data.ID"]),"celltype"] = "Glia"
write.csv(repeats, quote=F,file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/UMR_LMR_repeats_Overlap.csv")

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/CREs/figures/UMR_LMR_overlap_with_repeats.pdf")
x = repeats[which(repeats$repeats %in% c("SINE","LINE","Simple_repeat","DNA","LTR","Low_complexity","No repeats")),]
x$repeats = droplevels(x$repeats)
ggplot(x[which(x$mr=="UMR"),], aes(x = celltype, y = perc, fill = repeats)) + 
  theme_classic() + geom_bar(position = "fill",stat = "identity") +
  labs(fill="") +
  ylab("Proportion") + 
  xlab("") +
  ggtitle("Proportion Overlapping Repeats: UMR") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20), legend.title=element_blank()) + theme(legend.position="bottom")
ggplot(x[which(x$mr=="LMR"),], aes(x = celltype, y = perc, fill = repeats)) + 
  theme_classic() + geom_bar(position = "fill",stat = "identity") +
  labs(fill="") +
  ylab("Proportion") + 
  xlab("") +
  ggtitle("Proportion Overlapping Repeats: LMR") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20), legend.title=element_blank()) + theme(legend.position="bottom")
dev.off()


## assign genomic features

annotation = lapply(dt, function(x) lapply(x, function(y) y[,length(unique(regionID)), by = "annotation"]))
annotation = lapply(annotation, function(x) lapply(x, function(y) data.frame(y, perc = round(y$V1/sum(y$V1)*100,2))))
annotation = lapply(annotation, function(x) do.call(rbind, Map(cbind, x, id = as.list(names(x)))))
annotation = do.call(rbind, Map(cbind, annotation, mr = as.list(names(annotation))))
annotation$annotation = factor(annotation$annotation, levels = c("Promoter", "5'UTR", "CDS", "3'UTR", "Intron", "Intergenic"))
annotation[!(annotation$id %in% postnatalpd$Data.ID),"celltype"] = "Prenatal"
annotation[which(annotation$id %in% postnatalpd[which(postnatalpd$Cell.Type=="Neuron"),"Data.ID"]),"celltype"] = "Neuron"
annotation[which(annotation$id %in% postnatalpd[which(postnatalpd$Cell.Type=="Glia"),"Data.ID"]),"celltype"] = "Glia"
write.csv(annotation, quote=F,file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/UMR_LMR_annotation_Overlap.csv")

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/CREs/figures/UMR_LMR_overlap_with_annotation.pdf", width = 14)
ggplot(annotation, aes(x = celltype, y = perc, fill = annotation)) + 
  theme_classic() + geom_bar(position = "fill",stat = "identity") +
  labs(fill="") + facet_grid(. ~ mr) +
  ylab("Proportion") + 
  xlab("") +
  ggtitle("Proportion Overlapping Features") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20), legend.title=element_blank()) + theme(legend.position="bottom")
dev.off()


## Overlap with DMRs

dmr = lapply(dt, function(x) lapply(x, function(y) y[,length(unique(regionID)), by = "dmr"]))
dmr = lapply(dmr, function(x) lapply(x, function(y) data.frame(y, perc = round(y$V1/sum(y$V1)*100,2))))
models = lapply(dmr, function(x) lapply(x, function(y) data.frame(CT = sum(data.frame(y)[grep("CT", y$dmr),"perc"]), 
                                                                  Age = sum(data.frame(y)[grep("Age", y$dmr),"perc"]),
                                                                  Int = sum(data.frame(y)[grep("Int", y$dmr),"perc"]))))
dmr = lapply(dmr, function(x) do.call(rbind, Map(cbind, x, id = as.list(names(x)))))
models = lapply(models, function(x) do.call(rbind, Map(cbind, x, id = as.list(names(x)))))

dmr = do.call(rbind, Map(cbind, dmr, mr = as.list(names(dmr))))
models = do.call(rbind, Map(cbind, models, mr = as.list(names(models))))
dmr$dmr = factor(dmr$dmr, levels = c("no:no:no","CT:no:no","no:no:Int","CT:no:Int","no:Age:no","CT:Age:Int","no:Age:Int","CT:Age:no"))
dmr[!(dmr$id %in% postnatalpd$Data.ID),"celltype"] = "Prenatal"
dmr[which(dmr$id %in% postnatalpd[which(postnatalpd$Cell.Type=="Neuron"),"Data.ID"]),"celltype"] = "Neuron"
dmr[which(dmr$id %in% postnatalpd[which(postnatalpd$Cell.Type=="Glia"),"Data.ID"]),"celltype"] = "Glia"
models[!(models$id %in% postnatalpd$Data.ID),"celltype"] = "Prenatal"
models[which(models$id %in% postnatalpd[which(postnatalpd$Cell.Type=="Neuron"),"Data.ID"]),"celltype"] = "Neuron"
models[which(models$id %in% postnatalpd[which(postnatalpd$Cell.Type=="Glia"),"Data.ID"]),"celltype"] = "Glia"

write.csv(dmr, quote=F,file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/UMR_LMR_dmr_Overlap.csv")
write.csv(models, quote=F,file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/UMR_LMR_dmr_Overlap_byModel.csv")


pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/CREs/figures/UMR_LMR_overlap_with_dmr.pdf", width = 14)
ggplot(dmr, aes(x = celltype, y = perc, fill = dmr)) + 
  theme_classic() + geom_bar(position = "fill",stat = "identity") +
  labs(fill="") + facet_grid(. ~ mr) +
  ylab("Proportion") + 
  xlab("") +
  ggtitle("Proportion Overlapping DMRs") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20), legend.title=element_blank()) + theme(legend.position="bottom")
ggplot(dmr[which(dmr$dmr=="no:no:no"),], aes(x = celltype, y = perc, fill = dmr)) + 
  theme_classic() + geom_boxplot() +
  labs(fill="") + facet_grid(. ~ mr) +
  ylab("Percent") + 
  xlab("") +
  ggtitle("Percent Not Overlapping DMRs") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20), legend.title=element_blank()) + theme(legend.position="bottom")
ggplot(models, aes(x = celltype, y = CT)) + 
  theme_classic() + geom_boxplot() +
  labs(fill="") + facet_grid(. ~ mr) +
  ylab("Percent") + 
  xlab("") +
  ggtitle("Percent Overlapping Cell Type DMRs") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20), legend.title=element_blank()) + theme(legend.position="bottom")
ggplot(models, aes(x = celltype, y = Age)) + 
  theme_classic() + geom_boxplot() +
  labs(fill="") + facet_grid(. ~ mr) +
  ylab("Percent") + 
  xlab("") +
  ggtitle("Percent Overlapping Age DMRs") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20), legend.title=element_blank()) + theme(legend.position="bottom")
ggplot(models, aes(x = celltype, y = Int)) + 
  theme_classic() + geom_boxplot() +
  labs(fill="") + facet_grid(. ~ mr) +
  ylab("Percent") + 
  xlab("") +
  ggtitle("Percent Overlapping Interaction DMRs") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20), legend.title=element_blank()) + theme(legend.position="bottom")
dev.off()


## identify sequence present in all of a group

genes = list("All" = unique(unlist(lapply(dt, function(x) lapply(x, function(y) data.frame(y)[which(y$All=="All"),"EntrezID"])))),
             "Prenatal" = unique(unlist(lapply(dt, function(x) lapply(x, function(y) data.frame(y)[which(y$Prenatal=="Prenatal"),"EntrezID"])))),
             "Postnatal" = unique(unlist(lapply(dt, function(x) lapply(x, function(y) data.frame(y)[which(y$Postnatal=="Postnatal"),"EntrezID"])))),
             "Neurons" = unique(unlist(lapply(dt, function(x) lapply(x, function(y) data.frame(y)[which(y$Neurons=="Neurons"),"EntrezID"])))),
             "Glia" = unique(unlist(lapply(dt, function(x) lapply(x, function(y) data.frame(y)[which(y$Glia=="Glia"),"EntrezID"])))),
             "AllUMR" = unique(unlist(lapply(dt$UMR, function(x) data.frame(x)[which(x$All=="All"),"EntrezID"]))),
             "PrenatalUMR" = unique(unlist(lapply(dt$UMR, function(x) data.frame(x)[which(x$Prenatal=="Prenatal"),"EntrezID"]))),
             "PostnatalUMR" = unique(unlist(lapply(dt$UMR, function(x) data.frame(x)[which(x$Postnatal=="Postnatal"),"EntrezID"]))),
             "NeuronsUMR" = unique(unlist(lapply(dt$UMR, function(x) data.frame(x)[which(x$Neurons=="Neurons"),"EntrezID"]))),
             "GliaUMR" = unique(unlist(lapply(dt$UMR, function(x) data.frame(x)[which(x$Glia=="Glia"),"EntrezID"]))),
             "AllLMR" = unique(unlist(lapply(dt$LMR, function(x) data.frame(x)[which(x$All=="All"),"EntrezID"]))),
             "PrenatalLMR" = unique(unlist(lapply(dt$LMR, function(x) data.frame(x)[which(x$Prenatal=="Prenatal"),"EntrezID"]))),
             "PostnatalLMR" = unique(unlist(lapply(dt$LMR, function(x) data.frame(x)[which(x$Postnatal=="Postnatal"),"EntrezID"]))),
             "NeuronsLMR" = unique(unlist(lapply(dt$LMR, function(x) data.frame(x)[which(x$Neurons=="Neurons"),"EntrezID"]))),
             "GliaLMR" = unique(unlist(lapply(dt$LMR, function(x) data.frame(x)[which(x$Glia=="Glia"),"EntrezID"]))))
genes = lapply(genes, function(x) na.omit(as.character(x)))

# Compare the enriched terms between 7 groups
# KEGG
compareKegg = compareCluster(genes, fun="enrichKEGG", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
# Biological Process
compareBP = compareCluster(genes, fun="enrichGO", ont = "BP", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
# Molecular Function
compareMF = compareCluster(genes, fun="enrichGO",  ont = "MF", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
# Cellular Component
compareCC = compareCluster(genes, fun="enrichGO",  ont = "CC", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
# Disease Ontology
compareDO = compareCluster(genes, fun="enrichDO",  ont = "DO", qvalueCutoff = 0.05, pvalueCutoff = 0.05)

# save object
save(compareKegg, compareBP, compareMF, compareCC, compareDO,
     file="/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/UMR_LMR_KEGG_GO_DO_objects.rda")

# plot compared results
pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/CREs/figures/UMR_LMR_KEGG_GO_DO_plots.pdf", height = 160, width = 26)
plot(compareKegg, colorBy="p.adjust", showCategory = 1000, title= "KEGG Pathway Enrichment")
plot(compareMF, colorBy="p.adjust", showCategory = 1000, title= "Molecular Function GO Enrichment")
plot(compareCC, colorBy="p.adjust", showCategory = 1000, title= "Cellular Compartment GO Enrichment")
plot(compareDO, colorBy="p.adjust", showCategory = 1000, title= "Disease Ontology Enrichment")
dev.off()
pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/CREs/figures/UMR_LMR_BP_plot.pdf", height = 275, width = 26)
plot(compareBP, colorBy="p.adjust", showCategory = 1500, title= "Biological Process GO Enrichment")
dev.off()















