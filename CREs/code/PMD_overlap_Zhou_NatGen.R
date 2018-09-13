library(ggplot2)
library(GenomicRanges)


## Get PMDs

load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/PMDs_methylSeekR.rda")
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/PMDs_methylSeekR_prenatal.rda")
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/DMR/DMR_objects.rda")

pmds = c(lapply(PMDsegments.CG, function(x) x[x$type=="PMD"]), lapply(PMDsegments.CGpren, function(x) x[x$type=="PMD"]))
pmds = lapply(pmds, function(x) x[which(width(x)>=100000)])
centtelo = read.table("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/centromeres_telomeres_hg19.txt", header=T)
centtelo = makeGRangesFromDataFrame(centtelo, keep.extra.columns = T)
oo = lapply(pmds, function(x) findOverlaps(x, centtelo[which(centtelo$type != "heterochromatin")]))
pmds = mapply(function(p,o) p[-queryHits(o)], pmds, oo, SIMPLIFY = F)

## Identify overlap with common HMDs and PMDs

solo = read.table("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/PMD_coordinates_hg19.bed.gz",
                  col.names = c("chromosome","start","end","fraction","status","common"))
solo = unlist(lapply(split(solo, solo$common), function(x) split(x, x$status)[elementNROWS(split(x, x$status))>0]), recursive = F)
solo = lapply(solo, function(x) makeGRangesFromDataFrame(x, keep.extra.columns = T))
names(solo) = c("commonHMD", "commonPMD", "neither.HMD", "neither.PMD")

# Total bases
celltype = list(LIBDneuron = pmds[which(names(pmds) %in% pd[pd$Cell.Type=="Neuron","Data.ID"])],
                LIBDglia = pmds[which(names(pmds) %in% pd[pd$Cell.Type=="Glia","Data.ID"])],
                LIBDprenatal = pmds[-which(names(pmds) %in% pd[pd$Cell.Type %in% c("Neuron","Glia"),"Data.ID"])])
celltype = lapply(celltype, function(x) reduce(makeGRangesFromDataFrame(do.call(rbind, lapply(x,as.data.frame)), keep.extra.columns = T)))

overlap = lapply(celltype, function(cell) lapply(solo, function(com) c( OurOverlap = round(sum(width(Reduce(intersect, list(cell,com))))/sum(width(cell))*100,1),
                                                                        commonOverlap = round(sum(width(Reduce(intersect, list(cell,com))))/sum(width(com))*100,1))))
overlap = lapply(overlap, function(x) do.call(rbind, Map(cbind, lapply(x, function(y) data.frame(OurOverlap = y["OurOverlap"], commonOverlap = y["commonOverlap"])),
                                                         common = as.list(names(x)))))
overlap = do.call(rbind, Map(cbind, overlap, LIBD = as.list(names(celltype))))
overlap = reshape2::melt(overlap)
overlap$LIBD = gsub("LIBD", "", overlap$LIBD)
overlap$common = gsub("commonHMD", "common\nHMD", overlap$common)
overlap$common = gsub("commonPMD", "common\nPMD", overlap$common)
overlap$common = gsub("neither.", "other\n", overlap$common)
overlap$variable = gsub("Our", "LIBD\n", overlap$variable)
overlap$variable = gsub("common", "Zhou\n", overlap$variable)
colnames(overlap) = c("common", "LIBD", "Comparison", "Percent")
overlap$common = factor(overlap$common, levels = c("common\nPMD", "other\nPMD", "common\nHMD","other\nHMD"))

write.csv(overlap, quote=F,file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/PMD_Zhou_Overlap.csv")

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/CREs/figures/PMD_Zhou_Overlap.pdf", width = 10)
ggplot(overlap, aes(x = Comparison, y = Percent, fill = LIBD)) + 
  theme_classic() + geom_bar(stat = "identity") +
  labs(fill="") + facet_grid(common ~ LIBD) +
  ylab("Percent") + ylim(0,100) + xlab("") + 
  scale_fill_brewer(8, palette="Dark2") +
  ggtitle("Percent PMD Bases Replicated") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20), legend.title=element_blank()) + theme(legend.position="bottom")
dev.off()


## Individual sample overlap

overlap = lapply(pmds, function(cell) lapply(solo, function(com) c( OurOverlap = round(sum(width(Reduce(intersect, list(cell,com))))/sum(width(cell))*100,1),
                                                                        commonOverlap = round(sum(width(Reduce(intersect, list(cell,com))))/sum(width(com))*100,1))))
overlap = lapply(overlap, function(x) do.call(rbind, Map(cbind, lapply(x, function(y) data.frame(OurOverlap = y["OurOverlap"], commonOverlap = y["commonOverlap"])),
                                                         common = as.list(names(x)))))
overlap = do.call(rbind, Map(cbind, overlap, LIBD = as.list(names(pmds))))
overlap = reshape2::melt(overlap)
overlap$common = gsub("commonHMD", "common\nHMD", overlap$common)
overlap$common = gsub("commonPMD", "common\nPMD", overlap$common)
overlap$common = gsub("neither.", "other\n", overlap$common) 
overlap$variable = gsub("Our", "LIBD\n", overlap$variable)
overlap$variable = gsub("common", "Zhou\n", overlap$variable)
colnames(overlap) = c("common", "LIBD", "Comparison", "Percent")
overlap$common = factor(overlap$common, levels = c("common\nPMD", "other\nPMD", "common\nHMD","other\nHMD"))
overlap[which(overlap$LIBD %in% pd[pd$Cell.Type=="Neuron","Data.ID"]),"CellType"] = "Neuron"
overlap[which(overlap$LIBD %in% pd[pd$Cell.Type=="Glia","Data.ID"]),"CellType"] = "Glia"
overlap[-which(overlap$LIBD %in% pd[pd$Cell.Type %in% c("Neuron","Glia"),"Data.ID"]),"CellType"] = "Prenatal"

write.csv(overlap, quote=F,file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/PMD_bySample_Zhou_Overlap.csv")

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/CREs/figures/PMD_bySample_Zhou_Overlap.pdf", width = 8, height=4)
ggplot(overlap, aes(x = common, y = Percent, fill = CellType)) + 
  geom_boxplot() +
  labs(fill="") + facet_grid(. ~ Comparison) +
  ylab("Percent") + ylim(0,100) + xlab("") + 
  scale_fill_brewer(8, palette="Dark2") +
  ggtitle("Percent PMD Bases Replicated") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20), legend.title=element_blank()) + theme(legend.position="bottom")
dev.off()


