library(ggplot2)
library(GenomicRanges)


load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/DMR/DMR_objects.rda")
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/bumphunting/rda/limma_Neuron_CpGs_minCov_3_ageInfo_dmrs.Rdata")


## Get DMRs

dmrs = split(dmrs, dmrs$k6cluster_label)
names(dmrs) = c("Group 1 (G-N+)","Group 2 (G0N+)","Group 3 (G0N-)","Group 4 (G+N0)","Group 5 (G+N-)","Group 6 (G-N0)")
oo = lapply(dmrs, function(x) findOverlaps(x, makeGRangesFromDataFrame(DMR$Interaction)))
dmrs = lapply(oo, function(x) DMR$Interaction[subjectHits(x),])
DMRgr = lapply(c(DMR, dmrs), function(x) makeGRangesFromDataFrame(x[which(x$fwer<=0.05),], keep.extra.columns = T))
DMRgr = lapply(DMRgr, function(x) unique(granges(x)))


## Get single cell CpG DMRs

scDMRs = list()
for (i in 1:21) {
  scDMRs[[i]] = openxlsx::read.xlsx(xlsxFile = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/Luo_Ecker_scMethylSeq_TableS6.xlsx", sheet = i, startRow = 2)
}
names(scDMRs) = lapply(scDMRs, function(x) colnames(x)[1])
for (i in 1:21) { 
  colnames(scDMRs[[i]]) = c("chromosome","start","end")
  scDMRs[[i]] = scDMRs[[i]][-1,]
}
scDMRs = do.call(rbind, scDMRs)
scDMRs = makeGRangesFromDataFrame(scDMRs)


## Identify overlap with our DMRs and scDMRs

# Total bases

scDMRs = reduce(scDMRs)
DMRgr = lapply(DMRgr, reduce)

overlap = lapply(DMRgr, function(cell) c( OurOverlap = round(sum(width(Reduce(intersect, list(cell,scDMRs))))/sum(width(cell))*100,1),
                                          commonOverlap = round(sum(width(Reduce(intersect, list(cell,scDMRs))))/sum(width(scDMRs))*100,1)))

overlap = do.call(rbind, Map(cbind, Group = as.list(names(overlap)), lapply(overlap, function(y) data.frame(OurOverlap = y["OurOverlap"], scOverlap = y["commonOverlap"]))))

overlap = reshape2::melt(overlap)
overlap$variable = gsub("Our", "LIBD ", overlap$variable)
overlap$variable = gsub("sc", "Single Cell ", overlap$variable)
colnames(overlap) = c("Group","Comparison", "Percent")

write.csv(overlap, quote=F,file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/DMR/DMR_scMethyl_DMR_overlap.csv")

overlap$Group = gsub(" (", "\n(",overlap$Group, fixed=T)


pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/DMR/figures/DMR_scMethyl_DMR_overlap.pdf", width = 8, height = 6)
ggplot(overlap[grep("Group",overlap$Group)[grep("Group",overlap$Group) %in% grep("LIBD",overlap$Comparison)],] , aes(x = Group, y = Percent, fill = Group)) + 
  theme_classic() + geom_bar(stat = "identity") + labs(fill="") +
  ylab("Percent") + ylim(0,100) + xlab("") + 
  scale_fill_brewer(8, palette="Dark2") +
  ggtitle("Percent Neuronal Subtype-specific Bases") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20), legend.title=element_blank()) + theme(legend.position="bottom")
dev.off()





