library('RColorBrewer')
library(ggplot2)

load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/DMR/DMR_objects.rda")
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/DMR/Age/DMR_KEGG_GO_DO_objects_byAge.rda")
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/processed_beta_values_plusMap.rda")

x = goList_BP.dir[[4]]
pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/DMR/Age/figures/promoter_GO_BP_plot_byAge.pdf",width = 8, height = 8)
dotplot(x, showCategory =10)
dev.off()

select = dtage[sig=="FWER < 0.05" & nearestSymbol%in%c("SHH","TET3","SMAD9","VANGL2","DCHS1"),,]
coord = as.list(select$regionID)
names(coord) = c("VANGL2a","DCHS1","SHH","TET3","VANGL2b","SMAD9")
ov = lapply(coord, function(x) findOverlaps(gr,GRanges(x)))
means = c(lapply(ov[1:5], function(x) colMeans(meth[queryHits(x),])), list(SMAD9 = meth[queryHits(ov[[6]]),]))
names(means[[1]])
lapply(means, function(x) match(rownames(pd), names(x)))

means = lapply(means, function(x) data.frame(x, pd$Age, pd$Cell.Type))

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/DMR/Age/figures/interesting_DMRgenes_byAge_MethByAge_scatterplot.pdf")
for (i in 1:length(means)){
	colnames(means[[i]]) = c("Methylation","Age","CellType")
	g = ggplot(means[[i]], aes(x = Age, y = Methylation, color=CellType)) + geom_point() +
	geom_smooth(method=lm) + scale_color_brewer(palette = "Dark2") +
  	labs(fill="") +
  	ylab("Mean Methylation") + 
  	xlab("Age") +
  	ggtitle(paste0("Age DMR in ", names(means)[i],", r=", round(cor(x=means[[i]][,"Methylation"], y =means[[i]][,"Age"]), 2))) +
  	theme(title = element_text(size = 20)) +
  	theme(text = element_text(size = 20))
  	print(g)
}
dev.off()


# How many of the "Overall Age" DMRs overlap the "Age in a Cell Type" DMRs?

fwer = lapply(DMR, function(x) x[which(x$fwer<=0.05),])
fwer = lapply(fwer, function(x) makeGRangesFromDataFrame(x, keep.extra.columns=T))
fwer_df = lapply(fwer, as.data.frame)
hits = findOverlaps(fwer$Interaction, fwer$Age)
age = cbind(fwer_df$Age[subjectHits(hits),], fwer_df$Interaction[queryHits(hits),])
length(unique(subjectHits(hits))) # 59 of the 129 "Overall Age" DMRs overlap an "Interaction" DMR

overlaps <- pintersect(fwer$Interaction[queryHits(hits)], fwer$Age[subjectHits(hits)])
percentOverlap <- width(overlaps) / width(fwer$Age[subjectHits(hits)])
c(mean = mean(percentOverlap), sd = sd(percentOverlap), min = min(percentOverlap), max = max(percentOverlap))
#      mean         sd        min        max 
# 0.62571268 0.28743417 0.09381592 1.00000000