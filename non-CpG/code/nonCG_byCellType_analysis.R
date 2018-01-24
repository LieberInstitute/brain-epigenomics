library(data.table)
library(ggplot2)
library(plyr)

load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/non-CpG/CH_object.rda")
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/non-CpG/CHneurons_object.rda")
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/methylSeekR_objects.rda")

load("/media/Backup1_/amanda/CREs/methylSeekR_objects.rda")
load("/media/Backup1_/amanda/non-CpG/CH_object.rda")


table(CH$CT.dir=="pos")
# FALSE     TRUE 
# 7630202 33188540 
# 18.7% of CH is greater methylated in glia

table(CH[which(CH$padj.CellType<=0.05),"CT.dir"]=="pos")
#  FALSE    TRUE 
#  79239 7602836
# 99% of all significantly regulated sites by cell type are more methylated in neurons

table(CH$Age.dir=="pos")
#FALSE     TRUE 
#9843075 30975667
#75.9% are increasing over age

table(CH[which(CH$padj.Age<=0.05),"CT.dir"]=="pos")
#FALSE    TRUE 
#24967 3169651 
# 99.2% are increasing over age in significantly dmCH sites.


## Subset by trinucleotids context

CHdt = data.table(CH)
CHneuronsdt = data.table(CHneurons)
res = list("By Cell Type" = CHdt[,length(unique(regionID)), by=c("trinucleotide_context", "CT.sig")],
		   "By Age" = CHdt[,length(unique(regionID)), by=c("trinucleotide_context", "Age.sig")],
		   "By Age in a Cell Type" = CHdt[,length(unique(regionID)), by=c("trinucleotide_context", "Int.sig")],
		   "By Age in Neurons" = CHneuronsdt[,length(unique(regionID)), by=c("trinucleotide_context", "sig")],
		   "By Cell Type" = CHdt[CT.sig=="FDR < 0.05",length(unique(regionID)), by=c("trinucleotide_context", "CT.dir")],
		   "By Age" = CHdt[Age.sig=="FDR < 0.05",length(unique(regionID)), by=c("trinucleotide_context", "Age.dir")],
		   "By Age in a Cell Type" = CHdt[Int.sig=="FDR < 0.05",length(unique(regionID)), by=c("trinucleotide_context", "Int.dir")],
		   "By Age in Neurons" = CHneuronsdt[padj<=0.05,length(unique(regionID)), by=c("trinucleotide_context", "Dir")])
context = c("CAG", "CAC", "CTG", "CAT", "CAA", "CTC", "CCA", "CTT", "CTA", "CCT", "CCC", "CCG")
for (i in 1:length(res)) {
	for (j in 1:length(context)) {
		res[[i]][which(res[[i]][,"trinucleotide_context"] == context[j]),"perc"] = round(res[[i]][which(res[[i]][,"trinucleotide_context"] == context[j]),"V1"]/
		sum(res[[i]][which(res[[i]][,"trinucleotide_context"] == context[j]),"V1"])*100,1)
	}
}
for (i in 1:length(res)) {
	res[[i]]$trinucleotide_context = factor(res[[i]]$trinucleotide_context, levels = c("CAG", "CAC", "CTG", "CAT", "CAA", "CTC", "CCA", "CTT", "CTA", "CCT", "CCC", "CCG"))
}
res = lapply(res, data.frame)
res[[5]]$CT.dir = ifelse(res[[5]]$CT.dir=="pos", "Hypermethylated\nIn Neurons", "Hypermethylated\nIn Glia")
res[[6]]$Age.dir = ifelse(res[[6]]$Age.dir=="pos", "Increasingly\nMethylated", "Decreasingly\nMethylated")
res[[8]]$Dir = ifelse(res[[8]]$Dir=="pos", "Increasingly\nMethylated", "Decreasingly\nMethylated")
res <- lapply(res, function(x) ddply(x, .(trinucleotide_context), transform, pos = cumsum(V1) - (0.5 * V1)))
for (j in 1:length(context)) {
		res[[8]][which(res[[8]][,"trinucleotide_context"] == context[j]),"perc"] = round(res[[8]][which(res[[8]][,"trinucleotide_context"] == context[j]),"V1"]/
		sum(res[[8]][which(res[[8]][,"trinucleotide_context"] == context[j]),"V1"])*100,1)
}


pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/non-CpG/figures/nonCG_byContext.pdf", width = 10)
for (i in 1:length(res)) {
g = ggplot(res[[i]], aes(x = trinucleotide_context, y = V1)) + geom_bar(aes(fill=res[[i]][,2]),stat = "identity") +
  geom_text(aes(x = trinucleotide_context, y = pos, label = paste0(perc,"%"))) +
  labs(fill="") +
  ylab("Count") + 
  xlab("Context") +
  ggtitle(paste0("non-CpG Context: ", names(res)[i])) +
  theme(axis.text.x=element_text(angle=45,hjust=1)) +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20)) +
  theme(legend.position="bottom", legend.direction="horizontal", legend.title = element_blank())
  print(g)
 }
dev.off()

x = res[[4]]
x$perc = round(x$V1/sum(x[which(x$sig=="FDR < 0.05"),"V1"])*100, 1)
x$perc = paste0(x$perc,"%")
x[which(x$sig=="FDR > 0.05"),"perc"] = ""

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/non-CpG/figures/nonCG_byContext_byAge_inNeurons.pdf", width = 10)
  ggplot(x, aes(x = trinucleotide_context, y = V1), fill=sig) + geom_bar(stat = "identity") +
  geom_text(label = perc, vjust = -.5) +
  labs(fill="") +
  ylab("Count") + 
  xlab("Context") +
  ggtitle("non-CpG Context: By Age in Neurons") +
  theme(axis.text.x=element_text(angle=45,hjust=1)) +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20)) +
  theme(legend.position="bottom", legend.direction="horizontal", legend.title = element_blank())
dev.off()


pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/non-CpG/figures/nonCG_byContext_byAge_inNeurons.pdf", width = 10)
  ggplot(x, aes(x = trinucleotide_context, y = V1), fill=sig) + geom_bar(stat = "identity") +
  geom_text(label = perc, vjust = -.5) +
  labs(fill="") +
  ylab("Count") + 
  xlab("Context") +
  ggtitle("non-CpG Context: By Age in Neurons") +
  theme(axis.text.x=element_text(angle=45,hjust=1)) +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20)) +
  theme(legend.position="bottom", legend.direction="horizontal", legend.title = element_blank())
dev.off()


## Check coverage by CH context

covdf = lapply(CHlist, as.data.frame)
covdf = do.call(rbind, Map(cbind, covdf, ID = as.list(names(covdf))))
covdt = data.table(covdf)
covdt = covdt[,mean(T), by=c("ID", )]


