library(GenomicRanges)
library(GenomicFeatures)
library(data.table)
library(ggplot2)
library(clusterProfiler)
require(org.Hs.eg.db)

load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/non-CpG/CH_object.rda")


### Both Neurons and Glia ###

## How many sites are DM?

dim(CHres) # 40818742       19
dim(CHres[which(CHres$padj.CellType<=0.05),]) # 7682075      19 
dim(CHres[which(CHres$padj.Age<=0.05),]) # 3194618      19
dim(CHres[which(CHres$padj.Interaction<=0.05),]) # 76 19

dtCH = data.table(CH)

### Explore annotation of regions

## how many fall within CpG islands?

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/non-CpG/figures/nonCG_overlap_with_CpG_Islands.pdf")
x = dtCH[,length(unique(regionID)), by = "islands"]
x$perc = round(x$V1/sum(x$V1)*100,2)
ggplot(x, aes(x = islands, y = V1)) + geom_bar(stat = "identity") +
  geom_text(aes( label = paste0(perc,"%")), vjust = -.5) +
  labs(fill="") +
  ylab("Count") + 
  xlab("") +
  ggtitle("non-CpGs Overlapping CpG Islands:\nAll non-CpGs") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
x = dtCH[CT.sig=="FDR < 0.05",length(unique(regionID)), by = "islands"]
x$perc = round(x$V1/sum(x$V1)*100,2)
ggplot(x, aes(x = islands, y = V1)) + geom_bar(stat = "identity") +
  geom_text(aes( label = paste0(perc,"%")), vjust = -.5) +
  labs(fill="") +
  ylab("Count") + 
  xlab("") +
  ggtitle("non-CpGs Overlapping CpG Islands:\nFDR < 0.05 by Cell Type") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
x = dtCH[Age.sig=="FDR < 0.05",length(unique(regionID)), by = "islands"]
x$perc = round(x$V1/sum(x$V1)*100,2)
ggplot(x, aes(x = islands, y = V1)) + geom_bar(stat = "identity") +
  geom_text(aes( label = paste0(perc,"%")), vjust = -.5) +
  labs(fill="") +
  ylab("Count") + 
  xlab("") +
  ggtitle("non-CpGs Overlapping CpG Islands:\nFDR < 0.05 by Age") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
x = dtCH[Int.sig=="FDR < 0.05",length(unique(regionID)), by = "islands"]
x$perc = round(x$V1/sum(x$V1)*100,2)
ggplot(x, aes(x = islands, y = V1)) + geom_bar(stat = "identity") +
  geom_text(aes( label = paste0(perc,"%")), vjust = -.5) +
  labs(fill="") +
  ylab("Count") + 
  xlab("") +
  ggtitle("non-CpGs Overlapping CpG Islands:\nFDR < 0.05 by Cell Type and Age") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
x = dtCH[,length(unique(regionID)), by = c("islands", "CT.sig")]
x$perc = NA
x$perc[grep("<", x$CT.sig)] = unlist(c(round(x[grep("<", x$CT.sig),"V1"]/sum(x[grep("<", x$CT.sig),"V1"])*100,2)))
x$perc[grep(">", x$CT.sig)] = unlist(c(round(x[grep(">", x$CT.sig),"V1"]/sum(x[grep(">", x$CT.sig),"V1"])*100,2)))
ggplot(x, aes(x = islands, y = V1)) + geom_bar(stat = "identity") +
  geom_text(aes( label = paste0(perc,"%")), vjust = -.5) +
  facet_grid(. ~ CT.sig) +
  labs(fill="") +
  ylab("Count") + 
  xlab("") +
  ggtitle("non-CpGs Overlapping CpG Islands:\nCell Type") +
  theme(axis.text.x=element_text(angle=45,hjust=1)) +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
x = dtCH[,length(unique(regionID)), by = c("islands", "Age.sig")]
x$perc = NA
x$perc[grep("<", x$Age.sig)] = unlist(c(round(x[grep("<", x$Age.sig),"V1"]/sum(x[grep("<", x$Age.sig),"V1"])*100,2)))
x$perc[grep(">", x$Age.sig)] = unlist(c(round(x[grep(">", x$Age.sig),"V1"]/sum(x[grep(">", x$Age.sig),"V1"])*100,2)))
ggplot(x, aes(x = islands, y = V1)) + geom_bar(stat = "identity") +
  geom_text(aes( label = paste0(perc,"%")), vjust = -.5) +
  facet_grid(. ~ Age.sig) +
  labs(fill="") +
  ylab("Count") + 
  xlab("") +
  ggtitle("non-CpGs Overlapping CpG Islands:\nAge") +
  theme(axis.text.x=element_text(angle=45,hjust=1)) +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
x = dtCH[,length(unique(regionID)), by = c("islands", "Int.sig")]
x$perc = NA
x$perc[grep("<", x$Int.sig)] = unlist(c(round(x[grep("<", x$Int.sig),"V1"]/sum(x[grep("<", x$Int.sig),"V1"])*100,2)))
x$perc[grep(">", x$Int.sig)] = unlist(c(round(x[grep(">", x$Int.sig),"V1"]/sum(x[grep(">", x$Int.sig),"V1"])*100,2)))
ggplot(x, aes(x = islands, y = V1)) + geom_bar(stat = "identity") +
  geom_text(aes( label = paste0(perc,"%")), vjust = -.5) +
  facet_grid(. ~ Int.sig) +
  labs(fill="") +
  ylab("Count") + 
  xlab("") +
  ggtitle("non-CpGs Overlapping CpG Islands:\nInteraction") +
  theme(axis.text.x=element_text(angle=45,hjust=1)) +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
dev.off()

# Is there a relationship between being significantly DM and overlapping a CpG island?

fisher.test(data.frame(c(nrow(CH[which(CH$CT.sig=="FDR < 0.05" & CH$islands=="CpG_Island"),]),
                         nrow(CH[which(CH$CT.sig=="FDR > 0.05" & CH$islands=="CpG_Island"),])),
                       c(nrow(CH[which(CH$CT.sig=="FDR < 0.05" & CH$islands=="non-Island"),]),
                         nrow(CH[which(CH$CT.sig=="FDR > 0.05" & CH$islands=="non-Island"),]))))
# CpG islands are depleted in dmCH by cell type
#p-value < 2.2e-16
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.150302 0.155692
#sample estimates:
#  odds ratio 
#0.1529949
fisher.test(data.frame(c(nrow(CH[which(CH$Age.sig=="FDR < 0.05" & CH$islands=="CpG_Island"),]),
                         nrow(CH[which(CH$Age.sig=="FDR > 0.05" & CH$islands=="CpG_Island"),])),
                       c(nrow(CH[which(CH$Age.sig=="FDR < 0.05" & CH$islands=="non-Island"),]),
                         nrow(CH[which(CH$Age.sig=="FDR > 0.05" & CH$islands=="non-Island"),]))))
# CpG islands are depleted in dmCH by Age
#p-value < 2.2e-16
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.1581654 0.1671672
#sample estimates:
#  odds ratio 
#0.1626362
fisher.test(data.frame(c(nrow(CH[which(CH$Int.sig=="FDR < 0.05" & CH$islands=="CpG_Island"),]),
                         nrow(CH[which(CH$Int.sig=="FDR > 0.05" & CH$islands=="CpG_Island"),])),
                       c(nrow(CH[which(CH$Int.sig=="FDR < 0.05" & CH$islands=="non-Island"),]),
                         nrow(CH[which(CH$Int.sig=="FDR > 0.05" & CH$islands=="non-Island"),]))))
# No relationship between location of CpG islands and dmCH by Age and Cell type
#p-value = 1
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.000000 5.404446
#sample estimates:
#  odds ratio 
#0


# assign genomic features

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/non-CpG/figures/nonCG_annotation.pdf")
x = dtCH[,length(unique(regionID)), by = "annotation"]
x$perc = round(x$V1/sum(x$V1)*100,2)
ggplot(x, aes(x = annotation, y = V1)) + geom_bar(stat = "identity") +
  geom_text(aes( label = paste0(perc,"%")), vjust = -.5) +
  labs(fill="") +
  ylab("Count") + 
  xlab("") +
  ggtitle("non-CpG annotation: All non-CpGs") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
x = dtCH[CT.sig=="FDR < 0.05",length(unique(regionID)), by = "annotation"]
x$perc = round(x$V1/sum(x$V1)*100,2)
ggplot(x, aes(x = annotation, y = V1)) + geom_bar(stat = "identity") +
  geom_text(aes( label = paste0(perc,"%")), vjust = -.5) +
  labs(fill="") +
  ylab("Count") + 
  xlab("") +
  ggtitle("non-CpG annotation:\nFDR < 0.05 by Cell Type") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
x = dtCH[Age.sig=="FDR < 0.05",length(unique(regionID)), by = "annotation"]
x$perc = round(x$V1/sum(x$V1)*100,2)
ggplot(x, aes(x = annotation, y = V1)) + geom_bar(stat = "identity") +
  geom_text(aes( label = paste0(perc,"%")), vjust = -.5) +
  labs(fill="") +
  ylab("Count") + 
  xlab("") +
  ggtitle("non-CpG annotation:\nFDR < 0.05 by Age") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
x = dtCH[Int.sig=="FDR < 0.05",length(unique(regionID)), by = "annotation"]
x$perc = round(x$V1/sum(x$V1)*100,2)
ggplot(x, aes(x = annotation, y = V1)) + geom_bar(stat = "identity") +
  geom_text(aes( label = paste0(perc,"%")), vjust = -.5) +
  labs(fill="") +
  ylab("Count") + 
  xlab("") +
  ggtitle("non-CpG annotation:\nFDR < 0.05 by Cell Type and Age") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
x = dtCH[,length(unique(regionID)), by = c("annotation", "CT.sig")]
x$perc = NA
x$perc[grep("<", x$CT.sig)] = unlist(c(round(x[grep("<", x$CT.sig),"V1"]/sum(x[grep("<", x$CT.sig),"V1"])*100,2)))
x$perc[grep(">", x$CT.sig)] = unlist(c(round(x[grep(">", x$CT.sig),"V1"]/sum(x[grep(">", x$CT.sig),"V1"])*100,2)))
ggplot(x, aes(x = annotation, y = V1)) + geom_bar(stat = "identity") +
  geom_text(aes( label = paste0(perc,"%")), vjust = -.5) +
  facet_grid(. ~ sig) +
  labs(fill="") +
  ylab("Count") + 
  theme(axis.text.x=element_text(angle=45,hjust=1)) +
  xlab("") +
  ggtitle("non-CpG annotation:\nCell Type") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
x = dtCH[,length(unique(regionID)), by = c("annotation", "Age.sig")]
x$perc = NA
x$perc[grep("<", x$Age.sig)] = unlist(c(round(x[grep("<", x$Age.sig),"V1"]/sum(x[grep("<", x$Age.sig),"V1"])*100,2)))
x$perc[grep(">", x$Age.sig)] = unlist(c(round(x[grep(">", x$Age.sig),"V1"]/sum(x[grep(">", x$Age.sig),"V1"])*100,2)))
ggplot(x, aes(x = annotation, y = V1)) + geom_bar(stat = "identity") +
  geom_text(aes( label = paste0(perc,"%")), vjust = -.5) +
  facet_grid(. ~ sig) +
  labs(fill="") +
  ylab("Count") + 
  theme(axis.text.x=element_text(angle=45,hjust=1)) +
  xlab("") +
  ggtitle("non-CpG annotation:\nAge") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
x = dtCH[,length(unique(regionID)), by = c("annotation", "Int.sig")]
x$perc = NA
x$perc[grep("<", x$Int.sig)] = unlist(c(round(x[grep("<", x$Int.sig),"V1"]/sum(x[grep("<", x$Int.sig),"V1"])*100,2)))
x$perc[grep(">", x$Int.sig)] = unlist(c(round(x[grep(">", x$Int.sig),"V1"]/sum(x[grep(">", x$Int.sig),"V1"])*100,2)))
ggplot(x, aes(x = annotation, y = V1)) + geom_bar(stat = "identity") +
  geom_text(aes( label = paste0(perc,"%")), vjust = -.5) +
  facet_grid(. ~ sig) +
  labs(fill="") +
  ylab("Count") + 
  theme(axis.text.x=element_text(angle=45,hjust=1)) +
  xlab("") +
  ggtitle("non-CpG annotation:\nInteraction") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
dev.off()

# Is there a relationship between being significantly DM and overlapping a gene?

fisher.test(data.frame(c(nrow(CH[which(CH$CT.sig=="FDR < 0.05" & CH$distToGene==0),]),
                         nrow(CH[which(CH$CT.sig=="FDR > 0.05" & CH$distToGene==0),])),
                       c(nrow(CH[which(CH$CT.sig=="FDR < 0.05" & CH$distToGene>0),]),
                         nrow(CH[which(CH$CT.sig=="FDR > 0.05" & CH$distToGene>0),]))))
# Genes are underrepresented in dmCH by Cell Type
#p-value < 2.2e-16
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.8413600 0.8442556
#sample estimates:
#  odds ratio 
#0.8428441
fisher.test(data.frame(c(nrow(CH[which(CH$Age.sig=="FDR < 0.05" & CH$distToGene==0),]),
                         nrow(CH[which(CH$Age.sig=="FDR > 0.05" & CH$distToGene==0),])),
                       c(nrow(CH[which(CH$Age.sig=="FDR < 0.05" & CH$distToGene>0),]),
                         nrow(CH[which(CH$Age.sig=="FDR > 0.05" & CH$distToGene>0),]))))
# Genes are underrepresented in dmCH by Age
# p-value < 2.2e-16
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.8194398 0.8236559
#sample estimates:
#  odds ratio 
#0.8215612
fisher.test(data.frame(c(nrow(CH[which(CH$Int.sig=="FDR < 0.05" & CH$distToGene==0),]),
                         nrow(CH[which(CH$Int.sig=="FDR > 0.05" & CH$distToGene==0),])),
                       c(nrow(CH[which(CH$Int.sig=="FDR < 0.05" & CH$distToGene>0),]),
                         nrow(CH[which(CH$Int.sig=="FDR > 0.05" & CH$distToGene>0),]))))
# No relationship between location of genes and dmCH by Age and Cell type
#p-value = 0.8061
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.6283264 1.7398532
#sample estimates:
#  odds ratio 
#1.05829

# Is there a relationship between being significantly DM and overlapping a promoter?

fisher.test(data.frame(c(nrow(CH[which(CH$CT.sig=="FDR < 0.05" & CH$annotation=="Promoter"),]),
                         nrow(CH[which(CH$CT.sig=="FDR > 0.05" & CH$annotation=="Promoter"),])),
                       c(nrow(CH[which(CH$CT.sig=="FDR < 0.05" & CH$annotation!="Promoter"),]),
                         nrow(CH[which(CH$CT.sig=="FDR > 0.05" & CH$annotation!="Promoter"),]))))
# promoters are underrepresented in dmCH by Cell Type
#p-value < 2.2e-16
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.8529202 0.8629646
#sample estimates:
#  odds ratio 
#0.8579478
fisher.test(data.frame(c(nrow(CH[which(CH$Age.sig=="FDR < 0.05" & CH$annotation=="Promoter"),]),
                         nrow(CH[which(CH$Age.sig=="FDR > 0.05" & CH$annotation=="Promoter"),])),
                       c(nrow(CH[which(CH$Age.sig=="FDR < 0.05" & CH$annotation!="Promoter"),]),
                         nrow(CH[which(CH$Age.sig=="FDR > 0.05" & CH$annotation!="Promoter"),]))))
# promoters are underrepresented in dmCH by Age
#p-value < 2.2e-16
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.8975606 0.9129173
#sample estimates:
#  odds ratio 
#0.9051968 
fisher.test(data.frame(c(nrow(CH[which(CH$Int.sig=="FDR < 0.05" & CH$annotation=="Promoter"),]),
                         nrow(CH[which(CH$Int.sig=="FDR > 0.05" & CH$annotation=="Promoter"),])),
                       c(nrow(CH[which(CH$Int.sig=="FDR < 0.05" & CH$annotation!="Promoter"),]),
                         nrow(CH[which(CH$Int.sig=="FDR > 0.05" & CH$annotation!="Promoter"),]))))
# No relationship between location of promoters and dmCH by Cell Type and Age
#p-value = 0.198
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.4026757 6.0727879
#sample estimates:
#  odds ratio 
#1.997711

# Is there a relationship between being significantly DM and overlapping a gene and/or promoter?

fisher.test(data.frame(c(nrow(CH[which(CH$CT.sig=="FDR < 0.05" & CH$annotation=="Other"),]),
                         nrow(CH[which(CH$CT.sig=="FDR > 0.05" & CH$annotation=="Other"),])),
                       c(nrow(CH[which(CH$CT.sig=="FDR < 0.05" & CH$annotation!="Other"),]),
                         nrow(CH[which(CH$CT.sig=="FDR > 0.05" & CH$annotation!="Other"),]))))
# Genes/promoters are underrepresented in dmCH by Cell Type
#p-value < 2.2e-16
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  1.191807 1.195925
#sample estimates:
#  odds ratio 
#1.193906
fisher.test(data.frame(c(nrow(CH[which(CH$Age.sig=="FDR < 0.05" & CH$annotation=="Other"),]),
                         nrow(CH[which(CH$Age.sig=="FDR > 0.05" & CH$annotation=="Other"),])),
                       c(nrow(CH[which(CH$Age.sig=="FDR < 0.05" & CH$annotation!="Other"),]),
                         nrow(CH[which(CH$Age.sig=="FDR > 0.05" & CH$annotation!="Other"),]))))
# Genes/promoters are underrepresented in dmCH by Age
#p-value < 2.2e-16
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  1.217999 1.224176
#sample estimates:
#  odds ratio 
#1.221127  
fisher.test(data.frame(c(nrow(CH[which(CH$Int.sig=="FDR < 0.05" & CH$annotation=="Other"),]),
                         nrow(CH[which(CH$Int.sig=="FDR > 0.05" & CH$annotation=="Other"),])),
                       c(nrow(CH[which(CH$Int.sig=="FDR < 0.05" & CH$annotation!="Other"),]),
                         nrow(CH[which(CH$Int.sig=="FDR > 0.05" & CH$annotation!="Other"),]))))
# No relationship between location of genes/promoters and dmCH by Cell Type and Age
# p-value = 0.4659
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.5174125 1.3930275
#sample estimates:
#  odds ratio 
#0.8415662

### Gene Ontology: Cell Type

CTentrez = list(All = dtCH[CT.sig=="FDR < 0.05",list(na.omit(EntrezID)),], 
              GenesPlusPromoters = dtCH[annotation != "Other" & CT.sig=="FDR < 0.05",list(na.omit(EntrezID)),],
              Genes = dtCH[distToGene==0 & CT.sig=="FDR < 0.05",list(na.omit(EntrezID)),],
              Promoters = dtCH[annotation == "Promoter" & CT.sig=="FDR < 0.05",list(na.omit(EntrezID)),])
CTentrez.dir = list(All.pos = dtCH[Tstat.CellType>0 & CT.sig=="FDR < 0.05",list(na.omit(EntrezID)),], 
                  GenesPlusPromoters.pos = dtCH[Tstat.CellType>0 & annotation != "Other" & CT.sig=="FDR < 0.05",list(na.omit(EntrezID)),],
                  Genes.pos = dtCH[Tstat.CellType>0 & distToGene==0 & CT.sig=="FDR < 0.05",list(na.omit(EntrezID)),],
                  Promoters.pos = dtCH[Tstat.CellType>0 & annotation == "Promoter" & CT.sig=="FDR < 0.05",list(na.omit(EntrezID)),],
                  All.neg = dtCH[Tstat.CellType<0 & CT.sig=="FDR < 0.05",list(na.omit(EntrezID)),], 
                  GenesPlusPromoters.neg = dtCH[Tstat.CellType<0 & annotation != "Other" & CT.sig=="FDR < 0.05",list(na.omit(EntrezID)),],
                  Genes.neg = dtCH[Tstat.CellType<0 & distToGene==0 & CT.sig=="FDR < 0.05",list(na.omit(EntrezID)),],
                  Promoters.neg = dtCH[Tstat.CellType<0 & annotation == "Promoter" & CT.sig=="FDR < 0.05",list(na.omit(EntrezID)),])
CTentrez = lapply(CTentrez, function(x) as.character(unique(x$V1)))              
CTentrez.dir = lapply(CTentrez.dir, function(x) as.character(unique(x$V1)))       

GeneUniverse = dtCH[,list(na.omit(EntrezID)),]
GeneUniverse = as.character(unique(GeneUniverse$V1))

# Find enriched Pathways via KEGG
elementNROWS(CTentrez)
CTkeggList = lapply(CTentrez, function(x) enrichKEGG(as.character(x), organism="human", universe= GeneUniverse, 
                                                   minGSSize=5, pAdjustMethod="BH", qvalueCutoff=1))
CTkeggListdf = lapply(CTkeggList, function(x) as.data.frame(x))
CTkeggList.dir = lapply(CTentrez.dir, function(x) enrichKEGG(as.character(x), organism="human", universe= GeneUniverse, 
                                                   minGSSize=5, pAdjustMethod="BH", qvalueCutoff=1))
CTkeggList.dir.df = lapply(CTkeggList.dir, function(x) as.data.frame(x))

# Enriched Molecular Function GOs
CTgoList_MF = lapply(CTentrez, function(x) enrichGO(as.character(x), ont = "MF", OrgDb = org.Hs.eg.db, 
                                                  universe= GeneUniverse, minGSSize=5, pAdjustMethod="BH",
                                                  qvalueCutoff=1))
CTgoListdf_MF = lapply(CTgoList_MF, function(x) as.data.frame(x))
CTgoList_MF.dir = lapply(CTentrez.dir, function(x) enrichGO(as.character(x), ont = "MF", OrgDb = org.Hs.eg.db, 
                                                universe= GeneUniverse, minGSSize=5, pAdjustMethod="BH",
                                                qvalueCutoff=1))
CTgoListdf_MF.dir = lapply(CTgoList_MF.dir, function(x) as.data.frame(x))

# Biological Process GO enrichment
CTgoList_BP = lapply(CTentrez, function(x) enrichGO(as.character(x), ont = "BP", OrgDb = org.Hs.eg.db, 
                                                universe= GeneUniverse, minGSSize=5, pAdjustMethod="BH",
                                                qvalueCutoff=1))
CTgoListdf_BP = lapply(CTgoList_BP, function(x) as.data.frame(x))
CTgoList_BP.dir = lapply(CTentrez.dir, function(x) enrichGO(as.character(x), ont = "BP", OrgDb = org.Hs.eg.db, 
                                                        universe= GeneUniverse, minGSSize=5, pAdjustMethod="BH",
                                                        qvalueCutoff=1))
CTgoListdf_BP.dir = lapply(CTgoList_BP.dir, function(x) as.data.frame(x))

# Cellular Compartment GO enrichment
CTgoList_CC = lapply(CTentrez, function(x) enrichGO(as.character(x), ont = "CC", OrgDb = org.Hs.eg.db, 
                                                universe= GeneUniverse, minGSSize=5, pAdjustMethod="BH",
                                                qvalueCutoff=1))
CTgoListdf_CC = lapply(CTgoList_CC, function(x) as.data.frame(x))
CTgoList_CC.dir = lapply(CTentrez.dir, function(x) enrichGO(as.character(x), ont = "CC", OrgDb = org.Hs.eg.db, 
                                                        universe= GeneUniverse, minGSSize=5, pAdjustMethod="BH",
                                                        qvalueCutoff=1))
CTgoListdf_CC.dir = lapply(CTgoList_CC.dir, function(x) as.data.frame(x))

# Disease Ontology
CTgoList_DO = lapply(CTentrez, function(x) enrichDO(as.character(x), ont = "DO", universe= GeneUniverse, 
                                                  minGSSize=5, pAdjustMethod="BH", qvalueCutoff=1, readable=TRUE))
CTgoListdf_DO = lapply(CTgoList_DO, function(x) as.data.frame(x))
CTgoList_DO.dir = lapply(CTentrez.dir, function(x) enrichDO(as.character(x), ont = "DO", universe= GeneUniverse, 
                                                minGSSize=5, pAdjustMethod="BH", qvalueCutoff=1, readable=TRUE))
CTgoListdf_DO.dir = lapply(CTgoList_DO.dir, function(x) as.data.frame(x))

# Compare the enriched terms between 7 groups
# KEGG
CTcompareKegg = compareCluster(CTentrez, fun="enrichKEGG", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
CTcompareKegg.dir = compareCluster(CTentrez.dir, fun="enrichKEGG", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
# Biological Process
CTcompareBP = compareCluster(CTentrez, fun="enrichGO", ont = "BP", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
CTcompareBP.dir = compareCluster(CTentrez.dir, fun="enrichGO",  ont = "BP", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
# Molecular Function
CTcompareMF = compareCluster(CTentrez, fun="enrichGO",  ont = "MF", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
CTcompareMF.dir = compareCluster(CTentrez.dir, fun="enrichGO",  ont = "MF", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
# Cellular Component
CTcompareCC = compareCluster(CTentrez, fun="enrichGO",  ont = "CC", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
CTcompareCC.dir = compareCluster(CTentrez.dir, fun="enrichGO",  ont = "CC", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
# Disease Ontology
CTcompareDO = compareCluster(CTentrez, fun="enrichDO",  ont = "DO", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
CTcompareDO.dir = compareCluster(CTentrez.dir, fun="enrichDO",  ont = "DO", qvalueCutoff = 0.05, pvalueCutoff = 0.05)

# save, write to csv
save(CTcompareKegg, CTcompareKegg.dir, CTcompareBP, CTcompareBP.dir, CTcompareMF, CTcompareMF.dir, CTcompareCC, CTcompareCC.dir, CTcompareDO, CTcompareDO.dir,
     CTkeggList, CTkeggList.dir, CTgoList_BP, CTgoList_BP.dir, CTgoList_MF, CTgoList_MF.dir, CTgoList_CC, CTgoList_CC.dir, CTgoList_DO, CTgoList_DO.dir,
     file="/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/non-CpG/nonCG_KEGG_GO_DO_objects_byCellType.rda")

# plot compared results
pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/non-CpG/figures/nonCpG_KEGG_GO_DO_plots_byCellType.pdf", height = 20, width = 20)
plot(CTcompareKegg, colorBy="p.adjust", showCategory = 45, title= "non-CpG KEGG Pathway Enrichment\nby Cell Type")
plot(CTcompareBP, colorBy="p.adjust", showCategory = 45, title= "non-CpG Biological Process GO Enrichment\nby Cell Type")
plot(CTcompareMF, colorBy="p.adjust", showCategory = 45, title= "non-CpG Molecular Function GO Enrichment\nby Cell Type")
plot(CTcompareCC, colorBy="p.adjust", showCategory = 45, title= "non-CpG Cellular Compartment GO Enrichment\nby Cell Type")
plot(CTcompareDO, colorBy="p.adjust", showCategory = 30, title= "non-CpG Disease Ontology Enrichment\nby Cell Type")
plot(CTcompareKegg.dir, colorBy="p.adjust", showCategory = 45, title= "non-CpG KEGG Pathway Enrichment\nby Cell Type")
plot(CTcompareBP.dir, colorBy="p.adjust", showCategory = 45, title= "non-CpG Biological Process GO Enrichment\nby Cell Type")
plot(CTcompareMF.dir, colorBy="p.adjust", showCategory = 45, title= "non-CpG Molecular Function GO Enrichment\nby Cell Type")
plot(CTcompareCC.dir, colorBy="p.adjust", showCategory = 45, title= "non-CpG Cellular Compartment GO Enrichment\nby Cell Type")
plot(CTcompareDO.dir, colorBy="p.adjust", showCategory = 30, title= "non-CpG Disease Ontology Enrichment\nby Cell Type")
dev.off()


### Gene Ontology: Age
Ageentrez = list(All = dtCH[Age.sig=="FDR < 0.05",list(na.omit(EntrezID)),], 
                 GenesPlusPromoters = dtCH[annotation != "Other" & Age.sig=="FDR < 0.05",list(na.omit(EntrezID)),],
                 Genes = dtCH[distToGene==0 & Age.sig=="FDR < 0.05",list(na.omit(EntrezID)),],
                 Promoters = dtCH[annotation == "Promoter" & Age.sig=="FDR < 0.05",list(na.omit(EntrezID)),])
Ageentrez.dir = list(All.pos = dtCH[Tstat.Age>0 & Age.sig=="FDR < 0.05",list(na.omit(EntrezID)),], 
                     GenesPlusPromoters.pos = dtCH[Tstat.Age>0 & annotation != "Other" & Age.sig=="FDR < 0.05",list(na.omit(EntrezID)),],
                     Genes.pos = dtCH[Tstat.Age>0 & distToGene==0 & Age.sig=="FDR < 0.05",list(na.omit(EntrezID)),],
                     Promoters.pos = dtCH[Tstat.Age>0 & annotation == "Promoter" & Age.sig=="FDR < 0.05",list(na.omit(EntrezID)),],
                     All.neg = dtCH[Tstat.Age<0 & Age.sig=="FDR < 0.05",list(na.omit(EntrezID)),], 
                     GenesPlusPromoters.neg = dtCH[Tstat.Age<0 & annotation != "Other" & Age.sig=="FDR < 0.05",list(na.omit(EntrezID)),],
                     Genes.neg = dtCH[Tstat.Age<0 & distToGene==0 & Age.sig=="FDR < 0.05",list(na.omit(EntrezID)),],
                     Promoters.neg = dtCH[Tstat.Age<0 & annotation == "Promoter" & Age.sig=="FDR < 0.05",list(na.omit(EntrezID)),])
Ageentrez = lapply(Ageentrez, function(x) as.character(unique(x$V1)))              
Ageentrez.dir = lapply(Ageentrez.dir, function(x) as.character(unique(x$V1)))       

GeneUniverse = dtCH[,list(na.omit(EntrezID)),]
GeneUniverse = as.character(unique(GeneUniverse$V1))

# Find enriched Pathways via KEGG
elementNROWS(Ageentrez)
AgekeggList = lapply(Ageentrez, function(x) enrichKEGG(as.character(x), organism="human", universe= GeneUniverse, 
                                                       minGSSize=5, pAdjustMethod="BH", qvalueCutoff=1))
AgekeggListdf = lapply(AgekeggList, function(x) as.data.frame(x))
AgekeggList.dir = lapply(Ageentrez.dir, function(x) enrichKEGG(as.character(x), organism="human", universe= GeneUniverse, 
                                                               minGSSize=5, pAdjustMethod="BH", qvalueCutoff=1))
AgekeggList.dir.df = lapply(AgekeggList.dir, function(x) as.data.frame(x))

# Enriched Molecular Function GOs
AgegoList_MF = lapply(Ageentrez, function(x) enrichGO(as.character(x), ont = "MF", OrgDb = org.Hs.eg.db, 
                                                      universe= GeneUniverse, minGSSize=5, pAdjustMethod="BH",
                                                      qvalueCutoff=1))
AgegoListdf_MF = lapply(AgegoList_MF, function(x) as.data.frame(x))
AgegoList_MF.dir = lapply(Ageentrez.dir, function(x) enrichGO(as.character(x), ont = "MF", OrgDb = org.Hs.eg.db, 
                                                              universe= GeneUniverse, minGSSize=5, pAdjustMethod="BH",
                                                              qvalueCutoff=1))
AgegoListdf_MF.dir = lapply(AgegoList_MF.dir, function(x) as.data.frame(x))

# Biological Process GO enrichment
AgegoList_BP = lapply(Ageentrez, function(x) enrichGO(as.character(x), ont = "BP", OrgDb = org.Hs.eg.db, 
                                                      universe= GeneUniverse, minGSSize=5, pAdjustMethod="BH",
                                                      qvalueCutoff=1))
AgegoListdf_BP = lapply(AgegoList_BP, function(x) as.data.frame(x))
AgegoList_BP.dir = lapply(Ageentrez.dir, function(x) enrichGO(as.character(x), ont = "BP", OrgDb = org.Hs.eg.db, 
                                                              universe= GeneUniverse, minGSSize=5, pAdjustMethod="BH",
                                                              qvalueCutoff=1))
AgegoListdf_BP.dir = lapply(AgegoList_BP.dir, function(x) as.data.frame(x))

# Cellular Compartment GO enrichment
AgegoList_CC = lapply(Ageentrez, function(x) enrichGO(as.character(x), ont = "CC", OrgDb = org.Hs.eg.db, 
                                                      universe= GeneUniverse, minGSSize=5, pAdjustMethod="BH",
                                                      qvalueCutoff=1))
AgegoListdf_CC = lapply(AgegoList_CC, function(x) as.data.frame(x))
AgegoList_CC.dir = lapply(Ageentrez.dir, function(x) enrichGO(as.character(x), ont = "CC", OrgDb = org.Hs.eg.db, 
                                                              universe= GeneUniverse, minGSSize=5, pAdjustMethod="BH",
                                                              qvalueCutoff=1))
AgegoListdf_CC.dir = lapply(AgegoList_CC.dir, function(x) as.data.frame(x))

# Disease Ontology
AgegoList_DO = lapply(Ageentrez, function(x) enrichDO(as.character(x), ont = "DO", universe= GeneUniverse, 
                                                      minGSSize=5, pAdjustMethod="BH", qvalueCutoff=1, readable=TRUE))
AgegoListdf_DO = lapply(AgegoList_DO, function(x) as.data.frame(x))
AgegoList_DO.dir = lapply(Ageentrez.dir, function(x) enrichDO(as.character(x), ont = "DO", universe= GeneUniverse, 
                                                              minGSSize=5, pAdjustMethod="BH", qvalueCutoff=1, readable=TRUE))
AgegoListdf_DO.dir = lapply(AgegoList_DO.dir, function(x) as.data.frame(x))

# Compare the enriched terms between 7 groups
# KEGG
AgecompareKegg = compareCluster(Ageentrez, fun="enrichKEGG", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
AgecompareKegg.dir = compareCluster(Ageentrez.dir, fun="enrichKEGG", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
# Biological Process
AgecompareBP = compareCluster(Ageentrez, fun="enrichGO", ont = "BP", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
AgecompareBP.dir = compareCluster(Ageentrez.dir, fun="enrichGO",  ont = "BP", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
# Molecular Function
AgecompareMF = compareCluster(Ageentrez, fun="enrichGO",  ont = "MF", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
AgecompareMF.dir = compareCluster(Ageentrez.dir, fun="enrichGO",  ont = "MF", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
# Cellular Component
AgecompareCC = compareCluster(Ageentrez, fun="enrichGO",  ont = "CC", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
AgecompareCC.dir = compareCluster(Ageentrez.dir, fun="enrichGO",  ont = "CC", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
# Disease Ontology
AgecompareDO = compareCluster(Ageentrez, fun="enrichDO",  ont = "DO", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
AgecompareDO.dir = compareCluster(Ageentrez.dir, fun="enrichDO",  ont = "DO", qvalueCutoff = 0.05, pvalueCutoff = 0.05)

# save, write to csv
save(AgecompareKegg, AgecompareKegg.dir, AgecompareBP, AgecompareBP.dir, AgecompareMF, AgecompareMF.dir, AgecompareCC, AgecompareCC.dir, AgecompareDO, AgecompareDO.dir,
     AgekeggList, AgekeggList.dir, AgegoList_BP, AgegoList_BP.dir, AgegoList_MF, AgegoList_MF.dir, AgegoList_CC, AgegoList_CC.dir, AgegoList_DO, AgegoList_DO.dir,
     file="/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/non-CpG/nonCG_KEGG_GO_DO_objects_byAge.rda")

# plot compared results
pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/non-CpG/figures/nonCpG_KEGG_GO_DO_plots_byAge.pdf", height = 20, width = 20)
plot(AgecompareKegg, colorBy="p.adjust", showCategory = 45, title= "non-CpG KEGG Pathway Enrichment\nby Age")
plot(AgecompareBP, colorBy="p.adjust", showCategory = 45, title= "non-CpG Biological Process GO Enrichment\nby Age")
plot(AgecompareMF, colorBy="p.adjust", showCategory = 45, title= "non-CpG Molecular Function GO Enrichment\nby Age")
plot(AgecompareCC, colorBy="p.adjust", showCategory = 45, title= "non-CpG Cellular Compartment GO Enrichment\nby Age")
plot(AgecompareDO, colorBy="p.adjust", showCategory = 30, title= "non-CpG Disease Ontology Enrichment\nby Age")
plot(AgecompareKegg.dir, colorBy="p.adjust", showCategory = 45, title= "non-CpG KEGG Pathway Enrichment\nby Age")
plot(AgecompareBP.dir, colorBy="p.adjust", showCategory = 45, title= "non-CpG Biological Process GO Enrichment\nby Age")
plot(AgecompareMF.dir, colorBy="p.adjust", showCategory = 45, title= "non-CpG Molecular Function GO Enrichment\nby Age")
plot(AgecompareCC.dir, colorBy="p.adjust", showCategory = 45, title= "non-CpG Cellular Compartment GO Enrichment\nby Age")
plot(AgecompareDO.dir, colorBy="p.adjust", showCategory = 30, title= "non-CpG Disease Ontology Enrichment\nby Age")
dev.off()


### Gene Ontology: Interaction

Intentrez = list(All = dtCH[Int.sig=="FDR < 0.05",list(na.omit(EntrezID)),], 
                 GenesPlusPromoters = dtCH[annotation != "Other" & Int.sig=="FDR < 0.05",list(na.omit(EntrezID)),],
                 Genes = dtCH[distToGene==0 & Int.sig=="FDR < 0.05",list(na.omit(EntrezID)),],
                 Promoters = dtCH[annotation == "Promoter" & Int.sig=="FDR < 0.05",list(na.omit(EntrezID)),])
Intentrez.dir = list(All.pos = dtCH[Tstat.Interaction>0 & Int.sig=="FDR < 0.05",list(na.omit(EntrezID)),], 
                     GenesPlusPromoters.pos = dtCH[Tstat.Interaction>0 & annotation != "Other" & Int.sig=="FDR < 0.05",list(na.omit(EntrezID)),],
                     Genes.pos = dtCH[Tstat.Interaction>0 & distToGene==0 & Int.sig=="FDR < 0.05",list(na.omit(EntrezID)),],
                     Promoters.pos = dtCH[Tstat.Interaction>0 & annotation == "Promoter" & Int.sig=="FDR < 0.05",list(na.omit(EntrezID)),],
                     All.neg = dtCH[Tstat.Interaction<0 & Int.sig=="FDR < 0.05",list(na.omit(EntrezID)),], 
                     GenesPlusPromoters.neg = dtCH[Tstat.Interaction<0 & annotation != "Other" & Int.sig=="FDR < 0.05",list(na.omit(EntrezID)),],
                     Genes.neg = dtCH[Tstat.Interaction<0 & distToGene==0 & Int.sig=="FDR < 0.05",list(na.omit(EntrezID)),],
                     Promoters.neg = dtCH[Tstat.Interaction<0 & annotation == "Promoter" & Int.sig=="FDR < 0.05",list(na.omit(EntrezID)),])
Intentrez = lapply(Intentrez, function(x) as.character(unique(x$V1)))              
Intentrez.dir = lapply(Intentrez.dir, function(x) as.character(unique(x$V1)))       

GeneUniverse = dtCH[,list(na.omit(EntrezID)),]
GeneUniverse = as.character(unique(GeneUniverse$V1))

# Find enriched Pathways via KEGG
elementNROWS(Intentrez)
IntkeggList = lapply(Intentrez, function(x) enrichKEGG(as.character(x), organism="human", universe= GeneUniverse, 
                                                       minGSSize=5, pAdjustMethod="BH", qvalueCutoff=1))
IntkeggListdf = lapply(IntkeggList, function(x) as.data.frame(x))
IntkeggList.dir = lapply(Intentrez.dir, function(x) enrichKEGG(as.character(x), organism="human", universe= GeneUniverse, 
                                                               minGSSize=5, pAdjustMethod="BH", qvalueCutoff=1))
IntkeggList.dir.df = lapply(IntkeggList.dir, function(x) as.data.frame(x))

# Enriched Molecular Function GOs
IntgoList_MF = lapply(Intentrez, function(x) enrichGO(as.character(x), ont = "MF", OrgDb = org.Hs.eg.db, 
                                                      universe= GeneUniverse, minGSSize=5, pAdjustMethod="BH",
                                                      qvalueCutoff=1))
IntgoListdf_MF = lapply(IntgoList_MF, function(x) as.data.frame(x))
IntgoList_MF.dir = lapply(Intentrez.dir, function(x) enrichGO(as.character(x), ont = "MF", OrgDb = org.Hs.eg.db, 
                                                              universe= GeneUniverse, minGSSize=5, pAdjustMethod="BH",
                                                              qvalueCutoff=1))
IntgoListdf_MF.dir = lapply(IntgoList_MF.dir, function(x) as.data.frame(x))

# Biological Process GO enrichment
IntgoList_BP = lapply(Intentrez, function(x) enrichGO(as.character(x), ont = "BP", OrgDb = org.Hs.eg.db, 
                                                      universe= GeneUniverse, minGSSize=5, pAdjustMethod="BH",
                                                      qvalueCutoff=1))
IntgoListdf_BP = lapply(IntgoList_BP, function(x) as.data.frame(x))
IntgoList_BP.dir = lapply(Intentrez.dir, function(x) enrichGO(as.character(x), ont = "BP", OrgDb = org.Hs.eg.db, 
                                                              universe= GeneUniverse, minGSSize=5, pAdjustMethod="BH",
                                                              qvalueCutoff=1))
IntgoListdf_BP.dir = lapply(IntgoList_BP.dir, function(x) as.data.frame(x))

# Cellular Compartment GO enrichment
IntgoList_CC = lapply(Intentrez, function(x) enrichGO(as.character(x), ont = "CC", OrgDb = org.Hs.eg.db, 
                                                      universe= GeneUniverse, minGSSize=5, pAdjustMethod="BH",
                                                      qvalueCutoff=1))
IntgoListdf_CC = lapply(IntgoList_CC, function(x) as.data.frame(x))
IntgoList_CC.dir = lapply(Intentrez.dir, function(x) enrichGO(as.character(x), ont = "CC", OrgDb = org.Hs.eg.db, 
                                                              universe= GeneUniverse, minGSSize=5, pAdjustMethod="BH",
                                                              qvalueCutoff=1))
IntgoListdf_CC.dir = lapply(IntgoList_CC.dir, function(x) as.data.frame(x))

# Disease Ontology
IntgoList_DO = lapply(Intentrez, function(x) enrichDO(as.character(x), ont = "DO", universe= GeneUniverse, 
                                                      minGSSize=5, pAdjustMethod="BH", qvalueCutoff=1, readable=TRUE))
IntgoListdf_DO = lapply(IntgoList_DO, function(x) as.data.frame(x))
IntgoList_DO.dir = lapply(Intentrez.dir, function(x) enrichDO(as.character(x), ont = "DO", universe= GeneUniverse, 
                                                              minGSSize=5, pAdjustMethod="BH", qvalueCutoff=1, readable=TRUE))
IntgoListdf_DO.dir = lapply(IntgoList_DO.dir, function(x) as.data.frame(x))

# Compare the enriched terms between 7 groups
# KEGG
IntcompareKegg = compareCluster(Intentrez, fun="enrichKEGG", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
IntcompareKegg.dir = compareCluster(Intentrez.dir, fun="enrichKEGG", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
# Biological Process
IntcompareBP = compareCluster(Intentrez, fun="enrichGO", ont = "BP", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
IntcompareBP.dir = compareCluster(Intentrez.dir, fun="enrichGO",  ont = "BP", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
# Molecular Function
IntcompareMF = compareCluster(Intentrez, fun="enrichGO",  ont = "MF", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
IntcompareMF.dir = compareCluster(Intentrez.dir, fun="enrichGO",  ont = "MF", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
# Cellular Component
IntcompareCC = compareCluster(Intentrez, fun="enrichGO",  ont = "CC", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
IntcompareCC.dir = compareCluster(Intentrez.dir, fun="enrichGO",  ont = "CC", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
# Disease Ontology
IntcompareDO = compareCluster(Intentrez, fun="enrichDO",  ont = "DO", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
IntcompareDO.dir = compareCluster(Intentrez.dir, fun="enrichDO",  ont = "DO", qvalueCutoff = 0.05, pvalueCutoff = 0.05)

# save, write to csv
save(IntcompareKegg.dir, IntcompareBP.dir, IntcompareDO.dir,
     IntkeggList, IntkeggList.dir, IntgoList_BP, IntgoList_BP.dir, IntgoList_MF, IntgoList_MF.dir, IntgoList_CC, IntgoList_CC.dir, IntgoList_DO, IntgoList_DO.dir,
     file="/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/non-CpG/nonCG_KEGG_GO_DO_objects_byInteraction.rda")


# plot compared results
pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/non-CpG/figures/nonCpG_KEGG_GO_DO_plots_byInteraction.pdf", height = 20, width = 20)
plot(IntcompareKegg.dir, colorBy="p.adjust", showCategory = 45, title= "non-CpG KEGG Pathway Enrichment\nby Interaction")
plot(IntcompareBP.dir, colorBy="p.adjust", showCategory = 45, title= "non-CpG Biological Process GO Enrichment\nby Interaction")
plot(IntcompareDO.dir, colorBy="p.adjust", showCategory = 30, title= "non-CpG Disease Ontology Enrichment\nby Interaction")
dev.off()

