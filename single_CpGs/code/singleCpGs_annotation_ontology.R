library(GenomicRanges)
library(GenomicFeatures)
library(data.table)
library(ggplot2)
library(clusterProfiler)
require(org.Hs.eg.db)

load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/single_CpGs/CpG_object.rda")


### Both Neurons and Glia ###

## How many sites are DM?

nrow(CpG) # 18664892
nrow(CpG[which(CpG$padj.CellType<=0.05),]) # 4824804 
nrow(CpG[which(CpG$padj.Age<=0.05),]) # 536164
nrow(CpG[which(CpG$padj.Interaction<=0.05),]) # 90227

dtCpG = data.table(CpG)

### Explore annotation of regions

## how many fall within CpG islands?

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/single_CpGs/figures/single_CpGs_overlap_with_CpG Islands.pdf")
x = dtCpG[,length(unique(regionID)), by = "islands"]
x$perc = round(x$V1/sum(x$V1)*100,2)
ggplot(x, aes(x = islands, y = V1)) + geom_bar(stat = "identity") +
  geom_text(aes( label = paste0(perc,"%")), vjust = -.5) +
  labs(fill="") +
  ylab("Count") + 
  xlab("") +
  ggtitle("Single CpGs Overlapping CpG Islands:\nAll CpGs") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
x = dtCpG[CT.sig=="FDR < 0.05",length(unique(regionID)), by = "islands"]
x$perc = round(x$V1/sum(x$V1)*100,2)
ggplot(x, aes(x = islands, y = V1)) + geom_bar(stat = "identity") +
  geom_text(aes( label = paste0(perc,"%")), vjust = -.5) +
  labs(fill="") +
  ylab("Count") + 
  xlab("") +
  ggtitle("Single CpGs Overlapping CpG Islands:\nFDR < 0.05 by Cell Type") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
x = dtCpG[Age.sig=="FDR < 0.05",length(unique(regionID)), by = "islands"]
x$perc = round(x$V1/sum(x$V1)*100,2)
ggplot(x, aes(x = islands, y = V1)) + geom_bar(stat = "identity") +
  geom_text(aes( label = paste0(perc,"%")), vjust = -.5) +
  labs(fill="") +
  ylab("Count") + 
  xlab("") +
  ggtitle("Single CpGs Overlapping CpG Islands:\nFDR < 0.05 by Age") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
x = dtCpG[Int.sig=="FDR < 0.05",length(unique(regionID)), by = "islands"]
x$perc = round(x$V1/sum(x$V1)*100,2)
ggplot(x, aes(x = islands, y = V1)) + geom_bar(stat = "identity") +
  geom_text(aes( label = paste0(perc,"%")), vjust = -.5) +
  labs(fill="") +
  ylab("Count") + 
  xlab("") +
  ggtitle("Single CpGs Overlapping CpG Islands:\nFDR < 0.05 by Cell Type and Age") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
x = dtCpG[,length(unique(regionID)), by = c("islands", "CT.sig")]
x$perc = NA
x$perc[grep("<", x$CT.sig)] = unlist(c(round(x[grep("<", x$CT.sig),"V1"]/sum(x[grep("<", x$CT.sig),"V1"])*100,2)))
x$perc[grep("", x$CT.sig)] = unlist(c(round(x[grep("", x$CT.sig),"V1"]/sum(x[grep("", x$CT.sig),"V1"])*100,2)))
ggplot(x, aes(x = islands, y = V1)) + geom_bar(stat = "identity") +
   geom_text(aes( label = paste0(perc,"%")), vjust = -.5) +
   facet_grid(. ~ CT.sig) +
   labs(fill="") +
   ylab("Count") + 
   xlab("") +
   ggtitle("Single CpGs Overlapping CpG Islands:\nCell Type") +
   theme(axis.text.x=element_text(angle=45,hjust=1)) + 
   theme(title = element_text(size = 20)) + 
   theme(text = element_text(size = 20))
x = dtCpG[,length(unique(regionID)), by = c("islands", "Age.sig")]
x$perc = NA
x$perc[grep("<", x$Age.sig)] = unlist(c(round(x[grep("<", x$Age.sig),"V1"]/sum(x[grep("<", x$Age.sig),"V1"])*100,2)))
x$perc[grep("", x$Age.sig)] = unlist(c(round(x[grep("", x$Age.sig),"V1"]/sum(x[grep("", x$Age.sig),"V1"])*100,2)))
ggplot(x, aes(x = islands, y = V1)) + geom_bar(stat = "identity") + 
   geom_text(aes( label = paste0(perc,"%")), vjust = -.5) +
   facet_grid(. ~ Age.sig) +
   labs(fill="") +
   ylab("Count") +
   xlab("") +
   ggtitle("Single CpGs Overlapping CpG Islands:\nAge") +
   theme(axis.text.x=element_text(angle=45,hjust=1)) +
   theme(title = element_text(size = 20)) +
   theme(text = element_text(size = 20))
x = dtCpG[,length(unique(regionID)), by = c("islands", "Int.sig")]
x$perc = NA
x$perc[grep("<", x$Int.sig)] = unlist(c(round(x[grep("<", x$Int.sig),"V1"]/sum(x[grep("<", x$Int.sig),"V1"])*100,2)))
x$perc[grep("", x$Int.sig)] = unlist(c(round(x[grep("", x$Int.sig),"V1"]/sum(x[grep("", x$Int.sig),"V1"])*100,2)))
ggplot(x, aes(x = islands, y = V1)) + geom_bar(stat = "identity") +
   geom_text(aes( label = paste0(perc,"%")), vjust = -.5) +
   facet_grid(. ~ Int.sig) +
   labs(fill="") +
   ylab("Count") +
   xlab("") +
   ggtitle("Single CpGs Overlapping CpG Islands:\nInteraction") +
   theme(axis.text.x=element_text(angle=45,hjust=1)) +
   theme(title = element_text(size = 20)) +
   theme(text = element_text(size = 20)) 
 dev.off()


# Is there a relationship between being significantly DM and overlapping a CpG island?

islands = list(CT = data.frame(YesIsland = c(nrow(CpG[which(CpG$CT.sig=="FDR < 0.05" & CpG$islands=="CpG Island"),]),
                         					 nrow(CpG[which(CpG$CT.sig=="FDR > 0.05" & CpG$islands=="CpG Island"),])),
                       		   NoIsland = c(nrow(CpG[which(CpG$CT.sig=="FDR < 0.05" & CpG$islands=="non-Island"),]),
                         					nrow(CpG[which(CpG$CT.sig=="FDR > 0.05" & CpG$islands=="non-Island"),])), row.names=c("Sig","NS")),
               Age = data.frame(YesIsland = c(nrow(CpG[which(CpG$Age.sig=="FDR < 0.05" & CpG$islands=="CpG Island"),]),
                         					  nrow(CpG[which(CpG$Age.sig=="FDR > 0.05" & CpG$islands=="CpG Island"),])),
                       			NoIsland = c(nrow(CpG[which(CpG$Age.sig=="FDR < 0.05" & CpG$islands=="non-Island"),]),
                         					 nrow(CpG[which(CpG$Age.sig=="FDR > 0.05" & CpG$islands=="non-Island"),])), row.names=c("Sig","NS")),
               Int = data.frame(YesIsland = c(nrow(CpG[which(CpG$Int.sig=="FDR < 0.05" & CpG$islands=="CpG Island"),]),
                         					  nrow(CpG[which(CpG$Int.sig=="FDR > 0.05" & CpG$islands=="CpG Island"),])),
                       			NoIsland = c(nrow(CpG[which(CpG$Int.sig=="FDR < 0.05" & CpG$islands=="non-Island"),]),
                         					 nrow(CpG[which(CpG$Int.sig=="FDR > 0.05" & CpG$islands=="non-Island"),])), row.names=c("Sig","NS")))

fisher = lapply(islands, fisher.test)
df = t(data.frame(pval = unlist(lapply(fisher, function(x) x$p.value)), OR = unlist(lapply(fisher, function(x) x$estimate))))
#            CT       Age       Int
#pval 0.0000000 0.0000000 0.0000000
#OR   0.2754005 0.6021545 0.2396739


# assign genomic features

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/single_CpGs/figures/singleCpG_annotation.pdf")
x = dtCpG[,length(unique(regionID)), by = "annotation"]
x$perc = round(x$V1/sum(x$V1)*100,2)
ggplot(x, aes(x = annotation, y = V1)) + geom_bar(stat = "identity") +
  geom_text(aes( label = paste0(perc,"%")), vjust = -.5) +
  labs(fill="") +
  ylab("Count") + 
  xlab("") +
  ggtitle("Single CpG Annotation: All CpGs") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
x = dtCpG[CT.sig=="FDR < 0.05",length(unique(regionID)), by = "annotation"]
x$perc = round(x$V1/sum(x$V1)*100,2)
ggplot(x, aes(x = annotation, y = V1)) + geom_bar(stat = "identity") +
  geom_text(aes( label = paste0(perc,"%")), vjust = -.5) +
  labs(fill="") +
  ylab("Count") + 
  xlab("") +
  ggtitle("Single CpG Annotation:\nFDR < 0.05 by Cell Type") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
x = dtCpG[Age.sig=="FDR < 0.05",length(unique(regionID)), by = "annotation"]
x$perc = round(x$V1/sum(x$V1)*100,2)
ggplot(x, aes(x = annotation, y = V1)) + geom_bar(stat = "identity") +
  geom_text(aes( label = paste0(perc,"%")), vjust = -.5) +
  labs(fill="") +
  ylab("Count") + 
  xlab("") +
  ggtitle("Single CpG Annotation:\nFDR < 0.05 by Age") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
x = dtCpG[Int.sig=="FDR < 0.05",length(unique(regionID)), by = "annotation"]
x$perc = round(x$V1/sum(x$V1)*100,2)
ggplot(x, aes(x = annotation, y = V1)) + geom_bar(stat = "identity") +
  geom_text(aes( label = paste0(perc,"%")), vjust = -.5) +
  labs(fill="") +
  ylab("Count") + 
  xlab("") +
  ggtitle("Single CpG Annotation:\nFDR < 0.05 by Cell Type and Age") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
x = dtCpG[,length(unique(regionID)), by = c("annotation", "CT.sig")]
x$perc = NA
x$perc[grep("<", x$CT.sig)] = unlist(c(round(x[grep("<", x$CT.sig),"V1"]/sum(x[grep("<", x$CT.sig),"V1"])*100,2)))
x$perc[grep(">", x$CT.sig)] = unlist(c(round(x[grep(">", x$CT.sig),"V1"]/sum(x[grep(">", x$CT.sig),"V1"])*100,2)))
ggplot(x, aes(x = annotation, y = V1)) + geom_bar(stat = "identity") +
  geom_text(aes( label = paste0(perc,"%")), vjust = -.5) +
  facet_grid(. ~ CT.sig) +
  labs(fill="") +
  ylab("Count") + 
  theme(axis.text.x=element_text(angle=45,hjust=1)) +
  xlab("") +
  ggtitle("Single CpG Annotation:\nCell Type") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
x = dtCpG[,length(unique(regionID)), by = c("annotation", "Age.sig")]
x$perc = NA
x$perc[grep("<", x$Age.sig)] = unlist(c(round(x[grep("<", x$Age.sig),"V1"]/sum(x[grep("<", x$Age.sig),"V1"])*100,2)))
x$perc[grep(">", x$Age.sig)] = unlist(c(round(x[grep(">", x$Age.sig),"V1"]/sum(x[grep(">", x$Age.sig),"V1"])*100,2)))
ggplot(x, aes(x = annotation, y = V1)) + geom_bar(stat = "identity") +
  geom_text(aes( label = paste0(perc,"%")), vjust = -.5) +
  facet_grid(. ~ Age.sig) +
  labs(fill="") +
  ylab("Count") + 
  theme(axis.text.x=element_text(angle=45,hjust=1)) +
  xlab("") +
  ggtitle("Single CpG Annotation:\nAge") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
x = dtCpG[,length(unique(regionID)), by = c("annotation", "Int.sig")]
x$perc = NA
x$perc[grep("<", x$Int.sig)] = unlist(c(round(x[grep("<", x$Int.sig),"V1"]/sum(x[grep("<", x$Int.sig),"V1"])*100,2)))
x$perc[grep(">", x$Int.sig)] = unlist(c(round(x[grep(">", x$Int.sig),"V1"]/sum(x[grep(">", x$Int.sig),"V1"])*100,2)))
ggplot(x, aes(x = annotation, y = V1)) + geom_bar(stat = "identity") +
  geom_text(aes( label = paste0(perc,"%")), vjust = -.5) +
  facet_grid(. ~ Int.sig) +
  labs(fill="") +
  ylab("Count") + 
  theme(axis.text.x=element_text(angle=45,hjust=1)) +
  xlab("") +
  ggtitle("Single CpG Annotation:\nInteraction") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
dev.off()

# Is there a relationship between being significantly DM and overlapping a gene?

gene = list(CT = data.frame(In = c(nrow(CpG[which(CpG$CT.sig=="FDR < 0.05" & CpG$distToGene==0),]),
                         		   nrow(CpG[which(CpG$CT.sig=="FDR > 0.05" & CpG$distToGene==0),])),
                       		Out = c(nrow(CpG[which(CpG$CT.sig=="FDR < 0.05" & CpG$distToGene>0),]),
                         			nrow(CpG[which(CpG$CT.sig=="FDR > 0.05" & CpG$distToGene>0),])), row.names=c("Sig","NS")),
			Age = data.frame(In = c(nrow(CpG[which(CpG$Age.sig=="FDR < 0.05" & CpG$distToGene==0),]),
                        			nrow(CpG[which(CpG$Age.sig=="FDR > 0.05" & CpG$distToGene==0),])),
            			     Out = c(nrow(CpG[which(CpG$Age.sig=="FDR < 0.05" & CpG$distToGene>0),]),
                         nrow(CpG[which(CpG$Age.sig=="FDR > 0.05" & CpG$distToGene>0),])), row.names=c("Sig","NS")),
			Int = data.frame(In = c(nrow(CpG[which(CpG$Int.sig=="FDR < 0.05" & CpG$distToGene==0),]),
            			            nrow(CpG[which(CpG$Int.sig=="FDR > 0.05" & CpG$distToGene==0),])),
                 			 Out = c(nrow(CpG[which(CpG$Int.sig=="FDR < 0.05" & CpG$distToGene>0),]),
                         			 nrow(CpG[which(CpG$Int.sig=="FDR > 0.05" & CpG$distToGene>0),])), row.names=c("Sig","NS")))
fisher = lapply(gene, fisher.test)
df = t(data.frame(pval = unlist(lapply(fisher, function(x) x$p.value)), OR = unlist(lapply(fisher, function(x) x$estimate))))
#            CT       Age       Int
#pval 0.0000000 0.0000000 0.0008033991
#OR   0.9054789 0.8959855 0.9771336452


# Is there a relationship between being significantly DM and overlapping a promoter?

promoter = list(CT = data.frame(In = c(nrow(CpG[which(CpG$CT.sig=="FDR < 0.05" & CpG$annotation=="Promoter"),]),
                         		   nrow(CpG[which(CpG$CT.sig=="FDR > 0.05" & CpG$annotation=="Promoter"),])),
                       		Out = c(nrow(CpG[which(CpG$CT.sig=="FDR < 0.05" & CpG$annotation!="Promoter"),]),
                         			nrow(CpG[which(CpG$CT.sig=="FDR > 0.05" & CpG$annotation!="Promoter"),])), row.names=c("Sig","NS")),
			Age = data.frame(In = c(nrow(CpG[which(CpG$Age.sig=="FDR < 0.05" & CpG$annotation=="Promoter"),]),
                        			nrow(CpG[which(CpG$Age.sig=="FDR > 0.05" & CpG$annotation=="Promoter"),])),
            			     Out = c(nrow(CpG[which(CpG$Age.sig=="FDR < 0.05" & CpG$annotation!="Promoter"),]),
                         			 nrow(CpG[which(CpG$Age.sig=="FDR > 0.05" & CpG$annotation!="Promoter"),])), row.names=c("Sig","NS")),
			Int = data.frame(In = c(nrow(CpG[which(CpG$Int.sig=="FDR < 0.05" & CpG$annotation=="Promoter"),]),
            			            nrow(CpG[which(CpG$Int.sig=="FDR > 0.05" & CpG$annotation=="Promoter"),])),
                 			 Out = c(nrow(CpG[which(CpG$Int.sig=="FDR < 0.05" & CpG$annotation!="Promoter"),]),
                         			 nrow(CpG[which(CpG$Int.sig=="FDR > 0.05" & CpG$annotation!="Promoter"),])), row.names=c("Sig","NS")))
fisher = lapply(promoter, fisher.test)
df = t(data.frame(pval = unlist(lapply(fisher, function(x) x$p.value)), OR = unlist(lapply(fisher, function(x) x$estimate))))
#            CT          Age           Int
#pval 0.0000000 0.0001444202 6.136361e-109
#OR   0.7013446 0.9743144996  6.635192e-01


# Is there a relationship between being significantly DM and overlapping a gene and/or promoter?

Intergenic = list(CT = data.frame(In = c(nrow(CpG[which(CpG$CT.sig=="FDR < 0.05" & CpG$annotation=="Intergenic"),]),
                         		   	     nrow(CpG[which(CpG$CT.sig=="FDR > 0.05" & CpG$annotation=="Intergenic"),])),
                       			  Out = c(nrow(CpG[which(CpG$CT.sig=="FDR < 0.05" & CpG$annotation!="Intergenic"),]),
                         			      nrow(CpG[which(CpG$CT.sig=="FDR > 0.05" & CpG$annotation!="Intergenic"),])), row.names=c("Sig","NS")),
				  Age = data.frame(In = c(nrow(CpG[which(CpG$Age.sig=="FDR < 0.05" & CpG$annotation=="Intergenic"),]),
                        			      nrow(CpG[which(CpG$Age.sig=="FDR > 0.05" & CpG$annotation=="Intergenic"),])),
            			     	   Out = c(nrow(CpG[which(CpG$Age.sig=="FDR < 0.05" & CpG$annotation!="Intergenic"),]),
                         				   nrow(CpG[which(CpG$Age.sig=="FDR > 0.05" & CpG$annotation!="Intergenic"),])), row.names=c("Sig","NS")),
				  Int = data.frame(In = c(nrow(CpG[which(CpG$Int.sig=="FDR < 0.05" & CpG$annotation=="Intergenic"),]),
            			            	  nrow(CpG[which(CpG$Int.sig=="FDR > 0.05" & CpG$annotation=="Intergenic"),])),
                 			 Out = c(nrow(CpG[which(CpG$Int.sig=="FDR < 0.05" & CpG$annotation!="Intergenic"),]),
                         			 nrow(CpG[which(CpG$Int.sig=="FDR > 0.05" & CpG$annotation!="Intergenic"),])), row.names=c("Sig","NS")))
fisher = lapply(Intergenic, fisher.test)
df = t(data.frame(pval = unlist(lapply(fisher, function(x) x$p.value)), OR = unlist(lapply(fisher, function(x) x$estimate))))
#           CT      Age          Int
#pval 0.000000 0.000000 6.977052e-32
#OR   1.168861 1.127211 1.085610e+00



### Gene Ontology

GeneUniverse = as.character(na.omit(unique(CpG$EntrezID)))

##Cell Type

CTentrez = list(All = CpG[CpG$CT.sig=="FDR < 0.05","EntrezID"], 
              	GenesPlusPromoters = CpG[CpG$annotation != "Intergenic" & CpG$CT.sig=="FDR < 0.05","EntrezID"],
              	Genes = CpG[CpG$distToGene==0 & CpG$CT.sig=="FDR < 0.05","EntrezID"],
              	Promoters = CpG[CpG$annotation == "Promoter" & CpG$CT.sig=="FDR < 0.05","EntrezID"])
CTentrez.dir = list(All.pos = CpG[CpG$Tstat.CellType>0 & CpG$CT.sig=="FDR < 0.05","EntrezID"], 
                  GenesPlusPromoters.pos = CpG[CpG$Tstat.CellType>0 & CpG$annotation != "Intergenic" & CpG$CT.sig=="FDR < 0.05","EntrezID"],
                  Genes.pos = CpG[CpG$Tstat.CellType>0 & CpG$distToGene==0 & CpG$CT.sig=="FDR < 0.05","EntrezID"],
                  Promoters.pos = CpG[CpG$Tstat.CellType>0 & CpG$annotation == "Promoter" & CpG$CT.sig=="FDR < 0.05","EntrezID"],
                  All.neg = CpG[CpG$Tstat.CellType<0 & CpG$CT.sig=="FDR < 0.05","EntrezID"], 
                  GenesPlusPromoters.neg = CpG[CpG$Tstat.CellType<0 & CpG$annotation != "Intergenic" & CpG$CT.sig=="FDR < 0.05","EntrezID"],
                  Genes.neg = CpG[CpG$Tstat.CellType<0 & CpG$distToGene==0 & CpG$CT.sig=="FDR < 0.05","EntrezID"],
                  Promoters.neg = CpG[CpG$Tstat.CellType<0 & CpG$annotation == "Promoter" & CpG$CT.sig=="FDR < 0.05","EntrezID"])
CTentrez = lapply(CTentrez, function(x) as.character(na.omit(unique(x))))              
CTentrez.dir = lapply(CTentrez.dir, function(x) as.character(na.omit(unique(x))))       


# Find enriched Pathways via KEGG
elementNROWS(CTentrez)
CTkeggList = lapply(CTentrez, function(x) enrichKEGG(as.character(x), organism="human", universe= GeneUniverse, 
                                                   minGSSize=5, pAdjustMethod="BH", qvalueCutoff=1))
CTkeggList.dir = lapply(CTentrez.dir, function(x) enrichKEGG(as.character(x), organism="human", universe= GeneUniverse, 
                                                   minGSSize=5, pAdjustMethod="BH", qvalueCutoff=1))
# enriched Molecular Function GOs
CTgoList_MF = lapply(CTentrez, function(x) enrichGO(as.character(x), ont = "MF", OrgDb = org.Hs.eg.db, 
                                                  universe= GeneUniverse, minGSSize=5, pAdjustMethod="BH",qvalueCutoff=1))
CTgoList_MF.dir = lapply(CTentrez.dir, function(x) enrichGO(as.character(x), ont = "MF", OrgDb = org.Hs.eg.db, 
                                                universe= GeneUniverse, minGSSize=5, pAdjustMethod="BH",qvalueCutoff=1))
# Biological Process GO enrichment
CTgoList_BP = lapply(CTentrez, function(x) enrichGO(as.character(x), ont = "BP", OrgDb = org.Hs.eg.db, 
                                                universe= GeneUniverse, minGSSize=5, pAdjustMethod="BH",qvalueCutoff=1))
CTgoList_BP.dir = lapply(CTentrez.dir, function(x) enrichGO(as.character(x), ont = "BP", OrgDb = org.Hs.eg.db, 
                                                        universe= GeneUniverse, minGSSize=5, pAdjustMethod="BH",qvalueCutoff=1))
# Cellular Compartment GO enrichment
CTgoList_CC = lapply(CTentrez, function(x) enrichGO(as.character(x), ont = "CC", OrgDb = org.Hs.eg.db, 
                                                universe= GeneUniverse, minGSSize=5, pAdjustMethod="BH",qvalueCutoff=1))
CTgoList_CC.dir = lapply(CTentrez.dir, function(x) enrichGO(as.character(x), ont = "CC", OrgDb = org.Hs.eg.db, 
                                                        universe= GeneUniverse, minGSSize=5, pAdjustMethod="BH",qvalueCutoff=1))
# Disease Ontology
CTgoList_DO = lapply(CTentrez, function(x) enrichDO(as.character(x), ont = "DO", universe= GeneUniverse, 
                                                  minGSSize=5, pAdjustMethod="BH", qvalueCutoff=1, readable=TRUE))
CTgoList_DO.dir = lapply(CTentrez.dir, function(x) enrichDO(as.character(x), ont = "DO", universe= GeneUniverse, 
                                                minGSSize=5, pAdjustMethod="BH", qvalueCutoff=1, readable=TRUE))


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

# save
save(CTcompareKegg, CTcompareKegg.dir, CTcompareBP, CTcompareBP.dir, CTcompareMF, CTcompareMF.dir, CTcompareCC, CTcompareCC.dir, CTcompareDO, CTcompareDO.dir,
     CTkeggList, CTkeggList.dir, CTgoList_BP, CTgoList_BP.dir, CTgoList_MF, CTgoList_MF.dir, CTgoList_CC, CTgoList_CC.dir, CTgoList_DO, CTgoList_DO.dir,
     file="/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/single_CpGs/singleCpG_KEGG_GO_DO_objects_byCellType.rda")

# plot compared results
pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/single_CpGs/figures/single_CpGs_KEGG_GO_DO_plots_byCellType.pdf", height = 20, width = 20)
plot(CTcompareKegg, colorBy="p.adjust", showCategory = 45, title= "Single CpG KEGG Pathway enrichment\nby Cell Type")
plot(CTcompareBP, colorBy="p.adjust", showCategory = 45, title= "Single CpG Biological Process GO enrichment\nby Cell Type")
plot(CTcompareMF, colorBy="p.adjust", showCategory = 45, title= "Single CpG Molecular Function GO enrichment\nby Cell Type")
plot(CTcompareCC, colorBy="p.adjust", showCategory = 45, title= "Single CpG Cellular Compartment GO enrichment\nby Cell Type")
plot(CTcompareDO, colorBy="p.adjust", showCategory = 30, title= "Single CpG Disease Ontology enrichment\nby Cell Type")
plot(CTcompareKegg.dir, colorBy="p.adjust", showCategory = 45, title= "Single CpG KEGG Pathway enrichment\nby Cell Type")
plot(CTcompareBP.dir, colorBy="p.adjust", showCategory = 45, title= "Single CpG Biological Process GO enrichment\nby Cell Type")
plot(CTcompareMF.dir, colorBy="p.adjust", showCategory = 45, title= "Single CpG Molecular Function GO enrichment\nby Cell Type")
plot(CTcompareCC.dir, colorBy="p.adjust", showCategory = 45, title= "Single CpG Cellular Compartment GO enrichment\nby Cell Type")
plot(CTcompareDO.dir, colorBy="p.adjust", showCategory = 30, title= "Single CpG Disease Ontology enrichment\nby Cell Type")
dev.off()


### Gene Ontology: Age
Ageentrez = list(All = CpG[CpG$Age.sig=="FDR < 0.05","EntrezID"], 
              	GenesPlusPromoters = CpG[CpG$annotation != "Intergenic" & CpG$Age.sig=="FDR < 0.05","EntrezID"],
              	Genes = CpG[CpG$distToGene==0 & CpG$Age.sig=="FDR < 0.05","EntrezID"],
              	Promoters = CpG[CpG$annotation == "Promoter" & CpG$Age.sig=="FDR < 0.05","EntrezID"])
Ageentrez.dir = list(All.pos = CpG[CpG$Tstat.Age>0 & CpG$Age.sig=="FDR < 0.05","EntrezID"], 
                  GenesPlusPromoters.pos = CpG[CpG$Tstat.Age>0 & CpG$annotation != "Intergenic" & CpG$Age.sig=="FDR < 0.05","EntrezID"],
                  Genes.pos = CpG[CpG$Tstat.Age>0 & CpG$distToGene==0 & CpG$Age.sig=="FDR < 0.05","EntrezID"],
                  Promoters.pos = CpG[CpG$Tstat.Age>0 & CpG$annotation == "Promoter" & CpG$Age.sig=="FDR < 0.05","EntrezID"],
                  All.neg = CpG[CpG$Tstat.Age<0 & CpG$Age.sig=="FDR < 0.05","EntrezID"], 
                  GenesPlusPromoters.neg = CpG[CpG$Tstat.Age<0 & CpG$annotation != "Intergenic" & CpG$Age.sig=="FDR < 0.05","EntrezID"],
                  Genes.neg = CpG[CpG$Tstat.Age<0 & CpG$distToGene==0 & CpG$Age.sig=="FDR < 0.05","EntrezID"],
                  Promoters.neg = CpG[CpG$Tstat.Age<0 & CpG$annotation == "Promoter" & CpG$Age.sig=="FDR < 0.05","EntrezID"])
Ageentrez = lapply(Ageentrez, function(x) as.character(na.omit(unique(x))))              
Ageentrez.dir = lapply(Ageentrez.dir, function(x) as.character(na.omit(unique(x))))       


# Find enriched Pathways via KEGG
elementNROWS(Ageentrez)
AgekeggList = lapply(Ageentrez, function(x) enrichKEGG(as.character(x), organism="human", universe= GeneUniverse, 
                                                       minGSSize=5, pAdjustMethod="BH", qvalueCutoff=1))
AgekeggList.dir = lapply(Ageentrez.dir, function(x) enrichKEGG(as.character(x), organism="human", universe= GeneUniverse, 
                                                               minGSSize=5, pAdjustMethod="BH", qvalueCutoff=1))
# enriched Molecular Function GOs
AgegoList_MF = lapply(Ageentrez, function(x) enrichGO(as.character(x), ont = "MF", OrgDb = org.Hs.eg.db, 
                                                      universe= GeneUniverse, minGSSize=5, pAdjustMethod="BH", qvalueCutoff=1))
AgegoList_MF.dir = lapply(Ageentrez.dir, function(x) enrichGO(as.character(x), ont = "MF", OrgDb = org.Hs.eg.db, 
                                                              universe= GeneUniverse, minGSSize=5, pAdjustMethod="BH", qvalueCutoff=1))
# Biological Process GO enrichment
AgegoList_BP = lapply(Ageentrez, function(x) enrichGO(as.character(x), ont = "BP", OrgDb = org.Hs.eg.db, 
                                                      universe= GeneUniverse, minGSSize=5, pAdjustMethod="BH", qvalueCutoff=1))
AgegoList_BP.dir = lapply(Ageentrez.dir, function(x) enrichGO(as.character(x), ont = "BP", OrgDb = org.Hs.eg.db, 
                                                              universe= GeneUniverse, minGSSize=5, pAdjustMethod="BH", qvalueCutoff=1))
# Cellular Compartment GO enrichment
AgegoList_CC = lapply(Ageentrez, function(x) enrichGO(as.character(x), ont = "CC", OrgDb = org.Hs.eg.db, 
                                                      universe= GeneUniverse, minGSSize=5, pAdjustMethod="BH", qvalueCutoff=1))
AgegoList_CC.dir = lapply(Ageentrez.dir, function(x) enrichGO(as.character(x), ont = "CC", OrgDb = org.Hs.eg.db, 
                                                              universe= GeneUniverse, minGSSize=5, pAdjustMethod="BH",
                                                              qvalueCutoff=1))
# Disease Ontology
AgegoList_DO = lapply(Ageentrez, function(x) enrichDO(as.character(x), ont = "DO", universe= GeneUniverse, 
                                                      minGSSize=5, pAdjustMethod="BH", qvalueCutoff=1, readable=TRUE))
AgegoList_DO.dir = lapply(Ageentrez.dir, function(x) enrichDO(as.character(x), ont = "DO", universe= GeneUniverse, 
                                                              minGSSize=5, pAdjustMethod="BH", qvalueCutoff=1, readable=TRUE))
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
     file="/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/single_CpGs/singleCpG_KEGG_GO_DO_objects_byAge.rda")

# plot compared results
pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/single_CpGs/figures/single_CpGs_KEGG_GO_DO_plots_byAge.pdf", height = 20, width = 20)
plot(AgecompareKegg, colorBy="p.adjust", showCategory = 45, title= "Single CpG KEGG Pathway enrichment\nby Age")
plot(AgecompareBP, colorBy="p.adjust", showCategory = 45, title= "Single CpG Biological Process GO enrichment\nby Age")
plot(AgecompareMF, colorBy="p.adjust", showCategory = 45, title= "Single CpG Molecular Function GO enrichment\nby Age")
plot(AgecompareCC, colorBy="p.adjust", showCategory = 45, title= "Single CpG Cellular Compartment GO enrichment\nby Age")
plot(AgecompareDO, colorBy="p.adjust", showCategory = 30, title= "Single CpG Disease Ontology enrichment\nby Age")
plot(AgecompareKegg.dir, colorBy="p.adjust", showCategory = 45, title= "Single CpG KEGG Pathway enrichment\nby Age")
plot(AgecompareBP.dir, colorBy="p.adjust", showCategory = 45, title= "Single CpG Biological Process GO enrichment\nby Age")
plot(AgecompareMF.dir, colorBy="p.adjust", showCategory = 45, title= "Single CpG Molecular Function GO enrichment\nby Age")
plot(AgecompareCC.dir, colorBy="p.adjust", showCategory = 45, title= "Single CpG Cellular Compartment GO enrichment\nby Age")
plot(AgecompareDO.dir, colorBy="p.adjust", showCategory = 30, title= "Single CpG Disease Ontology enrichment\nby Age")
dev.off()


### Gene Ontology: Interaction

Intentrez = list(All = CpG[CpG$Int.sig=="FDR < 0.05","EntrezID"], 
              	GenesPlusPromoters = CpG[CpG$annotation != "Intergenic" & CpG$Int.sig=="FDR < 0.05","EntrezID"],
              	Genes = CpG[CpG$distToGene==0 & CpG$Int.sig=="FDR < 0.05","EntrezID"],
              	Promoters = CpG[CpG$annotation == "Promoter" & CpG$Int.sig=="FDR < 0.05","EntrezID"])
Intentrez.dir = list(All.pos = CpG[CpG$Tstat.Interaction>0 & CpG$Int.sig=="FDR < 0.05","EntrezID"], 
                  GenesPlusPromoters.pos = CpG[CpG$Tstat.Interaction>0 & CpG$annotation != "Intergenic" & CpG$Int.sig=="FDR < 0.05","EntrezID"],
                  Genes.pos = CpG[CpG$Tstat.Interaction>0 & CpG$distToGene==0 & CpG$Int.sig=="FDR < 0.05","EntrezID"],
                  Promoters.pos = CpG[CpG$Tstat.Interaction>0 & CpG$annotation == "Promoter" & CpG$Int.sig=="FDR < 0.05","EntrezID"],
                  All.neg = CpG[CpG$Tstat.Interaction<0 & CpG$Int.sig=="FDR < 0.05","EntrezID"], 
                  GenesPlusPromoters.neg = CpG[CpG$Tstat.Interaction<0 & CpG$annotation != "Intergenic" & CpG$Int.sig=="FDR < 0.05","EntrezID"],
                  Genes.neg = CpG[CpG$Tstat.Interaction<0 & CpG$distToGene==0 & CpG$Int.sig=="FDR < 0.05","EntrezID"],
                  Promoters.neg = CpG[CpG$Tstat.Interaction<0 & CpG$annotation == "Promoter" & CpG$.Interactionsig=="FDR < 0.05","EntrezID"])
Intentrez = lapply(Intentrez, function(x) as.character(na.omit(unique(x))))              
Intentrez.dir = lapply(Intentrez.dir, function(x) as.character(na.omit(unique(x))))       


# Find enriched Pathways via KEGG
elementNROWS(Intentrez)
IntkeggList = lapply(Intentrez, function(x) enrichKEGG(as.character(x), organism="human", universe= GeneUniverse, 
                                                       minGSSize=5, pAdjustMethod="BH", qvalueCutoff=1))
IntkeggList.dir = lapply(Intentrez.dir, function(x) enrichKEGG(as.character(x), organism="human", universe= GeneUniverse, 
                                                               minGSSize=5, pAdjustMethod="BH", qvalueCutoff=1))
# enriched Molecular Function GOs
IntgoList_MF = lapply(Intentrez, function(x) enrichGO(as.character(x), ont = "MF", OrgDb = org.Hs.eg.db, 
                                                      universe= GeneUniverse, minGSSize=5, pAdjustMethod="BH",
                                                      qvalueCutoff=1))
IntgoList_MF.dir = lapply(Intentrez.dir, function(x) enrichGO(as.character(x), ont = "MF", OrgDb = org.Hs.eg.db, 
                                                              universe= GeneUniverse, minGSSize=5, pAdjustMethod="BH",
                                                              qvalueCutoff=1))
# Biological Process GO enrichment
IntgoList_BP = lapply(Intentrez, function(x) enrichGO(as.character(x), ont = "BP", OrgDb = org.Hs.eg.db, 
                                                      universe= GeneUniverse, minGSSize=5, pAdjustMethod="BH",
                                                      qvalueCutoff=1))
IntgoList_BP.dir = lapply(Intentrez.dir, function(x) enrichGO(as.character(x), ont = "BP", OrgDb = org.Hs.eg.db, 
                                                              universe= GeneUniverse, minGSSize=5, pAdjustMethod="BH",
                                                              qvalueCutoff=1))
# Cellular Compartment GO enrichment
IntgoList_CC = lapply(Intentrez, function(x) enrichGO(as.character(x), ont = "CC", OrgDb = org.Hs.eg.db, 
                                                      universe= GeneUniverse, minGSSize=5, pAdjustMethod="BH",
                                                      qvalueCutoff=1))
IntgoList_CC.dir = lapply(Intentrez.dir, function(x) enrichGO(as.character(x), ont = "CC", OrgDb = org.Hs.eg.db, 
                                                              universe= GeneUniverse, minGSSize=5, pAdjustMethod="BH",
                                                              qvalueCutoff=1))
# Disease Ontology
IntgoList_DO = lapply(Intentrez, function(x) enrichDO(as.character(x), ont = "DO", universe= GeneUniverse, 
                                                      minGSSize=5, pAdjustMethod="BH", qvalueCutoff=1, readable=TRUE))
IntgoList_DO.dir = lapply(Intentrez.dir, function(x) enrichDO(as.character(x), ont = "DO", universe= GeneUniverse, 
                                                              minGSSize=5, pAdjustMethod="BH", qvalueCutoff=1, readable=TRUE))


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
save(IntcompareKegg, IntcompareKegg.dir, IntcompareBP, IntcompareBP.dir, IntcompareMF, IntcompareMF.dir, IntcompareCC, IntcompareCC.dir, IntcompareDO, IntcompareDO.dir,
     IntkeggList, IntkeggList.dir, IntgoList_BP, IntgoList_BP.dir, IntgoList_MF, IntgoList_MF.dir, IntgoList_CC, IntgoList_CC.dir, IntgoList_DO, IntgoList_DO.dir,
     file="/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/single_CpGs/singleCpG_KEGG_GO_DO_objects_byInteraction.rda")


# plot compared results
pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/single_CpGs/figures/single_CpGs_KEGG_GO_DO_plots_byInteraction.pdf", height = 20, width = 20)
plot(IntcompareKegg, colorBy="p.adjust", showCategory = 45, title= "Single CpG KEGG Pathway enrichment\nby Interaction")
plot(IntcompareBP, colorBy="p.adjust", showCategory = 45, title= "Single CpG Biological Process GO enrichment\nby Interaction")
plot(IntcompareMF, colorBy="p.adjust", showCategory = 45, title= "Single CpG Molecular Function GO enrichment\nby Interaction")
plot(IntcompareCC, colorBy="p.adjust", showCategory = 45, title= "Single CpG Cellular Compartment GO enrichment\nby Interaction")
plot(IntcompareDO, colorBy="p.adjust", showCategory = 30, title= "Single CpG Disease Ontology enrichment\nby Interaction")
plot(IntcompareKegg.dir, colorBy="p.adjust", showCategory = 45, title= "Single CpG KEGG Pathway enrichment\nby Interaction")
plot(IntcompareBP.dir, colorBy="p.adjust", showCategory = 45, title= "Single CpG Biological Process GO enrichment\nby Interaction")
plot(IntcompareMF.dir, colorBy="p.adjust", showCategory = 45, title= "Single CpG Molecular Function GO enrichment\nby Interaction")
plot(IntcompareCC.dir, colorBy="p.adjust", showCategory = 45, title= "Single CpG Cellular Compartment GO enrichment\nby Interaction")
plot(IntcompareDO.dir, colorBy="p.adjust", showCategory = 30, title= "Single CpG Disease Ontology enrichment\nby Interaction")
dev.off()

