library(ggplot2)


#### explore expression patterns of the three example cdDMRs in paper

ids = geneMap[geneMap$Symbol %in% c("HDAC4", "CACNA1B", "AKT3"),"gencodeID"]

## Get sorted RNAseq data

load("/dcl01/lieber/ajaffe/CellSorting/RNAseq_pipeline/rawCounts_CellSorting_July5_n12.rda", verbose=T)

# format phenotype table
metrics[grep("12", metrics$SampleID),"Age"] = "Neonate"
metrics[grep("13", metrics$SampleID),"Age"] = "Teenager"
metrics[grep("14", metrics$SampleID),"Age"] = "Toddler"
metrics$CellType = ifelse(metrics$NeuN=="NeuN_Minus", "Glia", "Neuron")
polya = metrics[which(metrics$Prep=="PolyA"),]
ribo = metrics[which(metrics$Prep=="Ribo"),]

# format counts
polyaCounts = geneCounts[,grep("PolyA", colnames(geneCounts))]
riboCounts = geneCounts[,grep("Ribo", colnames(geneCounts))]
combinedCounts = polyaCounts + riboCounts
match(rownames(polya), colnames(combinedCounts))

bgE = matrix(rep(colSums(combinedCounts)), nc = nrow(polya), nr = nrow(combinedCounts),    byrow=TRUE)
widE = matrix(rep(geneMap$Length), nr = nrow(combinedCounts), nc = nrow(polya),    byrow=FALSE)
sortedRpkm = combinedCounts/(widE/1000)/(bgE/1e6)

sortedRpkm = sortedRpkm[rownames(sortedRpkm) %in% ids,]
sortedRpkm = reshape::melt(sortedRpkm)
colnames(sortedRpkm) = c("gencodeID","sampleID", "RPKM")
sortedRpkm[grep("Minus", sortedRpkm$sampleID),"CellType"] = "Glia"
sortedRpkm[grep("Plus", sortedRpkm$sampleID),"CellType"] = "Neuron"
sortedRpkm[grep("12", sortedRpkm$sampleID),"Age"] = 0.36
sortedRpkm[grep("13", sortedRpkm$sampleID),"Age"] = 14.01
sortedRpkm[grep("14", sortedRpkm$sampleID),"Age"] = 2.71
sortedRpkm$Symbol = geneMap[match(sortedRpkm$gencodeID, geneMap$gencodeID),"Symbol"]

## get homogenate RPKM

load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/brainseq_pipeline/polyA_unstranded/rse_gene_polyA_dlpfc_n41.Rdata", verbose = T)

homogenate = getRPKM(rse_gene)
homogenate = homogenate[rownames(homogenate) %in% ids,]
homogenate = reshape::melt(homogenate)
colnames(homogenate) = c("gencodeID","sampleID", "RPKM")
homogenate$CellType = "Homogenate"
homogenate$Age = data.frame(colData(rse_gene))[match(homogenate$sampleID, rownames(data.frame(colData(rse_gene)))),"Age"]
homogenate$Symbol = geneMap[match(homogenate$gencodeID, geneMap$gencodeID),"Symbol"]
homogenate = homogenate[homogenate$Age>0,]


## plot the trajectories

df = rbind(sortedRpkm, homogenate)
df$Sorted = ifelse(df$CellType=="Homogenate", "Homogenate", "Sorted")
df$Sorted = factor(df$Sorted, levels = c("Sorted", "Homogenate"))
df$CellType = factor(df$CellType, levels = c("Glia","Neuron","Homogenate"))


pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/sorted_nuclear_RNA/figures/example_expression_forPaper.pdf", height = 4, width = 5)
ggplot(df[df$Symbol=="AKT3",], aes(x = Age, y = log(RPKM+1), colour = CellType)) + facet_grid(. ~ Sorted) +
  geom_point() + geom_smooth() +
  theme_classic() + scale_colour_brewer(8, palette="Dark2") +
  ylab("log(RPKM+1)") + xlab("Age") + xlim(0,23) + ylim(0,4) + 
  ggtitle("Gene Expression: AKT3") + 
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20), legend.title=element_blank()) + theme(legend.position="bottom")
ggplot(df[df$Symbol=="HDAC4",], aes(x = Age, y = log(RPKM+1), colour = CellType)) + facet_grid(. ~ Sorted) +
  geom_point() + geom_smooth() +
  theme_classic() + scale_colour_brewer(8, palette="Dark2") +
  ylab("log(RPKM+1)") + xlab("Age") + xlim(0,23) + ylim(0,4) + 
  ggtitle("Gene Expression: HDAC4") + 
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20), legend.title=element_blank()) + theme(legend.position="bottom")
ggplot(df[df$Symbol=="CACNA1B",], aes(x = Age, y = log(RPKM+1), colour = CellType)) + facet_grid(. ~ Sorted) +
  geom_point() + geom_smooth() +
  theme_classic() + scale_colour_brewer(8, palette="Dark2") +
  ylab("log(RPKM+1)") + xlab("Age") + xlim(0,23) + ylim(0,4) + 
  ggtitle("Gene Expression: CACNA1B") + 
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20), legend.title=element_blank()) + theme(legend.position="bottom")
dev.off()

