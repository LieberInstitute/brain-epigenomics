library(GenomicRanges)
library(ggplot2)
library(readxl)
library(jaffelab)
library(minfi)
library(bsseq)

load('/dcl01/lieber/ajaffe/lab/brain-epigenomics/bumphunting/BSobj_bsseqSmooth_Neuron_minCov_3.Rdata')


## find CpGs in our data that are subtype markers
tab = read_excel("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/Supplementary-Table-S6.xlsx", skip=5)
tab = as.data.frame(tab)

### six classes
classIndex = list(gliaLow = which(tab$GLIA.meanCpG < 0.15 & tab$GABA.meanCpG > 0.85 &  tab$GLU.meanCpG > 0.85 ),
			gliaHi = which(tab$GLIA.meanCpG > 0.85 & tab$GABA.meanCpG < 0.2 &  tab$GLU.meanCpG < 0.2 ),
			gluLow = which(tab$GLIA.meanCpG > 0.8 & tab$GABA.meanCpG > 0.8 &  tab$GLU.meanCpG < 0.15 ),
			gluHi = which(tab$GLIA.meanCpG < 0.3 & tab$GABA.meanCpG < 0.3 &  tab$GLU.meanCpG > 0.7 ),
			gabaLow = which(tab$GLIA.meanCpG > 0.8 & tab$GABA.meanCpG < 0.2 &  tab$GLU.meanCpG > 0.8 ),
			gabaHi = which(tab$GLIA.meanCpG < 0.2 & tab$GABA.meanCpG > 0.8 &  tab$GLU.meanCpG < 0.2 ))
tabSub = tab[unique(unlist(classIndex)),]


## add filtering stats
big = makeGRangesFromDataFrame(tabSub, seqnames.field = "CHR", start.field = "POSITION", end.field = "POSITION", keep.extra.columns = T)
oo = findOverlaps(big, BSobj)
BSobj_sub = BSobj[subjectHits(oo),]

## filter by mean coverage
keepIndex = which(rowMeans(getCoverage(BSobj_sub)) > 10)
BSobj_sub = BSobj_sub[keepIndex,]
tabSub = tabSub[keepIndex,]

## do deconvolution
coefs = as.matrix(tabSub[,c("GLIA.meanCpG", "GABA.meanCpG", "GLU.meanCpG")])
colnames(coefs) = ss(colnames(coefs), "\\.")
counts = minfi:::projectCellType(getMeth(BSobj_sub), coefs)
counts = as.data.frame(counts)

boxplot(counts$GLIA ~ colData(BSobj)$Cell.Type)
boxplot(counts$GLU ~ colData(BSobj)$Cell.Type)
boxplot(counts$GABA ~ colData(BSobj)$Cell.Type)

## plots
pd = colData(BSobj)
pd = cbind(pd, as.data.frame(counts[match(rownames(pd), rownames(counts)),]) )    
pd = as.data.frame(pd)

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/single_CpGs/figures/estimated_NeuronalSubtype_proportions_array.pdf", height=6, width=5)
ggplot(pd, aes(x = Age, y = GLU, colour = Cell.Type)) + geom_point() +
  theme_classic() + scale_colour_brewer(8, palette="Dark2") +
  ylab("Cell type proportion") + xlab("Age (years)") +
  ggtitle("Glutamatergic") + ylim(0,1) +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20), legend.title=element_blank()) + theme(legend.position="bottom")
ggplot(pd, aes(x = Age, y = GABA, colour = Cell.Type)) + geom_point() +
  theme_classic() + scale_colour_brewer(8, palette="Dark2") +
  ylab("Cell type proportion") + xlab("Age (years)") +
  ggtitle("GABAergic") + ylim(0,1) +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20), legend.title=element_blank()) + theme(legend.position="bottom")
ggplot(pd, aes(x = Age, y = GLIA, colour = Cell.Type)) + geom_point() +
  theme_classic() + scale_colour_brewer(8, palette="Dark2") +
  ylab("Cell type proportion") + xlab("Age (years)") +
  ggtitle("Glia") + ylim(0,1) +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20), legend.title=element_blank()) + theme(legend.position="bottom")
dev.off()

## regression
summary(lm(GABA ~ Age, data=pd, subset=Cell.Type == "Neuron"))
summary(lm(GLU ~ Age, data=pd, subset=Cell.Type == "Neuron"))

#######################
## using full data ####

load("/dcl02/lieber/ajaffe/psychENCODE_Data/Dracheva/cellType_weights_houseman.rda")

## filter for coverage
oo2 = findOverlaps(gr, BSobj)
BSobj_sub2 = BSobj[subjectHits(oo2),]
coefEsts_sub2 = coefEsts[queryHits(oo2),]

## filter by mean coverage
keepIndex = which(rowMeans(getCoverage(BSobj_sub2)) > 10)
BSobj_sub2 = BSobj_sub2[keepIndex,]
coefEsts_sub2 = coefEsts_sub2[keepIndex,]

## do deconvolution
counts2 = minfi:::projectCellType(getMeth(BSobj_sub2), coefEsts_sub2)
counts2 = as.data.frame(counts2)

boxplot(counts2$GLIA ~ colData(BSobj)$Cell.Type)
boxplot(counts2$GLU ~ colData(BSobj)$Cell.Type)
boxplot(counts2$GABA ~ colData(BSobj)$Cell.Type)

## plots
pd2 = colData(BSobj_sub2)
pd2 = cbind(pd2, as.data.frame(counts2[match(rownames(pd2), rownames(counts2)),]) )    
pd2 = as.data.frame(pd2)

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/single_CpGs/figures/estimated_NeuronalSubtype_proportions_array_fullData.pdf", height=6, width=5)
ggplot(pd2, aes(x = Age, y = GLU, colour = Cell.Type)) + geom_point() +
  theme_classic() + scale_colour_brewer(8, palette="Dark2") +
  ylab("Cell type proportion") + xlab("Age (years)") +
  ggtitle("Glutamatergic") + ylim(0,1) +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20), legend.title=element_blank()) + theme(legend.position="bottom")
ggplot(pd2, aes(x = Age, y = GABA, colour = Cell.Type)) + geom_point() +
  theme_classic() + scale_colour_brewer(8, palette="Dark2") +
  ylab("Cell type proportion") + xlab("Age (years)") +
  ggtitle("GABAergic") + ylim(0,1) +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20), legend.title=element_blank()) + theme(legend.position="bottom")
ggplot(pd2, aes(x = Age, y = GLIA, colour = Cell.Type)) + geom_point() +
  theme_classic() + scale_colour_brewer(8, palette="Dark2") +
  ylab("Cell type proportion") + xlab("Age (years)") +
  ggtitle("Glia") + ylim(0,1) +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20), legend.title=element_blank()) + theme(legend.position="bottom")
dev.off()

## regression
summary(lm(GABA ~ Age, data=pd2, subset=Cell.Type == "Neuron"))
summary(lm(GLU ~ Age, data=pd2, subset=Cell.Type == "Neuron"))


