library(GenomicRanges)
library(ggplot2)
library(readxl)
library(jaffelab)
library(minfi)
library(bsseq)

## read in wgbs data
load('/dcl01/lieber/ajaffe/lab/brain-epigenomics/bumphunting/BSobj_bsseqSmooth_Neuron_minCov_3.Rdata')

## coverage filter
keepIndex = which(rowMeans(getCoverage(BSobj)) > 10)
BSobj2 = BSobj[keepIndex,]

## read in array data, RGset, to do SQN together
load("/dcl02/lieber/ajaffe/psychENCODE_Data/Dracheva/Dracheva_GluGabaGlia_RGset.rda")
load("/dcl02/lieber/ajaffe/psychENCODE_Data/Dracheva/Dracheva_GluGabaGlia_Mset.rda") # for pheno
Mset = mapToGenome(RGset, mergeManifest=TRUE) # get raw M and U

## find overlapping CpGs
oo = findOverlaps(Mset, BSobj2)
BSobj2 = BSobj2[subjectHits(oo),]
Mset2 = Mset[queryHits(oo),]

## jointly normalize
M_seq = getCoverage(BSobj2, type = "M")
U_seq = getCoverage(BSobj2, type = "Cov") - M_seq

M_array = minfi::getMeth(Mset2)
U_array = minfi::getUnmeth(Mset2)

M_both = cbind(M_array, M_seq)
U_both = cbind(U_array, U_seq)

Mset_joint = GenomicMethylSet(gr = rowRanges(Mset2), 
	Meth = M_both,Unmeth = U_both,
	annotation = annotation(Mset2))
sex = stringr::str_trim(c(Mset_SQN$predictedSex, BSobj$Sex))

Mset_joint_SQN =  preprocessQuantile(Mset_joint, 
	fixOutliers = FALSE, #Quantile normalization
	removeBadSamples = FALSE, 
	quantileNormalize = TRUE, 
	stratified = FALSE, 
	mergeManifest = TRUE, 
	sex = sex)
	
Mset_array_joint = Mset_joint_SQN[,1:15]
Mset_seq_joint = Mset_joint_SQN[,-(1:15)]

## get cell type probes
colData(Mset_array_joint) = colData(Mset_SQN)
compData <- minfi:::pickCompProbes(mSet = Mset_array_joint, 
	cellTypes = unique(Mset_array_joint$CellType),
    compositeCellType = "Brain", probeSelect = "any")
coefEsts = compData$coefEsts
save(coefEsts, file = "joint_normalized_coefEsts.rda")

## get counts
meth_seq_joint = minfi::getBeta(Mset_seq_joint[rownames(coefEsts),])
counts = minfi:::projectCellType(meth_seq_joint, coefEsts,lessThanOne = TRUE)
counts = as.data.frame(counts)


## plots
pd = colData(BSobj)
pd = cbind(pd, as.data.frame(counts[match(rownames(pd), rownames(counts)),]) )    
pd = as.data.frame(pd)

## 
library(lattice)
theSeq = seq(0,1,by=0.001)                         
my.col <- colorRampPalette(c("blue", "white","red"))(length(theSeq))
toplot = meth_seq_joint
colnames(toplot) = pd$Cell.Type
colInd = order(pd$Cell.Type, pd$Age)
rowInd = order(coefEsts[,"GLU"])
pdf("levelplot_decon.pdf", w= 12, h=7)
print(levelplot(as.matrix(cbind(coefEsts, toplot[,colInd]))[rowInd,], 
    aspect = "fill", at = theSeq,pretty=TRUE,
    panel = panel.levelplot.raster, col.regions = my.col,
        ylab = "Sample", xlab = "CpG",
	scales=list(x=list(rot=90, cex=0.7))))
dev.off()
		
		
pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/single_CpGs/figures/estimated_NeuronalSubtype_proportions_array_QN.pdf", height=6, width=5)
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

