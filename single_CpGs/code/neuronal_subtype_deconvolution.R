library(bsseq)
library(GenomicRanges)
library(ggplot2)


load('/dcl01/lieber/ajaffe/lab/brain-epigenomics/bumphunting/BSobj_bsseqSmooth_Neuron_minCov_3.Rdata')

#### Estimate proportion of GABA and Glut neurons within our sorted samples

## Get our data

gr = granges(BSobj)
cov = getCoverage(BSobj)
meth = getMeth(BSobj)
pd = pData(BSobj)


## find CpGs in our data that are subtype markers

tab = openxlsx::read.xlsx("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/Supplementary-Table-S8.xlsx", startRow = 5)
head(tab)
tabgr = makeGRangesFromDataFrame(tab, seqnames.field = "CHR", start.field = "POSITION", end.field = "POSITION", keep.extra.columns = T)


# ones with >0.8 change

big = tabgr[abs(tabgr$deltaGABA.GLU)>=0.8]
oo = findOverlaps(big, gr)

biggr = gr[subjectHits(oo)]
bigcov = cov[subjectHits(oo),]

meanCov = rowMeans(bigcov)

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/single_CpGs/figures/coverage_dist_neuronalSubtypeCpGs.pdf")
hist(meanCov)
dev.off()


# ones with mean coverage of >15

biggr.HC = biggr[which(meanCov>15)]
oo = findOverlaps(tabgr, biggr.HC)
big.HC = tabgr[queryHits(oo)]

table(big.HC$DiffMeth=="DM.GLU")
#FALSE  TRUE 
#   53  1529


# Filter to the training set

DM.GLU = big.HC[big.HC$DiffMeth=="DM.GLU"]
DM.GLU = DM.GLU[order(DM.GLU$deltaGABA.GLU, decreasing = T)]
DM.GLU = DM.GLU[1:100]

subtypeCs = c(DM.GLU, big.HC[big.HC$DiffMeth=="DM.GABA"])

oo = findOverlaps(subtypeCs, gr)
subB = meth[subjectHits(oo),]
coefs = data.frame(GABAergic = subtypeCs$GABA.meanCpG, Glutamateric = subtypeCs$GLU.meanCpG)


# Estimate proportions

counts = minfi:::projectCellType(subB, as.matrix(coefs))    
pd = cbind(pd, as.data.frame(counts[match(rownames(pd), rownames(counts)),]) )    
pd = as.data.frame(pd)
write.csv(pd, file="/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/estimated_NeuronalSubtype_proportions.csv", quote = F)


### Plot results

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/single_CpGs/figures/estimated_NeuronalSubtype_proportions.pdf", height=6, width=5)

ggplot(pd, aes(x = Age, y = Glutamateric, colour = Cell.Type)) + geom_point() +
  theme_classic() + scale_colour_brewer(8, palette="Dark2") +
  ylab("Cell type proportion") + xlab("Age (years)") +
  ggtitle("Glutamatergic") + ylim(0,1) +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20), legend.title=element_blank()) + theme(legend.position="bottom")
ggplot(pd, aes(x = Age, y = GABAergic, colour = Cell.Type)) + geom_point() +
  theme_classic() + scale_colour_brewer(8, palette="Dark2") +
  ylab("Cell type proportion") + xlab("Age (years)") +
  ggtitle("GABAergic") + ylim(0,1) +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20), legend.title=element_blank()) + theme(legend.position="bottom")

dev.off()








