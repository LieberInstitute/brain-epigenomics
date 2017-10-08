library('bsseq')
library('devtools')
library('limma')
library('jaffelab')
library('shinycsv') # must download
library('RColorBrewer')

## Load data
load('/dcl01/lieber/ajaffe/lab/brain-epigenomics/bsseq/bsobj_by_chr/allChrs_postNatal_cleaned_nonCG_noHomogenate_highCov.Rdata')

## Size of the data
dim(BSobj) # 58109566       32

## extract pheno
pd <- pData(BSobj)

## Fix pheno data
pd$Race[pd$Race == "CAUC "] <- 'CAUC'
pd$Sex[pd$Sex == " M"] <- 'M'
pd$RIN <- as.numeric(gsub(" ", "", pd$RIN))
pd$pg.DNA.nuclei.input <- as.numeric(pd$pg.DNA.nuclei.input)
pd$Reads <- as.numeric(pd$Reads)
pd$Percent.GreaterThan.Q30 <- as.numeric(pd$Percent.GreaterThan.Q30)

## Get methylation info
gr <- granges(BSobj)
meth <- getMeth(BSobj, type = 'raw')
meth = meth[,which(colnames(meth) %in% rownames(pd))]
meth.g0 <- meth > 0
meth.filt <- rowSums(meth.g0) >= 5
meth.tab <- table(meth.filt)
meth.tab
round(meth.tab / sum(meth.tab) * 100, 2)
meth = meth[meth.filt, ]
gr = gr[meth.filt]
rm(meth.g0, meth.filt, BSobj)

# what proportion have non-zero DNAm per cell type?

neurons = meth[,which(colnames(meth) %in% rownames(pd[which(pd$Cell.Type=="Neuron"),]))]
glia = meth[,which(colnames(meth) %in% rownames(pd[which(pd$Cell.Type=="Glia"),]))]
neurons.g0 <- neurons > 0
glia.g0 <- glia > 0
neurons.filt <- rowSums(neurons.g0) >0
glia.filt <- rowSums(glia.g0) >0
neurons.tab <- table(neurons.filt)
neurons.tab
#FALSE     TRUE 
#471 40818271 
round(neurons.tab / sum(neurons.tab) * 100, 2)
#FALSE  TRUE 
#0   100 

glia.tab <- table(glia.filt)
glia.tab
#FALSE     TRUE
#8364639 32454103
round(glia.tab / sum(glia.tab) * 100, 2)
#FALSE  TRUE
#20.49 79.51


# compare mean DNAm within each cell type

meanNeuron = rowMeans(neurons)
meanGlia = rowMeans(glia)
sdNeuron = apply(neurons, 1, sd)
sdGlia = apply(glia, 1, sd)

df = data.frame(mean = c(meanNeuron, meanGlia), sd = c(sdNeuron, sdGlia),
                CellType = c(rep.int("Neuron", length(meanNeuron)), rep.int("Glia", length(meanGlia))))

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/non-CpG/figures/meanMethylation_byCellType.pdf")
ggplot(df, aes(x=CellType, y = mean)) + geom_boxplot() + 
  ylab("Mean Methylation") + 
  xlab("") +
  ggtitle("Mean nonCpG Methylation By Cell Type") + 
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
dev.off()

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/non-CpG/figures/sdMethylation_byCellType.pdf")
ggplot(df, aes(x=CellType, y = sd)) + geom_boxplot() + 
  ylab("SD") + 
  xlab("") +
  ggtitle("nonCpG Methylation Standard Deviation By Cell Type") + 
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
dev.off()

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/non-CpG/figures/meanMethylation_byCellType_density.pdf")
ggplot(df, aes(x=mean)) + geom_density(aes(group=CellType, colour=CellType)) + 
  ylab("Mean Methylation") + 
  xlab("") +
  ggtitle("Mean nonCpG Methylation By Cell Type") + 
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
dev.off()

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/non-CpG/figures/sdMethylation_byCellType_density.pdf")
ggplot(df, aes(x=sd)) + geom_density(aes(group=CellType, colour=CellType)) + 
  ylab("SD") + 
  xlab("") +
  ggtitle("nonCpG Methylation Standard Deviation By Cell Type") + 
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
dev.off()


t.test(x=meanNeuron, y=meanGlia)
#data:  meanNeuron and meanGlia
#t = 2851.7, df = 50404000, p-value < 2.2e-16
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  0.06564673 0.06573703
#sample estimates:
#  mean of x  mean of y 
#0.09234615 0.02665427
