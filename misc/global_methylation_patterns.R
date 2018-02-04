library("bsseq")

load('/dcl01/lieber/ajaffe/lab/brain-epigenomics/bumphunting/BSobj_bsseqSmooth_Neuron_minCov_3.Rdata')


# postnatal pd table

pd = pData(BSobj)
pd$Race[pd$Race== "CAUC "] <- 'CAUC'
pd$Sex[pd$Sex == " M"] <- 'M'
pd$RIN <- as.numeric(gsub(" ", "", pd$RIN))
pd$pg.DNA.nuclei.input <- as.numeric(pd$pg.DNA.nuclei.input)
pd$Reads <- as.numeric(pd$Reads)
pd$Percent.GreaterThan.Q30 <- as.numeric(pd$Percent.GreaterThan.Q30)

# load methylation matrix for filtered CpGs

Gmeth = getMeth(BSobj, type = 'raw')
dim(Gmeth) # 18664892 CpGs measured
neurGmeth = Gmeth[,which(colnames(Gmeth) %in% pd$Data.ID[which(pd$Cell.Type=="Neuron")])]
gliaGmeth = Gmeth[,which(colnames(Gmeth) %in% pd$Data.ID[which(pd$Cell.Type=="Glia")])]
neurGmeth$Mean = rowMeans(neurGmeth)
gliaGmeth$Mean = rowMeans(gliaGmeth)

## Get stats: postnatal CpGs

round(table(neurGmeth$Mean<0.2)/length(neurGmeth$Mean)*100,2)
#FALSE  TRUE 
#89.84 10.16 
round(table(neurGmeth$Mean>=0.2 & neurGmeth$Mean<=0.8)/length(neurGmeth$Mean)*100,2)
#FALSE  TRUE 
#86.15 13.85 
round(table(neurGmeth$Mean>0.8)/length(neurGmeth$Mean)*100,2)
#FALSE  TRUE 
#24.01 75.99 
round(table(gliaGmeth$Mean<0.2)/length(neurGmeth$Mean)*100,2)
#FALSE  TRUE 
#88.93 11.07 
round(table(gliaGmeth$Mean>=0.2 & gliaGmeth$Mean<=0.8)/length(neurGmeth$Mean)*100,2)
#FALSE  TRUE 
#82.24 17.76 
round(table(gliaGmeth$Mean>0.8)/length(neurGmeth$Mean)*100,2)
#FALSE  TRUE 
#28.83 71.17 


# load methylation matrix for non-CpGs

load('/dcl01/lieber/ajaffe/lab/brain-epigenomics/bsseq/bsobj_by_chr/allChrs_postNatal_cleaned_nonCG_noHomogenate_highCov.Rdata')

Hmeth <- getMeth(BSobj, type = 'raw')
dim(Hmeth) # 58109566 CpHs measured
neurHmeth = Hmeth[,which(colnames(Hmeth) %in% pd$Data.ID[which(pd$Cell.Type=="Neuron")])]
gliaHmeth = Hmeth[,which(colnames(Hmeth) %in% pd$Data.ID[which(pd$Cell.Type=="Glia")])]
neurHmean = rowMeans(neurHmeth)
gliaHmean = rowMeans(gliaHmeth)

## Get stats: postnatal CpHs

round(table(neurHmean<0.2)/length(neurHmean)*100,1)
#FALSE  TRUE 
#8.4  91.6 
round(table(neurHmean>=0.2 & neurHmean<=0.8)/length(neurHmean)*100,1)
#FALSE  TRUE 
#92     8 
round(table(neurHmean>0.8)/length(neurHmean)*100,1)
#FALSE  TRUE 
#99.6   0.4 
round(table(gliaHmean<0.2)/length(gliaHmean)*100,1)
#FALSE  TRUE 
#0.6  99.4 
round(table(gliaHmean>=0.2 & gliaHmean<=0.8)/length(gliaHmean)*100,1)
#FALSE  TRUE 
#99.5   0.5 
round(table(gliaHmean>0.8)/length(gliaHmean)*100,1)
#FALSE  TRUE 
#99.9   0.1


# load methylation matrix for prenatal CpGs

load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/bumphunting/BSobj_bsseqSmooth_Neuron_minCov_3_prenatal.Rdata")

fGmeth <- getMeth(BSobj, type = 'raw')
dim(fGmeth) # 58109566 CpHs measured
mean1 = rowMeans(fGmeth[1:9332446,], na.rm = T)
mean2 = rowMeans(fGmeth[9332447:nrow(fGmeth),], na.rm = T)
Mean = c(mean1,mean2)
length(na.omit(Mean)) # 18664892

## Get stats: prenatal CpGs

round(table(Mean<0.2)/length(Mean)*100,1)
#FALSE  TRUE 
#87.9  12.1 
round(table(Mean>=0.2 & Mean<=0.8)/length(Mean)*100,1)
#FALSE  TRUE 
#82.7  17.3 
round(table(Mean>0.8)/length(Mean)*100,1)
#FALSE  TRUE 
#29.4  70.6


## Get numbers of significant bases

load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/single_CpGs/singleCpG_object.rda")

nrow(CpG[which(CpG$padj.CellType<=0.05),]) # 4824804
round(nrow(CpG[which(CpG$padj.CellType<=0.05),])/nrow(CpG)*100,1) # 25.8%
nrow(CpG[which(CpG$padj.Age<=0.05),]) # 536164
round(nrow(CpG[which(CpG$padj.Age<=0.05),])/nrow(CpG)*100,1) # 2.9%
nrow(CpG[which(CpG$padj.Interaction<=0.05),]) # 90227
round(nrow(CpG[which(CpG$padj.Interaction<=0.05),])/nrow(CpG)*100,1) # 0.5%

load(file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/non-CpG/CH_object.rda")

nrow(CH[which(CH$padj.CellType<=0.05),]) # 7682075
round(nrow(CH[which(CH$padj.CellType<=0.05),])/nrow(CH)*100,1) # 18.8%
nrow(CH[which(CH$padj.Age<=0.05),]) # 536164
round(nrow(CH[which(CH$padj.Age<=0.05),])/nrow(CH)*100,1) # 7.8%
nrow(CH[which(CH$padj.Interaction<=0.05),]) # 76
round(nrow(CH[which(CH$padj.Interaction<=0.05),])/nrow(CH)*100,1) # 0%

load(file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/non-CpG/CHneurons_object.rda")

nrow(CHneurons[which(CHneurons$padj<=0.05),]) # 4020371
round(nrow(CHneurons[which(CHneurons$padj<=0.05),])/nrow(CH)*100,1) # 9.8%


## What is the global proportion of methylated C's, and how does that change by cell type and age?

# At CpGs
Gmeth = getMeth(BSobj, type = 'raw')
neurGmeth = Gmeth[,which(colnames(Gmeth) %in% pd$Data.ID[which(pd$Cell.Type=="Neuron")])]
gliaGmeth = Gmeth[,which(colnames(Gmeth) %in% pd$Data.ID[which(pd$Cell.Type=="Glia")])]

methCG = c()
for(i in 1:ncol(Gmeth)) {
  methCG[i] = nrow(Gmeth[which(Gmeth[,i]>0),])
}
methCG = data.frame(prop = round(methCG/nrow(Gmeth)*100,1), 
                    CellType = ifelse(colnames(Gmeth) %in% pd[which(pd$Cell.Type=="Neuron"),"Data.ID"],"Neuron","Glia"),
                    Age.Bin = pd[match(colnames(Gmeth), pd$Data.ID),"Age.Bin"],
                    Age = pd[match(colnames(Gmeth), pd$Data.ID),"Age"])

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/singleCpG/figures/global_methylation.pdf")
ggplot(df, aes(x=mean)) + geom_density(aes(group=CellType, colour=CellType), size=2) + scale_color_brewer(palette = "Dark2") +
  ylab("Mean Methylation") + 
  xlab("") +
  ggtitle("Mean nonCpG Methylation By Cell Type") + 
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
dev.off()






