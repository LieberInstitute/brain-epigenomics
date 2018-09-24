library("bsseq")
library(ggplot2)
library(data.table)

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

pd = pData(BSobj)
pd$Race[pd$Race== "CAUC "] <- 'CAUC'
pd$Sex[pd$Sex == " M"] <- 'M'
write.csv(pd, quote=F, file="/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/prenatal.WGBS.pheno.info.csv")

fGmeth <- getMeth(BSobj, type = 'raw')
dim(fGmeth) # 
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


# load methylation matrix for prenatal non-CpGs

load('/dcl01/lieber/ajaffe/lab/brain-epigenomics/bsseq/prenatal/allChrs_cleaned_CX_Prenatal_CpH.Rdata')

pd <- read.csv("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/prenatal.WGBS.pheno.info.csv")
BSobj <- updateObject(BSobj, verbose=TRUE)
Pmeth <- getMeth(BSobj, type = 'raw')
rm(BSobj)
dim(Pmeth) # 58109566 CpHs measured
Pmeth = as.matrix(Pmeth)

mean1 = rowMeans(Pmeth[1:9684929,], na.rm = T)
mean2 = rowMeans(Pmeth[9684930:19369856,], na.rm = T)
mean3 = rowMeans(Pmeth[19369857:29054784,], na.rm = T)
mean4 = rowMeans(Pmeth[29054785:38739712,], na.rm = T)
mean5 = rowMeans(Pmeth[38739713:48424640,], na.rm = T)
mean6 = rowMeans(Pmeth[48424641:nrow(Pmeth),], na.rm = T)

Mean = c(mean1,mean2,mean3,mean4,mean5,mean6)
length(na.omit(Mean)) # 58109566


## Get stats: prenatal CpHs

round(table(Mean<0.2)/length(Mean)*100,2)
#FALSE  TRUE 
# 0.25 99.75 
round(table(Mean>=0.2 & Mean<=0.8)/length(Mean)*100,2)
#FALSE  TRUE 
#99.83  0.17 
round(table(Mean>0.8)/length(Mean)*100,2)
#FALSE  TRUE 
#99.93  0.07 


# load methylation matrix for homogenate postnatal CpGs

load("/dcl01/lieber/WGBS/LIBD_Data/bsseqObj/bsseqObj_postNatal_cleaned_CpGonly.rda")

pd = pData(BSobj)
pd$Race[pd$Race== "CAUC "] <- 'CAUC'
pd$Sex[pd$Sex == " M"] <- 'M'
write.csv(pd, quote=F, file="/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/postnatal.homogenate.WGBS.pheno.info.csv")

hGmeth <- getMeth(BSobj, type = 'raw')
dim(hGmeth) #  28217448 CpGs measured
Mean = rowMeans(hGmeth, na.rm = T)
length(na.omit(Mean)) # 27983221

## Get stats: homogenate CpGs

round(table(Mean<0.2)/length(Mean)*100,1)
#FALSE  TRUE 
# 90.7   8.5
round(table(Mean>=0.2 & Mean<=0.8)/length(Mean)*100,1)
#FALSE  TRUE 
# 82.2  16.9
round(table(Mean>0.8)/length(Mean)*100,1)
#FALSE  TRUE 
# 25.4  73.7 


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
methCG$Age.Bin = factor(methCG$Age.Bin, levels = c("Neonate","Toddler","Child","Early.Teen","Teen","Young.Adult"))

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/single_CpGs/figures/global_methylation.pdf")
ggplot(methCG, aes(x = CellType, y = prop)) + geom_boxplot() +
  ylab("mCpG / CpG") + xlab("") + theme_classic() +
  ggtitle("Global Proportion of Methylated CpGs") + 
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
ggplot(methCG, aes(x = Age.Bin, y = prop)) + geom_boxplot() +
  ylab("mCpG / CpG") + xlab("") + theme_classic() +
  ggtitle("Global Proportion of Methylated CpGs") + 
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
ggplot(methCG, aes(x = Age.Bin, y = prop)) + geom_boxplot() +
  ylab("mCpG / CpG") + xlab("") + facet_grid(. ~ CellType) +
  ggtitle("Global Proportion of Methylated CpGs") + theme_classic() + 
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
ggplot(methCG, aes(x = Age, y = prop, colour = CellType)) + geom_point() +
  geom_smooth(method = "lm") +  
  theme_classic() + scale_colour_brewer(8, palette="Dark2") +
  ylab("mCpG / CpG") + xlab("") +
  ggtitle("Global Proportion of Methylated CpGs") + 
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20), legend.title=element_blank()) + theme(legend.position="bottom")
ggplot(methCG, aes(x = Age, y = prop, colour = CellType)) + geom_point() +
  geom_path() +  
  theme_classic() + scale_colour_brewer(8, palette="Dark2") +
  ylab("mCpG / CpG") + xlab("") +
  ggtitle("Global Proportion of Methylated CpGs") + 
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20), legend.title=element_blank()) + theme(legend.position="bottom")
dev.off()


# At nonCpGs

Hmeth = getMeth(BSobj, type = 'raw')
Hgr = granges(BSobj)
Hdf = data.frame(Hgr)
rownames(Hmeth) = paste0(Hdf$seqnames, ":", Hdf$start, "-", Hdf$end)

methCH = c()
for(i in 1:ncol(Hmeth)) {
  methCH[i] = nrow(Hmeth[which(Hmeth[,i]>0),])
}
methCH10 = c()
for(i in 1:ncol(Hmeth)) {
  methCH10[i] = nrow(Hmeth[which(Hmeth[,i]>0.1),])
}

methCH = data.frame(prop = round(methCH/nrow(Hmeth)*100,1), prop10 = round(methCH10/nrow(Hmeth)*100,1),
                    CellType = ifelse(colnames(Hmeth) %in% pd[which(pd$Cell.Type=="Neuron"),"Data.ID"],"Neuron","Glia"),
                    Age.Bin = pd[match(colnames(Hmeth), pd$Data.ID),"Age.Bin"],
                    Age = pd[match(colnames(Hmeth), pd$Data.ID),"Age"])
methCH$Age.Bin = factor(methCH$Age.Bin, levels = c("Neonate","Toddler","Child","Early.Teen","Teen","Young.Adult"))
write.csv(methCH, quote = F, file="/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/non-CpG/proportion_mCH_byCellType_byAge.csv")

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/non-CpG/figures/nonCG_global_methylation.pdf")
ggplot(methCH, aes(x = CellType, y = prop)) + geom_boxplot() +
  theme_classic() +
  ylab("mCpH / CpH") + xlab("") + ylim(0,50) +
  ggtitle("Global Proportion of Methylated CpHs") + 
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20), legend.title=element_blank()) + theme(legend.position="bottom")
ggplot(methCH, aes(x = Age.Bin, y = prop)) + geom_boxplot() +
  theme_classic() +
  ylab("mCpH / CpH") + xlab("") + ylim(0,50) +
  ggtitle("Global Proportion of Methylated CpHs") + 
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20), legend.title=element_blank()) + theme(legend.position="bottom")
ggplot(methCH, aes(x = Age.Bin, y = prop)) + geom_boxplot() +
  theme_classic() +  ylim(0,50) +
  ylab("mCpH / CpH") + xlab("") + facet_grid(. ~ CellType) +
  theme(axis.text.x=element_text(angle = 45, hjust = 1)) +
  ggtitle("Global Proportion of Methylated CpHs") + 
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20), legend.title=element_blank()) + theme(legend.position="bottom")
ggplot(methCH, aes(x = Age, y = prop, colour = CellType)) + geom_point() +
  geom_smooth(method = "lm") +  ylim(0,50) +
  theme_classic() + scale_colour_brewer(8, palette="Dark2") +
  ylab("mCpH / CpH") + xlab("") +
  ggtitle("Global Proportion of Methylated CpHs") + 
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20), legend.title=element_blank()) + theme(legend.position="bottom")
ggplot(methCH, aes(x = Age, y = prop, colour = CellType)) + geom_point() +
  geom_path() +  ylim(0,50) +
  theme_classic() + scale_colour_brewer(8, palette="Dark2") +
  ylab("mCpH / CpH") + xlab("") +
  ggtitle("Global Proportion of Methylated CpHs") + 
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20), legend.title=element_blank()) + theme(legend.position="bottom")

ggplot(methCH, aes(x = CellType, y = prop10)) + geom_boxplot() +
  theme_classic() +  ylim(0,30) +
  ylab("mCpH / CpH") + xlab("") +
  ggtitle("Global Proportion of\nCpHs Methylated >10%") + 
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20), legend.title=element_blank()) + theme(legend.position="bottom")
ggplot(methCH, aes(x = Age.Bin, y = prop10)) + geom_boxplot() +
  theme_classic() +  ylim(0,30) +
  ylab("mCpH / CpH") + xlab("") +
  ggtitle("Global Proportion of\nCpHs Methylated >10%") + 
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20), legend.title=element_blank()) + theme(legend.position="bottom")
ggplot(methCH, aes(x = Age.Bin, y = prop10)) + geom_boxplot() +
  theme_classic() +  ylim(0,30) +
  ylab("mCpH / CpH") + xlab("") + facet_grid(. ~ CellType) +
  theme(axis.text.x=element_text(angle = 45, hjust = 1)) +
  ggtitle("Global Proportion of\nCpHs Methylated >10%") + 
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20), legend.title=element_blank()) + theme(legend.position="bottom")
ggplot(methCH, aes(x = Age, y = prop10, colour = CellType)) + geom_point() +
  geom_smooth(method = "lm") +  ylim(0,30) +
  theme_classic() + scale_colour_brewer(8, palette="Dark2") +
  ylab("mCpH / CpH") + xlab("") +
  ggtitle("Global Proportion of\nCpHs Methylated >10%") + 
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20), legend.title=element_blank()) + theme(legend.position="bottom")
ggplot(methCH, aes(x = Age, y = prop10, colour = CellType)) + geom_point() +
  geom_path() +  ylim(0,30) +
  theme_classic() + scale_colour_brewer(8, palette="Dark2") +
  ylab("mCpH / CpH") + xlab("") +
  ggtitle("Global Proportion of\nCpHs Methylated >10%") + 
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20), legend.title=element_blank()) + theme(legend.position="bottom")
dev.off()


## Separate into trinucleotide contexts

freq = table(as.character(Hdf$trinucleotide_context))
x = cbind(as.data.frame(Hmeth), ID = rownames(Hmeth), trinucleotide_context = Hdf$trinucleotide_context)
pd$Group = paste(pd$Cell.Type, pd$Age.Bin, sep="\n")
colnames(x) = c(paste0(colnames(x)[1:32], ":", pd[match(colnames(x)[1:32], as.character(pd$Data.ID)),"Group"]), "ID", "trinucleotide_context")
Hdt = reshape2::melt(x)
rm(x)
Hdt$Group = gsub('.*\\:', "", Hdt$variable)
Hdt$variable = gsub("\\:.*", "", Hdt$variable)
Hdt = data.table(Hdt)

tri = Hdt[value>0,length(unique(ID)), by= c("variable","trinucleotide_context")]
tri10 = Hdt[value>0.10,length(unique(ID)), by= c("variable","trinucleotide_context")]
nrowHmeth = nrow(Hmeth)
rm(Hmeth)
save(tri,tri10, Hdt, file="/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/non-CpG/number_mCH_bytrinucleotideContext.rda")

methCH = data.frame(rbind(data.frame(tri), data.frame(tri10)), prop = c(round(tri$V1/nrowHmeth,3), round(tri10$V1/nrowHmeth,3)),
                    perc = c(rep.int("nonzero",nrow(tri)), rep.int("ten",nrow(tri10))),
                    triFreq = freq[match(c(as.character(tri$trinucleotide_context),as.character(tri10$trinucleotide_context)), names(freq))],
                    CellType = pd[match(c(as.character(tri$variable),as.character(tri10$variable)), pd$Data.ID),"Cell.Type"],
                    Age.Bin = pd[match(c(as.character(tri$variable),as.character(tri10$variable)), pd$Data.ID),"Age.Bin"],
                    Age = pd[match(c(as.character(tri$variable),as.character(tri10$variable)), pd$Data.ID),"Age"])
methCH$triprop = round(methCH$V1 / methCH$triFreq.Freq,3)
methCH$Age.Bin  = factor(methCH$Age.Bin, levels = c("Neonate","Toddler","Child","Early.Teen","Teen","Young.Adult"))
methCH$trinucleotide_context = factor(methCH$trinucleotide_context, levels = c("CAG","CAC","CAT","CAA","CTG","CTC","CCA","CTT","CTA","CCT","CCC","CCG"))
write.csv(methCH, quote = F, file="/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/non-CpG/proportion_mCH_byCellType_byAge_byTrinucleotideContext.csv")

ageprop = rbind(cbind(Hdt[value>0, length(unique(ID)), by= c("Group","trinucleotide_context")], perc = "nonzero"),
                cbind(Hdt[value>0.10,length(unique(ID)), by= c("Group","trinucleotide_context")], perc = "ten"))
ageprop = cbind(ageprop, triFreq = freq[match(as.character(ageprop$trinucleotide_context), names(freq))])
ageprop$triprop = round(ageprop$V1 / ageprop$triFreq.N,3)
write.csv(ageprop, quote = F, file="/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/non-CpG/proportion_mCH_byCellType_byAgeGroup_byTrinucleotideContext.csv")
ageprop = read.csv(file="/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/non-CpG/proportion_mCH_byCellType_byAgeGroup_byTrinucleotideContext.csv")
ageprop = data.frame(CellType = ageprop[seq.int(from=1,to=576,by=2),"Group"], Age.Bin = ageprop[-seq.int(from=1,to=576,by=2),"X"],
               trinucleotide_context = ageprop[-seq.int(from=1,to=576,by=2),"Group"], V1 = ageprop[-seq.int(from=1,to=576,by=2),"trinucleotide_context"],
               perc = ageprop[-seq.int(from=1,to=576,by=2),"V1"], triFreq = ageprop[-seq.int(from=1,to=576,by=2),"triFreq.V1"])
ageprop$triprop = round(ageprop$V1 / ageprop$triFreq,3)
ageprop$Age.Bin  = factor(ageprop$Age.Bin, levels = c("Neonate","Toddler","Child","Early.Teen","Teen","Young.Adult"))
ageprop$trinucleotide_context = factor(ageprop$trinucleotide_context, levels = c("CAG","CAC","CAT","CAA","CTG","CTC","CCA","CTT","CTA","CCT","CCC","CCG"))
write.csv(ageprop, quote = F, file="/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/non-CpG/proportion_mCH_byCellType_byAgeGroup_byTrinucleotideContext.csv")


pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/non-CpG/figures/nonCG_global_methylation_byTrinucleotideContext.pdf", width =10)
ggplot(ageprop[ageprop$perc=="nonzero",], aes(x = Age.Bin, y = triprop, fill=trinucleotide_context), color=trinucleotide_context) + 
  geom_bar(position = "fill",stat = "identity", width=0.75) + facet_grid(. ~ CellType) +
  ylab("Proportion mCpH / CpH") + 
  theme_classic() + xlab("") +
  ggtitle("Global Methylation of CpHs") + 
  theme(axis.text.x=element_text(angle = 45, hjust = 1)) +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20), legend.title=element_blank()) + theme(legend.position="bottom")
ggplot(ageprop[ageprop$perc=="ten",], aes(x = Age.Bin, y = triprop, fill=trinucleotide_context), color=trinucleotide_context) + 
  geom_bar(position = "fill",stat = "identity", width=0.75) + facet_grid(. ~ CellType) +
  ylab("Proportion mCpH / CpH") + 
  theme_classic() + xlab("") +
  ggtitle("Global Methylation of CpHs (>10%)") +
  theme(axis.text.x=element_text(angle = 45, hjust = 1)) +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20), legend.title=element_blank()) + theme(legend.position="bottom")
ggplot(methCH[methCH$perc=="nonzero",], aes(x = trinucleotide_context, y = triprop, fill=CellType)) + geom_boxplot() +
  theme_classic() + scale_fill_brewer(8, palette="Dark2") +
  ylab("mCpH / CpH") + xlab("") +
  ggtitle("Global Proportion of Methylated CpHs\nContext Proportion") + 
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20), legend.title=element_blank()) + theme(legend.position="bottom")
ggplot(methCH[methCH$perc=="nonzero",], aes(x = trinucleotide_context, y = triprop, fill=Age.Bin)) + geom_boxplot() +
  theme_classic() +
  ylab("mCpH / CpH") + xlab("") + 
  ggtitle("Global Proportion of Methylated CpHs\nContext Proportion") + 
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20), legend.title=element_blank()) + theme(legend.position="bottom")
ggplot(methCH[methCH$perc=="nonzero",], aes(x = trinucleotide_context, y = triprop, fill= Age.Bin)) + geom_boxplot() +
  theme_classic() +  
  ylab("mCpH / CpH") + xlab("") + facet_grid(. ~ CellType) +
  theme(axis.text.x=element_text(angle = 45, hjust = 1)) +
  ggtitle("Global Proportion of Methylated CpHs\nContext Proportion") + 
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20), legend.title=element_blank()) + theme(legend.position="bottom")
ggplot(methCH[which(methCH$perc=="nonzero" & methCH$trinucleotide_context %in% c("CAG","CAC")),], 
       aes(x = Age, y = triprop, colour = CellType)) + geom_point() + facet_grid(. ~ trinucleotide_context) +
  geom_smooth(method = "lm") +  
  theme_classic() + scale_colour_brewer(8, palette="Dark2") +
  ylab("mCpH / CpH") + xlab("") +
  ggtitle("Global Proportion of Methylated CpHs\nContext Proportion") + 
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20), legend.title=element_blank()) + theme(legend.position="bottom")
ggplot(methCH[which(methCH$perc=="nonzero" & methCH$trinucleotide_context %in% c("CAG","CAC")),], 
       aes(x = Age, y = prop, colour = CellType)) + geom_point() + facet_grid(. ~ trinucleotide_context) +
  geom_smooth(method = "lm") +  
  theme_classic() + scale_colour_brewer(8, palette="Dark2") +
  ylab("mCpH / CpH") + xlab("") +
  ggtitle("Global Proportion of Methylated CpHs\nTotal Proportion") + 
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20), legend.title=element_blank()) + theme(legend.position="bottom")

ggplot(methCH[methCH$perc=="ten",], aes(x = trinucleotide_context, y = triprop, fill=CellType)) + geom_boxplot() +
  theme_classic() +   scale_fill_brewer(8, palette="Dark2") +
  ylab("mCpH / CpH") + xlab("") +
  ggtitle("Global Proportion of CpHs Methylated >10%\nContext Proportion") + 
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20), legend.title=element_blank()) + theme(legend.position="bottom")
ggplot(methCH[methCH$perc=="ten",], aes(x = trinucleotide_context, y = triprop, fill= Age.Bin)) + geom_boxplot() +
  theme_classic() +  
  ylab("mCpH / CpH") + xlab("") +
  ggtitle("Global Proportion of CpHs Methylated >10%\nContext Proportion") + 
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20), legend.title=element_blank()) + theme(legend.position="bottom")
ggplot(methCH[methCH$perc=="ten",], aes(x = trinucleotide_context, y = triprop, fill= Age.Bin)) + geom_boxplot() +
  theme_classic() +  
  ylab("mCpH / CpH") + xlab("") + facet_grid(. ~ CellType) +
  theme(axis.text.x=element_text(angle = 45, hjust = 1)) +
  ggtitle("Global Proportion of CpHs Methylated >10%\nContext Proportion") + 
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20), legend.title=element_blank()) + theme(legend.position="bottom")
ggplot(methCH[which(methCH$perc=="ten" & methCH$trinucleotide_context %in% c("CAG","CAC")),], 
       aes(x = Age, y = triprop, colour = CellType)) + geom_point() + facet_grid(. ~ trinucleotide_context) +
  geom_smooth(method = "lm") +  
  theme_classic() + scale_colour_brewer(8, palette="Dark2") +
  ylab("mCpH / CpH") + xlab("") +
  ggtitle("Global Proportion of CpHs Methylated >10%\nContext Proportion") + 
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20), legend.title=element_blank()) + theme(legend.position="bottom")
ggplot(methCH[which(methCH$perc=="ten" & methCH$trinucleotide_context %in% c("CAG","CAC")),], 
       aes(x = Age, y = prop, colour = CellType)) + geom_point() + facet_grid(. ~ trinucleotide_context) +
  geom_smooth(method = "lm") +  
  theme_classic() + scale_colour_brewer(8, palette="Dark2") +
  ylab("mCpH / CpH") + xlab("") +
  ggtitle("Global Proportion of CpHs Methylated >10%\nTotal Proportion") + 
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20), legend.title=element_blank()) + theme(legend.position="bottom")
dev.off()


## Correlations of different context proportions by age

methCH$c_context = ifelse(methCH$trinucleotide_context %in% c("CTG", "CAG", "CCG"), "CHG", "CHH")
methCH$context_total = ifelse(methCH$c_context=="CHG", sum(freq[which(names(freq) %in% c("CTG", "CAG", "CCG"))]),sum(freq[-which(names(freq) %in% c("CTG", "CAG", "CCG"))]))
cont = data.table(methCH)[,sum(V1), by=c("c_context","variable", "perc","CellType","Age.Bin","Age","context_total")]
cont$prop = round(cont$V1 / cont$context_total,3)

cor = list(nonzero.glia.CHG = cor.test(cont[perc=="nonzero" & CellType=="Glia" & c_context=="CHG",,]$Age, cont[perc=="nonzero" & CellType=="Glia" & c_context=="CHG",,]$prop),
           nonzero.glia.CHH = cor.test(cont[perc=="nonzero" & CellType=="Glia" & c_context=="CHH",,]$Age, cont[perc=="nonzero" & CellType=="Glia" & c_context=="CHH",,]$prop),
           nonzero.neuron.CHG = cor.test(cont[perc=="nonzero" & CellType=="Neuron" & c_context=="CHG",,]$Age, cont[perc=="nonzero" & CellType=="Neuron" & c_context=="CHG",,]$prop),
           nonzero.neuron.CHH = cor.test(cont[perc=="nonzero" & CellType=="Neuron" & c_context=="CHH",,]$Age, cont[perc=="nonzero" & CellType=="Neuron" & c_context=="CHH",,]$prop),
           ten.glia.CHG = cor.test(cont[perc=="ten" & CellType=="Glia" & c_context=="CHG",,]$Age, cont[perc=="ten" & CellType=="Glia" & c_context=="CHG",,]$prop),
           ten.glia.CHH = cor.test(cont[perc=="ten" & CellType=="Glia" & c_context=="CHH",,]$Age, cont[perc=="ten" & CellType=="Glia" & c_context=="CHH",,]$prop),
           ten.neuron.CHG = cor.test(cont[perc=="ten" & CellType=="Neuron" & c_context=="CHG",,]$Age, cont[perc=="ten" & CellType=="Neuron" & c_context=="CHG",,]$prop),
           ten.neuron.CHH = cor.test(cont[perc=="ten" & CellType=="Neuron" & c_context=="CHH",,]$Age, cont[perc=="ten" & CellType=="Neuron" & c_context=="CHH",,]$prop))
write.csv(data.frame(comparison = names(cor), tstat=unlist(lapply(cor, function(x) x$statistic)), 
          pval = unlist(lapply(cor, function(x) x$p.value)), rho = unlist(lapply(cor, function(x) x$estimate))), quote=F,
          file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/non-CpG/proportion_mCH_byCellType_byAge_byTrinucleotideContext_correlation.csv")