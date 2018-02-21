library(bsseq)
library(GEOquery)
library(recount)
library(jaffelab)
library('devtools')
library('limma')
library(GenomicRanges)


load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/single_CpGs/singleCpG_object.rda")
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/bsseq/lister/BSobj_matched_lister_minCov_3.Rdata")


## load data

pheno = pData(getGEO("GSE47966")[[1]])
pd = pheno[colnames(BSobj),c(1,2,8:13)]

for(i in 1:ncol(pd)) pd[,i] = as.character(pd[,i])
colnames(pd)[5:8] = ss(as.character(pd[1,5:8]), ": ")
for(i in 5:8) pd[,i] = ss(pd[,i], ": ", 2)

## fix last row
pd$gender[13] = "female"
pd$age[13] = NA
names(pd)[6] = "cellType"

## redo age
pd$ageNumeric = ss(pd$age, " ")
pd$ageNumeric = as.numeric(pd$ageNumeric)
pd$ageNumeric[grep("week", pd$age)] = -(40-pd$ageNumeric[grep("week", pd$age)])/52
pd$ageNumeric[grep("day", pd$age)] = pd$ageNumeric[grep("day", pd$age)]/365

pd$cellType = ss(pd$cellType, " ")
save(pd, file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/lister_phenotype_corrected.rda")

## Get methylation and phenotype info

cov <- getCoverage(BSobj, type = 'Cov')
cov = cov[,which(colnames(cov) %in% pd$geo_accession)]
cov.ge1 <- cov >= 1
cov.filt <- rowSums(cov.ge1) == ncol(cov)
print("Number of bases filtered")
table(cov.filt)
#FALSE     TRUE 
#18184422   480470 

BSobj <- BSobj[cov.filt, ]
rm(cov, cov.ge1, cov.filt)

meth = getMeth(BSobj, type = 'raw')


## Get chr coordinates and methylation values

gr = granges(BSobj)
df = as.data.frame(gr)

dim(meth) # 480470       13
rownames(meth) = paste0(df$seqnames, ":", df$start)

meth = meth[,which(colnames(meth) %in% pd$geo_accession[which(pd$ageNumeric>0 & pd$cellType=="all")])]

dim(meth[which(rowSums(meth)>0),]) # 477344 8
meth = meth[which(rowSums(meth)>0),]

## Explore with limma
model = with(pd[which(pd$cellType=="all" & pd$ageNumeric>0),], model.matrix(~ ageNumeric))
fit = lmFit(meth, design = model)
eb = eBayes(fit)

homLister = data.frame(ID = rownames(meth), coef = fit$coef[,"ageNumeric"], tstat = eb$t[,"ageNumeric"], pval = eb$p.value[,"ageNumeric"])
homLister$padj = p.adjust(homLister$pval, method = "fdr")

save(pd, homLister, fit, eb, model, 
     file = '/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/single_CpGs/Lister_homogenate_age_limma.rda')

## Get comparable CpG values from our sorted data

lgr = gr
load('/dcl01/lieber/ajaffe/lab/brain-epigenomics/bumphunting/BSobj_bsseqSmooth_Neuron_minCov_3.Rdata')

# postnatal pd table

spd = pData(BSobj)
spd$Race[spd$Race== "CAUC "] <- 'CAUC'
spd$Sex[spd$Sex == " M"] <- 'M'
spd$RIN <- as.numeric(gsub(" ", "", spd$RIN))
spd$pg.DNA.nuclei.input <- as.numeric(spd$pg.DNA.nuclei.input)
spd$Reads <- as.numeric(spd$Reads)
spd$Percent.GreaterThan.Q30 <- as.numeric(spd$Percent.GreaterThan.Q30)

smeth = getMeth(BSobj, type = 'raw')
sgr = granges(BSobj)
sdf = as.data.frame(sgr)
rownames(smeth) = paste0(sdf$seqnames, ":", sdf$start)
oo = findOverlaps(sgr, lgr)
smeth = smeth[queryHits(oo),]

fitCell = lmFit(smeth, model.matrix(~Age*Cell.Type, data=spd))
ebCell = ebayes(fitCell)

colnames(homLister) = c("ID", "dm_age_hom", "t_age_hom", "pval_age_hom", "qval_age_hom")
stats = data.frame(sID = rownames(smeth),
                   dm_age_cell = fitCell$coef[,2],
                   dm_int_cell = fitCell$coef[,4],
                   dm_type_cell = fitCell$coef[,3],
                   t_age_cell = ebCell$t[,2],
                   t_int_cell = ebCell$t[,4],
                   t_type_cell = ebCell$t[,3],
                   pval_age_cell = ebCell$p[,2],
                   pval_int_cell = ebCell$p[,4],
                   pval_type_cell = ebCell$p[,3])
stats$qval_age_cell = p.adjust(stats$pval_age_cell, "fdr")
stats$qval_int_cell = p.adjust(stats$pval_int_cell, "fdr")
stats$qval_type_cell = p.adjust(stats$pval_type_cell, "fdr")
stats = data.frame(homLister, stats[match(homLister$ID, stats$sID),])

colSums(stats[,grep("qval", colnames(stats))] < 0.05)


## filter

sigIndex = which(rowSums(stats[,grep("pval", colnames(stats))] < 1e-4) > 0)
statsSig = stats[sigIndex,]

## metrics
colSums(statsSig[,grep("pval", colnames(statsSig))] < 1e-4)
colSums(statsSig[,grep("pval", colnames(statsSig))] < 1e-6)
colSums(statsSig[,grep("qval", colnames(statsSig))] < 0.05)


############ Make Plots

# cell type?
pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/single_CpGs/figures/dm_cellType_vs_homAge_usingListerHomogenate.pdf")
par(mar=c(5,6,2,2),cex.axis=2,cex.lab=2)
plot(dm_type_cell ~ dm_age_hom,data=statsSig,
     subset = pval_age_hom < 1e-4,pch=21,bg="grey",
     xlim=c(-0.05,0.05),ylim=c(-0.8,0.8),cex=1.4,
     xlab = "Age (Homogenate)",ylab="Cell Type")
abline(h=0,v=0,lty=3)
dev.off()

# combined
pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/single_CpGs/figures/dm_totalAge_vs_homAge_usingListerHomogenate.pdf")
par(mar=c(5,6,2,2),cex.axis=2,cex.lab=2)
b =(statsSig$dm_age_cell + statsSig$dm_int_cell)
plot(b ~ dm_age_hom,data=statsSig,
     subset = pval_age_hom < 1e-4,pch=21,bg="grey",
     xlim=c(-0.05,0.05),ylim=c(-0.05,0.05),
     xlab = "Age (Homogenate)",ylab="Total Cell Type-Specific")
abline(h=0,v=0,lty=3)
dev.off()

## separate
pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/single_CpGs/figures/dm_separateAge_vs_homAge_usingListerHomogenate.pdf")
par(mar=c(5,6,2,2),cex.axis=2,cex.lab=2)
plot(dm_age_cell ~ dm_age_hom,data=statsSig,
     subset = pval_age_hom < 1e-4,pch=21,bg="grey",
     xlim=c(-0.05,0.05),ylim=c(-0.05,0.05),
     xlab = "Age (Homogenate)",ylab="Overall Age (Cell)")
abline(h=0,v=0,lty=3)

# interaction
plot(dm_int_cell ~ dm_age_hom,data=statsSig,
     subset = pval_age_hom < 1e-4,pch=21,bg="grey",
     xlim=c(-0.05,0.05),ylim=c(-0.05,0.05),
     xlab = "Age (Homogenate)",ylab="Specific Age (Cell)")
abline(h=0,v=0,lty=3)
dev.off()