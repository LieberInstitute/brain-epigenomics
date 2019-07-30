library(GenomicFeatures)
library(ggplot2)
library(GenomicRanges)
library(bumphunter)
library(data.table)
library(VennDiagram)
library("limma")
library("edgeR")


load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/DMV_gene_comps.rda")
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/sorted_nuclear_RNA/DE_limma_results_objects.rda", verbose=T)
load("/dcl01/lieber/ajaffe/CellSorting/RNAseq_pipeline/rawCounts_CellSorting_July5_n12.rda")
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/homogenate_RNA/DE_limma_results_homogenateRNAseq.rda", verbose=T)
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/brainseq_pipeline/polyA_unstranded/rse_gene_polyA_dlpfc_n41.Rdata", verbose = T)


## Format pd objects

hompd = colData(rse_gene)
metrics[grep("12", metrics$SampleID),"Age"] = 0.36
metrics[grep("13", metrics$SampleID),"Age"] = 14.01
metrics[grep("14", metrics$SampleID),"Age"] = 2.71
metrics$CellType = ifelse(metrics$NeuN=="NeuN_Minus", "Glia", "Neuron")

# format counts

prenCounts = assay(rse_gene)
prenMap = rowData(rse_gene)
identical(rownames(prenCounts), rownames(geneCounts))

polyaCounts = geneCounts[,grep("PolyA", colnames(geneCounts))]
riboCounts = geneCounts[,grep("Ribo", colnames(geneCounts))]
combinedCounts = polyaCounts + riboCounts

ctCounts = cbind(combinedCounts, prenCounts[,colnames(prenCounts) %in% rownames(hompd[hompd$Age<0,])]) 

newpd = rbind(data.frame(sampleID = rownames(metrics[grep("PolyA",rownames(metrics)),]), CellType = metrics[grep("PolyA",rownames(metrics)),"CellType"], 
                         Age = metrics[grep("PolyA",rownames(metrics)),"Age"]),
              data.frame(sampleID = rownames(hompd[hompd$Age<0,]), CellType = "Prenatal", Age = hompd[hompd$Age<0,"Age"]))
newpd = newpd[order(newpd$sampleID),]
newpd$CellType = factor(newpd$CellType, levels = c("Prenatal","Glia","Neuron"))
ctCounts = ctCounts[,order(colnames(ctCounts))]
rownames(newpd) = newpd$sampleID
ctCounts = ctCounts[rowSums(ctCounts)>0,]

hompd = hompd[hompd$Age>0,]
hompd[hompd$Age<1,"Age.Bin"] = "Infant"
hompd[hompd$Age>=1 & hompd$Age<=12,"Age.Bin"] = "Child"
hompd[hompd$Age>12 & hompd$Age<=17,"Age.Bin"] = "Teen"
hompd[hompd$Age>17,"Age.Bin"] = "Adult"
hompd$Age.Bin = factor(hompd$Age.Bin, levels = c("Infant","Child","Teen","Adult"))
geneCounts = prenCounts[,colnames(prenCounts) %in% rownames(hompd)]
geneCounts = geneCounts[rowSums(geneCounts)>0,]


## Differential expression

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/CREs/figures/DMV_gene_expr.pdf")

newpd$CellType = factor(newpd$CellType, levels = c("Prenatal","Glia","Neuron"))
design <- model.matrix(~ newpd$CellType)
dge <- DGEList(counts = ctCounts)
dge <- calcNormFactors(dge)
dat <- voom(dge, design, plot = TRUE)
fit_gene_ct <- lmFit(dat, design)
fit_gene_ct <- eBayes(fit_gene_ct)
newpd$CellType = factor(newpd$CellType, levels = c("Neuron","Prenatal","Glia"))
design <- model.matrix(~ newpd$CellType)
dge <- DGEList(counts = ctCounts)
dge <- calcNormFactors(dge)
dat <- voom(dge, design, plot = TRUE)
fit_gene_ct2 <- lmFit(dat, design)
fit_gene_ct2 <- eBayes(fit_gene_ct2)
hompd$Age.Bin = factor(hompd$Age.Bin, levels = c("Infant","Child","Teen","Adult"))
design <- model.matrix(~ hompd$Age.Bin)
dge <- DGEList(counts = geneCounts)
dge <- calcNormFactors(dge)
dat <- voom(dge, design, plot = TRUE)
fit_gene_age <- lmFit(dat, design)
fit_gene_age <- eBayes(fit_gene_age)
hompd$Age.Bin = factor(hompd$Age.Bin, levels = c("Adult","Infant","Child","Teen"))
design <- model.matrix(~ hompd$Age.Bin)
dge <- DGEList(counts = geneCounts)
dge <- calcNormFactors(dge)
dat <- voom(dge, design, plot = TRUE)
fit_gene_age2 <- lmFit(dat, design)
fit_gene_age2 <- eBayes(fit_gene_age2)
hompd$Age.Bin = factor(hompd$Age.Bin, levels = c("Child","Adult","Infant","Teen"))
design <- model.matrix(~ hompd$Age.Bin)
dge <- DGEList(counts = geneCounts)
dge <- calcNormFactors(dge)
dat <- voom(dge, design, plot = TRUE)
fit_gene_age3 <- lmFit(dat, design)
fit_gene_age3 <- eBayes(fit_gene_age3)
dev.off()


## Make results objects

ageRNAres = as.data.frame(cbind(fit_gene_age$coefficients, fit_gene_age$t, fit_gene_age$p.value,
                                fit_gene_age2$coefficients, fit_gene_age2$t, fit_gene_age2$p.value,
                                fit_gene_age3$coefficients, fit_gene_age3$t, fit_gene_age3$p.value))
ageRNAres = ageRNAres[,-grep("Intercept", colnames(ageRNAres))]
colnames(ageRNAres) = c("Coeff.Infant.Child", "Coeff.Infant.Teen", "Coeff.Infant.Adult", "Tstat.Infant.Child", "Tstat.Infant.Teen", "Tstat.Infant.Adult",
                        "pval.Infant.Child",  "pval.Infant.Teen",   "pval.Infant.Adult", "Coeff.Adult.Infant",   "Coeff.Adult.Child",  "Coeff.Adult.Teen",
                        "Tstat.Adult.Infant", "Tstat.Adult.Child",  "Tstat.Adult.Teen", "pval.Adult.Infant", "pval.Adult.Child",  "pval.Adult.Teen",
                        "Coeff.Child.Adult",  "Coeff.Child.Infant", "Coeff.Child.Teen", "Tstat.Child.Adult",  "Tstat.Child.Infant", "Tstat.Child.Teen",
                        "pval.Child.Adult",  "pval.Child.Infant", "pval.Child.Teen")
ageRNAres = ageRNAres[,colnames(ageRNAres) %in% c("Coeff.Infant.Child", "Coeff.Infant.Teen", "Coeff.Infant.Adult", "Tstat.Infant.Child", 
                                                  "Tstat.Infant.Teen", "Tstat.Infant.Adult","pval.Infant.Child",  "pval.Infant.Teen",   
                                                  "pval.Infant.Adult", "Coeff.Adult.Child", "Coeff.Adult.Teen",
                                                  "Tstat.Adult.Child", "Tstat.Adult.Teen", "pval.Adult.Child", "pval.Adult.Teen",
                                                  "Coeff.Child.Teen", "Tstat.Child.Teen", "pval.Child.Teen")]
ageRNAres$FDR.Infant.Child = p.adjust(ageRNAres$pval.Infant.Child, method = "fdr")
ageRNAres$FDR.Infant.Teen = p.adjust(ageRNAres$pval.Infant.Teen, method = "fdr")
ageRNAres$FDR.Infant.Adult = p.adjust(ageRNAres$pval.Infant.Adult, method = "fdr")
ageRNAres$FDR.Adult.Child = p.adjust(ageRNAres$pval.Adult.Child, method = "fdr")
ageRNAres$FDR.Adult.Teen = p.adjust(ageRNAres$pval.Adult.Teen, method = "fdr")
ageRNAres$FDR.Child.Teen = p.adjust(ageRNAres$pval.Child.Teen, method = "fdr")
ageRNAres$gencodeID = rownames(ageRNAres)

ctRNAres = as.data.frame(cbind(fit_gene_ct$coefficients, fit_gene_ct$t, fit_gene_ct$p.value,
                               fit_gene_ct2$coefficients, fit_gene_ct2$t, fit_gene_ct2$p.value))
ctRNAres = ctRNAres[,-grep("Intercept", colnames(ctRNAres))]
colnames(ctRNAres) = c("Coeff.Prenatal.Glia","Coeff.Prenatal.Neuron", "Tstat.Prenatal.Glia", "Tstat.Prenatal.Neuron",
                        "pval.Prenatal.Glia",  "pval.Prenatal.Neuron", "Coeff.Neuron.Prenatal", "Coeff.Neuron.Glia",
                        "Tstat.Neuron.Prenatal", "Tstat.Neuron.Glia", "pval.Neuron.Prenatal", "pval.Neuron.Glia")
ctRNAres = ctRNAres[,colnames(ctRNAres) %in% c("Coeff.Prenatal.Glia","Coeff.Prenatal.Neuron", "Tstat.Prenatal.Glia", "Tstat.Prenatal.Neuron",
                                               "pval.Prenatal.Glia",  "pval.Prenatal.Neuron", "Coeff.Neuron.Glia",
                                               "Tstat.Neuron.Glia", "pval.Neuron.Glia")]
ctRNAres$FDR.Prenatal.Glia = p.adjust(ctRNAres$pval.Prenatal.Glia, method = "fdr")
ctRNAres$FDR.Prenatal.Neuron = p.adjust(ctRNAres$pval.Prenatal.Neuron, method = "fdr")
ctRNAres$FDR.Neuron.Glia = p.adjust(ctRNAres$pval.Neuron.Glia, method = "fdr")
ctRNAres$gencodeID = rownames(ctRNAres)


## Associate with gene expression

DMV.CTcomps = Map(cbind, lapply(DMV.CTcomps, function(x) x[,1:13]), lapply(DMV.CTcomps, function(x) ctRNAres[match(x$gencodeID, ctRNAres$gencodeID),]))
elementNROWS(DMV.CTcomps)
CT = do.call(rbind, list(PnotN = data.frame(Comp = "PnotN", tstat = DMV.CTcomps$PnotN$Tstat.Prenatal.Neuron, fdr = DMV.CTcomps$PnotN$pval.Prenatal.Neuron),
                         NnotP = data.frame(Comp = "NnotP", tstat = -(DMV.CTcomps$NnotP$Tstat.Prenatal.Neuron), fdr = DMV.CTcomps$NnotP$pval.Prenatal.Neuron),
                         PnotG = data.frame(Comp = "PnotG", tstat = DMV.CTcomps$PnotG$Tstat.Prenatal.Glia, fdr = DMV.CTcomps$PnotG$pval.Prenatal.Glia),
                         GnotP = data.frame(Comp = "GnotP", tstat = -(DMV.CTcomps$GnotP$Tstat.Prenatal.Glia), fdr = DMV.CTcomps$GnotP$pval.Prenatal.Glia),
                         GnotN = data.frame(Comp = "GnotN", tstat = -(DMV.CTcomps$GnotN$Tstat.Neuron.Glia), fdr = DMV.CTcomps$GnotN$pval.Neuron.Glia),
                         NnotG = data.frame(Comp = "NnotG", tstat = DMV.CTcomps$NnotG$Tstat.Neuron.Glia, fdr = DMV.CTcomps$NnotG$pval.Neuron.Glia)))
CT$Sig = ifelse(CT$fdr<=0.05,"FDR<0.05","FDR>0.05")

# Genes that are escaping the DMV state (likely accumulating DNAm) are higher expressed in the cell type in which the gene escapes

DMV.Agecomps = Map(cbind, lapply(DMV.Agecomps, function(x) x[,1:13]), lapply(DMV.Agecomps, function(x) ageRNAres[match(x$gencodeID, ageRNAres$gencodeID),]))
elementNROWS(DMV.Agecomps)
AgeN = do.call(rbind, list(InotC = data.frame(Comp = "InotC", tstat = DMV.Agecomps$InotC$Tstat.Infant.Child, fdr = DMV.Agecomps$InotC$pval.Infant.Child),
                           CnotI = data.frame(Comp = "CnotI", tstat = -(DMV.Agecomps$CnotI$Tstat.Infant.Child), fdr = DMV.Agecomps$CnotI$pval.Infant.Child),
                           InotT = data.frame(Comp = "InotT", tstat = DMV.Agecomps$InotT$Tstat.Infant.Teen, fdr = DMV.Agecomps$InotT$pval.Infant.Teen),
                           TnotI = data.frame(Comp = "TnotI", tstat = -(DMV.Agecomps$TnotI$Tstat.Infant.Teen), fdr = DMV.Agecomps$TnotI$pval.Infant.Teen),
                           InotA = data.frame(Comp = "InotA", tstat = DMV.Agecomps$InotA$Tstat.Infant.Adult, fdr = DMV.Agecomps$InotA$pval.Infant.Adult),
                           AnotI = data.frame(Comp = "AnotI", tstat = -(DMV.Agecomps$AnotI$Tstat.Infant.Adult), fdr = DMV.Agecomps$AnotI$pval.Infant.Adult),
                           CnotT = data.frame(Comp = "CnotT", tstat = DMV.Agecomps$CnotT$Tstat.Child.Teen, fdr = DMV.Agecomps$CnotT$pval.Child.Teen),
                           TnotC = data.frame(Comp = "TnotC", tstat = -(DMV.Agecomps$TnotC$Tstat.Child.Teen), fdr = DMV.Agecomps$TnotC$pval.Child.Teen),
                           CnotA = data.frame(Comp = "CnotA", tstat = -(DMV.Agecomps$CnotA$Tstat.Adult.Child), fdr = DMV.Agecomps$CnotA$pval.Adult.Child),
                           AnotC = data.frame(Comp = "AnotC", tstat = DMV.Agecomps$AnotC$Tstat.Adult.Child, fdr = DMV.Agecomps$AnotC$pval.Adult.Child),
                           TnotA = data.frame(Comp = "TnotA", tstat = -(DMV.Agecomps$TnotA$Tstat.Adult.Teen), fdr = DMV.Agecomps$TnotA$pval.Adult.Teen),
                           AnotT = data.frame(Comp = "AnotT", tstat = DMV.Agecomps$AnotT$Tstat.Adult.Teen, fdr = DMV.Agecomps$AnotT$pval.Adult.Teen)))
AgeN$Sig = ifelse(AgeN$fdr<=0.05,"FDR<0.05","FDR>0.05")

DMV.AgecompsG = Map(cbind, lapply(DMV.AgecompsG, function(x) x[,1:13]), lapply(DMV.AgecompsG, function(x) ageRNAres[match(x$gencodeID, ageRNAres$gencodeID),]))
elementNROWS(DMV.AgecompsG)
AgeG = do.call(rbind, list(InotC = data.frame(Comp = "InotC", tstat = DMV.AgecompsG$InotC$Tstat.Infant.Child, fdr = DMV.AgecompsG$InotC$pval.Infant.Child),
                           CnotI = data.frame(Comp = "CnotI", tstat = -(DMV.AgecompsG$CnotI$Tstat.Infant.Child), fdr = DMV.AgecompsG$CnotI$pval.Infant.Child),
                           InotT = data.frame(Comp = "InotT", tstat = DMV.AgecompsG$InotT$Tstat.Infant.Teen, fdr = DMV.AgecompsG$InotT$pval.Infant.Teen),
                           TnotI = data.frame(Comp = "TnotI", tstat = -(DMV.AgecompsG$TnotI$Tstat.Infant.Teen), fdr = DMV.AgecompsG$TnotI$pval.Infant.Teen),
                           InotA = data.frame(Comp = "InotA", tstat = DMV.AgecompsG$InotA$Tstat.Infant.Adult, fdr = DMV.AgecompsG$InotA$pval.Infant.Adult),
                           AnotI = data.frame(Comp = "AnotI", tstat = -(DMV.AgecompsG$AnotI$Tstat.Infant.Adult), fdr = DMV.AgecompsG$AnotI$pval.Infant.Adult),
                           CnotT = data.frame(Comp = "CnotT", tstat = DMV.AgecompsG$CnotT$Tstat.Child.Teen, fdr = DMV.AgecompsG$CnotT$pval.Child.Teen),
                           TnotC = data.frame(Comp = "TnotC", tstat = -(DMV.AgecompsG$TnotC$Tstat.Child.Teen), fdr = DMV.AgecompsG$TnotC$pval.Child.Teen),
                           CnotA = data.frame(Comp = "CnotA", tstat = -(DMV.AgecompsG$CnotA$Tstat.Adult.Child), fdr = DMV.AgecompsG$CnotA$pval.Adult.Child),
                           AnotC = data.frame(Comp = "AnotC", tstat = DMV.AgecompsG$AnotC$Tstat.Adult.Child, fdr = DMV.AgecompsG$AnotC$pval.Adult.Child),
                           TnotA = data.frame(Comp = "TnotA", tstat = -(DMV.AgecompsG$TnotA$Tstat.Adult.Teen), fdr = DMV.AgecompsG$TnotA$pval.Adult.Teen),
                           AnotT = data.frame(Comp = "AnotT", tstat = DMV.AgecompsG$AnotT$Tstat.Adult.Teen, fdr = DMV.AgecompsG$AnotT$pval.Adult.Teen)))
AgeG$Sig = ifelse(AgeG$fdr<=0.05,"FDR<0.05","FDR>0.05")

expr = rbind(cbind(CT, Model = "Cell Type"), cbind(AgeN, Model = "Neurons"), cbind(AgeG, Model = "Glia"))
write.csv(expr, quote=F, file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/DMV_gene_expression.csv")

expr = read.csv("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/DMV_gene_expression.csv")
expr$Model = factor(expr$Model, levels = c("Cell Type", "Glia", "Neurons"))

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/CREs/figures/DMV_gene_expression.pdf", 
    height=4.75, width = 6)
ggplot(expr[which(expr$Model=="Cell Type" & expr$fdr<=0.05),], 
       aes(x = Comp, y = tstat)) + geom_boxplot() +
  scale_fill_brewer(8, palette="Dark2") +
  labs(fill="") + geom_hline(yintercept=0, linetype="dotted") +
  ylab("T Statistic") + xlab("") +
  ggtitle("DMV Gene Expression Enrichment") +
  theme(title = element_text(size = 20),
        text = element_text(size = 20))

ggplot(expr[which(expr$Model!="Cell Type" & expr$fdr<=0.05),], 
       aes(x = Comp, y = tstat, fill = Model)) + geom_boxplot() +
  facet_grid(Model ~ .) + guides(fill=FALSE) +
  labs(fill="") + geom_hline(yintercept=0, linetype="dotted") +
  ylab("T Statistic") + xlab("") + scale_fill_brewer(8, palette="Dark2") +
  ggtitle("DMV Gene Expression Enrichment") +
  theme(title = element_text(size = 20), text = element_text(size = 20),
        axis.text.x = element_text(angle = 25, hjust = 1))
dev.off()

table(expr[which(expr$Model=="Cell Type"),"fdr"]<=0.05)
#FALSE  TRUE 
# 8484 19290 

table(expr[which(expr$Model=="Neurons"),"fdr"]<=0.05)
#FALSE  TRUE 
# 3530  1086 
1086/(1086+3530) # 0.2352686

table(expr[which(expr$Model=="Glia"),"fdr"]<=0.05)
#FALSE  TRUE 
# 3779  1043 
1043/(1043+3779) # 0.2163003


## How does this relate to DNAm at these genes?

load('/dcl01/lieber/ajaffe/lab/brain-epigenomics/bsseq/bsobj_by_chr/allChrs_postNatal_cleaned_nonCG_noHomogenate_highCov.Rdata')
BSobj_ch = BSobj
meth_ch = getMeth(BSobj_ch, type = 'raw')
methMap_ch = granges(BSobj_ch)

load('/dcl01/lieber/ajaffe/lab/brain-epigenomics/bumphunting/BSobj_bsseqSmooth_Neuron_minCov_3.Rdata')
postpd = pData(BSobj)
postpd$Race[postpd$Race== "CAUC "] <- 'CAUC'
postpd$Sex[postpd$Sex == " M"] <- 'M'
postpd$RIN <- as.numeric(gsub(" ", "", postpd$RIN))
postpd$pg.DNA.nuclei.input <- as.numeric(postpd$pg.DNA.nuclei.input)
postpd$Reads <- as.numeric(postpd$Reads)
postpd$Percent.GreaterThan.Q30 <- as.numeric(postpd$Percent.GreaterThan.Q30)
meth =getMeth(BSobj, type = 'raw')
methMap = granges(BSobj)


## CpG

ids = lapply(c(DMV.CTcomps,DMV.Agecomps,DMV.AgecompsG), function(x) paste0(x$Chr,":",x$Start,"-",x$End,":",x$Strand))
ids = lapply(ids, function(x) GRanges(x[which(x != "NA:NA-NA:NA")]))
oo = lapply(ids, function(x) findOverlaps(x, methMap))

meanMeth = lapply(oo, function(x) do.call("rbind", lapply(split(subjectHits(x), factor(queryHits(x))), function(ii) colMeans(t(t(meth[ii,]))))))
meanMeth = lapply(meanMeth, reshape2::melt)

tstat.CT = list(PnotN.NnotP = t.test(meanMeth$GnotN[meanMeth$GnotN$Var2 %in% pd[pd$Cell.Type=="Neuron","Data.ID"],"value"],
                                     meanMeth$NnotG[meanMeth$NnotG$Var2 %in% pd[pd$Cell.Type=="Glia","Data.ID"],"value"]),
                PnotG.GnotP = t.test(meanMeth$GnotN[meanMeth$GnotN$Var2 %in% pd[pd$Cell.Type=="Glia","Data.ID"],"value"],
                                     meanMeth$NnotG[meanMeth$NnotG$Var2 %in% pd[pd$Cell.Type=="Glia","Data.ID"],"value"]),
                GnotN.NnotG = t.test(meanMeth$GnotN[meanMeth$GnotN$Var2 %in% pd[pd$Cell.Type=="Glia","Data.ID"],"value"],
                                     meanMeth$NnotG[meanMeth$NnotG$Var2 %in% pd[pd$Cell.Type=="Neuron","Data.ID"],"value"]))
ct = data.frame(Tstat = unlist(lapply(tstat.CT, function(x) x$statistic)), mean1 = unlist(lapply(tstat.CT, function(x) x$estimate[1])),
                mean2 = unlist(lapply(tstat.CT, function(x) x$estimate[2])), pval = unlist(lapply(tstat.CT, function(x) x$p.value)), comps = names(tstat.CT), row.names = NULL)

names(meanMeth)[29:46] = paste0(names(meanMeth)[29:46],".G")
tstat.N = list(InotC.CnotI = t.test(meanMeth$InotC[meanMeth$InotC$Var2 %in% pd[pd$Age<=1 & pd$Cell.Type=="Neuron","Data.ID"],"value"],
                                    meanMeth$CnotI[meanMeth$CnotI$Var2 %in% pd[pd$Age>1 & pd$Age<=12 & pd$Cell.Type=="Neuron","Data.ID"],"value"]),
               InotT.TnotI = t.test(meanMeth$InotT[meanMeth$InotT$Var2 %in% pd[pd$Age<=1 & pd$Cell.Type=="Neuron","Data.ID"],"value"],
                                    meanMeth$TnotI[meanMeth$TnotI$Var2 %in% pd[pd$Age>12 & pd$Age<=17 & pd$Cell.Type=="Neuron","Data.ID"],"value"]),
               InotA.AnotI = t.test(meanMeth$InotA[meanMeth$InotA$Var2 %in% pd[pd$Age<=1 & pd$Cell.Type=="Neuron","Data.ID"],"value"],
                                    meanMeth$AnotI[meanMeth$AnotI$Var2 %in% pd[pd$Age>17 & pd$Cell.Type=="Neuron","Data.ID"],"value"]),
               CnotT.TnotC = t.test(meanMeth$CnotT[meanMeth$CnotT$Var2 %in% pd[pd$Age>1 & pd$Age<=12 & pd$Cell.Type=="Neuron","Data.ID"],"value"],
                                    meanMeth$TnotC[meanMeth$TnotC$Var2 %in% pd[pd$Age>12 & pd$Age<=17 & pd$Cell.Type=="Neuron","Data.ID"],"value"]),
               CnotA.AnotC = t.test(meanMeth$CnotA[meanMeth$CnotA$Var2 %in% pd[pd$Age>1 & pd$Age<=12 & pd$Cell.Type=="Neuron","Data.ID"],"value"],
                                    meanMeth$AnotC[meanMeth$AnotC$Var2 %in% pd[pd$Age>17 & pd$Cell.Type=="Neuron","Data.ID"],"value"]),
               TnotA.AnotT = t.test(meanMeth$TnotA[meanMeth$TnotA$Var2 %in% pd[pd$Age>12 & pd$Age<=17 & pd$Cell.Type=="Neuron","Data.ID"],"value"],
                                    meanMeth$AnotT[meanMeth$AnotT$Var2 %in% pd[pd$Age>17 & pd$Cell.Type=="Neuron","Data.ID"],"value"]))
neuro = data.frame(Tstat = unlist(lapply(tstat.N, function(x) x$statistic)), mean1 = unlist(lapply(tstat.N, function(x) x$estimate[1])),
                   mean2 = unlist(lapply(tstat.N, function(x) x$estimate[2])), pval = unlist(lapply(tstat.N, function(x) x$p.value)), comps = names(tstat.N), row.names = NULL)

tstat.G = list(InotC.CnotI = t.test(meanMeth$InotC.G[meanMeth$InotC.G$Var2 %in% pd[pd$Age<=1 & pd$Cell.Type=="Glia","Data.ID"],"value"],
                                    meanMeth$CnotI.G[meanMeth$CnotI.G$Var2 %in% pd[pd$Age>1 & pd$Age<=12 & pd$Cell.Type=="Glia","Data.ID"],"value"]),
               InotT.TnotI = t.test(meanMeth$InotT.G[meanMeth$InotT.G$Var2 %in% pd[pd$Age<=1 & pd$Cell.Type=="Glia","Data.ID"],"value"],
                                    meanMeth$TnotI.G[meanMeth$TnotI.G$Var2 %in% pd[pd$Age>12 & pd$Age<=17 & pd$Cell.Type=="Glia","Data.ID"],"value"]),
               InotA.AnotI = t.test(meanMeth$InotA.G[meanMeth$InotA.G$Var2 %in% pd[pd$Age<=1 & pd$Cell.Type=="Glia","Data.ID"],"value"],
                                    meanMeth$AnotI.G[meanMeth$AnotI.G$Var2 %in% pd[pd$Age>17 & pd$Cell.Type=="Glia","Data.ID"],"value"]),
               CnotT.TnotC = t.test(meanMeth$CnotT.G[meanMeth$CnotT.G$Var2 %in% pd[pd$Age>1 & pd$Age<=12 & pd$Cell.Type=="Glia","Data.ID"],"value"],
                                    meanMeth$TnotC.G[meanMeth$TnotC.G$Var2 %in% pd[pd$Age>12 & pd$Age<=17 & pd$Cell.Type=="Glia","Data.ID"],"value"]),
               CnotA.AnotC = t.test(meanMeth$CnotA.G[meanMeth$CnotA.G$Var2 %in% pd[pd$Age>1 & pd$Age<=12 & pd$Cell.Type=="Glia","Data.ID"],"value"],
                                    meanMeth$AnotC.G[meanMeth$AnotC.G$Var2 %in% pd[pd$Age>17 & pd$Cell.Type=="Glia","Data.ID"],"value"]),
               TnotA.AnotT = t.test(meanMeth$TnotA.G[meanMeth$TnotA.G$Var2 %in% pd[pd$Age>12 & pd$Age<=17 & pd$Cell.Type=="Glia","Data.ID"],"value"],
                                    meanMeth$AnotT.G[meanMeth$AnotT.G$Var2 %in% pd[pd$Age>17 & pd$Cell.Type=="Glia","Data.ID"],"value"]))
gli = data.frame(Tstat = unlist(lapply(tstat.G, function(x) x$statistic)), mean1 = unlist(lapply(tstat.G, function(x) x$estimate[1])),
                 mean2 = unlist(lapply(tstat.G, function(x) x$estimate[2])), pval = unlist(lapply(tstat.G, function(x) x$p.value)), comps = names(tstat.G), row.names = NULL)

dfCG = rbind(cbind(neuro, Comp = "Neurons Age"), cbind(gli, Comp = "Glia Age"), cbind(ct, Comp = "Cell Type"))



## CpH


ids = lapply(c(DMV.CTcomps,DMV.Agecomps,DMV.AgecompsG), function(x) paste0(x$Chr,":",x$Start,"-",x$End,":",x$Strand))
ids = lapply(ids, function(x) GRanges(x[which(x != "NA:NA-NA:NA")]))
oo_ch = lapply(ids, function(x) findOverlaps(x, methMap_ch))
meanMeth_ch = lapply(oo_ch, function(x) do.call("rbind", lapply(split(subjectHits(x), factor(queryHits(x), levels=1:length(unique(queryHits(x))))), 
                                                                function(ii) colMeans(t(t(meth_ch[ii,]))))))
meanMeth_ch = lapply(meanMeth_ch, reshape2::melt)


GnotN.NnotG = t.test(meanMeth_ch$GnotN[meanMeth_ch$GnotN$Var2 %in% pd[pd$Cell.Type=="Glia","Data.ID"],"value"],
                     meanMeth_ch$NnotG[meanMeth_ch$NnotG$Var2 %in% pd[pd$Cell.Type=="Neuron","Data.ID"],"value"])
ct = data.frame(Tstat = GnotN.NnotG$statistic, mean1 = GnotN.NnotG$estimate[1], mean2 = GnotN.NnotG$estimate[2], pval = GnotN.NnotG$p.value, comps = "GnotN.NnotG")

names(meanMeth_ch)[29:46] = paste0(names(meanMeth_ch)[29:46],".G")
tstat.N = list(InotC.CnotI = t.test(meanMeth_ch$InotC[meanMeth_ch$InotC$Var2 %in% pd[pd$Age<=1 & pd$Cell.Type=="Neuron","Data.ID"],"value"],
                                    meanMeth_ch$CnotI[meanMeth_ch$CnotI$Var2 %in% pd[pd$Age>1 & pd$Age<=12 & pd$Cell.Type=="Neuron","Data.ID"],"value"]),
               InotT.TnotI = t.test(meanMeth_ch$InotT[meanMeth_ch$InotT$Var2 %in% pd[pd$Age<=1 & pd$Cell.Type=="Neuron","Data.ID"],"value"],
                                    meanMeth_ch$TnotI[meanMeth_ch$TnotI$Var2 %in% pd[pd$Age>12 & pd$Age<=17 & pd$Cell.Type=="Neuron","Data.ID"],"value"]),
               InotA.AnotI = t.test(meanMeth_ch$InotA[meanMeth_ch$InotA$Var2 %in% pd[pd$Age<=1 & pd$Cell.Type=="Neuron","Data.ID"],"value"],
                                    meanMeth_ch$AnotI[meanMeth_ch$AnotI$Var2 %in% pd[pd$Age>17 & pd$Cell.Type=="Neuron","Data.ID"],"value"]),
               CnotT.TnotC = t.test(meanMeth_ch$CnotT[meanMeth_ch$CnotT$Var2 %in% pd[pd$Age>1 & pd$Age<=12 & pd$Cell.Type=="Neuron","Data.ID"],"value"],
                                    meanMeth_ch$TnotC[meanMeth_ch$TnotC$Var2 %in% pd[pd$Age>12 & pd$Age<=17 & pd$Cell.Type=="Neuron","Data.ID"],"value"]),
               CnotA.AnotC = t.test(meanMeth_ch$CnotA[meanMeth_ch$CnotA$Var2 %in% pd[pd$Age>1 & pd$Age<=12 & pd$Cell.Type=="Neuron","Data.ID"],"value"],
                                    meanMeth_ch$AnotC[meanMeth_ch$AnotC$Var2 %in% pd[pd$Age>17 & pd$Cell.Type=="Neuron","Data.ID"],"value"]),
               TnotA.AnotT = t.test(meanMeth_ch$TnotA[meanMeth_ch$TnotA$Var2 %in% pd[pd$Age>12 & pd$Age<=17 & pd$Cell.Type=="Neuron","Data.ID"],"value"],
                                    meanMeth_ch$AnotT[meanMeth_ch$AnotT$Var2 %in% pd[pd$Age>17 & pd$Cell.Type=="Neuron","Data.ID"],"value"]))
neuro = data.frame(Tstat = unlist(lapply(tstat.N, function(x) x$statistic)), mean1 = unlist(lapply(tstat.N, function(x) x$estimate[1])),
                   mean2 = unlist(lapply(tstat.N, function(x) x$estimate[2])), pval = unlist(lapply(tstat.N, function(x) x$p.value)), comps = names(tstat.N), row.names = NULL)

tstat.G = list(InotC.CnotI = t.test(meanMeth_ch$InotC.G[meanMeth_ch$InotC.G$Var2 %in% pd[pd$Age<=1 & pd$Cell.Type=="Glia","Data.ID"],"value"],
                                    meanMeth_ch$CnotI.G[meanMeth_ch$CnotI.G$Var2 %in% pd[pd$Age>1 & pd$Age<=12 & pd$Cell.Type=="Glia","Data.ID"],"value"]),
               InotT.TnotI = t.test(meanMeth_ch$InotT.G[meanMeth_ch$InotT.G$Var2 %in% pd[pd$Age<=1 & pd$Cell.Type=="Glia","Data.ID"],"value"],
                                    meanMeth_ch$TnotI.G[meanMeth_ch$TnotI.G$Var2 %in% pd[pd$Age>12 & pd$Age<=17 & pd$Cell.Type=="Glia","Data.ID"],"value"]),
               InotA.AnotI = t.test(meanMeth_ch$InotA.G[meanMeth_ch$InotA.G$Var2 %in% pd[pd$Age<=1 & pd$Cell.Type=="Glia","Data.ID"],"value"],
                                    meanMeth_ch$AnotI.G[meanMeth_ch$AnotI.G$Var2 %in% pd[pd$Age>17 & pd$Cell.Type=="Glia","Data.ID"],"value"]),
               CnotT.TnotC = t.test(meanMeth_ch$CnotT.G[meanMeth_ch$CnotT.G$Var2 %in% pd[pd$Age>1 & pd$Age<=12 & pd$Cell.Type=="Glia","Data.ID"],"value"],
                                    meanMeth_ch$TnotC.G[meanMeth_ch$TnotC.G$Var2 %in% pd[pd$Age>12 & pd$Age<=17 & pd$Cell.Type=="Glia","Data.ID"],"value"]),
               CnotA.AnotC = t.test(meanMeth_ch$CnotA.G[meanMeth_ch$CnotA.G$Var2 %in% pd[pd$Age>1 & pd$Age<=12 & pd$Cell.Type=="Glia","Data.ID"],"value"],
                                    meanMeth_ch$AnotC.G[meanMeth_ch$AnotC.G$Var2 %in% pd[pd$Age>17 & pd$Cell.Type=="Glia","Data.ID"],"value"]),
               TnotA.AnotT = t.test(meanMeth_ch$TnotA.G[meanMeth_ch$TnotA.G$Var2 %in% pd[pd$Age>12 & pd$Age<=17 & pd$Cell.Type=="Glia","Data.ID"],"value"],
                                    meanMeth_ch$AnotT.G[meanMeth_ch$AnotT.G$Var2 %in% pd[pd$Age>17 & pd$Cell.Type=="Glia","Data.ID"],"value"]))
gli = data.frame(Tstat = unlist(lapply(tstat.G, function(x) x$statistic)), mean1 = unlist(lapply(tstat.G, function(x) x$estimate[1])),
                 mean2 = unlist(lapply(tstat.G, function(x) x$estimate[2])), pval = unlist(lapply(tstat.G, function(x) x$p.value)), comps = names(tstat.G), row.names = NULL)


dfmeth = rbind(cbind(dfCG[,-1], Context = "CpG"), cbind(rbind(cbind(neuro, Comp = "Neurons Age"), cbind(gli, Comp = "Glia Age"), cbind(ct, Comp = "Cell Type")), Context = "CpH"))
dfmeth$FDR = p.adjust(dfmeth$pval, method = "fdr")

write.csv(dfmeth, quote=F, file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/t.test_DMV_gene_methylation.csv")



library(PWMEnrich.Hsapiens.background)
data(MotifDb.Hsap)
tfs = unlist(lapply(MotifDb.Hsap, function(x) x$name))
tfs = tfs[-grep("UW",tfs)]

not = tfs[!tfs %in% geneMap$Symbol]
not = na.omit(not[-grep("::", not, fixed=T)])

nomatch = c("SMCR7L","DUX4","MYF","MZF1_1-4","MZF1_5-13","RORA_1","RORA_2","MIZF","EWSR1-FLI1","ZNF238",
            "ZNF306","POU5F1P1","BHLHB2","BHLHB3","CART1","RAXL1","TRP53","TRP73","ZNF435","HNF1B")
nomatch[which(nomatch %in% geneMap$Symbol)]      
nomatch[which(nomatch %in% good)]      

tfs = unique(c(tfs[tfs %in% geneMap$Symbol], c("HINFP1","RBM8B"), unique(unlist(strsplit(not[grep("::", not, fixed=T)], "::", fixed=T)))))
tfs = geneMap[which(geneMap$Symbol %in% tfs),] # 830 genes
