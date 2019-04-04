library(GenomicRanges)
library(bumphunter)
library(RColorBrewer)
library(jaffelab)


# qrsh -l mem_free=150G,h_vmem=150G,h_fsize=200G

load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/non-CpG/CH_object.rda")
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/non-CpG/CHneurons_object.rda")
load("/users/ajaffe/Lieber/Projects/RNAseq/n36/finalCode/gwasResults_lifted.rda")

pgc = read.delim("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/tableS3_Nature2014_pgc2_loci.txt",as.is=TRUE)
pgc$chr = ss(pgc$Position..hg19.,":")
tmp = ss(pgc$Position..hg19.,":",2)
pgc$start = as.numeric(ss(tmp, "-"))
pgc$end = as.numeric(ss(tmp, "-",2))
pgcGR = GRanges(pgc$chr, IRanges(pgc$start,pgc$end))
pgcGR$Dx = "SCZ"

gwas = c(as.list(split(gwasLift, gwasLift$Dx)), list(SCZ = pgcGR))


# Find overlaps with CpHs in all three models

CHgr = makeGRangesFromDataFrame(CHdt)

oo = lapply(gwas, function(x) findOverlaps(CHgr, x))
CHdt$Alz  = ifelse(CHdt$rnum %in% queryHits(oo$Alz), "yes","no")
CHdt$Park = ifelse(CHdt$rnum %in% queryHits(oo$Park), "yes","no")
CHdt$T2D  = ifelse(CHdt$rnum %in% queryHits(oo$T2D), "yes","no")
CHdt$SCZ  = ifelse(CHdt$rnum %in% queryHits(oo$SCZ), "yes","no")

CHneuronsgr = makeGRangesFromDataFrame(CHneuronsdt)

oo = lapply(gwas, function(x) findOverlaps(CHneuronsgr, x))
CHneuronsdt$Alz  = ifelse(CHneuronsdt$rnum %in% queryHits(oo$Alz), "yes","no")
CHneuronsdt$Park = ifelse(CHneuronsdt$rnum %in% queryHits(oo$Park), "yes","no")
CHneuronsdt$T2D  = ifelse(CHneuronsdt$rnum %in% queryHits(oo$T2D), "yes","no")
CHneuronsdt$SCZ  = ifelse(CHneuronsdt$rnum %in% queryHits(oo$SCZ), "yes","no")


## make contingency tables

CH = data.frame(CHdt)
tests = c("CT.sig","Age.sig","Int.sig")
tables = list(list(),list(),list(),list())
for (j in 1:length(gwas)) {
  for (i in 1:length(tests)) {
    tables[[j]][[i]] = data.frame(YesLocus = c(nrow(CH[CH[,colnames(CH)==names(gwas)[j]]=="yes" & CH[,colnames(CH)==tests[i]]=="FDR < 0.05",]),
                                               nrow(CH[CH[,colnames(CH)==names(gwas)[j]]=="yes" & CH[,colnames(CH)==tests[i]]=="FDR > 0.05",])),
                                  NoLocus = c(nrow(CH[CH[,colnames(CH)==names(gwas)[j]]=="no" & CH[,colnames(CH)==tests[i]]=="FDR < 0.05",]),
                                              nrow(CH[CH[,colnames(CH)==names(gwas)[j]]=="no" & CH[,colnames(CH)==tests[i]]=="FDR > 0.05",])), 
                                  row.names = c("YesSig","NoSig"))
  }
  names(tables[[j]]) = c("Cell Type","Age","Interaction")
}
names(tables) = c("Alzheimer's Disease","Parkinson's Disease","Type II Diabetes","Schizophrenia")

CHneurons = data.frame(CHneuronsdt)
for (j in 1:length(gwas)) {
  tables[[j]][["Age in Neurons"]] = data.frame(YesLocus = c(nrow(CHneurons[CHneurons[,colnames(CHneurons)==names(gwas)[j]]=="yes" & CHneurons$sig=="FDR < 0.05",]),
                                                            nrow(CHneurons[CHneurons[,colnames(CHneurons)==names(gwas)[j]]=="yes" & CHneurons$sig=="FDR > 0.05",])),
                                               NoLocus = c(nrow(CHneurons[CHneurons[,colnames(CHneurons)==names(gwas)[j]]=="no" & CHneurons$sig=="FDR < 0.05",]),
                                                           nrow(CHneurons[CHneurons[,colnames(CHneurons)==names(gwas)[j]]=="no" & CHneurons$sig=="FDR > 0.05",])), 
                                               row.names = c("YesSig","NoSig"))
}
                                                
                                      
fisher = lapply(tables, function(x) lapply(x, fisher.test))
df = do.call(rbind, Map(cbind, GWAS = as.list(names(fisher)), 
                        Map(cbind, Model = lapply(fisher, names), lapply(fisher, function(x) 
                          do.call(rbind, lapply(x, function(y) data.frame(OR = y$estimate,
                                                                          lower = y$conf.int[1],
                                                                          upper = y$conf.int[2],
                                                                          pval = y$p.value)))))))
df$fdr = p.adjust(df$pval, method= "fdr")

ta = do.call(rbind, lapply(tables, function(x) do.call(rbind, lapply(x, function(y) data.frame(YesGWAS.YesDMR = y[1,1], 
                                                                                               NoGWAS.YesDMR = y[1,2],
                                                                                               YesGWAS.NoDMR = y[2,1], 
                                                                                               NoGWAS.NoDMR = y[2,2])))))
df = cbind(df, ta)
rownames(df) = NULL


write.csv(df, quote=F, file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/DMR/GWAS_fisher_CpH_results.csv")


## Combine with DMR results into one table

df2 = read.csv("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/DMR/GWAS_fisher_DMR_results.csv")
df2$Model = df2$Extended
df2 = df2[,!colnames(df2) %in% c("X","Extended")]
df = rbind(cbind("C Context" = "CpG-derived DMR", df2), cbind("C Context" = "CpH", df))
df$fdr = p.adjust(df$pval, method= "fdr")

write.csv(df, quote=F, file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/DMR/GWAS_fisher_results.csv")



## Test the genes within these loci

geneMapGR = makeGRangesFromDataFrame(geneMap, keep.extra.columns = T)
genes = lapply(gwas, function(x) findOverlaps(geneMapGR, x))
genes = lapply(genes, function(x) geneMapGR[queryHits(x)])


## Enrichment in mCpH by cell type, age and interaction

geneuniverse = na.omit(unique(geneMap$EntrezID))
gwasGenes = lapply(genes, function(x) na.omit(unique(as.character(x$EntrezID))))
sig = list("Cell Type" = unique(CH[CH$distToGene==0 & CH$CT.sig=="FDR < 0.05","EntrezID"]),
           "Age" = unique(CH[CH$distToGene==0 & CH$Age.sig=="FDR < 0.05","EntrezID"]),
           "Interaction" = unique(CH[CH$distToGene==0 & CH$Int.sig=="FDR < 0.05","EntrezID"]), 
           "Age in Neurons" = unique(CHneurons[CHneurons$distToGene==0 & CHneurons$sig=="FDR < 0.05","EntrezID"]))
sig = lapply(sig, na.omit)
notsig = lapply(sig, function(y) geneuniverse[!(geneuniverse %in% y)])

DMRenrich = mapply(function(sig,notsig) lapply(gwasGenes, function(x) {
  DE_OVERLAP = c( sum( sig %in% x),sum(!(sig %in% x)))
  NOT_DE_OVERLAP= c(sum(notsig %in% x), sum(!(notsig %in% x)))
  enrich_table = cbind(DE_OVERLAP, NOT_DE_OVERLAP)
  res = fisher.test(enrich_table)
  dat=c(res$p.value, res$estimate)
  names(dat) <- c("P.Value","Odds.Ratio")
  return(dat)
}), sig, notsig, SIMPLIFY =F) 

for (i in 1:length(DMRenrich)) {
  names(DMRenrich[[i]]) = c("Alzheimers Disease","Parkinsons Disease","Type II Diabetes","Schizophrenia")
}


DMRenrich = do.call(rbind, Map(cbind, Model = as.list(names(DMRenrich)), lapply(DMRenrich, function(x) 
  do.call(rbind, Map(cbind, GWAS = as.list(names(x)), lapply(x, function(y) data.frame(Odds.Ratio = y["Odds.Ratio"], P.Value = y["P.Value"])))))))
rownames(DMRenrich) = NULL
DMRenrich$FDR = p.adjust(DMRenrich$P.Value, method="fdr")
DMRenrich = DMRenrich[order(DMRenrich$GWAS),]
write.csv(DMRenrich, quote=F, file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/DMR/GWASgenes_fisher_CpH_results.csv")