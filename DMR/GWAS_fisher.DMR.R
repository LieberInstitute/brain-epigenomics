library(GenomicRanges)
library(bumphunter)
library(jaffelab)
library(scales)
library(gplots)
library(cowplot)
library(RColorBrewer)



# qrsh -l mem_free=150G,h_vmem=150G,h_fsize=200G

load("/users/ajaffe/Lieber/Projects/RNAseq/n36/finalCode/gwasResults_lifted.rda")
load('/dcl01/lieber/ajaffe/lab/brain-epigenomics/bumphunting/BSobj_bsseqSmooth_Neuron_minCov_3.Rdata')

pgc = read.delim("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/tableS3_Nature2014_pgc2_loci.txt",as.is=TRUE)
pgc$chr = ss(pgc$Position..hg19.,":")
tmp = ss(pgc$Position..hg19.,":",2)
pgc$start = as.numeric(ss(tmp, "-"))
pgc$end = as.numeric(ss(tmp, "-",2))
pgcGR = GRanges(pgc$chr, IRanges(pgc$start,pgc$end))
pgcGR$Dx = "SCZ"

gwas = c(as.list(split(gwasLift, gwasLift$Dx)), list(SCZ = pgcGR))
names(gwas) = c("Alzheimer's Disease", "Parkinson's Disease", "Type II Diabetes", "Schizophrenia")


# Identify all CpG clusters in the genome

gr <- granges(BSobj)
cl = clusterMaker( chr = as.character(seqnames(gr)), 
                   pos = start(gr),   maxGap = 1000)
gr.clusters = split(gr, cl)
gr.clusters = unlist(range(gr.clusters))
df.clusters = as.data.frame(gr.clusters)


# Find overlaps with DMRS in all three models

brain_categories <- readRDS("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/ldsc/categories.rds")

brain_categories_df <- data.frame(Category = names(brain_categories)[4:13],	
                                  Extended = c("Cell Type (Glia > Neuron)", "Cell Type (Neuron > Glia)", 
                                               "Age (Younger > Older)", "Age (Older > Younger)", 
                                               "Group 1 (Decreasing Glial; Increasing Neuronal)", 
                                               "Group 2 (Static Glial; Increasing Neuronal)", 
                                               "Group 3 (Static Glial; Decreasing Neuronal)", 
                                               "Group 4 (Increasing Glial; Static Neuronal)", 
                                               "Group 5 (Increasing Glial; Decreasing Neuronal)", 
                                               "Group 6 (Decreasing Glial; Static Neuronal)"))
brain_categories = brain_categories[4:13]

oo = lapply(brain_categories, function(x) findOverlaps(gr.clusters, x))
lapply(oo, function(x) length(unique(queryHits(x))))

Overlap = lapply(gwas, function(x) findOverlaps(gr.clusters, x))

df.clusters$regionID = paste0(df.clusters$seqnames,":",df.clusters$start,"-",df.clusters$end)
df.clusters$rnum = 1:length(gr.clusters)

for (i in 1:length(brain_categories)) {
  df.clusters[,names(brain_categories)[i]] = ifelse(df.clusters$rnum %in% queryHits(oo[[i]]), "yes","no")
}

for (i in 1:length(gwas)) {
  df.clusters[,names(Overlap)[i]] = ifelse(df.clusters$rnum %in% queryHits(Overlap[[i]]), "yes","no")
}


## make contingency tables

tables = list(list(),list(),list(),list())
for (j in 1:length(gwas)) {
  for (i in 1:length(names(brain_categories))) {
    tables[[j]][[i]] = data.frame(YesLocus = c(nrow(df.clusters[df.clusters[,colnames(df.clusters)==names(gwas)[j]]=="yes" & 
                                                                  df.clusters[,colnames(df.clusters)==names(brain_categories)[i]]=="yes",]),
                                      nrow(df.clusters[df.clusters[,colnames(df.clusters)==names(gwas)[j]]=="yes" & 
                                                         df.clusters[,colnames(df.clusters)==names(brain_categories)[i]]=="no",])),
                           NoLocus = c(nrow(df.clusters[df.clusters[,colnames(df.clusters)==names(gwas)[j]]=="no" & 
                                                          df.clusters[,colnames(df.clusters)==names(brain_categories)[i]]=="yes",]),
                                     nrow(df.clusters[df.clusters[,colnames(df.clusters)==names(gwas)[j]]=="no" & 
                                                        df.clusters[,colnames(df.clusters)==names(brain_categories)[i]]=="no",])), 
                           row.names = c("YesDMR","NoDMR"))
    }
  names(tables[[j]]) = names(brain_categories)
  }
names(tables) = names(gwas)

fisher = lapply(tables, function(x) lapply(x, fisher.test))
df = do.call(rbind, Map(cbind, GWAS = as.list(names(fisher)), 
                        Map(cbind, Model = lapply(fisher, names), 
                            lapply(fisher, function(x) do.call(rbind, lapply(x, function(y) data.frame(OR = y$estimate, 
                                                                                                       lower = y$conf.int[1],
                                                                                                       upper = y$conf.int[2], 
                                                                                                       pval = y$p.value)))))))
df$fdr = p.adjust(df$pval, method= "fdr")
rownames(df) = NULL
df$Extended = brain_categories_df[match(df$Model, brain_categories_df$Category),"Extended"]

ta = do.call(rbind, lapply(tables, function(x) do.call(rbind, lapply(x, function(y) data.frame(YesGWAS.YesDMR = y[1,1], 
                                                                                               NoGWAS.YesDMR = y[1,2],
                                                                                               YesGWAS.NoDMR = y[2,1], 
                                                                                               NoGWAS.NoDMR = y[2,2])))))
rownames(ta) = NULL
df = cbind(df, ta)

write.csv(df, quote=F, file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/DMR/GWAS_fisher_DMR_results.csv")

df[df$fdr<=0.05,c("GWAS","Model","OR","fdr")]
#   GWAS               Model        OR          fdr
#1   Alz Celltype_HypoNeuron  3.606411 3.362091e-03
#11 Park Celltype_HypoNeuron  2.984641 4.139130e-03
#12 Park   Celltype_HypoGlia  2.654812 3.854868e-03
#13 Park      Age_Decreasing 39.117463 4.936101e-03
#17 Park          Gr3_cdDMRs 18.489768 1.364264e-05
#20 Park          Gr6_cdDMRs  8.619130 2.366612e-03
#21  T2D Celltype_HypoNeuron  2.537779 2.214559e-02
#22  T2D   Celltype_HypoGlia  4.438351 7.070127e-08
#28  T2D          Gr4_cdDMRs 32.255689 1.081316e-03
#31  SCZ Celltype_HypoNeuron  2.908269 4.465282e-14
#32  SCZ   Celltype_HypoGlia  1.407507 4.142161e-02
#33  SCZ      Age_Decreasing 11.864517 2.877558e-03
#37  SCZ          Gr3_cdDMRs  7.867389 5.683037e-09
#39  SCZ          Gr5_cdDMRs  4.951127 5.786899e-03
#40  SCZ          Gr6_cdDMRs  2.770898 8.835915e-03


## Plot results

df$sig = ifelse(df$fdr<=0.05, "Significant", "Not Significant")
df$sig = factor(df$sig)
df$Extended = factor(df$Extended, levels = c("Cell Type (Glia > Neuron)", "Cell Type (Neuron > Glia)",
                                             "Age (Older > Younger)","Age (Younger > Older)",                          
                                             "Group 1 (Decreasing Glial; Increasing Neuronal)",
                                             "Group 2 (Static Glial; Increasing Neuronal)",    
                                             "Group 3 (Static Glial; Decreasing Neuronal)",    
                                             "Group 4 (Increasing Glial; Static Neuronal)",    
                                             "Group 5 (Increasing Glial; Decreasing Neuronal)",
                                             "Group 6 (Decreasing Glial; Static Neuronal)"))


pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/DMR/figures/disease_GWAS_loci_enrichment_dotplot.pdf",
    height = 2.5, width = 6)

ggplot(data = df, aes(x = Extended, y = OR, col = Extended, shape = sig)) +
  geom_point() + geom_pointrange(aes(ymin = lower, ymax = upper)) +
  facet_grid(. ~ GWAS) +
  theme_bw() + xlab("Feature") + ylab("Odds Ratio") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  geom_hline(yintercept = 0, lty = 2) +
  scale_shape_manual(values = c(1, 16)) +
  scale_size_manual(values = c(1, 2)) +
  scale_colour_manual(values = c("#FC8D62","#66C2A5","#377EB8","#E41A1C",
                                 "#1B9E77","#D95F02","#7570B3","#E7298A","#66A61E","#E6AB02")) +
  guides(col = FALSE, shape = FALSE)
ggplot(data = df, aes(x = Extended, y = log(OR), col = Extended, shape = sig)) +
  geom_point() + geom_pointrange(aes(ymin = log(lower), ymax = log(upper))) +
  facet_grid(. ~ GWAS) +
  theme_bw() + xlab("Feature") + ylab("log(Odds Ratio)") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  geom_hline(yintercept = 0, lty = 2) +
  scale_shape_manual(values = c(1, 16)) +
  scale_size_manual(values = c(1, 2)) +
  scale_colour_manual(values = c("#FC8D62","#66C2A5","#377EB8","#E41A1C",
                                 "#1B9E77","#D95F02","#7570B3","#E7298A","#66A61E","#E6AB02")) +
  guides(col = FALSE, shape = FALSE)
dev.off()



## Test the genes within these loci

geneMapGR = makeGRangesFromDataFrame(geneMap, keep.extra.columns = T)
genes = lapply(gwas, function(x) findOverlaps(geneMapGR, x))
genes = lapply(genes, function(x) geneMapGR[queryHits(x)])


## Enrichment in DMRs by cell type, age and interaction

geneuniverse = na.omit(unique(geneMap$EntrezID))
gwasGenes = lapply(genes, function(x) unique(as.character(x$EntrezID)))
sig = lapply(DMRgr, function(x) unique(as.character(x[x$distToGene==0]$EntrezID)))
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


DMRenrich = do.call(rbind, Map(cbind, Model = as.list(names(DMRenrich)), lapply(DMRenrich, function(x) 
            do.call(rbind, Map(cbind, GWAS = as.list(names(x)), lapply(x, function(y) data.frame(Odds.Ratio = y["Odds.Ratio"], P.Value = y["P.Value"])))))))
rownames(DMRenrich) = NULL
DMRenrich$FDR = p.adjust(DMRenrich$P.Value, method="fdr")
DMRenrich = DMRenrich[order(DMRenrich$GWAS),]
write.csv(DMRenrich, quote=F, file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/DMR/GWASgenes_fisher_DMR_results.csv")


DMRenrich[DMRenrich$FDR<=0.05,]

#                       Model GWAS Odds.Ratio      P.Value          FDR
#2  Hypomethylated in Neurons Park   2.214807 1.539965e-02 4.738353e-02
#14               Interaction Park   2.744962 1.188121e-02 4.320441e-02
#26           Group 3 (G0,N-) Park   6.010057 5.264312e-03 2.339694e-02
#3  Hypomethylated in Neurons  T2D   4.106108 1.217342e-07 2.434684e-06
#7     Hypomethylated in Glia  T2D   3.133085 9.052546e-06 9.052546e-05
#15               Interaction  T2D   2.740540 7.972642e-03 3.189057e-02
#31           Group 4 (G+,N0)  T2D  21.701037 5.633151e-06 7.510868e-05
#4  Hypomethylated in Neurons  SCZ   3.077028 5.005744e-17 2.002298e-15
#8     Hypomethylated in Glia  SCZ   1.624898 6.919659e-04 3.459829e-03
#12                       Age  SCZ   3.717664 1.326620e-02 4.422065e-02
#16               Interaction  SCZ   2.290865 2.198304e-05 1.758643e-04
#20           Group 1 (G-,N+)  SCZ   3.414708 1.236026e-04 7.063006e-04
#28           Group 3 (G0,N-)  SCZ   3.941949 2.864395e-05 1.909597e-04
