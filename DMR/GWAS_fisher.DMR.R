library(GenomicRanges)
library(bumphunter)
library(RColorBrewer)
library(jaffelab)


# qrsh -l mem_free=150G,h_vmem=150G,h_fsize=200G

load("/users/ajaffe/Lieber/Projects/RNAseq/n36/finalCode/gwasResults_lifted.rda")
load('/dcl01/lieber/ajaffe/lab/brain-epigenomics/bumphunting/BSobj_bsseqSmooth_Neuron_minCov_3.Rdata')
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/DMR/DMR_objects.rda")
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/bumphunting/rda/limma_Neuron_CpGs_minCov_3_ageInfo_dmrs.Rdata")

pgc = read.delim("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/tableS3_Nature2014_pgc2_loci.txt",as.is=TRUE)
pgc$chr = ss(pgc$Position..hg19.,":")
tmp = ss(pgc$Position..hg19.,":",2)
pgc$start = as.numeric(ss(tmp, "-"))
pgc$end = as.numeric(ss(tmp, "-",2))
pgcGR = GRanges(pgc$chr, IRanges(pgc$start,pgc$end))
pgcGR$Dx = "SCZ"

gwas = c(as.list(split(gwasLift, gwasLift$Dx)), list(SCZ = pgcGR))


# Identify all CpG clusters in the genome

gr <- granges(BSobj)
cl = clusterMaker( chr = as.character(seqnames(gr)), 
                   pos = start(gr),   maxGap = 1000)
gr.clusters = split(gr, cl)
gr.clusters = unlist(range(gr.clusters))
df.clusters = as.data.frame(gr.clusters)

# Find overlaps with DMRS in all three models

dmrs = split(dmrs, dmrs$k6cluster_label)
names(dmrs) = c("Group 1 (G-N+)","Group 2 (G0N+)","Group 3 (G0N-)","Group 4 (G+N0)","Group 5 (G+N-)","Group 6 (G-N0)")
oo = lapply(dmrs, function(x) findOverlaps(x, makeGRangesFromDataFrame(DMR$Interaction)))
dmrs = lapply(oo, function(x) DMR$Interaction[subjectHits(x),])
CT = split(DMR$CellType, DMR$CellType$Dir)
names(CT) = c("Hypomethylated in Neurons", "Hypomethylated in Glia")
DMRgr = lapply(c(CT, DMR[names(DMR)!="CellType"], dmrs), function(x) makeGRangesFromDataFrame(x[which(x$fwer<=0.05),], keep.extra.columns = T))

oo = lapply(DMRgr, function(x) findOverlaps(gr.clusters, x))
lapply(oo, function(x) length(unique(queryHits(x))))

Overlap = lapply(gwas, function(x) findOverlaps(x, gr.clusters))

df.clusters$regionID = paste0(df.clusters$seqnames,":",df.clusters$start,"-",df.clusters$end)
df.clusters$rnum = 1:length(gr.clusters)
df.clusters$"Hypomethylated in Neurons" = ifelse(df.clusters$rnum %in% queryHits(oo$"Hypomethylated in Neurons"), "yes","no")
df.clusters$"Hypomethylated in Glia" = ifelse(df.clusters$rnum %in% queryHits(oo$"Hypomethylated in Glia"), "yes","no")
df.clusters$Age = ifelse(df.clusters$rnum %in% queryHits(oo$Age), "yes","no")
df.clusters$Interaction = ifelse(df.clusters$rnum %in% queryHits(oo$Interaction), "yes","no")
df.clusters$`Group 1 (G-N+)` = ifelse(df.clusters$rnum %in% queryHits(oo$`Group 1 (G-N+)`), "yes","no")
df.clusters$`Group 2 (G0N+)` = ifelse(df.clusters$rnum %in% queryHits(oo$`Group 2 (G0N+)`), "yes","no")
df.clusters$`Group 3 (G0N-)` = ifelse(df.clusters$rnum %in% queryHits(oo$`Group 3 (G0N-)`), "yes","no")
df.clusters$`Group 4 (G+N0)` = ifelse(df.clusters$rnum %in% queryHits(oo$`Group 4 (G+N0)`), "yes","no")
df.clusters$`Group 5 (G+N-)` = ifelse(df.clusters$rnum %in% queryHits(oo$`Group 5 (G+N-)`), "yes","no")
df.clusters$`Group 6 (G-N0)` = ifelse(df.clusters$rnum %in% queryHits(oo$`Group 6 (G-N0)`), "yes","no")

df.clusters$Alz  = ifelse(df.clusters$rnum %in% subjectHits(Overlap$Alz), "yes","no")
df.clusters$Park = ifelse(df.clusters$rnum %in% subjectHits(Overlap$Park), "yes","no")
df.clusters$T2D  = ifelse(df.clusters$rnum %in% subjectHits(Overlap$T2D), "yes","no")
df.clusters$SCZ  = ifelse(df.clusters$rnum %in% subjectHits(Overlap$SCZ), "yes","no")


## make contingency tables

tables = list(list(),list(),list(),list())
for (j in 1:length(gwas)) {
  for (i in 1:length(names(DMRgr))) {
    tables[[j]][[i]] = data.frame(YesLocus = c(nrow(df.clusters[df.clusters[,colnames(df.clusters)==names(gwas)[j]]=="yes" & df.clusters[,colnames(df.clusters)==names(DMRgr)[i]]=="yes",]),
                                      nrow(df.clusters[df.clusters[,colnames(df.clusters)==names(gwas)[j]]=="yes" & df.clusters[,colnames(df.clusters)==names(DMRgr)[i]]=="no",])),
                           NoLocus = c(nrow(df.clusters[df.clusters[,colnames(df.clusters)==names(gwas)[j]]=="no" & df.clusters[,colnames(df.clusters)==names(DMRgr)[i]]=="yes",]),
                                     nrow(df.clusters[df.clusters[,colnames(df.clusters)==names(gwas)[j]]=="no" & df.clusters[,colnames(df.clusters)==names(DMRgr)[i]]=="no",])), 
                           row.names = c("YesDMR","NoDMR"))
    }
  names(tables[[j]]) = names(DMRgr)
  }
names(tables) = names(gwas)

fisher = lapply(tables, function(x) lapply(x, fisher.test))
df = do.call(rbind, Map(cbind, GWAS = as.list(names(fisher)), Map(cbind, Model = lapply(fisher, names), lapply(fisher, function(x) do.call(rbind, lapply(x, function(y) data.frame(OR = y$estimate, pval = y$p.value)))))))
df$fdr = p.adjust(df$pval, method= "fdr")
rownames(df) = NULL

write.csv(df, quote=F, file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/DMR/GWAS_fisher_DMR_results.csv")

df[df$fdr<=0.05,]
#   GWAS                     Model        OR         pval          fdr
#1   Alz Hypomethylated in Neurons  3.606411 6.724181e-04 2.988525e-03
#11 Park Hypomethylated in Neurons  2.984641 1.034783e-03 3.595757e-03
#12 Park    Hypomethylated in Glia  2.654812 8.673454e-04 3.469382e-03
#13 Park                       Age 14.327692 9.131028e-03 2.282757e-02
#14 Park               Interaction  5.468792 9.786045e-06 6.524030e-05
#17 Park           Group 3 (G0,N-) 18.489768 1.364264e-06 1.091411e-05
#20 Park           Group 6 (G-,N0)  8.619130 3.549918e-04 1.774959e-03
#21  T2D Hypomethylated in Neurons  2.537779 7.750956e-03 2.066922e-02
#22  T2D    Hypomethylated in Glia  4.438351 5.302595e-09 7.070127e-08
#24  T2D               Interaction  4.045316 1.078727e-03 3.595757e-03
#28  T2D           Group 4 (G+,N0) 32.255689 1.351645e-04 7.723688e-04
#31  SCZ Hypomethylated in Neurons  2.908269 1.116321e-15 4.465282e-14
#32  SCZ    Hypomethylated in Glia  1.407507 1.553310e-02 3.654848e-02
#33  SCZ                       Age  4.222380 1.686004e-02 3.746675e-02
#34  SCZ               Interaction  2.684525 2.125739e-07 2.125739e-06
#37  SCZ           Group 3 (G0,N-)  7.867389 2.841518e-10 5.683037e-09
#39  SCZ           Group 5 (G+,N-)  4.951127 1.736070e-03 5.341753e-03
#40  SCZ           Group 6 (G-,N0)  2.770898 2.871673e-03 8.204779e-03

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
