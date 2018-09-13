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
names(dmrs) = c("Group 1 (G-,N+)","Group 2 (G0,N+)","Group 3 (G0,N-)","Group 4 (G+,N0)","Group 5 (G+,N-)","Group 6 (G-,N0)")
oo = lapply(dmrs, function(x) findOverlaps(x, makeGRangesFromDataFrame(DMR$Interaction)))
dmrs = lapply(oo, function(x) DMR$Interaction[subjectHits(x),])

DMRgr = lapply(c(DMR, dmrs), function(x) makeGRangesFromDataFrame(x[which(x$fwer<=0.05),], keep.extra.columns=T))


oo = lapply(DMRgr, function(x) findOverlaps(gr.clusters, x))
lapply(oo, function(x) length(unique(queryHits(x))))

Overlap = lapply(gwas, function(x) findOverlaps(x, gr.clusters))

df.clusters$regionID = paste0(df.clusters$seqnames,":",df.clusters$start,"-",df.clusters$end)
df.clusters$rnum = 1:length(gr.clusters)
df.clusters$CellType = ifelse(df.clusters$rnum %in% queryHits(oo$CellType), "yes","no")
df.clusters$Age = ifelse(df.clusters$rnum %in% queryHits(oo$Age), "yes","no")
df.clusters$Interaction = ifelse(df.clusters$rnum %in% queryHits(oo$Interaction), "yes","no")
df.clusters$`Group 1 (G-,N+)` = ifelse(df.clusters$rnum %in% queryHits(oo$`Group 1 (G-,N+)`), "yes","no")
df.clusters$`Group 2 (G0,N+)` = ifelse(df.clusters$rnum %in% queryHits(oo$`Group 2 (G0,N+)`), "yes","no")
df.clusters$`Group 3 (G0,N-)` = ifelse(df.clusters$rnum %in% queryHits(oo$`Group 3 (G0,N-)`), "yes","no")
df.clusters$`Group 4 (G+,N0)` = ifelse(df.clusters$rnum %in% queryHits(oo$`Group 4 (G+,N0)`), "yes","no")
df.clusters$`Group 5 (G+,N-)` = ifelse(df.clusters$rnum %in% queryHits(oo$`Group 5 (G+,N-)`), "yes","no")
df.clusters$`Group 6 (G-,N0)` = ifelse(df.clusters$rnum %in% queryHits(oo$`Group 6 (G-,N0)`), "yes","no")

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
#   GWAS           Model        OR         pval          fdr
#10 Park        CellType  2.264214 1.087124e-03 3.913646e-03
#11 Park             Age 14.327692 9.131028e-03 2.528592e-02
#12 Park     Interaction  5.468792 9.786045e-06 5.871627e-05
#15 Park Group 3 (G0,N-) 18.489768 1.364264e-06 9.822698e-06
#18 Park Group 6 (G-,N0)  8.619130 3.549918e-04 1.597463e-03

#19  T2D        CellType  3.555008 1.134421e-08 1.913923e-07
#21  T2D     Interaction  4.045316 1.078727e-03 3.913646e-03
#25  T2D Group 4 (G+,N0) 32.255689 1.351645e-04 6.951319e-04

#28  SCZ        CellType  1.817179 1.594936e-08 1.913923e-07
#29  SCZ             Age  4.222380 1.686004e-02 4.335438e-02
#30  SCZ     Interaction  2.684525 2.125739e-07 1.913165e-06
#33  SCZ Group 3 (G0,N-)  7.867389 2.841518e-10 1.022947e-08
#35  SCZ Group 5 (G+,N-)  4.951127 1.736070e-03 5.681683e-03
#36  SCZ Group 6 (G-,N0)  2.770898 2.871673e-03 8.615018e-03

## Test the genes within these loci

geneMapGR = makeGRangesFromDataFrame(geneMap, keep.extra.columns = T)
genes = lapply(gwas, function(x) findOverlaps(geneMapGR, x))
genes = lapply(genes, function(x) geneMapGR[queryHits(x)])
DMRgenes = lapply(DMRgr, function(x) unique(as.character(x[x$distToGene==0]$nearestID)))


## Enrichment in DMRs by cell type, age and interaction

geneuniverse = na.omit(unique(geneMap$gencodeID))
gwasGenes = lapply(genes, function(x) unique(as.character(x[which(x$gencodeID %in% geneuniverse)]$gencodeID))) # drop genes that are not present in the test set
sig = lapply(DMRgr, function(x) unique(as.character(x[x$distToGene==0]$nearestID)))
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

#             Model GWAS Odds.Ratio      P.Value          FDR
#10     Interaction Park   2.744962 1.188121e-02 4.752485e-02
#22 Group 3 (G0,N-) Park   6.010057 5.264312e-03 2.707360e-02
#3         CellType  T2D   3.568781 2.225193e-09 4.005347e-08
#11     Interaction  T2D   2.740540 7.972642e-03 3.587689e-02
#27 Group 4 (G+,N0)  T2D  21.701037 5.633151e-06 6.759781e-05
#4         CellType  SCZ   2.179325 1.027636e-12 3.699490e-11
#8              Age  SCZ   3.717664 1.326620e-02 4.775831e-02
#12     Interaction  SCZ   2.290865 2.198304e-05 1.978474e-04
#16 Group 1 (G-,N+)  SCZ   3.414708 1.236026e-04 7.416157e-04
#24 Group 3 (G0,N-)  SCZ   3.941949 2.864395e-05 2.062364e-04






