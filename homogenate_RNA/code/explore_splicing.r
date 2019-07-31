library(GenomicRanges)
library(bumphunter)
library(SGSeq)


# load meQTLs

load('/dcl01/lieber/ajaffe/lab/brain-epigenomics/meth_vs_expr/rda/meqtl_mres_gene_using_near_meth11_proteincoding.Rdata', verbose = TRUE)
class(mres)
names(mres)
names(mres$nonCpG)
mres$nonCpG$expr
mres$nonCpG$eqtls
gmres = mres
load('/dcl01/lieber/ajaffe/lab/brain-epigenomics/meth_vs_expr/rda/meqtl_mres_exon_using_near_meth11_proteincoding.Rdata', verbose = TRUE)
emres = mres
load('/dcl01/lieber/ajaffe/lab/brain-epigenomics/meth_vs_expr/rda/meqtl_mres_psi_using_near_meth11_proteincoding.Rdata', verbose = TRUE)
pmres = mres

genH = list(gene = unique(gmres$nonCpG$eqtls$gene), exon = unique(emres$nonCpG$eqtls$gene), psi = unique(pmres$nonCpG$eqtls$gene))
elementNROWS(genH)
# gene exon  psi 
# 2413 4413  593 
genG = list(gene = unique(gmres$CpG$eqtls$gene), exon = unique(emres$CpG$eqtls$gene), psi = unique(pmres$CpG$eqtls$gene))
elementNROWS(genG)
#gene exon  psi 
#3370 4806 1885 
gen = mapply(function(h,g) unique(c(h,g)), genH,genG, SIMPLIFY = F)
elementNROWS(gen)
#gene exon  psi 
#5023 8817 2113 

1-(elementNROWS(gen)/(elementNROWS(genH)+elementNROWS(genG)))
# gene      exon       psi 
# 0.1314197 0.0436056 0.1472962 
# small number of features identified in both contexts with association


# How many of the sub-gene features are also gene features in a context?

genH = list(gene = unique(gmres$nonCpG$eqtls$gene), exon = unique(rowData(emres$nonCpG$expr)$gencodeID), psi = unique(pmres$nonCpG$eqtls$gene))
elementNROWS(genH)
# gene exon  psi 
# 2413 1355  593 
genG = list(gene = unique(gmres$CpG$eqtls$gene), exon = unique(rowData(emres$CpG$expr)$gencodeID), psi = unique(pmres$CpG$eqtls$gene))
elementNROWS(genG)
#gene exon  psi 
#3370 1646 1885 
gen = mapply(function(h,g) unique(c(h,g)), genH,genG, SIMPLIFY = F)


length(genG$exon[genG$exon %in% genG$gene])/length(genG$exon)*100 # 52.55164
length(genG$psi[genG$psi %in% genG$gene])/length(genG$psi)*100 # 35.91512
length(genH$exon[genH$exon %in% genH$gene])/length(genH$exon)*100 # 65.97786
length(genH$psi[genH$psi %in% genH$gene])/length(genH$psi)*100 # 37.94266


# Format DMR lists

spl = list(Exon = list(CpH = list(Pos = unique(rowData(emres$nonCpG$expr)$EntrezID[which(emres$nonCpG$eqtl$statistic>0)]),
                                  Neg = unique(rowData(emres$nonCpG$expr)$EntrezID[which(emres$nonCpG$eqtl$statistic<0)])),
                       CpG = list(Pos = unique(rowData(emres$CpG$expr)$EntrezID[which(emres$CpG$eqtl$statistic>0)]),
                                  Neg = unique(rowData(emres$CpG$expr)$EntrezID[which(emres$CpG$eqtl$statistic<0)]))),
           PSI = list(CpH = list(Pos = unique(pmres$nonCpG$eqtl$gene[which(pmres$nonCpG$eqtl$statistic>0)]),
                                 Neg = unique(pmres$nonCpG$eqtl$gene[which(pmres$nonCpG$eqtl$statistic<0)])),
                      CpG = list(Pos = unique(pmres$CpG$eqtl$gene[which(pmres$CpG$eqtl$statistic>0)]),
                                 Neg = unique(pmres$CpG$eqtl$gene[which(pmres$CpG$eqtl$statistic<0)]))),
           Gene = list(CpH = list(Pos = unique(rowData(gmres$nonCpG$expr)$EntrezID[which(gmres$nonCpG$eqtl$statistic>0)]),
                                  Neg = unique(rowData(gmres$nonCpG$expr)$EntrezID[which(gmres$nonCpG$eqtl$statistic<0)])),
                       CpG = list(Pos = unique(rowData(gmres$CpG$expr)$EntrezID[which(gmres$CpG$eqtl$statistic>0)]),
                                  Neg = unique(rowData(gmres$CpG$expr)$EntrezID[which(gmres$CpG$eqtl$statistic<0)]))))
spl$PSI = lapply(spl$PSI, function(x) lapply(x, function(y) unique(geneMap[match(y, geneMap$gencodeID),"EntrezID"])))


aej_sets = openxlsx::read.xlsx('/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/Birnbaum_2013_AJP_Supplementary_table.xlsx')

## Enrichment in eqtls

geneuniverse = na.omit(unique(geneMap$EntrezID))
aej_sets_expressed = aej_sets[which(aej_sets$EntrezGene.ID %in% geneuniverse), ] # drop genes that are not present in the test set
splitSets = split(aej_sets_expressed[!aej_sets_expressed$Gene.Set %in% c("SCZ PGC GWAS","BPAD GWAS","SCZ Meta-analysis"),], 
                  aej_sets_expressed$Gene.Set[!aej_sets_expressed$Gene.Set %in% c("SCZ PGC GWAS","BPAD GWAS","SCZ Meta-analysis")])
notspl = lapply(spl, function(y) lapply(y, function(z) lapply(z, function(x) geneuniverse[!(geneuniverse %in% x)])))

DMRenrich = list()
for (h in 1:length(spl)) {
  DMRenrich[[h]] = list()
  for (i in 1:length(spl[[h]])) {
    DMRenrich[[h]][[i]] = list()
    for (j in 1:length(spl[[h]][[i]])) {
      DMRenrich[[h]][[i]][[j]] = lapply(splitSets, function(x) {
        DE_OVERLAP = c( sum( spl[[h]][[i]][[j]] %in% x$EntrezGene.ID),sum(!(spl[[h]][[i]][[j]] %in% x$EntrezGene.ID)))
        NOT_DE_OVERLAP= c(sum(notspl[[h]][[i]][[j]] %in% x$EntrezGene.ID), sum(!(notspl[[h]][[i]][[j]] %in% x$EntrezGene.ID)))
        enrich_table = cbind(DE_OVERLAP, NOT_DE_OVERLAP)
        res = fisher.test(enrich_table)
        dat=c(res$p.value, res$estimate, res$conf.int[1], res$conf.int[2])
        names(dat) <- c("pval","OR", "lower","upper")
        return(dat)
      })
    }
    names(DMRenrich[[h]][[i]]) = names(spl[[h]][[i]])
  }
  names(DMRenrich[[h]]) = names(spl[[h]])
}
names(DMRenrich) = names(spl)

enrich_table = list()
for (h in 1:length(spl)) {
  enrich_table[[h]] = list()
  for (i in 1:length(spl[[h]])) {
    enrich_table[[h]][[i]] = list()
    for (j in 1:length(spl[[h]][[i]])) {
      enrich_table[[h]][[i]][[j]] = lapply(splitSets, function(x) {
        DE_OVERLAP = c( sum( spl[[h]][[i]][[j]] %in% x$EntrezGene.ID),sum(!(spl[[h]][[i]][[j]] %in% x$EntrezGene.ID)))
        NOT_DE_OVERLAP= c(sum(notspl[[h]][[i]][[j]] %in% x$EntrezGene.ID), sum(!(notspl[[h]][[i]][[j]] %in% x$EntrezGene.ID)))
        cbind(DE_OVERLAP, NOT_DE_OVERLAP)
      })
    }
    names(enrich_table[[h]][[i]]) = names(spl[[h]][[i]])
  }
  names(enrich_table[[h]]) = names(spl[[h]])
}
names(enrich_table) = names(spl)


df = lapply(DMRenrich, function(x) lapply(x, function(z) lapply(z, function(h) 
                do.call(rbind, Map(cbind,lapply(h, function(y) data.frame(P.Value = y["pval"],
                                                                          Odds.Ratio = y["OR"],
                                                                          Lower = y["lower"],
                                                                          Upper = y["upper"])), GeneSet = as.list(names(h)))))))
df = lapply(df, function(x) lapply(x, function(z) do.call(rbind, Map(cbind, z, Direction = as.list(names(z))))))
df = lapply(df, function(x) do.call(rbind, Map(cbind, x, Context = as.list(names(x)))))
df = do.call(rbind, Map(cbind, df, Feature = as.list(names(df))))

df$FDR = p.adjust(df$P.Value, method = "fdr")
rownames(df) = NULL

ta = lapply(enrich_table, function(x) lapply(x, function(z) lapply(z, function(h) 
              do.call(rbind, Map(cbind,lapply(h, function(y) data.frame(YesGeneSet.YesExpr = y[1,1],
                                                                        NoGeneSet.YesExpr = y[2,1],
                                                                        YesGeneSet.NoExpr = y[1,2],
                                                                        NoGeneSet.NoExpr = y[2,2])), GeneSet = as.list(names(h)))))))
ta = lapply(ta, function(x) lapply(x, function(z) do.call(rbind, Map(cbind, z, Direction = as.list(names(z))))))
ta = lapply(ta, function(x) do.call(rbind, Map(cbind, x, Context = as.list(names(x)))))
ta = do.call(rbind, Map(cbind, ta, Feature = as.list(names(ta))))

rownames(ta) = NULL
df = cbind(df, ta[,1:4])

write.csv(df,file="/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/Birnbaum_geneSet_enrichment_splicingResults_updated.csv",quote=F)

df[which(df$FDR<=0.05),colnames(df) %in% c("GeneSet","Odds.Ratio", "FDR","Direction","Context","Feature")]
#   Odds.Ratio      GeneSet Direction Context Feature          FDR
#4    8.306668          NDD       Pos     CpH    Exon 1.553603e-02
#16   3.764205 ASD DATABASE       Pos     CpG    Exon 5.311130e-03
#37   3.752510 ASD DATABASE       Neg     CpH     PSI 3.163048e-04
#42   2.928525      SCZ SNV       Neg     CpH     PSI 1.553603e-02
#44   3.082437 ASD DATABASE       Pos     CpG     PSI 3.196115e-07
#46   4.219748          NDD       Pos     CpG     PSI 2.087327e-02
#49   2.525660      SCZ SNV       Pos     CpG     PSI 3.163048e-04
#51   2.583172 ASD DATABASE       Neg     CpG     PSI 3.163048e-04
#53   4.386468          NDD       Neg     CpG     PSI 2.520312e-02
#56   2.116158      SCZ SNV       Neg     CpG     PSI 2.141714e-02
#58   2.484406 ASD DATABASE       Pos     CpH    Gene 9.791413e-03
#65   2.787712 ASD DATABASE       Neg     CpH    Gene 1.392849e-05
#72   4.330241 ASD DATABASE       Pos     CpG    Gene 1.362547e-14
#74   3.959069          NDD       Pos     CpG    Gene 2.514168e-02
#77   2.273289      SCZ SNV       Pos     CpG    Gene 2.568041e-03
#79   2.057475 ASD DATABASE       Neg     CpG    Gene 3.596949e-03

df = read.csv("./Desktop/BAMS/Birnbaum_geneSet_enrichment_splicingResults_updated.csv")
table(df[which(df$FDR<=0.05),colnames(df) %in% c("Direction","Context","Feature")])

## Enrichment in GWAS genes

load("/users/ajaffe/Lieber/Projects/RNAseq/n36/finalCode/gwasResults_lifted.rda")
pgc = read.delim("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/tableS3_Nature2014_pgc2_loci.txt",as.is=TRUE)
pgc$chr = ss(pgc$Position..hg19.,":")
tmp = ss(pgc$Position..hg19.,":",2)
pgc$start = as.numeric(ss(tmp, "-"))
pgc$end = as.numeric(ss(tmp, "-",2))
pgcGR = GRanges(pgc$chr, IRanges(pgc$start,pgc$end))
pgcGR$Dx = "SCZ"

gwas = c(as.list(split(gwasLift, gwasLift$Dx)), list(SCZ = pgcGR))

geneMapGR = makeGRangesFromDataFrame(geneMap, keep.extra.columns = T)
genes = lapply(gwas, function(x) findOverlaps(geneMapGR, x))
genes = lapply(genes, function(x) geneMapGR[queryHits(x)])

gwasGenes = lapply(genes, function(x) unique(as.character(x$EntrezID)))
notspl = lapply(spl, function(y) lapply(y, function(z) lapply(z, function(x) geneuniverse[!(geneuniverse %in% x)])))

enrich = list()
for (h in 1:length(spl)) {
  enrich[[h]] = list()
  for (i in 1:length(spl[[h]])) {
    enrich[[h]][[i]] = list()
    for (j in 1:length(spl[[h]][[i]])) {
      enrich[[h]][[i]][[j]] = lapply(gwasGenes, function(x) {
        DE_OVERLAP = c( sum( spl[[h]][[i]][[j]] %in% x),sum(!(spl[[h]][[i]][[j]] %in% x)))
        NOT_DE_OVERLAP= c(sum(notspl[[h]][[i]][[j]] %in% x), sum(!(notspl[[h]][[i]][[j]] %in% x)))
        enrich_table = cbind(DE_OVERLAP, NOT_DE_OVERLAP)
        res = fisher.test(enrich_table)
        dat=c(res$p.value, res$estimate)
        names(dat) <- c("P.Value","Odds.Ratio")
        return(dat)
      })
    }
    names(enrich[[h]][[i]]) = names(spl[[h]][[i]])
  }
  names(enrich[[h]]) = names(spl[[h]])
}
names(enrich) = names(spl)


enrich = lapply(enrich, function(x) lapply(x, function(z) 
  lapply(z, function(h) do.call(rbind, Map(cbind,lapply(h, function(y) data.frame(OR = y[names(y)=="Odds.Ratio"], pval = y[names(y)=="P.Value"])), GeneSet = as.list(names(h)))))))
enrich = lapply(enrich, function(x) lapply(x, function(z) do.call(rbind, Map(cbind, z, Direction = as.list(names(z))))))
enrich = lapply(enrich, function(x) do.call(rbind, Map(cbind, x, Context = as.list(names(x)))))
enrich = do.call(rbind, Map(cbind, enrich, Feature = as.list(names(enrich))))

DMRenrich = rbind(DMRenrich, enrich)
DMRenrich$FDR = p.adjust(DMRenrich$pval, method="fdr")
rownames(DMRenrich) = NULL
DMRenrich = data.frame(Feature = DMRenrich$Feature, Context = DMRenrich$Context, Direction = DMRenrich$Direction, GeneSet = DMRenrich$GeneSet,
                       P.value = DMRenrich$pval, Odds.Ratio = DMRenrich$OR, FDR = DMRenrich$FDR)

write.csv(DMRenrich,file="/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/Birnbaum_geneSet_enrichment_splicingResults.csv",quote=F)

DMRenrich[which(DMRenrich$FDR<=0.05),]
#    Feature Context Direction      GeneSet      P.value Odds.Ratio          FDR
#5      Exon     CpH       Pos          NDD 2.089300e-03   8.306668 1.917240e-02
#20     Exon     CpG       Pos ASD DATABASE 5.690496e-04   3.764205 8.877174e-03
#47      PSI     CpH       Neg ASD DATABASE 1.765525e-05   3.752510 5.874233e-04
#48      PSI     CpH       Neg    BPAD GWAS 2.612524e-03   3.603506 2.145020e-02
#54      PSI     CpH       Neg      SCZ SNV 2.219433e-03   2.928525 1.923509e-02
#56      PSI     CpG       Pos ASD DATABASE 7.609797e-09   3.082437 5.935641e-07
#59      PSI     CpG       Pos          NDD 3.230387e-03   4.219748 2.399716e-02
#63      PSI     CpG       Pos      SCZ SNV 2.259320e-05   2.525660 5.874233e-04
#65      PSI     CpG       Neg ASD DATABASE 2.008817e-05   2.583172 5.874233e-04
#68      PSI     CpG       Neg          NDD 4.800595e-03   4.386468 2.880357e-02
#72      PSI     CpG       Neg      SCZ SNV 3.569523e-03   2.116158 2.423865e-02
#74     Gene     CpH       Pos ASD DATABASE 1.165644e-03   2.484406 1.236333e-02
#83     Gene     CpH       Neg ASD DATABASE 4.974462e-07   2.787712 2.586720e-05
#92     Gene     CpG       Pos ASD DATABASE 1.622079e-16   4.330241 2.530444e-14
#95     Gene     CpG       Pos          NDD 4.489586e-03   3.959069 2.801502e-02
#99     Gene     CpG       Pos      SCZ SNV 2.140034e-04   2.273289 4.769218e-03
#101    Gene     CpG       Neg ASD DATABASE 3.425665e-04   2.057475 5.937820e-03
#115    Exon     CpH       Neg          T2D 1.539480e-03   3.654442 1.500993e-02
#127     PSI     CpH       Pos          T2D 1.188782e-03   5.505959 1.236333e-02
#131     PSI     CpH       Neg          T2D 3.573647e-03   4.386497 2.423865e-02
#134     PSI     CpG       Pos         Park 4.404323e-03   2.775663 2.801502e-02
#135     PSI     CpG       Pos          T2D 1.097990e-03   3.010241 1.236333e-02
#136     PSI     CpG       Pos          SCZ 2.925908e-03   1.732568 2.282208e-02
#141    Gene     CpH       Pos          Alz 7.912097e-04   4.520591 1.077322e-02
#147    Gene     CpH       Neg          T2D 2.572376e-04   3.384374 5.016133e-03
#152    Gene     CpG       Pos          SCZ 8.287094e-04   1.828156 1.077322e-02



### What are some examples?

examples = list()
for (h in 1:length(spl)) {
  examples[[h]] = list()
  for (i in 1:length(spl[[h]])) {
    examples[[h]][[i]] = list()
    for (j in 1:length(spl[[h]][[i]])) {
      examples[[h]][[i]][[j]] = lapply(splitSets, function(x) spl[[h]][[i]][[j]][which(spl[[h]][[i]][[j]] %in% x$EntrezGene.ID)])
    }
    names(examples[[h]][[i]]) = names(spl[[h]][[i]])
  }
  names(examples[[h]]) = names(spl[[h]])
}
names(examples) = names(spl)


psi = lapply(examples$PSI, function(x) lapply(x, function(z) lapply(z, function(y) unique(geneMap[match(y, geneMap$EntrezID),"gencodeID"]))))

res = list(Exon = list(CpH = list(Pos = lapply(examples$Exon$CpH$Pos, function(x) emres$nonCpG$eqtl[which(emres$nonCpG$eqtl$statistic>0 & rowData(emres$nonCpG$expr)$EntrezID %in% x),]),
                                  Neg = lapply(examples$Exon$CpH$Neg, function(x) emres$nonCpG$eqtl[which(emres$nonCpG$eqtl$statistic<0 & rowData(emres$nonCpG$expr)$EntrezID %in% x),])),
                       CpG = list(Pos = lapply(examples$Exon$CpG$Pos, function(x) emres$CpG$eqtl[which(emres$CpG$eqtl$statistic>0 & rowData(emres$CpG$expr)$EntrezID %in% x),]),
                                  Neg = lapply(examples$Exon$CpG$Neg, function(x) emres$CpG$eqtl[which(emres$CpG$eqtl$statistic<0 & rowData(emres$CpG$expr)$EntrezID %in% x),]))),
           PSI = list(CpH = list(Pos = lapply(psi$CpH$Pos, function(x) pmres$nonCpG$eqtl[which(pmres$nonCpG$eqtl$statistic>0 & pmres$nonCpG$eqtl$gene %in% x),]),
                                 Neg = lapply(psi$CpH$Neg, function(x) pmres$nonCpG$eqtl[which(pmres$nonCpG$eqtl$statistic<0 & pmres$nonCpG$eqtl$gene %in% x),])),
                      CpG = list(Pos = lapply(psi$CpG$Pos, function(x) pmres$CpG$eqtl[which(pmres$CpG$eqtl$statistic>0 & pmres$CpG$eqtl$gene %in% x),]),
                                 Neg = lapply(psi$CpG$Neg, function(x) pmres$CpG$eqtl[which(pmres$CpG$eqtl$statistic<0 & pmres$CpG$eqtl$gene %in% x),]))),
           Gene = list(CpH = list(Pos = lapply(examples$Gene$CpH$Pos, function(x) gmres$nonCpG$eqtl[which(gmres$nonCpG$eqtl$statistic>0 & rowData(gmres$nonCpG$expr)$EntrezID %in% x),]),
                                  Neg = lapply(examples$Gene$CpH$Neg, function(x) gmres$nonCpG$eqtl[which(gmres$nonCpG$eqtl$statistic<0 & rowData(gmres$nonCpG$expr)$EntrezID %in% x),])),
                       CpG = list(Pos = lapply(examples$Gene$CpG$Pos, function(x) gmres$CpG$eqtl[which(gmres$CpG$eqtl$statistic>0 & rowData(gmres$CpG$expr)$EntrezID %in% x),]),
                                  Neg = lapply(examples$Gene$CpG$Neg, function(x) gmres$CpG$eqtl[which(gmres$CpG$eqtl$statistic<0 & rowData(gmres$CpG$expr)$EntrezID %in% x),]))))
res = lapply(res, function(x) lapply(x, function(z) lapply(z, function(h) do.call(rbind, lapply(h, as.data.frame)))))
res = lapply(res, function(x) lapply(x, function(z) Map(cbind, z, GeneSets = lapply(z, function(a) gsub('\\..*', '', rownames(a))))))
res = lapply(res, function(x) lapply(x, function(z) do.call(rbind, Map(cbind, z, Direction = as.list(names(z))))))
res = lapply(res, function(x) do.call(rbind, Map(cbind, x, Context = as.list(names(x)))))
res = do.call(rbind, Map(cbind, res, Feature = as.list(names(res))))
rownames(res) = NULL


save(res, file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/disease_associated_splicing_events.rda")

x = res[res$Feature=="PSI" & res$Context=="CpH" & res$FDR<=0.01,]
x$Symbol = geneMap[match(x$gene, geneMap$gencodeID),"Symbol"]
x = x[-which(x$gene %in% res[res$Feature=="Gene","gene"]),]
head(x[order(x$FDR),])

# For PSI, choose an AE event in DOCK1, a CpH negative direction SCZ SNV associate event, seems to overlap a CTCF binding site in UCSC

splitSets$"SCZ SNV"[splitSets$"SCZ SNV"$Gene.Symbol=="DOCK1",]
#Gene.Set Gene.Symbol EntrezGene.ID                Description
#430  SCZ SNV       DOCK1          1793 dedicator of cytokinesis 1
#Location Genetic.Evidence Beta.Coefficient.("Fetal.Effect")
#430 10q26.13-q26.3 Exome sequencing                         -1.362662
#t-statistic      p-value Adjusted.R2  Probe.ID Ref.Seq.Accession
#430   -14.11782 2.976505e-33    0.385762 hHC001930         NM_001380


x = res[res$Feature=="Exon" & res$Context=="CpH" & res$FDR<=0.01,]
Exon = unique(rbind(cbind(rbind(cbind(rowData(emres$nonCpG$expr)[which(emres$nonCpG$eqtl$statistic>0),], Direction = "Pos"),
                         cbind(rowData(emres$nonCpG$expr)[which(emres$nonCpG$eqtl$statistic<0),], Direction = "Neg")), Context = "CpH"),
             cbind(rbind(cbind(rowData(emres$CpG$expr)[which(emres$CpG$eqtl$statistic>0),], Direction = "Pos"),
                         cbind(rowData(emres$CpG$expr)[which(emres$CpG$eqtl$statistic<0),], Direction = "Neg")), Context = "CpG")))
x$gencode = Exon[match(as.character(x$gene), Exon$exon_libdID),"gencodeID"]
x$Symbol = geneMap[match(x$gencode, geneMap$gencodeID),"Symbol"]
x = x[-which(x$gencode %in% res[res$Feature=="Gene","gene"]),]

head(x[order(x$FDR),])
# For PSI, choose an AE event in DOCK1, a CpH negative direction SCZ SNV associate event, seems to overlap a CTCF binding site in UCSC

splitSets$"ASD DATABASE"[splitSets$"ASD DATABASE"$Gene.Symbol=="TTN",]
#Gene.Set Gene.Symbol EntrezGene.ID Description Location
#763 ASD DATABASE         TTN          7273       titin     2q31
#Genetic.Evidence Beta.Coefficient.("Fetal.Effect") t-statistic
#763 Rare Single Gene variant                         0.4744086     3.87321
#p-value Adjusted.R2  Probe.ID Ref.Seq.Accession
#763 0.0001390042 -0.06348465 hHA035073         NM_133378


## Do the Cs fall in a DMR?

load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/DMR/DMR_objects.rda")
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/bumphunting/rda/limma_Neuron_CpGs_minCov_3_ageInfo_dmrs.Rdata")

dmrs = split(dmrs, dmrs$k6cluster_label)
names(dmrs) = c("Group 1 (G-N+)","Group 2 (G0N+)","Group 3 (G0N-)","Group 4 (G+N0)","Group 5 (G+N-)","Group 6 (G-N0)")
oo = lapply(dmrs, function(x) findOverlaps(x, makeGRangesFromDataFrame(DMR$Interaction)))
dmrs = lapply(oo, function(x) DMR$Interaction[subjectHits(x),])
CT = split(DMR$CellType, DMR$CellType$Dir)
names(CT) = c("Hypomethylated in Neurons", "Hypomethylated in Glia")
DMRgr = lapply(c(CT, DMR[names(DMR)!="CellType"], dmrs), function(x) makeGRangesFromDataFrame(x[which(x$fwer<=0.05),], keep.extra.columns = T))
DMRgr = lapply(DMRgr, function(x) unique(granges(x)))

spl = list(Exon = list(CpH = list(Pos = granges(emres$nonCpG$meth)[which(emres$nonCpG$eqtl$statistic>0)],
                                  Neg = granges(emres$nonCpG$meth)[which(emres$nonCpG$eqtl$statistic<0)]),
                       CpG = list(Pos = granges(emres$CpG$meth)[which(emres$CpG$eqtl$statistic>0)],
                                  Neg = granges(emres$CpG$meth)[which(emres$CpG$eqtl$statistic<0)])),
           PSI = list(CpH = list(Pos = granges(pmres$nonCpG$meth)[which(pmres$nonCpG$eqtl$statistic>0)],
                                 Neg = granges(pmres$nonCpG$meth)[which(pmres$nonCpG$eqtl$statistic<0)]),
                      CpG = list(Pos = granges(pmres$CpG$meth)[which(pmres$CpG$eqtl$statistic>0)],
                                 Neg = granges(pmres$CpG$meth)[which(pmres$CpG$eqtl$statistic<0)])),
           Gene = list(CpH = list(Pos = granges(gmres$nonCpG$meth)[which(gmres$nonCpG$eqtl$statistic>0)],
                                  Neg = granges(gmres$nonCpG$meth)[which(gmres$nonCpG$eqtl$statistic<0)]),
                       CpG = list(Pos = granges(gmres$CpG$meth)[which(gmres$CpG$eqtl$statistic>0)],
                                  Neg = granges(gmres$CpG$meth)[which(gmres$CpG$eqtl$statistic<0)])))
spl = lapply(spl, function(x) lapply(x, function(y) lapply(y, function(z) unique(granges(z)))))


oo = lapply(DMRgr, function(x) lapply(spl, function(f) lapply(f, function(con) lapply(con, function(d) findOverlaps(d,x)))))

res = lapply(oo, function(x) lapply(x, function(y) lapply(y, function(a) lapply(a, function(z)
  data.frame(querylength = queryLength(z), subjectlength = subjectLength(z), 
             queryHits = length(unique(queryHits(z))), subjectHits = (length(unique(subjectHits(z)))))))))

res = do.call(rbind, Map(cbind, DMR = as.list(names(res)), 
                         lapply(res, function(x) do.call(rbind, Map(cbind, Feature = as.list(names(x)),
                                                                    lapply(x, function(y) do.call(rbind, Map(cbind, Context = as.list(names(y)),
                                                                                                             lapply(y, function(a) do.call(rbind, Map(cbind, Direction = as.list(names(a)), a)))))))))))
res$queryPerc = res$queryHits/res$querylength*100
res$subjectPerc = res$subjectHits/res$subjectlength*100

write.csv(res, quote = F, file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/splicingCs_DMR_overlap.csv")

# collapse 2 levels since it's so little overlap

spl = list(Exon = list(CpH = granges(emres$nonCpG$meth), CpG = granges(emres$CpG$meth)),
           PSI = list(CpH = granges(pmres$nonCpG$meth), CpG = granges(pmres$CpG$meth)),
           Gene = list(CpH = granges(gmres$nonCpG$meth), CpG = granges(gmres$CpG$meth)))
spl = lapply(spl, function(x) lapply(x, function(z) unique(granges(z))))


oo = lapply(spl, function(f) lapply(f, function(d) findOverlaps(d,reduce(c(DMRgr$"Hypomethylated in Neurons",DMRgr$"Hypomethylated in Glia",DMRgr$"Age", DMRgr$"Interaction")))))

res = lapply(oo, function(x) lapply(x, function(z)
  data.frame(querylength = queryLength(z), subjectlength = subjectLength(z), 
             queryHits = length(unique(queryHits(z))), subjectHits = (length(unique(subjectHits(z)))))))

res = do.call(rbind, Map(cbind, Feature = as.list(names(res)), lapply(res, function(x) do.call(rbind, Map(cbind, Context = as.list(names(x)), x)))))
res$queryPerc = res$queryHits/res$querylength*100
res$subjectPerc = res$subjectHits/res$subjectlength*100
range(res$queryPerc)




### Where is the C?

emres$CpG$eqtls$gencode = Exon[match(as.character(emres$CpG$eqtls$gene), Exon$exon_libdID),"gencodeID"]
x = emres$CpG$eqtls[-which(emres$CpG$eqtls$gencode %in% gmres$CpG$eqtls$gene),]
x[order(x$FDR),]
row150207
granges(emres$CpG$meth)[which(!emres$CpG$eqtls$gencode %in% gmres$CpG$eqtls$gene & emres$CpG$eqtls$FDR<=0.01)]
rowRanges(emres$CpG$expr)[which(!emres$CpG$eqtls$gencode %in% gmres$CpG$eqtls$gene & emres$CpG$eqtls$FDR<=0.01)]
emres$CpG$eqtls[which(!emres$CpG$eqtls$gencode %in% gmres$CpG$eqtls$gene)][emres$CpG$eqtls$FDR<=0.01]
emres$CpG$eqtls[emres$CpG$eqtls$gene=="e224400"]
chr3:52824923-52824923
chr3:52824923-52824923

chr22:37492219-37492219
chr22:37491946-37492130

chr16:57993770-57993770
chr16:57994404-57994434

emres$nonCpG$eqtls$gencode = Exon[match(as.character(emres$nonCpG$eqtls$gene), Exon$exon_libdID),"gencodeID"]
x = emres$nonCpG$eqtls[-which(emres$nonCpG$eqtls$gencode %in% gmres$nonCpG$eqtls$gene),]
x[order(x$FDR),]
row150207
granges(emres$nonCpG$meth)[which(!emres$nonCpG$eqtls$gencode %in% gmres$nonCpG$eqtls$gene & emres$nonCpG$eqtls$FDR<=0.01)]
rowRanges(emres$nonCpG$expr)[which(!emres$nonCpG$eqtls$gencode %in% gmres$nonCpG$eqtls$gene & emres$nonCpG$eqtls$FDR<=0.01)]
emres$nonCpG$eqtls[which(!emres$nonCpG$eqtls$gencode %in% gmres$nonCpG$eqtls$gene)][emres$nonCpG$eqtls$FDR<=0.01]
emres$nonCpG$eqtls[emres$nonCpG$eqtls$gene=="e224400"]

chr21:45879814-45882078
chr21:45879472-45879472

chr17:38632080-38633981
chr17:38631168-38631168

chr1:24416028-24416155
chr1:24416434-24416434

#They are all over the place


## load compare cluster outputs

load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/meth_vs_expr/rda/meqtl_venn_go_exon_using_near.Rdata")
library(clusterProfiler)

bpexon = go_cluster_comp[["BP"]]
mfexon = go_cluster_comp[["MF"]]
ccexon = go_cluster_comp[["CC"]]
mfexon = simplify(mfexon)
mfexon = data.frame(mfexon)
mfexon$Description = gsub(",", ";", mfexon$Description, fixed=TRUE)

write.csv(mfexon, quote=FALSE, row.names=FALSE, 
          file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/meth_vs_expr/rda/mfexon_table.csv")

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/meth_vs_expr/rda/mfexon_plot_full.pdf", height = 30, width = 21)
dotplot(mfexon, showCategory = 20, title= "Molecular Function GO: Exons")
dev.off()

pdf("./wgbs_development/Third submission/Figures/Figure S14 - molecular function gene ontology enrichment for spliced exons/mfexon_plot_top20.pdf", 
    height = 11, width = 21)
dotplot(mfexon, showCategory = 20, title= "Molecular Function GO: Exons")
dev.off()


bpexon@compareClusterResult = bpexon@compareClusterResult[which(bpexon@compareClusterResult$Cluster!="CpGmarg" & 
                                                                  bpexon@compareClusterResult$Cluster!="CpG:CpGmarg"),]
mfexon@compareClusterResult = mfexon@compareClusterResult[which(mfexon@compareClusterResult$Cluster!="CpGmarg" & 
                                                                  mfexon@compareClusterResult$Cluster!="CpG:CpGmarg"),]
ccexon@compareClusterResult = ccexon@compareClusterResult[which(ccexon@compareClusterResult$Cluster!="CpGmarg" & 
                                                                  ccexon@compareClusterResult$Cluster!="CpG:CpGmarg"),]

bpexon@compareClusterResult = bpexon@compareClusterResult[which(bpexon@compareClusterResult$p.adjust<=0.05),]
mfexon@compareClusterResult = mfexon@compareClusterResult[which(mfexon@compareClusterResult$p.adjust<=0.05),]
ccexon@compareClusterResult = ccexon@compareClusterResult[which(ccexon@compareClusterResult$p.adjust<=0.05),]

mfexon@compareClusterResult = mfexon@compareClusterResult[which(mfexon@compareClusterResult$Description %in% 
                                                                  c("extracellular matrix structural constituent","voltage−gated sodium channel activity",
                                                                    "voltage−gated ion channel activity involved in regulation of postsynaptic membrane potential",
                                                                    "secondary active transmembrane transporter activity","active transmembrane transporter activity",
                                                                    "neurotransmitter transporter activity","neurotransmitter:sodium symporter activity",
                                                                    "sodium ion transmembrane transporter activity","protein tyrosine kinase activity",
                                                                    "cation channel activity","active ion transmembrane transporter activity",
                                                                    "substrate−specific channel activity")),]
  
pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/meth_vs_expr/rda/mfexon_plot.pdf", height = 6, width = 14.5)
plot(mfexon, colorBy= "p.adjust",  showCategory = 8, title= "Molecular Function GO: Exons")
dev.off()


pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/meth_vs_expr/rda/bpexon_plot.pdf", height = 60, width = 20)
plot(bpexon, colorBy= "p.adjust",  showCategory = 700, title= paste0("Biological Process GO: Exons"))
plot(mfexon, colorBy= "p.adjust",  showCategory = 700, title= paste0("Molecular Function GO: Exons"))
plot(ccexon, colorBy= "p.adjust",  showCategory = 700, title= paste0("Cellular Compartment GO: Exons"))
dev.off()


# exon CpG and CpH both enriched for CC cell projection membrane


load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/meth_vs_expr/rda/meqtl_venn_go_psi_using_near.Rdata", verbose=T)

bppsi = go_cluster_comp[["BP"]]
mfpsi = go_cluster_comp[["MF"]]
ccpsi = go_cluster_comp[["CC"]]

bppsi@compareClusterResult = bppsi@compareClusterResult[which(bppsi@compareClusterResult$p.adjust<=0.05),]
mfpsi@compareClusterResult = mfpsi@compareClusterResult[which(mfpsi@compareClusterResult$p.adjust<=0.05),]
ccpsi@compareClusterResult = ccpsi@compareClusterResult[which(ccpsi@compareClusterResult$p.adjust<=0.05),]

ccpsi@compareClusterResult = ccpsi@compareClusterResult[which(ccpsi@compareClusterResult$Cluster!="CpGmarg"),]

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/meth_vs_expr/rda/ccpsi_plot.pdf", height = 6, width = 7)
plot(ccpsi, colorBy= "p.adjust",  showCategory = 8, title= paste0("Cellular Compartment GO: PSI"))
dev.off()


bppsi@compareClusterResult = bppsi@compareClusterResult[which(bppsi@compareClusterResult$Cluster!="CpGmarg" & 
                                                                  bppsi@compareClusterResult$Cluster!="CpG:CpGmarg"),]
mfpsi@compareClusterResult = mfpsi@compareClusterResult[which(mfpsi@compareClusterResult$Cluster!="CpGmarg" & 
                                                                  mfpsi@compareClusterResult$Cluster!="CpG:CpGmarg"),]


pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/meth_vs_expr/rda/bppsi_plot.pdf", height = 60, width = 20)
plot(bppsi, colorBy= "p.adjust",  showCategory = 700, title= paste0("Biological Process GO: psi"))
plot(mfpsi, colorBy= "p.adjust",  showCategory = 700, title= paste0("Molecular Function GO: psi"))
plot(ccpsi, colorBy= "p.adjust",  showCategory = 700, title= paste0("Cellular Compartment GO: psi"))
dev.off()






