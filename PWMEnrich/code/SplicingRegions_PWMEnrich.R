library(PWMEnrich)
library(GenomicRanges)
library(PWMEnrich.Hsapiens.background)
library(BSgenome.Hsapiens.UCSC.hg19)
library(clusterProfiler)
library(org.Hs.eg.db)
library('bsseq')
library('SGSeq')


## explore exon data
load('/dcl01/lieber/ajaffe/lab/brain-epigenomics/meth_vs_expr/rda/meqtl_mres_exon_using_near_meth11_proteincoding.Rdata', verbose = TRUE)
class(mres)
names(mres)
names(mres$nonCpG)
mres$nonCpG$expr
mres$nonCpG$eqtls
load('/dcl01/lieber/ajaffe/lab/brain-epigenomics/meth_vs_expr/rda/meqtl_c_by_gene_exon_using_near.Rdata', verbose = TRUE)
class(c_by_gene)
names(c_by_gene)
class(c_by_gene$nonCpG)
head(c_by_gene$nonCpG)
load('/dcl01/lieber/ajaffe/lab/brain-epigenomics/meth_vs_expr/rda/meqtl_data_venn_summ_exon_using_near.Rdata', verbose = TRUE)
class(data_venn_summ)
head(data_venn_summ)
load('/dcl01/lieber/ajaffe/lab/brain-epigenomics/meth_vs_expr/rda/meqtl_age_coef_exon_using_near.Rdata', verbose = TRUE)
class(age_coef)
names(age_coef)
class(age_coef$nonCpG)
head(age_coef$nonCpG)
load('/dcl01/lieber/ajaffe/lab/brain-epigenomics/meth_vs_expr/rda/meqtl_data_by_venn_exon_using_near.Rdata', verbose = TRUE)
class(data_by_venn)
dim(data_by_venn)
head(data_by_venn)


## How many of the mC's fall within the exon with which they associate?

gr = list(CpG.pos = rowRanges(mres$CpG$expr)[which(mres$CpG$eqtls$statistic>0),], nonCpG.pos = rowRanges(mres$nonCpG$expr)[which(mres$nonCpG$eqtls$statistic>0),],
          CpG.neg = rowRanges(mres$CpG$expr)[which(mres$CpG$eqtls$statistic<0),], nonCpG.neg = rowRanges(mres$nonCpG$expr)[which(mres$nonCpG$eqtls$statistic<0),])
Cgr = list(CpG.pos = granges(mres$CpG$meth)[which(mres$CpG$eqtls$statistic>0),], nonCpG.pos = granges(mres$nonCpG$meth)[which(mres$nonCpG$eqtls$statistic>0),],
           CpG.neg = granges(mres$CpG$meth)[which(mres$CpG$eqtls$statistic<0),], nonCpG.neg = granges(mres$nonCpG$meth)[which(mres$nonCpG$eqtls$statistic<0),])

oo = mapply(function(exon, mC) findOverlaps(exon,mC), gr, Cgr, SIMPLIFY = F)
oo = lapply(oo, as.data.frame)
do.call(rbind, mapply(function(dat, len) data.frame(total = len, SameExon = length(which(dat$queryHits==dat$subjectHits)), 
                                                    Percent = round(length(which(dat$queryHits==dat$subjectHits))/len*100,1)),oo,elementNROWS(Cgr), SIMPLIFY = F))
#           total SameExon Percent
#CpG.pos     1328      447    33.7
#nonCpG.pos  1288       78     6.1
#CpG.neg     9505     2311    24.3
#nonCpG.neg  7874      605     7.7


## Get sequences from different annotations in the correct format

seq = lapply(gr, function(x) getSeq(Hsapiens, x[which(width(x)>30)]))

### Run PWMEnrich

useBigMemoryPWMEnrich(TRUE)
registerCoresPWMEnrich(6)

# load the pre-compiled lognormal background computed using promoters

data(PWMLogn.hg19.MotifDb.Hsap)
splice.exon = lapply(seq, function(x) motifEnrichment(x, PWMLogn.hg19.MotifDb.Hsap, verbose=F))

## Test for Differential TF binding

pos = splice.exon[grep("pos", names(splice.exon))]
neg = splice.exon[grep("neg", names(splice.exon))]
seqpos = seq[grep("pos", names(seq))]
seqneg = seq[grep("neg", names(seq))]

exonDiff = list()
for (i in 1:length(pos)) {
  exonDiff[[i]] = motifDiffEnrichment(sequences1 = seqpos[[i]], sequences2 = seqneg[[i]],
                                    res1 = pos[[i]], res2 = neg[[i]], PWMLogn.hg19.MotifDb.Hsap, verbose=FALSE)
}
names(exonDiff) = c("CpG", "nonCpG")

exonDiff = c(exonDiff, list(CpG.nonCpG.pos = motifDiffEnrichment(sequences1 = seqpos$CpG.pos, sequences2 = seqpos$nonCpG.pos,
                                                                 res1 = pos$CpG.pos, res2 = pos$nonCpG.pos, PWMLogn.hg19.MotifDb.Hsap, verbose=FALSE),
                            CpG.nonCpG.neg = motifDiffEnrichment(sequences1 = seqneg$CpG.neg, sequences2 = seqneg$nonCpG.neg,
                                                                 res1 = neg$CpG.neg, res2 = neg$nonCpG.neg, PWMLogn.hg19.MotifDb.Hsap, verbose=FALSE)))

save(exonDiff, splice.exon, file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/PWMEnrich/spliceRegions_PWMEnrich_objects.rda")


## Save lists of differentially enriched TFs

load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/PWMEnrich/filtered_TF_list_expressionInfo.rda")

group = lapply(splice.exon, groupReport)
group = lapply(group, as.data.frame)
group = lapply(group, function(x) x[which(x$target %in% names(PostnataltargettogeneID)),])
group = Map(cbind, group, padj = lapply(group, function(x) p.adjust(x$p.value, method = "fdr")))
group = lapply(group, function(x) x[order(x$id),])

exonDiff = lapply(exonDiff, function(x) x$group.bg[which(names(x$group.bg) %in% group$CpG.pos$id)])
do.call(rbind, lapply(exonDiff, function(x) quantile(x, na.rm = T)))

exonDiffQ25 = mapply(function(tf, q) tf[which(tf<=q[2] & abs(tf)>=2)], exonDiff, lapply(exonDiff, function(x) quantile(x, na.rm = T)))
exonDiffQ75 = mapply(function(tf, q) tf[which(tf>=q[4] & abs(tf)>=2)], exonDiff, lapply(exonDiff, function(x) quantile(x, na.rm = T)))

comps = list(CpG = c("CpG.pos","CpG.neg"), nonCpG = c("nonCpG.pos", "nonCpG.neg"),
             CpG.nonCpG.pos = c("CpG.pos","nonCpG.pos"),
             CpG.nonCpG.neg = c("CpG.neg","nonCpG.neg"))
x = list()
for (i in 1:length(comps)) {
  x[[i]] = cbind(group[[comps[[i]][1]]], group[[comps[[i]][2]]][match(group[[comps[[i]][1]]][,"id"], group[[comps[[i]][2]]][,"id"]),])
  x[[i]] = x[[i]][,-grep("rank", colnames(x[[i]]))]
  colnames(x[[i]]) = gsub("padj", "FDR", colnames(x[[i]]))
}

exonDiffQ25 = mapply(function(d,x) data.frame(Quantile = "Q25","Difference Statistic" = d, x[match(names(d), x[,2]),]), exonDiffQ25, x, SIMPLIFY = F) 
exonDiffQ75 = mapply(function(d,x) data.frame(Quantile = "Q75","Difference Statistic" = d, x[match(names(d), x[,2]),]), exonDiffQ75, x, SIMPLIFY = F) 

exonDiffQ25 = lapply(exonDiffQ25, function(x) x[which(x[,14]<=0.01),])
exonDiffQ75 = lapply(exonDiffQ75, function(x) x[which(x[,8]<=0.01),])

for (i in 1:length(exonDiffQ25)) { 
  write.csv(rbind(exonDiffQ25[[i]], exonDiffQ75[[i]]), quote = F, row.names = F,
            file = paste0("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/PWMEnrich/differentially_enriched_TFs_", names(exonDiffQ75)[i], ".csv"))  }


## GO enrichment of these factors

exonDiffQ25 = lapply(exonDiffQ25, function(x) PostnataltargettogeneID[match(x[,3],names(PostnataltargettogeneID))])
exonDiffQ75 = lapply(exonDiffQ75, function(x) PostnataltargettogeneID[match(x[,3],names(PostnataltargettogeneID))])
exonDiffQ25 = lapply(exonDiffQ25, function(x) geneMap[which(geneMap$gencodeID %in% x),])
exonDiffQ75 = lapply(exonDiffQ75, function(x) geneMap[which(geneMap$gencodeID %in% x),])

entrez = mapply(function(x,y) list(Q25 = unique(na.omit(as.character(x$EntrezID))), Q75 = unique(na.omit(as.character(y$EntrezID)))),
                exonDiffQ25, exonDiffQ75, SIMPLIFY = F) 


## Assess enriched terms limiting the gene universe to the terms associated with the master list of TFs

GeneUniverse = unique(na.omit(as.character(geneMap[which(geneMap$gencodeID %in% PostnataltargettogeneID), "EntrezID"])))
length(GeneUniverse) # 616

# Find enriched pathways and processes

keggList = lapply(unlist(entrez, recursive=F), function(x) enrichKEGG(x, organism="human", universe= GeneUniverse, minGSSize=5, 
                                                                      pAdjustMethod="BH", qvalueCutoff=1))
goList_MF = lapply(unlist(entrez, recursive=F), function(x) enrichGO(x, ont = "MF", OrgDb = org.Hs.eg.db, universe= GeneUniverse, 
                                                                     minGSSize=5, pAdjustMethod="BH", qvalueCutoff=1))
goList_BP = lapply(unlist(entrez, recursive=F), function(x) enrichGO(x, ont = "BP", OrgDb = org.Hs.eg.db,universe= GeneUniverse, 
                                                                     minGSSize=5, pAdjustMethod="BH", qvalueCutoff=1))
goList_CC = lapply(unlist(entrez, recursive=F), function(x) enrichGO(x, ont = "CC", OrgDb = org.Hs.eg.db, universe= GeneUniverse, 
                                                                     minGSSize=5, pAdjustMethod="BH", qvalueCutoff=1))
goList_DO = lapply(unlist(entrez, recursive=F), function(x) enrichDO(x, ont = "DO", universe= GeneUniverse, minGSSize=5, 
                                                                     pAdjustMethod="BH", qvalueCutoff=1, readable=TRUE))


# Compare the enriched terms

compareKegg = lapply(entrez, function(x) compareCluster(x, fun="enrichKEGG", qvalueCutoff = 0.05, pvalueCutoff = 0.05))
compareBP = lapply(entrez, function(x) compareCluster(x, fun="enrichGO", ont = "BP", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05))
compareMF = lapply(entrez, function(x) compareCluster(x, fun="enrichGO",  ont = "MF", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05))
compareCC = lapply(entrez, function(x) compareCluster(x, fun="enrichGO",  ont = "CC", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05))
compareDO = lapply(entrez, function(x) compareCluster(x, fun="enrichDO",  ont = "DO", qvalueCutoff = 0.05, pvalueCutoff = 0.05))


save(keggList, goList_MF, goList_BP, goList_CC, goList_DO, compareBP, compareMF, compareCC, compareKegg,
     file="/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/PWMEnrich/spliceRegions_GO.objects.PWMEnrich.rda")


pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/PWMEnrich/figures/KEGG_MF_CC_DO_spliceRegions_TFs.pdf", width=12,height=60)
for (i in 1:length(compareKegg)) {
  print(plot(compareKegg[[i]], colorBy= "p.adjust",  showCategory = 700, title= paste0("KEGG Pathway Enrichment: ", names(compareKegg)[i])))
  print(plot(compareMF[[i]], colorBy= "p.adjust",  showCategory = 700, title= paste0("MF Pathway Enrichment: ", names(compareMF)[i])))
  print(plot(compareCC[[i]], colorBy= "p.adjust",  showCategory = 700, title= paste0("CC Pathway Enrichment: ", names(compareCC)[i])))
  print(plot(compareBP[[i]], colorBy= "p.adjust",  showCategory = 700, title= paste0("Biological Process GO: ", names(compareBP)[i])))
  print(plot(compareDO[[i]], colorBy= "p.adjust",  showCategory = 700, title= paste0("Disease Ontology: ", names(compareDO)[i])))
}
dev.off()



## Doing some manual checking of known splice factors...

do.call(rbind, lapply(group, function(x) x[which(x$target=="MECP2"),]))
#             rank target    id raw.score p.value top.motif.prop padj
#CpG.pos    1732.5  MECP2 MECP2 0.7442695       1     0.01466049    1
#nonCpG.pos 1493.0  MECP2 MECP2 0.9724046       1     0.02680747    1
#CpG.neg    1678.0  MECP2 MECP2 0.7607640       1     0.01512892    1
#nonCpG.neg 1721.5  MECP2 MECP2 0.8783678       1     0.02998462    1

exonDiff$CpG["MECP2"] # 0
exonDiff$nonCpG["MECP2"] # 0.002241363
exonDiff$CpG.nonCpG.pos["MECP2"] # -0.002241363 
exonDiff$CpG.nonCpG.neg["MECP2"] # 0

do.call(rbind, lapply(exonDiffQ25, function(x) x[which(x[,3]=="CTCF"),]))
#        Quantile Difference.Statistic target                      id raw.score
#CpG.255      Q25        -9.897927e+12   CTCF                    CTCF 23.492773
#CpG.301      Q25        -1.265043e+08   CTCF Hsapiens-jolma2013-CTCF  4.226071
#nonCpG       Q25        -3.643336e+11   CTCF                    CTCF 69.754491
#             p.value top.motif.prop          FDR target.1
#CpG.255 2.110034e-42     0.06018519 1.479581e-41     CTCF
#CpG.301 3.192906e-34     0.10185185 1.884792e-33     CTCF
#nonCpG  3.093087e-26     0.08610885 5.186524e-25     CTCF
#                           id.1 raw.score.1     p.value.1 top.motif.prop.1
#CpG.255                    CTCF    4.222023 4.890218e-197       0.07234179
#CpG.301 Hsapiens-jolma2013-CTCF    3.514215  5.511270e-78       0.06658854
#nonCpG                     CTCF    5.761351 1.922173e-156       0.07816504
#                FDR.1
#CpG.255 3.160884e-196
#CpG.301  2.214076e-77
#nonCpG  1.555212e-155

# Greater enriched in CpG than non-CpG in both expression direction effects, but does not pass the filtering steps
exonDiff$CpG.nonCpG.pos["CTCF"] # 765565
exonDiff$CpG.nonCpG.neg["CTCF"] # 9.533595e+12


do.call(rbind,lapply(group, function(x) x[which(x$target=="HP1BP3"),]))
#             rank target     id raw.score     p.value top.motif.prop
#CpG.pos    1732.5 HP1BP3 HP1BP3  1.308901 1.000000000     0.03009259
#nonCpG.pos  717.0 HP1BP3 HP1BP3  2.093420 0.001504384     0.09017059
#CpG.neg    1678.0 HP1BP3 HP1BP3  1.499761 1.000000000     0.03878116
#nonCpG.neg 1721.5 HP1BP3 HP1BP3  1.671395 1.000000000     0.03908252
#                  padj
#CpG.pos    1.000000000
#nonCpG.pos 0.005045138
#CpG.neg    1.000000000
#nonCpG.neg 1.000000000

unlist(lapply(exonDiff, function(x) x[names(x)=="HP1BP3"]))
#           CpG.HP1BP3         nonCpG.HP1BP3 CpG.nonCpG.pos.HP1BP3 CpG.nonCpG.neg.HP1BP3
#              0.00000              19.43043             -19.43043               0.00000  

do.call(rbind,lapply(group, function(x) x[which(x$target=="CBX3"),]))
#           rank target        id raw.score      p.value top.motif.prop
#CpG.pos     813   CBX3 LOC653972 0.8142757 8.246305e-01     0.02623457
#nonCpG.pos  745   CBX3 LOC653972 0.8090082 4.322330e-03     0.04386677
#CpG.neg     882   CBX3 LOC653972 0.8405795 1.718781e-02     0.04549329
#nonCpG.neg  781   CBX3 LOC653972 0.8325335 3.490931e-11     0.04472066
#                   padj
#CpG.pos    1.000000e+00
#nonCpG.pos 1.389149e-02
#CpG.neg    4.100266e-02
#nonCpG.neg 9.803415e-11

do.call(rbind,lapply(group, function(x) x[which(x$target=="CBX7"),]))
#           rank target   id raw.score      p.value top.motif.prop         padj
#CpG.pos    1085   CBX7 CBX7  1.093078 1.000000e+00     0.04475309 1.000000e+00
#nonCpG.pos  686   CBX7 CBX7  1.241742 4.146173e-04     0.07067425 1.458092e-03
#CpG.neg     706   CBX7 CBX7  1.262174 3.823963e-20     0.06211379 1.120082e-19
#nonCpG.neg  558   CBX7 CBX7  1.329684 5.331403e-52     0.07432086 2.234940e-51

unlist(lapply(exonDiff, function(x) x[names(x)=="CBX7"]))
#           CpG.CBX7         nonCpG.CBX7 CpG.nonCpG.pos.CBX7 CpG.nonCpG.neg.CBX7 
#      -9.118516e+03       -3.713625e+06       -2.829789e+01       -3.704535e+06 



### Explore binding motif score distribution

## Get PWMs

data(MotifDb.Hsap)
pwms = c(MotifDb.Hsap[which(names(MotifDb.Hsap)=="CTCF")], MotifDb.Hsap[grep("MECP2", names(MotifDb.Hsap))],
         MotifDb.Hsap[grep("CBX7", names(MotifDb.Hsap))])
# Can't find HP1 family PWMS (CBX1,CBX5,CBX3). Shame...

pwms = MotifDb.Hsap[which(names(MotifDb.Hsap) %in% group$CpG.pos$id)] # all the motifs included in my other analyses


## Get coordinates of expressed exons

load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/brainseq_pipeline/polyA_unstranded/rse_exon_polyA_dlpfc_n41.Rdata", verbose = T)

exonMap = rowRanges(rse_exon)
exonCounts = assays(rse_exon)$counts
hompd = colData(rse_exon)
exonCounts = exonCounts[,which(colnames(exonCounts) %in% rownames(hompd[which(hompd$Age>0),]))]
exonCounts = exonCounts[which(rowSums(exonCounts)>0),]
exonMap = exonMap[names(exonMap) %in% rownames(exonCounts)]
exonMap = exonMap[which(width(x)>30)]
exonMap = exonMap[which(as.character(seqnames(exonMap)) %in% paste0("chr",c(1:22,"X","Y")))]

sampleset = sample(exonMap, 25000, replace = FALSE)


## Make new log-normal background on exons

sampleseq = getSeq(Hsapiens, sampleset)
lognBcgnd = makePWMLognBackground(sampleseq, pwms, algorithm = "human", bg.len = 100, verbose=F)


## Run motifEnrichment on associated exons

splice.exonBack = lapply(seq, function(x) motifEnrichment(x, lognBcgnd, verbose=F))

# Make table of raw scores for distribution analysis

score = pval = list()
for (i in 1:length(splice.exonBack)) {
  score[[i]] = list(vector("list", length(splice.exonBack[[i]]$sequences)))
  pval[[i]] = list(vector("list", length(splice.exonBack[[i]]$sequences)))
  for (j in 1:length(splice.exonBack[[i]]$sequences)) {
    score[[i]][[j]] = sequenceReport(splice.exonBack[[i]], seq.id=j)
    pval[[i]][[j]] = sequenceReport(splice.exonBack[[i]], seq.id=j)
  }
}
names(score) = names(pval) = names(splice.exonBack)
score = lapply(score, function(x) lapply(x, function(y) as.data.frame(y)[order(y$id),]))
score = lapply(score, function(y) do.call(cbind, lapply(y, function(z) z$raw.score)))
for (i in 1:length(score)) { rownames(score[[i]]) = group$CpG.pos$id }
for (i in 1:length(score)) { colnames(score[[i]]) = names(splice.exonBack[[i]]$sequences) }
score = lapply(score, as.data.frame)
score = Map(cbind, score, Context = lapply(as.list(names(score)), function(x) gsub('\\..*', '', x)),
            Direction = lapply(as.list(names(score)), function(x) gsub('.*\\.', '', x)),
            id = lapply(score, rownames))
cols = lapply(score, colnames)
cols = lapply(cols, function(x) data.frame(exonID = x[1:(length(x)-3)], col = paste0("col",1:length(x[1:(length(x)-3)]))))
for (i in 1:length(score)) { colnames(score[[i]]) = c(as.character(cols[[i]][,"col"]),"Context","Direction","motifID") }
score = lapply(score, reshape2::melt)

pval = lapply(pval, function(x) lapply(x, function(y) as.data.frame(y)[order(y$id),]))
pval = lapply(pval, function(y) do.call(cbind, lapply(y, function(z) p.adjust(z$p.value, method = "fdr"))))
for (i in 1:length(pval)) { rownames(pval[[i]]) = group$CpG.pos$id }
for (i in 1:length(pval)) { colnames(pval[[i]]) = names(splice.exonBack[[i]]$sequences) }
pval = lapply(pval, as.data.frame)
pval = Map(cbind, pval, Context = lapply(as.list(names(pval)), function(x) gsub('\\..*', '', x)),
            Direction = lapply(as.list(names(pval)), function(x) gsub('.*\\.', '', x)),
            id = lapply(pval, rownames))
cols = lapply(pval, colnames)
cols = lapply(cols, function(x) data.frame(exonID = x[1:(length(x)-3)], col = paste0("col",1:length(x[1:(length(x)-3)]))))
for (i in 1:length(pval)) { colnames(pval[[i]]) = c(as.character(cols[[i]][,"col"]),"Context","Direction","motifID") }
pval = lapply(pval, reshape2::melt)

# does it overlap?
gr = list(CpG.pos = rowRanges(mres$CpG$expr)[which(mres$CpG$eqtls$statistic>0),], nonCpG.pos = rowRanges(mres$nonCpG$expr)[which(mres$nonCpG$eqtls$statistic>0),],
          CpG.neg = rowRanges(mres$CpG$expr)[which(mres$CpG$eqtls$statistic<0),], nonCpG.neg = rowRanges(mres$nonCpG$expr)[which(mres$nonCpG$eqtls$statistic<0),])
Cgr = list(CpG.pos = granges(mres$CpG$meth)[which(mres$CpG$eqtls$statistic>0),], nonCpG.pos = granges(mres$nonCpG$meth)[which(mres$nonCpG$eqtls$statistic>0),],
           CpG.neg = granges(mres$CpG$meth)[which(mres$CpG$eqtls$statistic<0),], nonCpG.neg = granges(mres$nonCpG$meth)[which(mres$nonCpG$eqtls$statistic<0),])

oo = mapply(function(exon, mC) findOverlaps(exon,mC), lapply(gr, function(x) x[width(x)>30]),
            mapply(function(c,g) c[width(g)>30], Cgr, gr, SIMPLIFY = F), SIMPLIFY = F)
overlapping = lapply(oo, function(x) x[which(queryHits(x)==subjectHits(x))])

colOverlaps = mapply(function(dat, ov) as.character(dat[queryHits(ov),"col"]), cols, overlapping, SIMPLIFY = F)

x = Map(cbind, score, C_in_exon = mapply(function(dat, ov) ifelse(dat$variable %in% ov, "Yes", "No"),
                                             score, colOverlaps, SIMPLIFY = F),
            Target = lapply(score, function(x) group$CpG.pos[match(x$motifID, group$CpG.pos$id), "target"]),
            exonID = mapply(function(x, co) co[match(x$variable, co$col), "exonID"], score, cols, SIMPLIFY = F))
score = do.call(rbind, x)

x = Map(cbind, pval, C_in_exon = mapply(function(dat, ov) ifelse(dat$variable %in% ov, "Yes", "No"),
                                         pval, colOverlaps, SIMPLIFY = F),
        Target = lapply(pval, function(x) group$CpG.pos[match(x$motifID, group$CpG.pos$id), "target"]),
        exonID = mapply(function(x, co) co[match(x$variable, co$col), "exonID"], pval, cols, SIMPLIFY = F))
pval = do.call(rbind, x)

score$Direction = ifelse(score$Direction=="pos", "Positive Association", "Negative Association")
score$C_in_exon = ifelse(score$C_in_exon=="Yes", "C in Exon", "C not in Exon")
pval$Direction = ifelse(pval$Direction=="pos", "Positive Association", "Negative Association")
pval$C_in_exon = ifelse(pval$C_in_exon=="Yes", "C in Exon", "C not in Exon")

save(score, splice.exonBack, lognBcgnd,
     file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/PWMEnrich/spliceRegions_PWMEnrich_objects_exonAdjustedBackgrounds.rda")


## Make a log-normal background for the canonical splice sequences

# PFMs from https://www.ncbi.nlm.nih.gov/pubmed/9536098?dopt=Abstract
splsitesPFM = list(highGC.3SS = rbind("A" = c(10,8,7,8,6,6,4,8,8,7,6,19,2,100,0,21,19),
                                      "C" = c(41,42,41,40,38,43,42,46,49,54,45,38,82,0,0,13,21),
                                      "G" = c(15,14,14,13,13,12,13,14,10,8,8,26,0,0,100,58,29),
                                      "T" = c(34,36,38,39,43,39,41,32,33,31,41,17,16,0,0,8,31)),
                   lowGC.3SS = rbind("A" = c(15,14,13,11,10,10,11,12,13,11,10,26,7,100,0,26,24),
                                     "C" = c(24,21,20,22,21,22,25,28,28,25,22,25,55,0,0,11,15),
                                     "G" = c(10,12,10,9,10,9,10,10,8,5,5,15,1,0,100,50,20),
                                     "T" = c(51,53,57,58,59,59,54,50,51,59,63,33,37,0,0,13,41)),
                   highGC.5SS = rbind("A" = c(32,56,8,0,0,38,70,5,13),
                                      "C" = c(38,15,4,0,0,4,9,5,21),
                                      "G" = c(19,15,80,100,0,56,14,86,25),
                                      "T" = c(11,14,8,0,100,2,7,4,41)),
                   lowGC.5SS = rbind("A" = c(38,62,12,0,0,71,73,11,21),
                                     "C" = c(31,10,4,0,0,2,6,6,10),
                                     "G" = c(18,12,77,100,0,24,8,75,14),
                                     "T" = c(13,16,7,0,100,3,13,8,55)))
genomic.acgt = getBackgroundFrequencies("hg19")
splsitesPWM = toPWM(splsitesPFM, prior=genomic.acgt)


exonMap = rowRanges(rse_exon)
exonMap = split(exonMap, exonMap$gencodeID)
ranges = lapply(exonMap, range)
intronMap = mapply(function(exon, r) setdiff(r,exon), exonMap, ranges, SIMPLIFY = F)

exonMap = do.call(getMethod(c, "GenomicRanges"), GRangesList(exonMap))
exonMap = exonMap[which(width(exonMap)>30)]
exonMap = exonMap[which(as.character(seqnames(exonMap)) %in% paste0("chr",c(1:22,"X","Y")))]
intronMap = do.call(getMethod(c, "GenomicRanges"), GRangesList(intronMap))
intronMap = intronMap[which(as.character(seqnames(intronMap)) %in% paste0("chr",c(1:22,"X","Y")))]
intronMap = intronMap[which(width(intronMap)>30)]

sampleset = c(sample(granges(exonMap), 12500, replace = FALSE),sample(intronMap, 12500, replace = FALSE))
sampleset = sortSeqlevels(sampleset)
sampleset = sort(sampleset)

## Make new log-normal background on mix of exons and introns for canonical splice donor and acceptor sequence

sampleseq = getSeq(Hsapiens, sampleset)
lognSpliceSeq = makePWMLognBackground(sampleseq, splsitesPWM, algorithm = "human", bg.len = 50, verbose=F)


## Run motifEnrichment on regions around associated Cs

Cgr = list(CpG.pos = granges(mres$CpG$meth)[which(mres$CpG$eqtls$statistic>0),], nonCpG.pos = granges(mres$nonCpG$meth)[which(mres$nonCpG$eqtls$statistic>0),],
           CpG.neg = granges(mres$CpG$meth)[which(mres$CpG$eqtls$statistic<0),], nonCpG.neg = granges(mres$nonCpG$meth)[which(mres$nonCpG$eqtls$statistic<0),])
for (i in 1:length(Cgr)) {
  names(Cgr[[i]]) = paste0(as.character(seqnames(Cgr[[i]])), ":",as.character(start(Cgr[[i]])), "-",as.character(end(Cgr[[i]])))
  start(Cgr[[i]]) = start(Cgr[[i]]) - 15
  end(Cgr[[i]]) = end(Cgr[[i]]) + 15
}

Cseq = lapply(Cgr, function(x) getSeq(Hsapiens, x))

C.allMotifsBack = lapply(Cseq, function(x) motifEnrichment(x, lognBcgnd, verbose=F))
C.canonicalspliceBack = lapply(Cseq, function(x) motifEnrichment(x, lognSpliceSeq, verbose=F))

save(score, splice.exonBack, lognBcgnd, lognSpliceSeq, C.allMotifsBack, C.canonicalspliceBack,
     file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/PWMEnrich/spliceRegions_PWMEnrich_objects_exonAdjustedBackgrounds.rda")


## make a score object for the C's

spliceMotif.group = lapply(C.canonicalspliceBack, groupReport)
spliceMotif.group = lapply(spliceMotif.group, as.data.frame)
spliceMotif.group = Map(cbind, spliceMotif.group, padj = lapply(spliceMotif.group, function(x) p.adjust(x$p.value, method = "fdr")))
spliceMotif.group = lapply(spliceMotif.group, function(x) x[order(x$id),])
spliceMotif.group = do.call(rbind, Map(cbind, spliceMotif.group, Context = lapply(as.list(names(spliceMotif.group)), function(x) gsub('\\..*', '', x)),
                                       Direction = lapply(as.list(names(spliceMotif.group)), function(x) gsub('.*\\.', '', x))))

scoreC = pvalC = list()
for (i in 1:length(C.allMotifsBack)) {
  scoreC[[i]] = pvalC[[i]] = list(vector("list", length(C.allMotifsBack[[i]]$sequences)))
  for (j in 1:length(C.allMotifsBack[[i]]$sequences)) {
    scoreC[[i]][[j]] = list(sequenceReport(C.allMotifsBack[[i]], seq.id=j),
                            sequenceReport(C.canonicalspliceBack[[i]], seq.id=j))
    pvalC[[i]][[j]] = list(sequenceReport(C.allMotifsBack[[i]], seq.id=j),
                           sequenceReport(C.canonicalspliceBack[[i]], seq.id=j))
  }
}
names(scoreC) = names(pvalC) = names(C.allMotifsBack)

scoreC = lapply(scoreC, function(x) lapply(x, function(z) lapply(z, function(y) as.data.frame(y)[order(y$id),])))
scoreC = list(lapply(scoreC, function(y) do.call(cbind, lapply(y, function(z) z[[1]]$raw.score))),
              lapply(scoreC, function(y) do.call(cbind, lapply(y, function(z) z[[2]]$raw.score))))
for (i in 1:length(scoreC[[1]])) { rownames(scoreC[[1]][[i]]) = group$CpG.pos$id }
for (i in 1:length(scoreC[[2]])) { rownames(scoreC[[2]][[i]]) = c("highGC.3SS","highGC.5SS","lowGC.3SS","lowGC.5SS") }
scoreC = mapply(function(all, spl) rbind(all, spl), scoreC[[1]], scoreC[[2]], SIMPLIFY = F)
for (i in 1:length(scoreC)) { colnames(scoreC[[i]]) = names(C.canonicalspliceBack[[i]]$sequences) }
scoreC = lapply(scoreC, as.data.frame)
scoreC = Map(cbind, scoreC, Context = lapply(as.list(names(scoreC)), function(x) gsub('\\..*', '', x)),
            Direction = lapply(as.list(names(scoreC)), function(x) gsub('.*\\.', '', x)), id = lapply(scoreC, rownames))
cols = lapply(scoreC, colnames)
cols = lapply(cols, function(x) data.frame(exonID = x[1:(length(x)-3)], col = paste0("col",1:length(x[1:(length(x)-3)]))))
for (i in 1:length(scoreC)) { colnames(scoreC[[i]]) = c(as.character(cols[[i]][,"col"]),"Context","Direction","motifID") }
scoreC = lapply(scoreC, reshape2::melt)

pvalC = lapply(pvalC, function(x) lapply(x, function(z) lapply(z, function(y) as.data.frame(y)[order(y$id),])))
pvalC = list(lapply(pvalC, function(y) do.call(cbind, lapply(y, function(z) p.adjust(z[[1]]$p.value, method = "fdr")))),
             lapply(pvalC, function(y) do.call(cbind, lapply(y, function(z) p.adjust(z[[2]]$p.value, method = "fdr")))))
for (i in 1:length(pvalC[[1]])) { rownames(pvalC[[1]][[i]]) = as.character(group$CpG.pos$id) }
for (i in 1:length(pvalC[[2]])) { rownames(pvalC[[2]][[i]]) = c("highGC.3SS","highGC.5SS","lowGC.3SS","lowGC.5SS") }
pvalC = mapply(function(all, spl) rbind(all, spl), pvalC[[1]], pvalC[[2]], SIMPLIFY = F)
for (i in 1:length(pvalC)) { colnames(pvalC[[i]]) = names(C.canonicalspliceBack[[i]]$sequences) }
pvalC = lapply(pvalC, as.data.frame)
pvalC = Map(cbind, pvalC, Context = lapply(as.list(names(pvalC)), function(x) gsub('\\..*', '', x)),
            Direction = lapply(as.list(names(pvalC)), function(x) gsub('.*\\.', '', x)), id = lapply(pvalC, rownames))
cols = lapply(pvalC, colnames)
cols = lapply(cols, function(x) data.frame(exonID = x[1:(length(x)-3)], col = paste0("col",1:length(x[1:(length(x)-3)]))))
for (i in 1:length(pvalC)) { colnames(pvalC[[i]]) = c(as.character(cols[[i]][,"col"]),"Context","Direction","motifID") }
pvalC = lapply(pvalC, reshape2::melt)



# does it overlap?
gr = list(CpG.pos = rowRanges(mres$CpG$expr)[which(mres$CpG$eqtls$statistic>0),], nonCpG.pos = rowRanges(mres$nonCpG$expr)[which(mres$nonCpG$eqtls$statistic>0),],
          CpG.neg = rowRanges(mres$CpG$expr)[which(mres$CpG$eqtls$statistic<0),], nonCpG.neg = rowRanges(mres$nonCpG$expr)[which(mres$nonCpG$eqtls$statistic<0),])
Cgr = list(CpG.pos = granges(mres$CpG$meth)[which(mres$CpG$eqtls$statistic>0),], nonCpG.pos = granges(mres$nonCpG$meth)[which(mres$nonCpG$eqtls$statistic>0),],
           CpG.neg = granges(mres$CpG$meth)[which(mres$CpG$eqtls$statistic<0),], nonCpG.neg = granges(mres$nonCpG$meth)[which(mres$nonCpG$eqtls$statistic<0),])

oo = mapply(function(exon, mC) findOverlaps(exon,mC), lapply(gr, function(x) x[width(x)>30]),
            mapply(function(c,g) c[width(g)>30], Cgr, gr, SIMPLIFY = F), SIMPLIFY = F)
overlapping = lapply(oo, function(x) x[which(queryHits(x)==subjectHits(x))])

colOverlaps = mapply(function(dat, ov) as.character(dat[queryHits(ov),"col"]), cols, overlapping, SIMPLIFY = F)

x = Map(cbind, scoreC, C_in_exon = mapply(function(dat, ov) ifelse(dat$variable %in% ov, "Yes", "No"),
                                         scoreC, colOverlaps, SIMPLIFY = F),
        Target = lapply(scoreC, function(x) group$CpG.pos[match(x$motifID, group$CpG.pos$id), "target"]),
        C_ID = mapply(function(x, co) co[match(x$variable, co$col), "exonID"], scoreC, cols, SIMPLIFY = F))
scoreC = do.call(rbind, x)

x = Map(cbind, pvalC, C_in_exon = mapply(function(dat, ov) ifelse(dat$variable %in% ov, "Yes", "No"),
                                         pvalC, colOverlaps, SIMPLIFY = F),
        Target = lapply(pvalC, function(x) group$CpG.pos[match(x$motifID, group$CpG.pos$id), "target"]),
        C_ID = mapply(function(x, co) co[match(x$variable, co$col), "exonID"], pvalC, cols, SIMPLIFY = F))
pvalC = do.call(rbind, x)

scoreC$Direction = ifelse(scoreC$Direction=="pos", "Positive Association", "Negative Association")
score$Direction = ifelse(score$Direction=="pos", "Positive Association", "Negative Association")
scoreC$C_in_exon = ifelse(scoreC$C_in_exon=="Yes", "C in Exon", "C not in Exon")
score$C_in_exon = ifelse(score$C_in_exon=="Yes", "C in Exon", "C not in Exon")
pvalC$Direction = ifelse(pvalC$Direction=="pos", "Positive Association", "Negative Association")
pvalC$C_in_exon = ifelse(pvalC$C_in_exon=="Yes", "C in Exon", "C not in Exon")
colnames(pvalC)[colnames(pvalC)=="value"] = "FDR"

save(score, pval, splice.exonBack, lognBcgnd, lognSpliceSeq, C.allMotifsBack, C.canonicalspliceBack, spliceMotif.group, scoreC, pvalC,
     file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/PWMEnrich/spliceRegions_PWMEnrich_objects_exonAdjustedBackgrounds.rda")


## Compare distribution of the raw score for CTCF

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/PWMEnrich/figures/CTCF_rawscore_distribution.pdf", width =10)
ggplot(scoreC[which(scoreC$motifID=="CTCF" & scoreC$Context=="CpG"),], aes(value, fill = C_in_exon)) + scale_fill_brewer(palette = "Set1") +
  geom_density(alpha = 1/2) + theme_classic() + 
  facet_grid(. ~ Direction) +
  ylab("Density") + 
  xlab("") + xlim(0,10) +
  ggtitle("CTCF Enrichent around C: CpG context") + 
  theme(title = element_text(size = 20)) + theme(text = element_text(size = 20), legend.title=element_blank()) + theme(legend.position="bottom")
ggplot(scoreC[which(scoreC$motifID=="CTCF" & scoreC$Context=="nonCpG"),], aes(value, fill = C_in_exon)) + scale_fill_brewer(palette = "Set1") +
  geom_density(alpha = 1/2) + theme_classic() + 
  facet_grid(. ~ Direction) +
  ylab("Density") + 
  xlab("") + xlim(0,10) +
  ggtitle("CTCF Enrichent around C: CpH context") + 
  theme(title = element_text(size = 20)) + theme(text = element_text(size = 20), legend.title=element_blank()) + theme(legend.position="bottom")
ggplot(score[which(score$motifID=="CTCF" & score$Context=="CpG"),], aes(value, fill = C_in_exon)) + scale_fill_brewer(palette = "Set1") +
  geom_density(alpha = 1/2) + theme_classic() + 
  facet_grid(. ~ Direction) +
  ylab("Density") + 
  xlab("") + xlim(0,10) +
  ggtitle("CTCF Enrichent around Exon: CpG context") + 
  theme(title = element_text(size = 20)) + theme(text = element_text(size = 20), legend.title=element_blank()) + theme(legend.position="bottom")
ggplot(score[which(score$motifID=="CTCF" & score$Context=="nonCpG"),], aes(value, fill = C_in_exon)) + scale_fill_brewer(palette = "Set1") +
  geom_density(alpha = 1/2) + theme_classic() + 
  facet_grid(. ~ Direction) +
  ylab("Density") + 
  xlab("") + xlim(0,10) +
  ggtitle("CTCF Enrichent around Exon: CpH context") + 
  theme(title = element_text(size = 20)) + theme(text = element_text(size = 20), legend.title=element_blank()) + theme(legend.position="bottom")
dev.off()

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/PWMEnrich/figures/MECP2_rawscore_distribution.pdf", width =10)
ggplot(scoreC[which(scoreC$motifID=="MECP2" & scoreC$Context=="CpG"),], aes(value, fill = C_in_exon)) + scale_fill_brewer(palette = "Set1") +
  geom_density(alpha = 1/2) + theme_classic() + 
  facet_grid(. ~ Direction) +
  ylab("Density") + 
  xlab("") + xlim(0,10) +
  ggtitle("MECP2 Enrichent around C: CpG context") + 
  theme(title = element_text(size = 20)) + theme(text = element_text(size = 20), legend.title=element_blank()) + theme(legend.position="bottom")
ggplot(scoreC[which(scoreC$motifID=="MECP2" & scoreC$Context=="nonCpG"),], aes(value, fill = C_in_exon)) + scale_fill_brewer(palette = "Set1") +
  geom_density(alpha = 1/2) + theme_classic() + 
  facet_grid(. ~ Direction) +
  ylab("Density") + 
  xlab("") + xlim(0,10) +
  ggtitle("MECP2 Enrichent around C: CpH context") + 
  theme(title = element_text(size = 20)) + theme(text = element_text(size = 20), legend.title=element_blank()) + theme(legend.position="bottom")
ggplot(score[which(score$motifID=="MECP2" & score$Context=="CpG"),], aes(value, fill = C_in_exon)) + scale_fill_brewer(palette = "Set1") +
  geom_density(alpha = 1/2) + theme_classic() + 
  facet_grid(. ~ Direction) +
  ylab("Density") + 
  xlab("") + xlim(0,10) +
  ggtitle("MECP2 Enrichent around Exon: CpG context") + 
  theme(title = element_text(size = 20)) + theme(text = element_text(size = 20), legend.title=element_blank()) + theme(legend.position="bottom")
ggplot(score[which(score$motifID=="MECP2" & score$Context=="nonCpG"),], aes(value, fill = C_in_exon)) + scale_fill_brewer(palette = "Set1") +
  geom_density(alpha = 1/2) + theme_classic() + 
  facet_grid(. ~ Direction) +
  ylab("Density") + 
  xlab("") + xlim(0,10) +
  ggtitle("MECP2 Enrichent around Exon: CpH context") + 
  theme(title = element_text(size = 20)) + theme(text = element_text(size = 20), legend.title=element_blank()) + theme(legend.position="bottom")
dev.off()


## CTCF and MECP2 enrichment in the C regions

splPvalC = split(pvalC, pvalC$Context)
df = lapply(splPvalC, function(x) list(CTCF.Sig.by.Dir = data.frame(Pos = c(nrow(x[which(x$motifID=="CTCF" & x$FDR<=0.05 & x$Direction=="Positive Association"),]),
                                                          nrow(x[which(x$motifID=="CTCF" & x$FDR>0.05 & x$Direction=="Positive Association"),])),
                                                  Neg = c(nrow(x[which(x$motifID=="CTCF" & x$FDR<=0.05 & x$Direction=="Negative Association"),]),
                                                          nrow(x[which(x$motifID=="CTCF" & x$FDR>0.05 & x$Direction=="Negative Association"),])), row.names = c("sigY","sigN")),
                                       CTCF.Sig.by.inExon.PosDir = data.frame(inExon = c(nrow(x[which(x$motifID=="CTCF" & x$FDR<=0.05 & x$Direction=="Positive Association" & x$C_in_exon=="C in Exon"),]),
                                                          nrow(x[which(x$motifID=="CTCF" & x$FDR>0.05 & x$Direction=="Positive Association" & x$C_in_exon=="C in Exon"),])),
                                                  not = c(nrow(x[which(x$motifID=="CTCF" & x$FDR<=0.05 & x$Direction=="Positive Association" & x$C_in_exon=="C not in Exon"),]),
                                                         nrow(x[which(x$motifID=="CTCF" & x$FDR>0.05 & x$Direction=="Positive Association" & x$C_in_exon=="C not in Exon"),])), row.names = c("sigY","sigN")),
                                       CTCF.Sig.by.inExon.NegDir = data.frame(inExon = c(nrow(x[which(x$motifID=="CTCF" & x$FDR<=0.05 & x$Direction=="Negative Association" & x$C_in_exon=="C in Exon"),]),
                                                          nrow(x[which(x$motifID=="CTCF" & x$FDR>0.05 & x$Direction=="Negative Association" & x$C_in_exon=="C in Exon"),])),
                                                  not = c(nrow(x[which(x$motifID=="CTCF" & x$FDR<=0.05 & x$Direction=="Negative Association" & x$C_in_exon=="C not in Exon"),]),
                                                         nrow(x[which(x$motifID=="CTCF" & x$FDR>0.05 & x$Direction=="Negative Association" & x$C_in_exon=="C not in Exon"),])), row.names = c("sigY","sigN")),
                                       MECP2.Sig.by.Dir = data.frame(Pos = c(nrow(x[which(x$motifID=="MECP2" & x$FDR<=0.05 & x$Direction=="Positive Association"),]),
                                                          nrow(x[which(x$motifID=="MECP2" & x$FDR>0.05 & x$Direction=="Positive Association"),])),
                                                  Neg = c(nrow(x[which(x$motifID=="MECP2" & x$FDR<=0.05 & x$Direction=="Negative Association"),]),
                                                          nrow(x[which(x$motifID=="MECP2" & x$FDR>0.05 & x$Direction=="Negative Association"),])), row.names = c("sigY","sigN")),
                                       MECP2.Sig.by.inExon.PosDir = data.frame(inExon = c(nrow(x[which(x$motifID=="MECP2" & x$FDR<=0.05 & x$Direction=="Positive Association" & x$C_in_exon=="C in Exon"),]),
                                                          nrow(x[which(x$motifID=="MECP2" & x$FDR>0.05 & x$Direction=="Positive Association" & x$C_in_exon=="C in Exon"),])),
                                                  not = c(nrow(x[which(x$motifID=="MECP2" & x$FDR<=0.05 & x$Direction=="Positive Association" & x$C_in_exon=="C not in Exon"),]),
                                                         nrow(x[which(x$motifID=="MECP2" & x$FDR>0.05 & x$Direction=="Positive Association" & x$C_in_exon=="C not in Exon"),])), row.names = c("sigY","sigN")),
                                       MECP2.Sig.by.inExon.NegDir = data.frame(inExon = c(nrow(x[which(x$motifID=="MECP2" & x$FDR<=0.05 & x$Direction=="Negative Association" & x$C_in_exon=="C in Exon"),]),
                                                          nrow(x[which(x$motifID=="MECP2" & x$FDR>0.05 & x$Direction=="Negative Association" & x$C_in_exon=="C in Exon"),])),
                                                  not = c(nrow(x[which(x$motifID=="MECP2" & x$FDR<=0.05 & x$Direction=="Negative Association" & x$C_in_exon=="C not in Exon"),]),
                                                         nrow(x[which(x$motifID=="MECP2" & x$FDR>0.05 & x$Direction=="Negative Association" & x$C_in_exon=="C not in Exon"),])), row.names = c("sigY","sigN"))))
fish = lapply(df, function(x) lapply(x, fisher.test))
unlist(lapply(fish, function(x) lapply(x, function(y) y$p.value)))

pvalC[which(pvalC$motifID=="CTCF" & pvalC$FDR<=0.05),]
#                   Context            Direction motifID variable        FDR
#CpG.neg.7020616        CpG Negative Association    CTCF  col6048 0.02023266
#nonCpG.neg.8816683  nonCpG Negative Association    CTCF  col7595 0.04258352
#C_in_exon Target                     C_ID
#CpG.neg.7020616    C not in Exon   CTCF     chr4:8307741-8307741
#nonCpG.neg.8816683 C not in Exon   CTCF chr2:233346659-233346659
pvalC[which(pvalC$motifID=="MECP2" & pvalC$FDR<=0.05),]
#none

splscoreC = split(scoreC, scoreC$Context)
stat = lapply(splscoreC, function(x) list(CTCF.byDir = wilcox.test(x[which(x$motifID=="CTCF" & x$Direction=="Positive Association"),"value"],
                                                          x[which(x$motifID=="CTCF" & x$Direction=="Negative Association"),"value"]),
                                           CTCF.CinExon.PosDir = wilcox.test(x[which(x$motifID=="CTCF" & x$Direction=="Positive Association" & x$C_in_exon=="C in Exon"),"value"],
                                                                   x[which(x$motifID=="CTCF" & x$Direction=="Positive Association" & x$C_in_exon=="C not in Exon"),"value"]),
                                           CTCF.CinExon.NegDir = wilcox.test(x[which(x$motifID=="CTCF" & x$Direction=="Negative Association" & x$C_in_exon=="C in Exon"),"value"],
                                                                   x[which(x$motifID=="CTCF" & x$Direction=="Negative Association" & x$C_in_exon=="C not in Exon"),"value"]),
                                           MECP2.byDir = wilcox.test(x[which(x$motifID=="MECP2" & x$Direction=="Positive Association"),"value"],
                                                          x[which(x$motifID=="MECP2" & x$Direction=="Negative Association"),"value"]),
                                           MECP2.CinExon.PosDir = wilcox.test(x[which(x$motifID=="MECP2" & x$Direction=="Positive Association" & x$C_in_exon=="C in Exon"),"value"],
                                                                   x[which(x$motifID=="MECP2" & x$Direction=="Positive Association" & x$C_in_exon=="C not in Exon"),"value"]),
                                           MECP2.CinExon.NegDir = wilcox.test(x[which(x$motifID=="MECP2" & x$Direction=="Negative Association" & x$C_in_exon=="C in Exon"),"value"],
                                                                   x[which(x$motifID=="MECP2" & x$Direction=="Negative Association" & x$C_in_exon=="C not in Exon"),"value"])))
df = data.frame(pval = unlist(lapply(stat, function(x) lapply(x, function(y) y$p.value))),
                stat = unlist(lapply(stat, function(x) lapply(x, function(y) y$statistic))),
                comp = gsub(".W","",names(unlist(lapply(stat, function(x) lapply(x, function(y) y$statistic))))))
## CpG.byDir: negatively regulating Cs are more likely to have a CTCF binding site 
#                                    pval       stat
#CpG.CTCF.byDir              8.730035e-05 -3.92504098
#CpG.CTCF.CinExon.PosDir     6.092987e-01  0.51131125
#CpG.CTCF.CinExon.NegDir     7.824377e-01 -0.27615464
#CpG.MECP2.byDir             6.797704e-02 -1.82632631
#CpG.MECP2.CinExon.PosDir    7.850724e-01 -0.27278958
#CpG.MECP2.CinExon.NegDir    6.349520e-01  0.47480731
#nonCpG.CTCF.byDir           1.154910e-01 -1.57417282
#nonCpG.CTCF.CinExon.PosDir  1.613027e-01 -1.40194357
#nonCpG.CTCF.CinExon.NegDir  9.889814e-01 -0.01381454
#nonCpG.MECP2.byDir          5.745768e-01  0.56142618
#nonCpG.MECP2.CinExon.PosDir 8.834245e-01 -0.14704552
#nonCpG.MECP2.CinExon.NegDir 4.409520e-01  0.77102946

splscore = split(score, scoreC$Context)
stat = lapply(splscore, function(x) list(CTCF.byDir = wilcox.test(x[which(x$motifID=="CTCF" & x$Direction=="Positive Association"),"value"],
                                                          x[which(x$motifID=="CTCF" & x$Direction=="Negative Association"),"value"]),
                                          CTCF.CinExon.PosDir = wilcox.test(x[which(x$motifID=="CTCF" & x$Direction=="Positive Association" & x$C_in_exon=="C in Exon"),"value"],
                                                                   x[which(x$motifID=="CTCF" & x$Direction=="Positive Association" & x$C_in_exon=="C not in Exon"),"value"]),
                                          CTCF.CinExon.NegDir = wilcox.test(x[which(x$motifID=="CTCF" & x$Direction=="Negative Association" & x$C_in_exon=="C in Exon"),"value"],
                                                                   x[which(x$motifID=="CTCF" & x$Direction=="Negative Association" & x$C_in_exon=="C not in Exon"),"value"]),
                                          MECP2.byDir = wilcox.test(x[which(x$motifID=="MECP2" & x$Direction=="Positive Association"),"value"],
                                                          x[which(x$motifID=="MECP2" & x$Direction=="Negative Association"),"value"]),
                                          MECP2.CinExon.PosDir = wilcox.test(x[which(x$motifID=="MECP2" & x$Direction=="Positive Association" & x$C_in_exon=="C in Exon"),"value"],
                                                                   x[which(x$motifID=="MECP2" & x$Direction=="Positive Association" & x$C_in_exon=="C not in Exon"),"value"]),
                                          MECP2.CinExon.NegDir = wilcox.test(x[which(x$motifID=="MECP2" & x$Direction=="Negative Association" & x$C_in_exon=="C in Exon"),"value"],
                                                                   x[which(x$motifID=="MECP2" & x$Direction=="Negative Association" & x$C_in_exon=="C not in Exon"),"value"])))
df = data.frame(pval = unlist(lapply(stat, function(x) lapply(x, function(y) y$p.value))),
                 stat = unlist(lapply(stat, function(x) lapply(x, function(y) y$statistic))),
                comp = gsub(".W","",names(unlist(lapply(stat, function(x) lapply(x, function(y) y$statistic))))))
#                                    pval      stat                        comp
#CpG.CTCF.byDir              8.787657e-03  2.6239977              CpG.CTCF.byDir
#CpG.CTCF.CinExon.PosDir     4.746548e-03 -2.8309108     CpG.CTCF.CinExon.PosDir
#CpG.CTCF.CinExon.NegDir     1.415374e-01  1.4704656     CpG.CTCF.CinExon.NegDir
#CpG.MECP2.byDir             4.949280e-01  0.6826348             CpG.MECP2.byDir
#CpG.MECP2.CinExon.PosDir    7.594054e-01 -0.3063272    CpG.MECP2.CinExon.PosDir
#CpG.MECP2.CinExon.NegDir    6.210428e-12  6.8884467    CpG.MECP2.CinExon.NegDir
#nonCpG.CTCF.byDir           2.389861e-02  2.2616365           nonCpG.CTCF.byDir
#nonCpG.CTCF.CinExon.PosDir  2.686990e-02 -2.2162465  nonCpG.CTCF.CinExon.PosDir
#nonCpG.CTCF.CinExon.NegDir  7.043395e-03 -2.6955590  nonCpG.CTCF.CinExon.NegDir
#nonCpG.MECP2.byDir          6.138577e-03  2.7440296          nonCpG.MECP2.byDir
#nonCpG.MECP2.CinExon.PosDir 1.772224e-03  3.1920322 nonCpG.MECP2.CinExon.PosDir
#nonCpG.MECP2.CinExon.NegDir 6.165861e-11  6.6184307 nonCpG.MECP2.CinExon.NegDir

df[df$pval<=0.05,]

score[which(score$motifID %in% c("MECP2","CTCF")),]
scoreC[which(scoreC$motifID %in% c("MECP2","CTCF")),]
pval[which(pval$motifID %in% c("MECP2","CTCF")),]
pvalC[which(pvalC$motifID %in% c("MECP2","CTCF")),]


## Try making 2x2 contingency tables about scores > or <= 2

stat = list(CTCF.byDir.CpG = data.frame(c(nrow(splscoreC$CpG[which(splscoreC$CpG$motifID=="CTCF" & splscoreC$CpG$Direction=="Positive Association" & splscoreC$CpG$value>2),]),
                                          nrow(splscoreC$CpG[which(splscoreC$CpG$motifID=="CTCF" & splscoreC$CpG$Direction=="Positive Association" & splscoreC$CpG$value<=2),])),
                                        c(nrow(splscoreC$CpG[which(splscoreC$CpG$motifID=="CTCF" & splscoreC$CpG$Direction=="Negative Association" & splscoreC$CpG$value>2),]),
                                          nrow(splscoreC$CpG[which(splscoreC$CpG$motifID=="CTCF" & splscoreC$CpG$Direction=="Negative Association" & splscoreC$CpG$value<=2),]))),
            CTCF.byDir.CpH = data.frame(c(nrow(splscoreC$nonCpG[which(splscoreC$nonCpG$motifID=="CTCF" & splscoreC$nonCpG$Direction=="Positive Association" & splscoreC$nonCpG$value>2),]),
                                          nrow(splscoreC$nonCpG[which(splscoreC$nonCpG$motifID=="CTCF" & splscoreC$nonCpG$Direction=="Positive Association" & splscoreC$nonCpG$value<=2),])),
                                        c(nrow(splscoreC$nonCpG[which(splscoreC$nonCpG$motifID=="CTCF" & splscoreC$nonCpG$Direction=="Negative Association" & splscoreC$nonCpG$value>2),]),
                                          nrow(splscoreC$nonCpG[which(splscoreC$nonCpG$motifID=="CTCF" & splscoreC$nonCpG$Direction=="Negative Association" & splscoreC$nonCpG$value<=2),]))),
            CTCF.CinExon.PosDir.CpG = data.frame(c(nrow(splscoreC$CpG[which(splscoreC$CpG$motifID=="CTCF" & splscoreC$CpG$Direction=="Positive Association" & splscoreC$CpG$C_in_exon=="C in Exon" & splscoreC$CpG$value>2),]),
                                                   nrow(splscoreC$CpG[which(splscoreC$CpG$motifID=="CTCF" & splscoreC$CpG$Direction=="Positive Association" & splscoreC$CpG$C_in_exon=="C in Exon" & splscoreC$CpG$value<=2),])),
                                                 c(nrow(splscoreC$CpG[which(splscoreC$CpG$motifID=="CTCF" & splscoreC$CpG$Direction=="Positive Association" & splscoreC$CpG$C_in_exon=="C not in Exon" & splscoreC$CpG$value>2),]),
                                                   nrow(splscoreC$CpG[which(splscoreC$CpG$motifID=="CTCF" & splscoreC$CpG$Direction=="Positive Association" & splscoreC$CpG$C_in_exon=="C not in Exon" & splscoreC$CpG$value<=2),]))),
            CTCF.CinExon.PosDir.CpH = data.frame(c(nrow(splscoreC$nonCpG[which(splscoreC$nonCpG$motifID=="CTCF" & splscoreC$nonCpG$Direction=="Positive Association" & splscoreC$nonCpG$C_in_exon=="C in Exon" & splscoreC$nonCpG$value>2),]),
                                                   nrow(splscoreC$nonCpG[which(splscoreC$nonCpG$motifID=="CTCF" & splscoreC$nonCpG$Direction=="Positive Association" & splscoreC$nonCpG$C_in_exon=="C in Exon" & splscoreC$nonCpG$value<=2),])),
                                                 c(nrow(splscoreC$nonCpG[which(splscoreC$nonCpG$motifID=="CTCF" & splscoreC$nonCpG$Direction=="Positive Association" & splscoreC$nonCpG$C_in_exon=="C not in Exon" & splscoreC$nonCpG$value>2),]),
                                                   nrow(splscoreC$nonCpG[which(splscoreC$nonCpG$motifID=="CTCF" & splscoreC$nonCpG$Direction=="Positive Association" & splscoreC$nonCpG$C_in_exon=="C not in Exon" & splscoreC$nonCpG$value<=2),]))),
            CTCF.CinExon.NegDir.CpG = data.frame(c(nrow(splscoreC$CpG[which(splscoreC$CpG$motifID=="CTCF" & splscoreC$CpG$Direction=="Negative Association" & splscoreC$CpG$C_in_exon=="C in Exon" & splscoreC$CpG$value>2),]),
                                                   nrow(splscoreC$CpG[which(splscoreC$CpG$motifID=="CTCF" & splscoreC$CpG$Direction=="Negative Association" & splscoreC$CpG$C_in_exon=="C in Exon" & splscoreC$CpG$value<=2),])),
                                                 c(nrow(splscoreC$CpG[which(splscoreC$CpG$motifID=="CTCF" & splscoreC$CpG$Direction=="Negative Association" & splscoreC$CpG$C_in_exon=="C not in Exon" & splscoreC$CpG$value>2),]),
                                                   nrow(splscoreC$CpG[which(splscoreC$CpG$motifID=="CTCF" & splscoreC$CpG$Direction=="Negative Association" & splscoreC$CpG$C_in_exon=="C not in Exon" & splscoreC$CpG$value<=2),]))),
            CTCF.CinExon.NegDir.CpH = data.frame(c(nrow(splscoreC$nonCpG[which(splscoreC$nonCpG$motifID=="CTCF" & splscoreC$nonCpG$Direction=="Negative Association" & splscoreC$nonCpG$C_in_exon=="C in Exon" & splscoreC$nonCpG$value>2),]),
                                                   nrow(splscoreC$nonCpG[which(splscoreC$nonCpG$motifID=="CTCF" & splscoreC$nonCpG$Direction=="Negative Association" & splscoreC$nonCpG$C_in_exon=="C in Exon" & splscoreC$nonCpG$value<=2),])),
                                                 c(nrow(splscoreC$nonCpG[which(splscoreC$nonCpG$motifID=="CTCF" & splscoreC$nonCpG$Direction=="Negative Association" & splscoreC$nonCpG$C_in_exon=="C not in Exon" & splscoreC$nonCpG$value>2),]),
                                                   nrow(splscoreC$nonCpG[which(splscoreC$nonCpG$motifID=="CTCF" & splscoreC$nonCpG$Direction=="Negative Association" & splscoreC$nonCpG$C_in_exon=="C not in Exon" & splscoreC$nonCpG$value<=2),]))),
            CTCF.CinExon.Dir.CpG = data.frame(c(nrow(splscoreC$CpG[which(splscoreC$CpG$motifID=="CTCF" & splscoreC$CpG$Direction=="Positive Association" & splscoreC$CpG$C_in_exon=="C in Exon" & splscoreC$CpG$value>2),]),
                                                   nrow(splscoreC$CpG[which(splscoreC$CpG$motifID=="CTCF" & splscoreC$CpG$Direction=="Positive Association" & splscoreC$CpG$C_in_exon=="C in Exon" & splscoreC$CpG$value<=2),])),
                                                 c(nrow(splscoreC$CpG[which(splscoreC$CpG$motifID=="CTCF" & splscoreC$CpG$Direction=="Negative Association" & splscoreC$CpG$C_in_exon=="C in Exon" & splscoreC$CpG$value>2),]),
                                                   nrow(splscoreC$CpG[which(splscoreC$CpG$motifID=="CTCF" & splscoreC$CpG$Direction=="Negative Association" & splscoreC$CpG$C_in_exon=="C in Exon" & splscoreC$CpG$value<=2),]))),
            CTCF.CinExon.Dir.CpH = data.frame(c(nrow(splscoreC$nonCpG[which(splscoreC$nonCpG$motifID=="CTCF" & splscoreC$nonCpG$Direction=="Positive Association" & splscoreC$nonCpG$C_in_exon=="C in Exon" & splscoreC$nonCpG$value>2),]),
                                                   nrow(splscoreC$nonCpG[which(splscoreC$nonCpG$motifID=="CTCF" & splscoreC$nonCpG$Direction=="Positive Association" & splscoreC$nonCpG$C_in_exon=="C in Exon" & splscoreC$nonCpG$value<=2),])),
                                                 c(nrow(splscoreC$nonCpG[which(splscoreC$nonCpG$motifID=="CTCF" & splscoreC$nonCpG$Direction=="Negative Association" & splscoreC$nonCpG$C_in_exon=="C in Exon" & splscoreC$nonCpG$value>2),]),
                                                   nrow(splscoreC$nonCpG[which(splscoreC$nonCpG$motifID=="CTCF" & splscoreC$nonCpG$Direction=="Negative Association" & splscoreC$nonCpG$C_in_exon=="C in Exon" & splscoreC$nonCpG$value<=2),]))),
            CTCF.CNOTinExon.Dir.CpG = data.frame(c(nrow(splscoreC$CpG[which(splscoreC$CpG$motifID=="CTCF" & splscoreC$CpG$Direction=="Positive Association" & splscoreC$CpG$C_in_exon=="C not in Exon" & splscoreC$CpG$value>2),]),
                                                   nrow(splscoreC$CpG[which(splscoreC$CpG$motifID=="CTCF" & splscoreC$CpG$Direction=="Positive Association" & splscoreC$CpG$C_in_exon=="C not in Exon" & splscoreC$CpG$value<=2),])),
                                                 c(nrow(splscoreC$CpG[which(splscoreC$CpG$motifID=="CTCF" & splscoreC$CpG$Direction=="Negative Association" & splscoreC$CpG$C_in_exon=="C not in Exon" & splscoreC$CpG$value>2),]),
                                                   nrow(splscoreC$CpG[which(splscoreC$CpG$motifID=="CTCF" & splscoreC$CpG$Direction=="Negative Association" & splscoreC$CpG$C_in_exon=="C not in Exon" & splscoreC$CpG$value<=2),]))),
            CTCF.CNOTinExon.Dir.CpH = data.frame(c(nrow(splscoreC$nonCpG[which(splscoreC$nonCpG$motifID=="CTCF" & splscoreC$nonCpG$Direction=="Positive Association" & splscoreC$nonCpG$C_in_exon=="C not in Exon" & splscoreC$nonCpG$value>2),]),
                                                   nrow(splscoreC$nonCpG[which(splscoreC$nonCpG$motifID=="CTCF" & splscoreC$nonCpG$Direction=="Positive Association" & splscoreC$nonCpG$C_in_exon=="C not in Exon" & splscoreC$nonCpG$value<=2),])),
                                                 c(nrow(splscoreC$nonCpG[which(splscoreC$nonCpG$motifID=="CTCF" & splscoreC$nonCpG$Direction=="Negative Association" & splscoreC$nonCpG$C_in_exon=="C not in Exon" & splscoreC$nonCpG$value>2),]),
                                                   nrow(splscoreC$nonCpG[which(splscoreC$nonCpG$motifID=="CTCF" & splscoreC$nonCpG$Direction=="Negative Association" & splscoreC$nonCpG$C_in_exon=="C not in Exon" & splscoreC$nonCpG$value<=2),]))),
            MECP2.byDir.CpG = data.frame(c(nrow(splscoreC$CpG[which(splscoreC$CpG$motifID=="MECP2" & splscoreC$CpG$Direction=="Positive Association" & splscoreC$CpG$value>2),]),
                                           nrow(splscoreC$CpG[which(splscoreC$CpG$motifID=="MECP2" & splscoreC$CpG$Direction=="Positive Association" & splscoreC$CpG$value<=2),])),
                                         c(nrow(splscoreC$CpG[which(splscoreC$CpG$motifID=="MECP2" & splscoreC$CpG$Direction=="Negative Association" & splscoreC$CpG$value>2),]),
                                           nrow(splscoreC$CpG[which(splscoreC$CpG$motifID=="MECP2" & splscoreC$CpG$Direction=="Negative Association" & splscoreC$CpG$value<=2),]))),
            MECP2.byDir.CpH = data.frame(c(nrow(splscoreC$nonCpG[which(splscoreC$nonCpG$motifID=="MECP2" & splscoreC$nonCpG$Direction=="Positive Association" & splscoreC$nonCpG$value>2),]),
                                           nrow(splscoreC$nonCpG[which(splscoreC$nonCpG$motifID=="MECP2" & splscoreC$nonCpG$Direction=="Positive Association" & splscoreC$nonCpG$value<=2),])),
                                         c(nrow(splscoreC$nonCpG[which(splscoreC$nonCpG$motifID=="MECP2" & splscoreC$nonCpG$Direction=="Negative Association" & splscoreC$nonCpG$value>2),]),
                                           nrow(splscoreC$nonCpG[which(splscoreC$nonCpG$motifID=="MECP2" & splscoreC$nonCpG$Direction=="Negative Association" & splscoreC$nonCpG$value<=2),]))),
            MECP2.CinExon.PosDir.CpG = data.frame(c(nrow(splscoreC$CpG[which(splscoreC$CpG$motifID=="MECP2" & splscoreC$CpG$Direction=="Positive Association" & splscoreC$CpG$C_in_exon=="C in Exon" & splscoreC$CpG$value>2),]),
                                                    nrow(splscoreC$CpG[which(splscoreC$CpG$motifID=="MECP2" & splscoreC$CpG$Direction=="Positive Association" & splscoreC$CpG$C_in_exon=="C in Exon" & splscoreC$CpG$value<=2),])),
                                                  c(nrow(splscoreC$CpG[which(splscoreC$CpG$motifID=="MECP2" & splscoreC$CpG$Direction=="Positive Association" & splscoreC$CpG$C_in_exon=="C not in Exon" & splscoreC$CpG$value>2),]),
                                                    nrow(splscoreC$CpG[which(splscoreC$CpG$motifID=="MECP2" & splscoreC$CpG$Direction=="Positive Association" & splscoreC$CpG$C_in_exon=="C not in Exon" & splscoreC$CpG$value<=2),]))),
            MECP2.CinExon.PosDir.CpH = data.frame(c(nrow(splscoreC$nonCpG[which(splscoreC$nonCpG$motifID=="MECP2" & splscoreC$nonCpG$Direction=="Positive Association" & splscoreC$nonCpG$C_in_exon=="C in Exon" & splscoreC$nonCpG$value>2),]),
                                                    nrow(splscoreC$nonCpG[which(splscoreC$nonCpG$motifID=="MECP2" & splscoreC$nonCpG$Direction=="Positive Association" & splscoreC$nonCpG$C_in_exon=="C in Exon" & splscoreC$nonCpG$value<=2),])),
                                                  c(nrow(splscoreC$nonCpG[which(splscoreC$nonCpG$motifID=="MECP2" & splscoreC$nonCpG$Direction=="Positive Association" & splscoreC$nonCpG$C_in_exon=="C not in Exon" & splscoreC$nonCpG$value>2),]),
                                                    nrow(splscoreC$nonCpG[which(splscoreC$nonCpG$motifID=="MECP2" & splscoreC$nonCpG$Direction=="Positive Association" & splscoreC$nonCpG$C_in_exon=="C not in Exon" & splscoreC$nonCpG$value<=2),]))),
            MECP2.CinExon.NegDir.CpG = data.frame(c(nrow(splscoreC$CpG[which(splscoreC$CpG$motifID=="MECP2" & splscoreC$CpG$Direction=="Negative Association" & splscoreC$CpG$C_in_exon=="C in Exon" & splscoreC$CpG$value>2),]),
                                                    nrow(splscoreC$CpG[which(splscoreC$CpG$motifID=="MECP2" & splscoreC$CpG$Direction=="Negative Association" & splscoreC$CpG$C_in_exon=="C in Exon" & splscoreC$CpG$value<=2),])),
                                                  c(nrow(splscoreC$CpG[which(splscoreC$CpG$motifID=="MECP2" & splscoreC$CpG$Direction=="Negative Association" & splscoreC$CpG$C_in_exon=="C not in Exon" & splscoreC$CpG$value>2),]),
                                                    nrow(splscoreC$CpG[which(splscoreC$CpG$motifID=="MECP2" & splscoreC$CpG$Direction=="Negative Association" & splscoreC$CpG$C_in_exon=="C not in Exon" & splscoreC$CpG$value<=2),]))),
            MECP2.CinExon.NegDir.CpH = data.frame(c(nrow(splscoreC$nonCpG[which(splscoreC$nonCpG$motifID=="MECP2" & splscoreC$nonCpG$Direction=="Negative Association" & splscoreC$nonCpG$C_in_exon=="C in Exon" & splscoreC$nonCpG$value>2),]),
                                                    nrow(splscoreC$nonCpG[which(splscoreC$nonCpG$motifID=="MECP2" & splscoreC$nonCpG$Direction=="Negative Association" & splscoreC$nonCpG$C_in_exon=="C in Exon" & splscoreC$nonCpG$value<=2),])),
                                                  c(nrow(splscoreC$nonCpG[which(splscoreC$nonCpG$motifID=="MECP2" & splscoreC$nonCpG$Direction=="Negative Association" & splscoreC$nonCpG$C_in_exon=="C not in Exon" & splscoreC$nonCpG$value>2),]),
                                                    nrow(splscoreC$nonCpG[which(splscoreC$nonCpG$motifID=="MECP2" & splscoreC$nonCpG$Direction=="Negative Association" & splscoreC$nonCpG$C_in_exon=="C not in Exon" & splscoreC$nonCpG$value<=2),]))),
            MECP2.CinExon.Dir.CpG = data.frame(c(nrow(splscoreC$CpG[which(splscoreC$CpG$motifID=="MECP2" & splscoreC$CpG$Direction=="Positive Association" & splscoreC$CpG$C_in_exon=="C in Exon" & splscoreC$CpG$value>2),]),
                                                nrow(splscoreC$CpG[which(splscoreC$CpG$motifID=="MECP2" & splscoreC$CpG$Direction=="Positive Association" & splscoreC$CpG$C_in_exon=="C in Exon" & splscoreC$CpG$value<=2),])),
                                              c(nrow(splscoreC$CpG[which(splscoreC$CpG$motifID=="MECP2" & splscoreC$CpG$Direction=="Negative Association" & splscoreC$CpG$C_in_exon=="C in Exon" & splscoreC$CpG$value>2),]),
                                                nrow(splscoreC$CpG[which(splscoreC$CpG$motifID=="MECP2" & splscoreC$CpG$Direction=="Negative Association" & splscoreC$CpG$C_in_exon=="C in Exon" & splscoreC$CpG$value<=2),]))),
            MECP2.CinExon.Dir.CpH = data.frame(c(nrow(splscoreC$nonCpG[which(splscoreC$nonCpG$motifID=="MECP2" & splscoreC$nonCpG$Direction=="Positive Association" & splscoreC$nonCpG$C_in_exon=="C in Exon" & splscoreC$nonCpG$value>2),]),
                                                nrow(splscoreC$nonCpG[which(splscoreC$nonCpG$motifID=="MECP2" & splscoreC$nonCpG$Direction=="Positive Association" & splscoreC$nonCpG$C_in_exon=="C in Exon" & splscoreC$nonCpG$value<=2),])),
                                              c(nrow(splscoreC$nonCpG[which(splscoreC$nonCpG$motifID=="MECP2" & splscoreC$nonCpG$Direction=="Negative Association" & splscoreC$nonCpG$C_in_exon=="C in Exon" & splscoreC$nonCpG$value>2),]),
                                                nrow(splscoreC$nonCpG[which(splscoreC$nonCpG$motifID=="MECP2" & splscoreC$nonCpG$Direction=="Negative Association" & splscoreC$nonCpG$C_in_exon=="C in Exon" & splscoreC$nonCpG$value<=2),]))),
            MECP2.CNOTinExon.Dir.CpG = data.frame(c(nrow(splscoreC$CpG[which(splscoreC$CpG$motifID=="MECP2" & splscoreC$CpG$Direction=="Positive Association" & splscoreC$CpG$C_in_exon=="C not in Exon" & splscoreC$CpG$value>2),]),
                                                   nrow(splscoreC$CpG[which(splscoreC$CpG$motifID=="MECP2" & splscoreC$CpG$Direction=="Positive Association" & splscoreC$CpG$C_in_exon=="C not in Exon" & splscoreC$CpG$value<=2),])),
                                                 c(nrow(splscoreC$CpG[which(splscoreC$CpG$motifID=="MECP2" & splscoreC$CpG$Direction=="Negative Association" & splscoreC$CpG$C_in_exon=="C not in Exon" & splscoreC$CpG$value>2),]),
                                                   nrow(splscoreC$CpG[which(splscoreC$CpG$motifID=="MECP2" & splscoreC$CpG$Direction=="Negative Association" & splscoreC$CpG$C_in_exon=="C not in Exon" & splscoreC$CpG$value<=2),]))),
            MECP2.CNOTinExon.Dir.CpH = data.frame(c(nrow(splscoreC$nonCpG[which(splscoreC$nonCpG$motifID=="MECP2" & splscoreC$nonCpG$Direction=="Positive Association" & splscoreC$nonCpG$C_in_exon=="C not in Exon" & splscoreC$nonCpG$value>2),]),
                                                   nrow(splscoreC$nonCpG[which(splscoreC$nonCpG$motifID=="MECP2" & splscoreC$nonCpG$Direction=="Positive Association" & splscoreC$nonCpG$C_in_exon=="C not in Exon" & splscoreC$nonCpG$value<=2),])),
                                                 c(nrow(splscoreC$nonCpG[which(splscoreC$nonCpG$motifID=="MECP2" & splscoreC$nonCpG$Direction=="Negative Association" & splscoreC$nonCpG$C_in_exon=="C not in Exon" & splscoreC$nonCpG$value>2),]),
                                                   nrow(splscoreC$nonCpG[which(splscoreC$nonCpG$motifID=="MECP2" & splscoreC$nonCpG$Direction=="Negative Association" & splscoreC$nonCpG$C_in_exon=="C not in Exon" & splscoreC$nonCpG$value<=2),]))))
fisher = lapply(stat, fisher.test)
df = data.frame(comp = names(fisher), Odds.Ratio = unlist(lapply(fisher, function(y) y$estimate)), pval = unlist(lapply(fisher, function(y) y$p.value)))
df$FDR = p.adjust(df$pval, method = "fdr")
rownames(df) = NULL


