library(PWMEnrich)
library(GenomicRanges)
library(PWMEnrich.Hsapiens.background)
library(BSgenome.Hsapiens.UCSC.hg19)
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
#CpG.pos     1349      454    33.7
#nonCpG.pos  1288       78     6.1
#CpG.neg     9674     2346    24.3
#nonCpG.neg  7874      605     7.7


## Get sequences from different annotations in the correct format

seq = lapply(gr, function(x) getSeq(Hsapiens, x))
seq = lapply(seq, function(x) x[which(width(x)>30)])


### Run PWMEnrich

useBigMemoryPWMEnrich(TRUE)
registerCoresPWMEnrich(5)

# load the pre-compiled lognormal background computed using promoters
data(PWMLogn.hg19.MotifDb.Hsap)

splice.exon = list()
splice.exon$CpG.pos = motifEnrichment(seq$CpG.pos, PWMLogn.hg19.MotifDb.Hsap, verbose=F)
splice.exon$nonCpG.pos = motifEnrichment(seq$nonCpG.pos, PWMLogn.hg19.MotifDb.Hsap, verbose=F)
splice.exon$CpG.neg = motifEnrichment(seq$CpG.neg, PWMLogn.hg19.MotifDb.Hsap, verbose=F)
splice.exon$nonCpG.neg = motifEnrichment(seq$nonCpG.neg, PWMLogn.hg19.MotifDb.Hsap, verbose=F)


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
group = Map(cbind, group, padj = lapply(group, function(x) p.adjust(x$p.value, method = "fdr")))
group = lapply(group, function(x) x[which(x$target %in% names(PostnataltargettogeneID)),])
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
#CpG.pos    1732.5  MECP2 MECP2 0.7565686       1     0.01595745    1
#nonCpG.pos 1493.0  MECP2 MECP2 0.9724046       1     0.02680747    1
#CpG.neg    1679.5  MECP2 MECP2 0.7646823       1     0.01571010    1
#nonCpG.neg 1721.5  MECP2 MECP2 0.8783678       1     0.02998462    1

exonDiff$CpG["MECP2"] # 0
exonDiff$nonCpG["MECP2"] # 0.002241363
exonDiff$CpG.nonCpG.pos["MECP2"] # -0.002241363 
exonDiff$CpG.nonCpG.neg["MECP2"] # 0

do.call(rbind, lapply(exonDiffQ25, function(x) x[which(x[,3]=="CTCF"),]))
#        Quantile Difference.Statistic pos.target                  pos.id
#CpG.254      Q25        -9.463450e+12       CTCF                    CTCF
#CpG.307      Q25        -1.035549e+08       CTCF Hsapiens-jolma2013-CTCF
#nonCpG       Q25        -3.643336e+11       CTCF                    CTCF
#        pos.raw.score  pos.p.value pos.top.motif.prop      pos.FDR neg.target
#CpG.254      23.04849 1.032534e-41         0.06155015 9.296870e-41       CTCF
#CpG.307       4.14569 1.085640e-32         0.10030395 8.087486e-32       CTCF
#nonCpG       69.75449 3.093087e-26         0.08610885 4.812170e-25       CTCF
#                         neg.id neg.raw.score   neg.p.value neg.top.motif.prop
#CpG.254                    CTCF      4.183613 1.874522e-196         0.07184751
#CpG.307 Hsapiens-jolma2013-CTCF      3.464304  2.291995e-76         0.06598240
#nonCpG                     CTCF      5.761351 1.922173e-156         0.07816504
#              neg.FDR
#CpG.254 1.283542e-195
#CpG.307  9.890176e-76
#nonCpG  1.441314e-155

# Greater enriched in CpG than non-CpG in both expression direction effects, but does not pass the filtering steps
exonDiff$CpG.nonCpG.pos["CTCF"] # 677137.7
exonDiff$CpG.nonCpG.neg["CTCF"] # 9.099117e+12


do.call(rbind,lapply(group, function(x) x[which(x$target=="HP1BP3"),]))
#             rank target     id raw.score     p.value top.motif.prop
#CpG.pos    1732.5 HP1BP3 HP1BP3  1.314786 1.000000000     0.02963526
#nonCpG.pos  717.0 HP1BP3 HP1BP3  2.093420 0.001504384     0.09017059
#CpG.neg    1679.5 HP1BP3 HP1BP3  1.504900 1.000000000     0.03917051
#nonCpG.neg 1721.5 HP1BP3 HP1BP3  1.671395 1.000000000     0.03908252
#                  padj
#CpG.pos    1.000000000
#nonCpG.pos 0.004798503
#CpG.neg    1.000000000
#nonCpG.neg 1.000000000

unlist(lapply(exonDiff, function(x) x[names(x)=="HP1BP3"]))
#CpG.HP1BP3         nonCpG.HP1BP3 CpG.nonCpG.pos.HP1BP3 CpG.nonCpG.neg.HP1BP3
#   0.00000              19.43043             -19.43043               0.00000
 
do.call(rbind,lapply(group, function(x) x[which(x$target=="CBX3"),]))
#           rank target        id raw.score      p.value top.motif.prop
#CpG.pos     815   CBX3 LOC653972 0.8138333 8.433663e-01     0.02583587
#nonCpG.pos  745   CBX3 LOC653972 0.8090082 4.322330e-03     0.04386677
#CpG.neg     882   CBX3 LOC653972 0.8402273 1.555749e-02     0.04566401
#nonCpG.neg  781   CBX3 LOC653972 0.8325335 3.490931e-11     0.04472066
#                   padj
#CpG.pos    1.000000e+00
#nonCpG.pos 1.326868e-02
#CpG.neg    4.034012e-02
#nonCpG.neg 1.022248e-10

do.call(rbind,lapply(group, function(x) x[which(x$target=="CBX7"),]))
#           rank target   id raw.score      p.value top.motif.prop         padj
#CpG.pos    1081   CBX7 CBX7  1.092068 1.000000e+00     0.04635258 1.000000e+00
#nonCpG.pos  686   CBX7 CBX7  1.241742 4.146173e-04     0.07067425 1.382259e-03
#CpG.neg     700   CBX7 CBX7  1.262932 3.120737e-21     0.06263092 1.019590e-20
#nonCpG.neg  558   CBX7 CBX7  1.329684 5.331403e-52     0.07432086 2.185111e-51

unlist(lapply(exonDiff, function(x) x[names(x)=="CBX7"]))
#   CpG.CBX7         nonCpG.CBX7 CpG.nonCpG.pos.CBX7 CpG.nonCpG.neg.CBX7 
#-11918.4571       -3713625.2774            -28.2978       -3701735.1181 
