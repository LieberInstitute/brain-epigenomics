library('bsseq')
library('devtools')
library('jaffelab')
library('RColorBrewer')
library(clusterProfiler)
require(org.Hs.eg.db)
library(GenomicRanges)
library(data.table)
library(VennDiagram)



## Load data
load('/dcl01/lieber/ajaffe/lab/brain-epigenomics/bsseq/bsobj_by_chr/allChrs_postNatal_cleaned_nonCG_noHomogenate_highCov.Rdata')
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/non-CpG/CH_object.rda")
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/non-CpG/CHneurons_object.rda")
load('/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/non-CpG/limma_exploration_mean_mCH_inGeneBody.Rdata')


## Explore gene-level mCH results

nrow(CHgene) # 40563
nrow(CHgene[which(CHgene$padj.CellType<=0.05),]) # 31920
nrow(CHgene[which(CHgene$padj.Age<=0.05),]) # 29658
nrow(CHgene[which(CHgene$padj.Interaction<=0.05),]) # 17022
nrow(CHgene[which(CHgene$padj.AgeNeuron<=0.05),]) # 28543

nrow(CHgene[which(CHgene$padj.CellType<=0.05 & CHgene$Tstat.CellType>0),]) # 31874/31920
nrow(CHgene[which(CHgene$padj.Age<=0.05 & CHgene$Tstat.Age>0),]) # 29622/29658
nrow(CHgene[which(CHgene$padj.Interaction<=0.05 & CHgene$Tstat.Interaction>0),]) # 16911/17022
nrow(CHgene[which(CHgene$padj.AgeNeuron<=0.05 & CHgene$Tstat.AgeNeuron>0),]) # 28487/28543


### Gene Ontology
entrez = list("Hypermethylated\nin Glia" = geneMap[match(CHgene[which(CHgene$padj.CellType<=0.05 & CHgene$Tstat.CellType<0),"geneID"], geneMap$gencodeID),"EntrezID"], 
              "Hypermethylated\nin Neurons" = geneMap[match(CHgene[which(CHgene$padj.CellType<=0.05 & CHgene$Tstat.CellType>0),"geneID"], geneMap$gencodeID),"EntrezID"],
              "Hypermethylated\nin Older Neurons" = geneMap[match(CHgene[which(CHgene$padj.AgeNeuron<=0.05 & CHgene$Tstat.AgeNeuron>0),"geneID"], geneMap$gencodeID),"EntrezID"],
              "Hypermethylated\nin Younger Neurons" = geneMap[match(CHgene[which(CHgene$padj.AgeNeuron<=0.05 & CHgene$Tstat.AgeNeuron<0),"geneID"], geneMap$gencodeID),"EntrezID"])
entrez = lapply(entrez, function(x) as.character(na.omit(unique(x))))              
GeneUniverse = as.character(na.omit(unique(geneMap[match(CHgene$geneID, geneMap$gencodeID),"EntrezID"])))

# Find enriched Pathways via KEGG
elementNROWS(entrez)
keggList = lapply(entrez, function(x) enrichKEGG(as.character(x), organism="human", universe= GeneUniverse, minGSSize=5, pAdjustMethod="BH", qvalueCutoff=1))
keggListdf = lapply(keggList, function(x) as.data.frame(x))
# Enriched Molecular Function GOs
goList_MF = lapply(entrez, function(x) enrichGO(as.character(x), ont = "MF", OrgDb = org.Hs.eg.db, universe= GeneUniverse, minGSSize=5, pAdjustMethod="BH", qvalueCutoff=1))
goListdf_MF = lapply(goList_MF, function(x) as.data.frame(x))
# Biological Process GO enrichment
goList_BP = lapply(entrez, function(x) enrichGO(as.character(x), ont = "BP", OrgDb = org.Hs.eg.db, universe= GeneUniverse, minGSSize=5, pAdjustMethod="BH", qvalueCutoff=1))
goListdf_BP = lapply(goList_BP, function(x) as.data.frame(x))
# Cellular Compartment GO enrichment
goList_CC = lapply(entrez, function(x) enrichGO(as.character(x), ont = "CC", OrgDb = org.Hs.eg.db, universe= GeneUniverse, minGSSize=5, pAdjustMethod="BH", qvalueCutoff=1))
goListdf_CC = lapply(goList_CC, function(x) as.data.frame(x))
# Disease Ontology
goList_DO = lapply(entrez, function(x) enrichDO(as.character(x), ont = "DO", universe= GeneUniverse, minGSSize=5, pAdjustMethod="BH", qvalueCutoff=1, readable=TRUE))
goListdf_DO = lapply(goList_DO, function(x) as.data.frame(x))

# Compare the enriched terms between 7 groups

compareKegg = compareCluster(entrez, fun="enrichKEGG", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
compareBP = compareCluster(entrez, fun="enrichGO", ont = "BP", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
compareMF = compareCluster(entrez, fun="enrichGO",  ont = "MF", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
compareCC = compareCluster(entrez, fun="enrichGO",  ont = "CC", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
compareDO = compareCluster(entrez, fun="enrichDO",  ont = "DO", qvalueCutoff = 0.05, pvalueCutoff = 0.05)

# plot compared results
pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/non-CpG/figures/geneBody_mCH_KEGG_GO_DO_plots.pdf", height = 180, width = 20)
plot(compareKegg, colorBy="p.adjust", showCategory = 1000, title= "KEGG Pathway Enrichment - Gene Body mCH")
plot(compareBP, colorBy="p.adjust", showCategory = 1000, title= "Biological Process GO Enrichment - Gene Body mCH")
plot(compareMF, colorBy="p.adjust", showCategory = 1000, title= "Molecular Function GO Enrichment - Gene Body mCH")
plot(compareCC, colorBy="p.adjust", showCategory = 1000, title= "Cellular Compartment GO Enrichment - Gene Body mCH")
plot(compareDO, colorBy="p.adjust", showCategory = 1000, title= "Disease Ontology Enrichment - Gene Body mCH")
dev.off()


## GSEA
tstat = list("Cell Type" = CHgene[order(CHgene$Tstat.CellType, decreasing = T),"Tstat.CellType"],
             "Age\nin Neurons" = CHgene[order(CHgene$Tstat.AgeNeuron, decreasing = T),"Tstat.AgeNeuron"])
names(tstat$"Cell Type") = geneMap[match(CHgene[order(CHgene$Tstat.CellType, decreasing = T),"geneID"],geneMap$gencodeID),"EntrezID"]
names(tstat$"Age\nin Neurons") = geneMap[match(CHgene[order(CHgene$Tstat.AgeNeuron, decreasing = T),"geneID"],geneMap$gencodeID),"EntrezID"]
tstat = lapply(tstat, function(x) x[as.character(na.omit(unique(names(x))))])
gsea = lapply(tstat, function(x) gseGO(geneList = x, OrgDb = org.Hs.eg.db, ont = "BP", nPerm = 1000, minGSSize = 100, maxGSSize = 500, pvalueCutoff = 0.05, verbose = FALSE))

#as.data.frame(gsea[[1]][gsea[[1]][,"Description"]=="stem cell differentiation",])
#ID               Description setSize enrichmentScore
#GO:0048863 GO:0048863 stem cell differentiation     222       0.3302079
#NES      pvalue   p.adjust   qvalues rank
#GO:0048863 1.208821 0.008991009 0.03333057 0.0215605 5081
#leading_edge
#GO:0048863 tags=32%, list=24%, signal=25%
#core_enrichment
#GO:0048863 5695/652/124540/80326/3084/11281/117854/2627/2263/6688/157506/10371/6899/8829/143471/1021/3720/4780/93649/3516/6422/5788/10512/6469/2099/8091/5693/1910/10413/860/655/399694/2626/29947/6886/1909/5076/6938/23603/23081/6774/5682/4005/7476/25/6862/4851/5979/6608/1908/58495/2103/8626/4602/81606/2066/51176/51755/2668/89780/340665/182/79727/4254/55553/7042/6929/590/3815/5698/8463/2625

stcells = unlist(strsplit(gsea[[1]]$core_enrichment[which(gsea[[1]]$Description=="stem cell differentiation")], "/", fixed = T))
stcells = inGene[which(rownames(inGene) %in% geneMap[which(as.character(geneMap$EntrezID) %in% stcells),"gencodeID"]),]
stcells = reshape2::melt(stcells)
stcells$Age = pd[match(stcells$Var2, pd$Data.ID),"Age"]
stcells$Cell.Type = pd[match(stcells$Var2, pd$Data.ID),"Cell.Type"]
stcellsdt = data.table::data.table(stcells)
corr = stcellsdt[,cor.test(value,Age), by = c("Var1","Cell.Type")]
corr$padj = p.adjust(corr$p.value, method="fdr")
write.csv(cbind(Symbol = geneMap[match(corr$Var1, geneMap$gencodeID),"Symbol"], as.data.frame(corr)), quote = F, file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/non-CpG/increasing_stemCell_geneBody_mCH_overAge_correlation_tstats.pdf")

# save
save(compareKegg, compareBP, compareMF, compareCC, compareDO, keggList, goList_BP, goList_MF, goList_CC, goList_DO, keggListdf, goListdf_BP, goListdf_MF, 
     goListdf_CC, goListdf_DO, gsea, file="/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/non-CpG/geneBody_mCH_KEGG_GO_DO_objects.rda")


### Explore trinucleotide context

## Subset by trinucleotide context

CHdt = data.table(CH)
CHneuronsdt = data.table(CHneurons)
res = list("By Cell Type" = CHdt[,length(unique(regionID)), by=c("trinucleotide_context", "CT.sig")],
           "By Age" = CHdt[,length(unique(regionID)), by=c("trinucleotide_context", "Age.sig")],
           "By Age in a Cell Type" = CHdt[,length(unique(regionID)), by=c("trinucleotide_context", "Int.sig")],
           "By Age in Neurons" = CHneuronsdt[,length(unique(regionID)), by=c("trinucleotide_context", "sig")],
           "By Cell Type" = CHdt[CT.sig=="FDR < 0.05",length(unique(regionID)), by=c("trinucleotide_context", "CT.dir")],
           "By Age" = CHdt[Age.sig=="FDR < 0.05",length(unique(regionID)), by=c("trinucleotide_context", "Age.dir")],
           "By Age in a Cell Type" = CHdt[Int.sig=="FDR < 0.05",length(unique(regionID)), by=c("trinucleotide_context", "Int.dir")],
           "By Age in Neurons" = CHneuronsdt[padj<=0.05,length(unique(regionID)), by=c("trinucleotide_context", "Dir")])
for (i in 1:length(res)) {
  for (j in 1:length(as.character(unique(data.frame(res[[i]])[,"trinucleotide_context"])))) {
    res[[i]][which(res[[i]][,"trinucleotide_context"] == as.character(unique(data.frame(res[[i]])[,"trinucleotide_context"]))[j]),"perc"] = 
      round(res[[i]][which(res[[i]][,"trinucleotide_context"] == as.character(unique(data.frame(res[[i]])[,"trinucleotide_context"]))[j]),"V1"]/ 
              sum(res[[i]][which(res[[i]][,"trinucleotide_context"] == as.character(unique(data.frame(res[[i]])[,"trinucleotide_context"]))[j]),"V1"])*100,1)
    res[[i]][which(res[[i]][,"trinucleotide_context"] == as.character(unique(data.frame(res[[i]])[,"trinucleotide_context"]))[j]),"pos"] = 
      sum(res[[i]][which(res[[i]][,"trinucleotide_context"] == as.character(unique(data.frame(res[[i]])[,"trinucleotide_context"]))[j]),"V1"])
  }
}
for (i in 1:length(res)) {
  res[[i]]$trinucleotide_context = factor(res[[i]]$trinucleotide_context, levels = c("CAG", "CAC", "CTG", "CAT", "CAA", "CTC", "CCA", "CTT", "CTA", "CCT", "CCC", "CCG"))
}
for (i in 1:4) { 
  res[[i]][,"perc"] = paste0(data.frame(res[[i]])[,"perc"],"%")
  colnames(res[[i]]) = c("trinucleotide_context","sig","V1","perc","pos")
  res[[i]][which(res[[i]][,"sig"]=="FDR > 0.05"),"perc"] = ""
  res[[i]][,"sig"] = factor(data.frame(res[[i]])[,"sig"])
}
for (i in 5:8) { 
  res[[i]][,"perc"] = paste0(data.frame(res[[i]])[,"perc"],"%")
  colnames(res[[i]]) = c("trinucleotide_context","dir","V1","perc","pos")
  res[[i]][which(res[[i]][,"dir"]=="pos"),"perc"] = ""
}
res[[5]]$dir = ifelse(res[[5]]$dir=="pos", "Hypomethylated\nIn Glia", "Hypomethylated\nIn Neurons")
res[[6]]$dir = ifelse(res[[6]]$dir=="pos", "Increasingly\nMethylated", "Decreasingly\nMethylated")
res[[8]]$dir = ifelse(res[[8]]$dir=="pos", "Increasingly\nMethylated", "Decreasingly\nMethylated")


pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/non-CpG/figures/nonCG_byContext.pdf", width = 10)
for (i in 1:4) {
  g = ggplot(res[[i]], aes(x = trinucleotide_context, y = V1/1000000)) + geom_bar(aes(fill=sig),stat = "identity") +
    geom_text(aes(x = trinucleotide_context, y = pos/1000000, label = perc), vjust = -0.75) +
    scale_fill_grey() + theme_classic() +
    labs(fill="") +
    ylab("Count (Million)") + xlab("Context") +
    ggtitle(paste0("non-CpG Context: ", names(res)[i])) +
    theme(axis.text.x=element_text(angle=45,hjust=1)) +
    theme(title = element_text(size = 20)) +
    theme(text = element_text(size = 20)) +
    theme(legend.position="bottom", legend.direction="horizontal", legend.title = element_blank())
  print(g)
}
for (i in 5:8) {
  g = ggplot(res[[i]], aes(x = trinucleotide_context, y = V1/1000000)) + geom_bar(aes(fill=dir),stat = "identity") +
    geom_text(aes(x = trinucleotide_context, y = pos/1000000, label = perc), vjust = -0.75) +
    scale_fill_grey() + theme_classic() +
    labs(fill="") +
    ylab("Count (Million)") + xlab("Context") +
    ggtitle(paste0("non-CpG Context: ", names(res)[i])) +
    theme(axis.text.x=element_text(angle=45,hjust=1)) +
    theme(title = element_text(size = 20)) +
    theme(text = element_text(size = 20)) +
    theme(legend.position="bottom", legend.direction="horizontal", legend.title = element_blank())
  print(g)
}
dev.off()


res = as.data.frame(CHneuronsdt[,length(unique(regionID)), by=c("trinucleotide_context", "sig")])
colnames(res)[colnames(res)=="trinucleotide_context"] = "TC"
res$TC = as.character(res$TC)

for (j in seq(length(unique(res$TC)))) {
  res[which(res$TC == unique(res$TC)[j]),"perc"] = 
    round(res[which(res$TC == unique(res$TC)[j]),"V1"]/
            sum(res[which(res$TC == unique(res$TC)[j]),"V1"])*100,1)
  
  res[which(res$TC == unique(res$TC)[j]),"pos"] = 
    sum(res[which(res$TC == unique(res$TC)[j]),"V1"])
}
res
res$TC = factor(res$TC, levels = c("CAG", "CAC", "CTG", "CAT", "CAA", "CTC", 
                                   "CCA", "CTT", "CTA", "CCT", "CCC", "CCG"))

res$perc = paste0(res$perc,"%")
res[which(res$sig=="FDR > 0.05"),"perc"] = ""
res$sig = factor(res$sig)
res$perc1 = c("", "", "", "", "", "", "", "", "12%", "", "", "", "",
              "", "", "", "", "", "", "", "", "", "", "")
res$perc2 = c("", "", "", "11.6%", "30.6%", "", "", "1.9%", "", "", "", "", "",
              "18.8%", "", "", "2%", "", "9%", "0.5%", "0.6%", "0.9%", "3.4%", "1.8%")


pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/non-CpG/figures/remake_CH_byContext.pdf", 
    width = 5.75, height = 4.25)
ggplot(res, aes(x = TC, y = V1/1000000)) + 
  geom_bar(aes(fill=sig),stat = "identity") +
  geom_text(aes(x = TC, y = pos/1000000, label = perc2), vjust = -0.75) +
  geom_text(aes(x = TC, y = pos/1000000, label = perc1), vjust = 1.6, color = "white") +
  scale_fill_grey() + theme_bw() +
  labs(fill="") +
  ylab("Count (Million)") + xlab("Context") +
  ggtitle("CpH Context: By Age in Neurons") +
  theme(axis.text.x = element_text(angle=25,hjust=0.75),
        title = element_text(size = 20),
        text = element_text(size = 20),
        legend.position=c(0.8, 0.9), 
        legend.title = element_blank(),
        legend.background = element_blank())

dev.off()


## Is there a difference in the proportion of CAG vs CAC in C's more methylated in glia vs neurons and vice versa?

table(CHdt$CT.sig,CHdt$c_context)
#                   CG      CHG      CHH
#  FDR < 0.05        0  2284532  5397543
#  FDR > 0.05        0 12744000 20392667
fisher.test(data.frame(c(2284532,12744000),c(5397543,20392667)))
#p-value < 2.2e-16
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.6761427 0.6784939
#sample estimates:
#  odds ratio 
#0.6772835 


table(data.frame(CHdt[CT.sig=="FDR < 0.05",,])[,"CT.dir"],data.frame(CHdt[CT.sig=="FDR < 0.05",,])[,"c_context"])
#         CG     CHG     CHH
#neg       0   36439   42800
#pos       0 2248093 5354743
fisher.test(data.frame(c(36439,2248093),c(42800,5354743)))
#p-value < 2.2e-16
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  1.999487 2.056607
#sample estimates:
#  odds ratio 
#2.02797 
# mCH that is enriched in glia is more likely to be CHG than CHH

table(data.frame(CHdt[CT.sig=="FDR < 0.05" & trinucleotide_context %in% c("CAG","CAC"),,])[,"CT.dir"],data.frame(CHdt[CT.sig=="FDR < 0.05" & trinucleotide_context %in% c("CAG","CAC"),,])[,"c_context"])
#         CG     CHG     CHH
#neg       0   12500    3454
#pos       0 2072603 2365368
fisher.test(data.frame(c(12500,2072603),c(3454,2365368)))
#p-value < 2.2e-16
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  3.976673 4.290161
#sample estimates:
#  odds ratio 
#4.130595 


## Are the GO terms for glial mCH-enriched genes affected by trinucleotide context?

CTentrez = list("Glial Genes\nCHG" = CHdt[CT.sig=="FDR < 0.05" & CT.dir=="neg" & c_context=="CHG",list(na.omit(EntrezID)),], 
                "Glial Genes\nCHH" = CHdt[CT.sig=="FDR < 0.05" & CT.dir=="neg" & c_context=="CHH",list(na.omit(EntrezID)),],
                "Neuronal Genes\nCHG" = CHdt[CT.sig=="FDR < 0.05" & CT.dir=="pos" & c_context=="CHG",list(na.omit(EntrezID)),], 
                "Neuronal Genes\nCHH" = CHdt[CT.sig=="FDR < 0.05" & CT.dir=="pos" & c_context=="CHH",list(na.omit(EntrezID)),],
                "Glial\nOverlapping Genes\nCHG" = CHdt[CT.sig=="FDR < 0.05" & CT.dir=="neg" & c_context=="CHG" & distToGene==0,list(na.omit(EntrezID)),], 
                "Glial\nOverlapping Genes\nCHH" = CHdt[CT.sig=="FDR < 0.05" & CT.dir=="neg" & c_context=="CHH" & distToGene==0,list(na.omit(EntrezID)),],
                "Neuronal\nOverlapping Genes\nCHG" = CHdt[CT.sig=="FDR < 0.05" & CT.dir=="pos" & c_context=="CHG" & distToGene==0,list(na.omit(EntrezID)),], 
                "Neuronal\nOverlapping Genes\nCHH" = CHdt[CT.sig=="FDR < 0.05" & CT.dir=="pos" & c_context=="CHH" & distToGene==0,list(na.omit(EntrezID)),])
CTentrez = lapply(CTentrez, function(x) as.character(unique(x$V1)))              

# Compare the enriched terms between 7 groups
CTcompareKegg = compareCluster(CTentrez, fun="enrichKEGG", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
CTcompareBP = compareCluster(CTentrez, fun="enrichGO", ont = "BP", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
CTcompareMF = compareCluster(CTentrez, fun="enrichGO",  ont = "MF", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
CTcompareCC = compareCluster(CTentrez, fun="enrichGO",  ont = "CC", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
CTcompareDO = compareCluster(CTentrez, fun="enrichDO",  ont = "DO", qvalueCutoff = 0.05, pvalueCutoff = 0.05)

# save, write to csv
save(CTcompareKegg, CTcompareMF, CTcompareCC, CTcompareDO, CTcompareBP,
     file="/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/non-CpG/nonCG_KEGG_GO_DO_objects_byCellType_byC_Context.rda")

# plot compared results
pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/non-CpG/figures/nonCG_KEGG_GO_DO_plots_byCellType_byC_Context.pdf", height = 250, width = 25)
plot(CTcompareKegg, colorBy="p.adjust", showCategory = 1100, title= "non-CpG KEGG Pathway Enrichment by Cell Type")
plot(CTcompareBP, colorBy="p.adjust", showCategory = 1100, title= "non-CpG Biological Process GO Enrichment by Cell Type")
plot(CTcompareMF, colorBy="p.adjust", showCategory = 1100, title= "non-CpG Molecular Function GO Enrichment by Cell Type")
plot(CTcompareCC, colorBy="p.adjust", showCategory = 1100, title= "non-CpG Cellular Compartment GO Enrichment by Cell Type")
plot(CTcompareDO, colorBy="p.adjust", showCategory = 1100, title= "non-CpG Disease Ontology Enrichment by Cell Type")
dev.off()


## Compare Cell Type-associated nonCGs by C context and overlap with genomic features

load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/UMRs_LMRs_annotated.rda")
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/dmvs_annotated.rda")
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/PMDs_annotated.rda")

## collate regions
CREs = list(LMR = lapply(lDMR, makeGRangesFromDataFrame), UMR = lapply(uDMR, makeGRangesFromDataFrame),
            DMV = lapply(dmvDMR, makeGRangesFromDataFrame), PMD = lapply(pmdDMR, makeGRangesFromDataFrame))
CREs = lapply(CREs, function(x) reduce(do.call(getMethod(c, "GenomicRanges"), GRangesList(x))))
rm(lDMR,uDMR,dmvDMR, pmdDMR)

oo = lapply(CREs, function(y) findOverlaps(makeGRangesFromDataFrame(CHdt), y))
CHdt$rnum = 1:nrow(CHdt)
CHdt$UMR = ifelse(CHdt$rnum %in% queryHits(oo$UMR), "UMR", "no")
CHdt$LMR = ifelse(CHdt$rnum %in% queryHits(oo$LMR), "LMR", "no")
CHdt$DMV = ifelse(CHdt$rnum %in% queryHits(oo$DMV), "DMV", "no")
CHdt$PMD = ifelse(CHdt$rnum %in% queryHits(oo$PMD), "PMD", "no")
save(CHdt, geneMap, pd, file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/non-CpG/CH_object.rda")


CTtables = list(islands = c(table(CHdt$CT.sig, CHdt$islands=="CpG_Island", CHdt$c_context), table(data.frame(CHdt[CT.sig=="FDR < 0.05",,])[,"CT.dir"],
                                data.frame(CHdt[CT.sig=="FDR < 0.05",,])[,"islands"]=="CpG_Island",data.frame(CHdt[CT.sig=="FDR < 0.05",,])[,"c_context"])),
                genes = c(table(CHdt$CT.sig, CHdt$distToGene==0,CHdt$c_context), table(data.frame(CHdt[CT.sig=="FDR < 0.05",,])[,"CT.dir"],
                              data.frame(CHdt[CT.sig=="FDR < 0.05",,])[,"distToGene"]==0,data.frame(CHdt[CT.sig=="FDR < 0.05",,])[,"c_context"])),
                UMR = c(table(CHdt$CT.sig, CHdt$UMR=="UMR",CHdt$c_context), table(data.frame(CHdt[CT.sig=="FDR < 0.05",,])[,"CT.dir"],
                              data.frame(CHdt[CT.sig=="FDR < 0.05",,])[,"UMR"]=="UMR",data.frame(CHdt[CT.sig=="FDR < 0.05",,])[,"c_context"])),
                LMR = c(table(CHdt$CT.sig, CHdt$LMR=="LMR", CHdt$c_context), table(data.frame(CHdt[CT.sig=="FDR < 0.05",,])[,"CT.dir"],
                              data.frame(CHdt[CT.sig=="FDR < 0.05",,])[,"LMR"]=="LMR",data.frame(CHdt[CT.sig=="FDR < 0.05",,])[,"c_context"])),
                DMV = c(table(CHdt$CT.sig, CHdt$DMV=="DMV", CHdt$c_context), table(data.frame(CHdt[CT.sig=="FDR < 0.05",,])[,"CT.dir"],
                              data.frame(CHdt[CT.sig=="FDR < 0.05",,])[,"DMV"]=="DMV",data.frame(CHdt[CT.sig=="FDR < 0.05",,])[,"c_context"])),
                PMD = c(table(CHdt$CT.sig, CHdt$PMD=="PMD", CHdt$c_context), table(data.frame(CHdt[CT.sig=="FDR < 0.05",,])[,"CT.dir"],
                              data.frame(CHdt[CT.sig=="FDR < 0.05",,])[,"PMD"]=="PMD",data.frame(CHdt[CT.sig=="FDR < 0.05",,])[,"c_context"])))

tab = unlist(lapply(CTtables, function(x) list(Sig_feature_CHG = data.frame(true = x[7:8], false = x[5:6], row.names = c("FDR < 0.05", "FDR > 0.05")), 
                                        Sig_feature_CHH = data.frame(true = x[11:12], false = x[9:10], row.names = c("FDR < 0.05", "FDR > 0.05")),
                                        forSig_Dir_feature_CHG = data.frame(true = x[19:20], false = x[17:18], row.names = c("neg","pos")), 
                                        forSig_Dir_feature_CHH = data.frame(true = x[23:24], false = x[21:22], row.names = c("neg","pos")))), recursive=F)
fisher = lapply(tab, fisher.test)
write.csv(data.frame(Test = names(unlist(lapply(fisher, function(x) x$p.value))), pval = unlist(lapply(fisher, function(x) x$p.value)), 
                     OR = unlist(lapply(fisher, function(x) x$estimate))), quote=F,
          file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/non-CpG/CHG_vs_CHH_feature_enrichment_fishertest.csv")


### Is there a difference in the proportion of CAG vs CAC in C's more methylated in older neurons vs younger neurons and vice versa?

table(CHneuronsdt$sig,CHneuronsdt$c_context)
#                   CG      CHG      CHH
#  FDR < 0.05        0  1135137  2885234
#  FDR > 0.05        0 12473649 19605564
fisher.test(data.frame(c(1135137,12473649),c(2885234,19605564)))
#p-value < 2.2e-16
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.6169567 0.6197734
#sample estimates:
#  odds ratio 
#0.6183777 


table(data.frame(CHneuronsdt[sig=="FDR < 0.05",,])[,"Dir"],data.frame(CHneuronsdt[sig=="FDR < 0.05",,])[,"c_context"])
#           CG     CHG     CHH
#  neg       0   14567   15800
#  pos       0 1120570 2869434
fisher.test(data.frame(c(14567,1120570),c(15800,2869434)))
#p-value < 2.2e-16
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  2.307906 2.415169
#sample estimates:
#  odds ratio 
#2.360879 
## mCH that is enriched in younger neurons is more likely to be CHG than CHH

table(data.frame(CHneuronsdt[sig=="FDR < 0.05" & trinucleotide_context %in% c("CAG","CAC"),,])[,"Dir"],
      data.frame(CHneuronsdt[sig=="FDR < 0.05" & trinucleotide_context %in% c("CAG","CAC"),,])[,"c_context"])
#         CG     CHG     CHH
#neg       0    5901    1338
#pos       0 1045294 1394793
fisher.test(data.frame(c(5901,1045294),c(1338,1394793)))
#p-value < 2.2e-16
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  5.545168 6.250295
#sample estimates:
#  odds ratio 
#5.88395 



## Are the GO terms for younger mCH-enriched genes affected by trinucleotide context?

Ageentrez = list("Younger Genes\nCHG" = CHneuronsdt[sig=="FDR < 0.05" & Dir=="neg" & c_context=="CHG",list(na.omit(EntrezID)),], 
                "Younger Genes\nCHH" = CHneuronsdt[sig=="FDR < 0.05" & Dir=="neg" & c_context=="CHH",list(na.omit(EntrezID)),],
                "Older Genes\nCHG" = CHneuronsdt[sig=="FDR < 0.05" & Dir=="pos" & c_context=="CHG",list(na.omit(EntrezID)),], 
                "Older Genes\nCHH" = CHneuronsdt[sig=="FDR < 0.05" & Dir=="pos" & c_context=="CHH",list(na.omit(EntrezID)),],
                "Younger\nOverlapping Genes\nCHG" = CHneuronsdt[sig=="FDR < 0.05" & Dir=="neg" & c_context=="CHG" & distToGene==0,list(na.omit(EntrezID)),], 
                "Younger\nOverlapping Genes\nCHH" = CHneuronsdt[sig=="FDR < 0.05" & Dir=="neg" & c_context=="CHH" & distToGene==0,list(na.omit(EntrezID)),],
                "Older\nOverlapping Genes\nCHG" = CHneuronsdt[sig=="FDR < 0.05" & Dir=="pos" & c_context=="CHG" & distToGene==0,list(na.omit(EntrezID)),], 
                "Older\nOverlapping Genes\nCHH" = CHneuronsdt[sig=="FDR < 0.05" & Dir=="pos" & c_context=="CHH" & distToGene==0,list(na.omit(EntrezID)),])
Ageentrez = lapply(Ageentrez, function(x) as.character(unique(x$V1)))              

# Compare the enriched terms between 7 groups
AgecompareKegg = compareCluster(Ageentrez, fun="enrichKEGG", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
AgecompareBP = compareCluster(Ageentrez, fun="enrichGO", ont = "BP", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
AgecompareMF = compareCluster(Ageentrez, fun="enrichGO",  ont = "MF", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
AgecompareCC = compareCluster(Ageentrez, fun="enrichGO",  ont = "CC", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
AgecompareDO = compareCluster(Ageentrez, fun="enrichDO",  ont = "DO", qvalueCutoff = 0.05, pvalueCutoff = 0.05)

# save
save(AgecompareKegg, AgecompareMF, AgecompareCC, AgecompareDO, AgecompareBP,
     file="/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/non-CpG/nonCG_KEGG_GO_DO_objects_byAgeinNeurons_byC_Context.rda")

# plot compared results
pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/non-CpG/figures/nonCG_KEGG_GO_DO_plots_byAgeinNeurons_byC_Context.pdf", height = 280, width = 25)
plot(AgecompareKegg, colorBy="p.adjust", showCategory = 1400, title= "non-CpG KEGG Pathway Enrichment by Age in Neurons")
plot(AgecompareBP, colorBy="p.adjust", showCategory = 1400, title= "non-CpG Biological Process GO Enrichment by Age in Neurons")
plot(AgecompareMF, colorBy="p.adjust", showCategory = 1400, title= "non-CpG Molecular Function GO Enrichment by Age in Neurons")
plot(AgecompareCC, colorBy="p.adjust", showCategory = 1400, title= "non-CpG Cellular Compartment GO Enrichment by Age in Neurons")
plot(AgecompareDO, colorBy="p.adjust", showCategory = 1400, title= "non-CpG Disease Ontology Enrichment by Age in Neurons")
dev.off()


comp = as.data.frame(AgecompareBP)
comp = comp[grep("Overlap", comp$Cluster),]
comp = split(comp$Description, comp$Cluster)
comp = comp[elementNROWS(comp)>0]
ov = calculate.overlap(comp)

plotExample = AgecompareBP # clusterProfiler output
plotExample@compareClusterResult = plotExample@compareClusterResult[which(plotExample@compareClusterResult$Description %in% 
                                                                            unlist(ov[names(ov) %in% c("a9","a14","a1","a3")])),]   
plotExample@compareClusterResult = plotExample@compareClusterResult[grep("Overlap",plotExample@compareClusterResult$Cluster),] 

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/non-CpG/figures/nonCG_BP_filtered_plots_byAgeinNeurons_byC_Context.pdf", height = 14, width = 12)
plot(plotExample, colorBy="p.adjust", showCategory = 10, title= "non-CpG BP Enrichment by Age in Neurons")
dev.off()


## Limit to CAG and CAC context only: GO terms for younger mCH-enriched genes affected by trinucleotide context?

Ageentrez = list("Decreasing\nmCAG" = CHneuronsdt[sig=="FDR < 0.05" & Dir=="neg" & trinucleotide_context=="CAG",list(na.omit(EntrezID)),], 
                 "Decreasing\nmCAC" = CHneuronsdt[sig=="FDR < 0.05" & Dir=="neg" & trinucleotide_context=="CAC",list(na.omit(EntrezID)),],
                 "Increasing\nmCAG" = CHneuronsdt[sig=="FDR < 0.05" & Dir=="pos" & trinucleotide_context=="CAG",list(na.omit(EntrezID)),], 
                 "Increasing\nmCAC" = CHneuronsdt[sig=="FDR < 0.05" & Dir=="pos" & trinucleotide_context=="CAC",list(na.omit(EntrezID)),],
                 "Decreasing\n(Overlapping Genes)\nmCAG" = CHneuronsdt[sig=="FDR < 0.05" & Dir=="neg" & trinucleotide_context=="CAG" & distToGene==0,list(na.omit(EntrezID)),], 
                 "Decreasing\n(Overlapping Genes)\nmCAC" = CHneuronsdt[sig=="FDR < 0.05" & Dir=="neg" & trinucleotide_context=="CAC" & distToGene==0,list(na.omit(EntrezID)),],
                 "Increasing\n(Overlapping Genes)\nmCAG" = CHneuronsdt[sig=="FDR < 0.05" & Dir=="pos" & trinucleotide_context=="CAG" & distToGene==0,list(na.omit(EntrezID)),], 
                 "Increasing\n(Overlapping Genes)\nmCAC" = CHneuronsdt[sig=="FDR < 0.05" & Dir=="pos" & trinucleotide_context=="CAC" & distToGene==0,list(na.omit(EntrezID)),])
Ageentrez = lapply(Ageentrez, function(x) as.character(unique(x$V1)))              

# Compare the enriched terms between 7 groups
AgecompareKegg = compareCluster(Ageentrez, fun="enrichKEGG", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
AgecompareBP = compareCluster(Ageentrez, fun="enrichGO", ont = "BP", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
AgecompareMF = compareCluster(Ageentrez, fun="enrichGO",  ont = "MF", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
AgecompareCC = compareCluster(Ageentrez, fun="enrichGO",  ont = "CC", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
AgecompareDO = compareCluster(Ageentrez, fun="enrichDO",  ont = "DO", qvalueCutoff = 0.05, pvalueCutoff = 0.05)

# save
save(AgecompareKegg, AgecompareMF, AgecompareCC, AgecompareDO, AgecompareBP,
     file="/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/non-CpG/nonCG_KEGG_GO_DO_objects_byAgeinNeurons_byCAG_orCAC.rda")

# plot compared results
pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/non-CpG/figures/nonCG_KEGG_GO_DO_plots_byAgeinNeurons_byCAG_orCAC.pdf", height = 80, width = 25)
dotplot(AgecompareKegg, showCategory = 1400, title= "non-CpG KEGG Pathway Enrichment by Age in Neurons")
dotplot(AgecompareMF, showCategory = 1400, title= "non-CpG Molecular Function GO Enrichment by Age in Neurons")
dotplot(AgecompareCC, showCategory = 1400, title= "non-CpG Cellular Compartment GO Enrichment by Age in Neurons")
dotplot(AgecompareDO, showCategory = 1400, title= "non-CpG Disease Ontology Enrichment by Age in Neurons")
dev.off()


comp = as.data.frame(AgecompareBP)
comp = comp[grep("Overlap", comp$Cluster),]
comp = split(comp$Description, comp$Cluster)
comp = comp[elementNROWS(comp)>0]
ov = calculate.overlap(comp)
elementNROWS(ov)

plotExample = AgecompareBP # clusterProfiler output
plotExample@compareClusterResult = plotExample@compareClusterResult[which(plotExample@compareClusterResult$Description %in% 
                                                                            unlist(ov[names(ov) %in% c("a9","a14","a1","a3")])),]  
plotExample@compareClusterResult = plotExample@compareClusterResult[grep("Overlap",plotExample@compareClusterResult$Cluster),]
plotExample@compareClusterResult$Description = gsub("positive regulation of adaptive immune response based on somatic recombination of immune receptors built from immunoglobulin superfamily domains",
                                                    "positive regulation of adaptive immune response based on\nsomatic recombination of immune receptors built from immunoglobulin superfamily domains",plotExample@compareClusterResult$Description)

plotExample2 = AgecompareBP
plotExample2@compareClusterResult = plotExample2@compareClusterResult[grep("Overlap",plotExample2@compareClusterResult$Cluster),]

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/non-CpG/figures/nonCG_BP_filtered_plots_byAgeinNeurons_byC_Context.pdf", height = 14, width = 14)
dotplot(plotExample, showCategory = 10, title= "non-CpG BP Enrichment by Age in Neurons")
dotplot(plotExample2, showCategory = 10, title= "non-CpG BP Enrichment by Age in Neurons")
dev.off()


## now in only the genes that are either CAG or CAC but not both

ent = list("Decreasing\nmCAG" = Ageentrez$"Decreasing\n(Overlapping Genes)\nmCAG"[-which(Ageentrez$"Decreasing\n(Overlapping Genes)\nmCAG" %in% 
                                                                                           Ageentrez$"Decreasing\n(Overlapping Genes)\nmCAC")], 
           "Decreasing\nmCAC" = Ageentrez$"Decreasing\n(Overlapping Genes)\nmCAC"[-which(Ageentrez$"Decreasing\n(Overlapping Genes)\nmCAC" %in% 
                                                                                           Ageentrez$"Decreasing\n(Overlapping Genes)\nmCAG")],
           "Increasing\nmCAG" = Ageentrez$"Increasing\n(Overlapping Genes)\nmCAG"[-which(Ageentrez$"Increasing\n(Overlapping Genes)\nmCAG" %in% 
                                                                                           Ageentrez$"Increasing\n(Overlapping Genes)\nmCAC")],
           "Increasing\nmCAC" = Ageentrez$"Increasing\n(Overlapping Genes)\nmCAC"[-which(Ageentrez$"Increasing\n(Overlapping Genes)\nmCAC" %in% 
                                                                                           Ageentrez$"Increasing\n(Overlapping Genes)\nmCAG")])
AgecompareBP = compareCluster(ent, fun="enrichGO", ont = "BP", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
AgecompareBP = simplify(AgecompareBP)

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/non-CpG/figures/nonCG_BP_plots_byAgeinNeurons_CAG_or_CAC_only.pdf", 
    height = 12, width = 10)
dotplot(AgecompareBP, showCategory = 20, 
        title= "CpH Biological Process Enrichment:\nBy Age in Neurons")
dev.off()

elementNROWS(ent)
# Decreasing\nmCAG Decreasing\nmCAC Increasing\nmCAG Increasing\nmCAC 
# 1406              338              870             2416 


## Compare Age-associated nonCGs in neurons by C context and overlap with genomic features

oo = lapply(CREs, function(y) findOverlaps(makeGRangesFromDataFrame(CHneuronsdt), y))
CHneuronsdt$rnum = 1:nrow(CHneuronsdt)
CHneuronsdt$UMR = ifelse(CHneuronsdt$rnum %in% queryHits(oo$UMR), "UMR", "no")
CHneuronsdt$LMR = ifelse(CHneuronsdt$rnum %in% queryHits(oo$LMR), "LMR", "no")
CHneuronsdt$DMV = ifelse(CHneuronsdt$rnum %in% queryHits(oo$DMV), "DMV", "no")
CHneuronsdt$PMD = ifelse(CHneuronsdt$rnum %in% queryHits(oo$PMD), "PMD", "no")
save(CHneuronsdt, geneMap, pd, file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/non-CpG/CHneurons_object.rda")


Agetables = list(islands = c(table(CHneuronsdt$sig, CHneuronsdt$islands=="CpG_Island", CHneuronsdt$c_context), table(data.frame(CHneuronsdt[sig=="FDR < 0.05",,])[,"Dir"],
                             data.frame(CHneuronsdt[sig=="FDR < 0.05",,])[,"islands"]=="CpG_Island",data.frame(CHneuronsdt[sig=="FDR < 0.05",,])[,"c_context"])),
                genes = c(table(CHneuronsdt$sig, CHneuronsdt$distToGene==0,CHneuronsdt$c_context), table(data.frame(CHneuronsdt[sig=="FDR < 0.05",,])[,"Dir"],
                          data.frame(CHneuronsdt[sig=="FDR < 0.05",,])[,"distToGene"]==0,data.frame(CHneuronsdt[sig=="FDR < 0.05",,])[,"c_context"])),
                UMR = c(table(CHneuronsdt$sig, CHneuronsdt$UMR=="UMR",CHneuronsdt$c_context), table(data.frame(CHneuronsdt[sig=="FDR < 0.05",,])[,"Dir"],
                        data.frame(CHneuronsdt[sig=="FDR < 0.05",,])[,"UMR"]=="UMR",data.frame(CHneuronsdt[sig=="FDR < 0.05",,])[,"c_context"])),
                LMR = c(table(CHneuronsdt$sig, CHneuronsdt$LMR=="LMR", CHneuronsdt$c_context), table(data.frame(CHneuronsdt[sig=="FDR < 0.05",,])[,"Dir"],
                        data.frame(CHneuronsdt[sig=="FDR < 0.05",,])[,"LMR"]=="LMR",data.frame(CHneuronsdt[sig=="FDR < 0.05",,])[,"c_context"])),
                DMV = c(table(CHneuronsdt$sig, CHneuronsdt$DMV=="DMV", CHneuronsdt$c_context), table(data.frame(CHneuronsdt[sig=="FDR < 0.05",,])[,"Dir"],
                        data.frame(CHneuronsdt[sig=="FDR < 0.05",,])[,"DMV"]=="DMV",data.frame(CHneuronsdt[sig=="FDR < 0.05",,])[,"c_context"])),
                PMD = c(table(CHneuronsdt$sig, CHneuronsdt$PMD=="PMD", CHneuronsdt$c_context), table(data.frame(CHneuronsdt[sig=="FDR < 0.05",,])[,"Dir"],
                        data.frame(CHneuronsdt[sig=="FDR < 0.05",,])[,"PMD"]=="PMD",data.frame(CHneuronsdt[sig=="FDR < 0.05",,])[,"c_context"])))

tab = unlist(lapply(Agetables, function(x) list(Sig_feature_CHG = data.frame(true = x[7:8], false = x[5:6], row.names = c("FDR < 0.05", "FDR > 0.05")), 
                                               Sig_feature_CHH = data.frame(true = x[11:12], false = x[9:10], row.names = c("FDR < 0.05", "FDR > 0.05")),
                                               forSig_Dir_feature_CHG = data.frame(true = x[19:20], false = x[17:18], row.names = c("neg","pos")), 
                                               forSig_Dir_feature_CHH = data.frame(true = x[23:24], false = x[21:22], row.names = c("neg","pos")))), recursive=F)
fisher = lapply(tab, fisher.test)
write.csv(data.frame(Test = names(unlist(lapply(fisher, function(x) x$p.value))), pval = unlist(lapply(fisher, function(x) x$p.value)), 
                     OR = unlist(lapply(fisher, function(x) x$estimate))), quote=F,
          file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/non-CpG/CHG_vs_CHH_feature_enrichment_fishertest_byAgeinNeurons.csv")

