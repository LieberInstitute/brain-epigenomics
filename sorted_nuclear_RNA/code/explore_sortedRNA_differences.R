library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(cowplot)


load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/sorted_nuclear_RNA/DE_limma_results_objects.rda")


## Check for appropriate enrichment of marker genes: by cell type

sigCT = nucRNAres[which(nucRNAres$padj.CellTypeNeuron<=0.05),]
dim(nucRNAres) # 46367
dim(sigCT) # 9994

genes = c("GFAP" = "Glia","ALDH1L1" = "Glia","OLIG2" = "Glia","AQP4" = "Glia",
          "SLC17A7" = "Neuron","SNAP25" = "Neuron","RELN" = "Neuron","TBR1" = "Neuron","GAD1" = "Neuron")
markers = nucRNAres[match(geneMap[which(geneMap$Symbol %in% names(genes)),"gencodeID"], nucRNAres$gencodeID),] 
markers$Sym = geneMap[match(markers$gencodeID,geneMap$gencodeID),"Symbol"]
markers$marker = genes[match(markers$Sym, names(genes))]
markers$Sym = factor(markers$Sym, levels = c("GFAP", "ALDH1L1", "OLIG2", "AQP4", "SLC17A7", "SNAP25", "TBR1", "GAD1", "RELN"))

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/sorted_nuclear_RNA/figures/cellType_markerGenes_tstats.pdf", width =4, height=3.75)
ggplot(markers, aes(x=Sym, y =Tstat.CellTypeNeuron)) + scale_fill_brewer(8, palette="Dark2") + 
  geom_bar(stat = 'identity', aes(fill = marker), position = 'dodge', col = 'transparent') +
  ylab("T Statistic") + xlab("") +
  ggtitle("Cell Type Marker Genes") + 
  theme(title = element_text(size = 20),
        text = element_text(size = 20), 
        legend.title=element_blank(),axis.text.x = element_text(angle = 30, hjust = 1),
        legend.position = c(0.6, 0.2))
dev.off()


# Compare the enriched terms

nucRNAres$entrez = geneMap[match(nucRNAres$gencodeID, geneMap$gencodeID), "EntrezID"]
entrez = list(Neuronal = nucRNAres$entrez[which(nucRNAres$Tstat.CellTypeNeuron>0 & nucRNAres$padj.CellTypeNeuron<=0.05)],
              Glial = nucRNAres$entrez[which(nucRNAres$Tstat.CellTypeNeuron<0 & nucRNAres$padj.CellTypeNeuron<=0.05)])
entrez = lapply(entrez, function(x) na.omit(as.character(unique(x))))

compareKegg = compareCluster(entrez, fun="enrichKEGG", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
compareBP = compareCluster(entrez, fun="enrichGO", ont = "BP", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
compareMF = compareCluster(entrez, fun="enrichGO",  ont = "MF", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
compareCC = compareCluster(entrez, fun="enrichGO",  ont = "CC", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
compareDO = compareCluster(entrez, fun="enrichDO",  ont = "DO", qvalueCutoff = 0.05, pvalueCutoff = 0.05)

# Save
save(compareBP, compareMF, compareCC, compareDO, compareKegg,
     file="/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/sorted_nuclear_RNA/GO.objects.sorted_RNA_byCellType.rda")

## plot Gene ontology enrichment

bp = as.data.frame(compareBP)
head(bp$Description,30)
toplot = rbind(bp[1:4,], bp[which(bp$Description %in% 
                                    c("gliogenesis","glial cell differentiation","regulation of I-kappaB kinase/NF-kappaB signaling","regulation of leukocyte mediated immunity")),])
toplot$log10 = -log10(toplot$p.adjust)
toplot = toplot[,colnames(toplot) !="geneID"]
toplot$Description = factor(c("Modulation of chemical synaptic transmission",
                              "Regulation of synaptic plasticity",
                              "Regulation of membrane potential",
                              "Regulation of neuron projection development",
                              "Gliogenesis","Glial cell differentiation",
                              "Regulation of\nI-kappaB kinase/NF-kappaB signaling",
                              "Regulation of leukocyte mediated immunity"), 
                            levels = c("Modulation of chemical synaptic transmission",
                                       "Regulation of synaptic plasticity",
                                       "Regulation of membrane potential", 
                                       "Regulation of neuron projection development",
                                       "Gliogenesis","Glial cell differentiation",
                                       "Regulation of\nI-kappaB kinase/NF-kappaB signaling", 
                                       "Regulation of leukocyte mediated immunity"))
toplot$Cluster = factor(toplot$Cluster, levels = c("Glial", "Neuronal"))


pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/sorted_nuclear_RNA/figures/BP_byCellType_sortedNuclearRNA.pdf", width=9.5,height=5.25)
ggplot(toplot, aes(x=Description, y =log10)) + scale_fill_brewer(8, palette="Dark2") + 
  geom_bar(stat = 'identity', aes(fill = Cluster), position = 'dodge', col = 'transparent') +
  geom_text(aes(label=GeneRatio), hjust=1.1, color="black", position = position_dodge(1), size=5) +
  coord_flip() +
  ylab("-log10(FDR)") + 
  xlab("") +
  ggtitle("Biological Process Enrichment:\nBy Cell Type in Nuclear RNA") + 
  theme(title = element_text(size = 20),
        text = element_text(size = 20), 
        legend.title=element_blank(),
        legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent", color = "transparent"),
        legend.position=c(0.8, 0.8))
dev.off()




