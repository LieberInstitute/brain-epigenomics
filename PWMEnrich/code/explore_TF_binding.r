library(GenomicRanges)
library(clusterProfiler)
library(org.Hs.eg.db)
library(BSgenome.Hsapiens.UCSC.hg19)
library(PWMEnrich.Hsapiens.background)
library("VennDiagram")
library(ggplot2)

load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/UMRs_LMRs_annotated.rda")
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/DMR/DMR_objects.rda")
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/PWMEnrich/DMR_PWMEnrich_objects.rda")
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/PWMEnrich/UMR_LMR_PWMEnrich_object.rda")
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/PWMEnrich/filtered_TF_list_expressionInfo.rda")

write.csv(pd, quote=F, file="/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/postnatal.sorted.WGBS.pheno.info.csv")
write.csv(hompd, quote=F, file="/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/homogenate.RNAseq.pheno.info.csv")


## Objects for analysis:
  # all_int: PWMEnrich result for all interaction DMRs
  # intergenic_int: PWMEnrich result for interaction DMRs that overlap intergenic regions
  # introns_int: PWMEnrich result for interaction DMRs that overlap introns
  # promoters_int: PWMEnrich result for interaction DMRs that overlap promoters
  # LMR_UMR_pwmenrich: PWMEnrich result for shared and not shared LMR and UMR sequence for All, Prenatal, Postnatal, Neurons, and Glia 
  # all_split: PWMEnrich result for all regions in three models split by direction
  # intergenic_split: PWMEnrich result for intergenic regions in three models split by direction
  # introns_split: PWMEnrich result for introns in three models split by direction
  # promoters_split: PWMEnrich result for promoters in three models split by direction
  # DMR, uDMR, lDMR: lists of annotated methylated regions 
  # geneMap
  # hompd, pd: phenotype tables
  # targettogeneID: map of target ID to gencodeID
  # TFhomRPKM, TFnucres: RNA-seq results for TFs in homogenate (hom) and nuclear RNA (nuc)
  # TFdiff: difference in test statistic between positive and negative beta values of DMRs (all, intergenic, introns, and promoters)
  # TFnucres: DE results by cell type for TF genes in nuclear RNA
  # ulTFdiff: difference in test statistic between UMRs and LMRs, and between cell types/ages in UMRs and LMRs


## collate TF group report results

TF = list(allTF = lapply(all_split, groupReport), promTF = lapply(promoters_split, groupReport),
          intronTF = lapply(introns_split, groupReport), intergenicTF = lapply(intergenic_split, groupReport))
TF = lapply(TF, function(x) lapply(x, as.data.frame))
TF = lapply(TF, function(t) Map(cbind, t, padj = lapply(t, function(x) p.adjust(x$p.value, method = "fdr"))))
TF = lapply(TF, function(t) lapply(t, function(x) x[which(x$target %in% names(targettogeneID)),]))
TF = lapply(TF, function(t) lapply(t, function(x) x[order(x$id),]))
TF = lapply(TF, function(t) Map(cbind, t, transf = lapply(t, function(x) -log10(x$padj))))
TF = lapply(TF, function(t) t[order(names(t))])

ulTF = lapply(LMR_UMR_pwmenrich, groupReport)
ulTF = lapply(ulTF, as.data.frame)
ulTF = Map(cbind, ulTF, padj = lapply(ulTF, function(x) p.adjust(x$p.value, method = "fdr")))
ulTF = lapply(ulTF, function(x) x[which(x$target %in% names(targettogeneID)),])
ulTF = Map(cbind, ulTF, transf = lapply(ulTF, function(x) -log10(x$padj)))
ulTF = lapply(ulTF, function(x) x[order(x$id),])
ulTF = ulTF[order(names(ulTF))]


### Are different genomic features enriched for different TFs?

## for DMRs by feature
# Correlate TF enrichment between features

df = list(Prom.Intron = mapply(function(prom,intn) if (length(prom$transf) > 0 & length(intn$transf) > 0) {
                               data.frame(Promoter = prom$transf, Intron = intn$transf) }, TF$promTF, TF$intronTF, SIMPLIFY = F),
          Prom.Intergenic = mapply(function(prom,intg) if (length(prom$transf) > 0 & length(intg$transf) > 0) {
                               data.frame(Promoter = prom$transf, Intergenic = intg$transf) }, TF$promTF[which(names(TF$promTF)!="Age.neg")], TF$intergenicTF, SIMPLIFY = F),
          Intron.Intergenic = mapply(function(intn,intg) if (length(intn$transf) > 0 & length(intg$transf) > 0) {
                               data.frame(Intron = intn$transf, Intergenic = intg$transf) }, TF$intronTF[which(names(TF$promTF)!="Age.neg")], TF$intergenicTF, SIMPLIFY = F))
df = lapply(df, function(d) lapply(d, function(x) x[is.finite(rowSums(x)),]))
cor = lapply(df, function(d) lapply(d, function(x) cor.test(x[,1], x[,2])))

x = data.frame(Estimate = unlist(lapply(cor, function(x) unlist(lapply(x, function(y) y$estimate)))),
                Stat = unlist(lapply(cor, function(x) unlist(lapply(x, function(y) y$statistic)))),
                pval = unlist(lapply(cor, function(x) unlist(lapply(x, function(y) y$p.value)))))
x[grep("Prom.Intron", rownames(x)),"Comparison"] = "Promoter:Intron"
x[grep("Prom.Intergenic", rownames(x)),"Comparison"] = "Promoter:Intergenic"
x[grep("Intron.Intergenic", rownames(x)),"Comparison"] = "Intron:Intergenic"
x[grep("Age.pos", rownames(x)),"Regions"] = "Age (+Beta)"
x[grep("CellType.pos", rownames(x)),"Regions"] = "Cell Type (+Beta)"
x[grep("Interaction.pos", rownames(x)),"Regions"] = "Interaction (+Beta)"
x[grep("Age.neg", rownames(x)),"Regions"] = "Age (-Beta)"
x[grep("CellType.neg", rownames(x)),"Regions"] = "Cell Type (-Beta)"
x[grep("Interaction.neg", rownames(x)),"Regions"] = "Interaction (-Beta)"
write.csv(x, quote=F, file="/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/PWMEnrich/Intron_vs_Promoter_vs_Intergenic_TF_correlation.csv")

names(df) = c("Promoter:Intron","Promoter:Intergenic","Intron:Intergenic")
pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/PWMEnrich/figures/Intron_vs_Promoter_vs_Intergenic_TF_enrichment.pdf", width = 10, height = 10)
for (i in 1:length(df)) {
  names(df[[i]]) = gsub(".neg", " (-Beta)", names(df[[i]]))
  names(df[[i]]) = gsub(".pos", " (+Beta)", names(df[[i]]))
  for (j in 1:length(df[[i]])) {
    g = ggplot(df[[i]][[j]], aes(x = df[[i]][[j]][,1], y = df[[i]][[j]][,2])) + geom_point() + geom_smooth(method = "lm") + theme_classic() + 
      ylab(colnames(df[[i]][[j]])[2]) + xlab(colnames(df[[i]][[j]])[1]) +
      ggtitle(paste0("-log10(FDR) TF Enrichment\nBetween Features: ",names(df)[i],"; ",names(df[[i]])[j], "\ncor=", 
                     round(cor[[i]][[j]]$estimate,2), ", p=", round(cor[[i]][[j]]$p.value,2))) + 
      theme(title = element_text(size = 20)) +
      theme(text = element_text(size = 20), legend.title=element_blank()) + theme(legend.position="bottom")
    print(g) }
}
dev.off()

# Venn diagram of significantly detected binding motifs

venn = list()
n = c("Age.pos", "CellType.neg", "CellType.pos", "Interaction.neg", "Interaction.pos")
for (i in 1:length(n)) {
  venn[[i]] = list(Promoter = unique(TF$promTF[[n[i]]][which(TF$promTF[[n[i]]][,"padj"]<=0.01),"id"]),
                   Intron = unique(TF$intronTF[[n[i]]][which(TF$intronTF[[n[i]]][,"padj"]<=0.01),"id"]),
                   Intergenic = unique(TF$intergenicTF[[n[i]]][which(TF$intergenicTF[[n[i]]][,"padj"]<=0.01),"target"]))
  venn.diagram(venn[[i]], paste0("/dcl01/lieber/ajaffe/lab/brain-epigenomics/PWMEnrich/figures/",n[i],
                                 "_Intron_Promoter_Intergenic_TFoverlap.jpeg"), 
                                 main=paste0("Intron-Promoter-Intergenic TF Overlap\n",n[i]),
                                 col = "transparent",
                                 fill = c("lightpink2","cornflowerblue", "olivedrab2"),
                                 alpha = 0.50,
                                 fontfamily = "Arial",
                                 fontface = "bold",
                                 cat.col = c("palevioletred4", "darkblue", "olivedrab4"),
                                 cat.fontfamily = "Arial", margin=0.2)
}


## for DMRs by beta direction
df = c(lapply(TF[names(TF)!="intergenicTF"], function(x) list(
                      "Cell Type" = data.frame("+Beta" = x$CellType.pos$transf, "-Beta" = x$CellType.neg$transf),
                      "Age" = data.frame("+Beta" = x$Age.pos$transf, "-Beta" = x$Age.neg$transf),
                      "Interaction" = data.frame("+Beta" = x$Interaction.pos$transf, "-Beta" = x$Interaction.neg$transf))),
       list(intergenicTF = list("Cell Type" = data.frame("+Beta" = TF$intergenicTF$CellType.pos$transf, "-Beta" = TF$intergenicTF$CellType.neg$transf),
                      "Interaction" = data.frame("+Beta" = TF$intergenicTF$Interaction.pos$transf, "-Beta" = TF$intergenicTF$Interaction.neg$transf))))
df = lapply(df, function(d) lapply(d, function(x) x[is.finite(rowSums(x)),]))
cor = lapply(df, function(d) lapply(d, function(x) cor.test(x[,1], x[,2])))

x = data.frame(Estimate = unlist(lapply(cor, function(x) unlist(lapply(x, function(y) y$estimate)))),
               Stat = unlist(lapply(cor, function(x) unlist(lapply(x, function(y) y$statistic)))),
               pval = unlist(lapply(cor, function(x) unlist(lapply(x, function(y) y$p.value)))))
x[grep("allTF", rownames(x)),"Feature"] = "All"
x[grep("promTF", rownames(x)),"Feature"] = "Promoters"
x[grep("intronTF", rownames(x)),"Feature"] = "Introns"
x[grep("intergenicTF", rownames(x)),"Feature"] = "Intergenic"
x[grep("Age", rownames(x)),"Comparison"] = "Age"
x[grep("Cell", rownames(x)),"Comparison"] = "Cell Type"
x[grep("Interaction", rownames(x)),"Comparison"] = "Interaction"
write.csv(x, quote=F, row.names = F, file="/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/PWMEnrich/DMR_byDir_byFeature_TF_correlation.csv")

names(df) = c("All","Promoters","Introns","Intergenic")
pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/PWMEnrich/figures/DMR_byDir_byFeature_TF_enrichment.pdf", width = 10, height = 10)
for (i in 1:length(df)) {
  for (j in 1:length(df[[i]])) {
    g = ggplot(df[[i]][[j]], aes(x = df[[i]][[j]][,1], y = df[[i]][[j]][,2])) + geom_point() + geom_smooth(method = "lm") + theme_classic() + 
      ylab(colnames(df[[i]][[j]])[2]) + xlab(colnames(df[[i]][[j]])[1]) +
      ggtitle(paste0("-log10(FDR) TF Enrichment\nBy Sign: ",names(df[[i]])[j], " (",names(df)[i], ")\ncor=", 
                     round(cor[[i]][[j]]$estimate,2), ", p=", round(cor[[i]][[j]]$p.value,5))) + 
      theme(title = element_text(size = 20)) +
      theme(text = element_text(size = 20), legend.title=element_blank()) + theme(legend.position="bottom")
    print(g) }
}
dev.off()

# Venn diagram of significantly detected binding motifs

venn = list()
comps = list(All_vs_Promoters = c("allTF", "promTF"), All_vs_Introns = c("allTF", "intronTF"), 
             All_vs_Intergenic = c("allTF", "intergenicTF"), Promoters_vs_Introns = c("promTF", "intronTF"), 
             Promoters_vs_Intergenic = c("promTF", "intergenicTF"), Introns_vs_Intergenic = c("intronTF","intergenicTF"))
TFsig = lapply(TF, function(y) lapply(y, function(x) x[which(x$padj<=0.01),]))
for (i in 1:length(comps)) {
  if ("intergenicTF" %in% comps[[i]]) {
  venn[[i]] = list("Cell Type" = list(unique(TFsig[[comps[[i]][1]]][["CellType.pos"]][,"target"]), 
                                      unique(TFsig[[comps[[i]][1]]][["CellType.neg"]][,"target"]),
                                      unique(TFsig[[comps[[i]][2]]][["CellType.pos"]][,"target"]), 
                                      unique(TFsig[[comps[[i]][2]]][["CellType.neg"]][,"target"])),
                   "Interaction" = list(unique(TFsig[[comps[[i]][1]]][["Interaction.pos"]][,"target"]), 
                                        unique(TFsig[[comps[[i]][1]]][["Interaction.neg"]][,"target"]),
                                        unique(TFsig[[comps[[i]][2]]][["Interaction.pos"]][,"target"]), 
                                        unique(TFsig[[comps[[i]][2]]][["Interaction.neg"]][,"target"]))) } 
  else {
    venn[[i]] = list("Cell Type" = list(unique(TFsig[[comps[[i]][1]]][["CellType.pos"]][,"target"]), 
                                        unique(TFsig[[comps[[i]][1]]][["CellType.neg"]][,"target"]),
                                        unique(TFsig[[comps[[i]][2]]][["CellType.pos"]][,"target"]), 
                                        unique(TFsig[[comps[[i]][2]]][["CellType.neg"]][,"target"])),
                     "Age" = list(unique(TFsig[[comps[[i]][1]]][["Age.pos"]][,"target"]), 
                                  unique(TFsig[[comps[[i]][1]]][["Age.neg"]][,"target"]),
                                  unique(TFsig[[comps[[i]][2]]][["Age.pos"]][,"target"]), 
                                  unique(TFsig[[comps[[i]][2]]][["Age.neg"]][,"target"])),
                     "Interaction" = list(unique(TFsig[[comps[[i]][1]]][["Interaction.pos"]][,"target"]), 
                                          unique(TFsig[[comps[[i]][1]]][["Interaction.neg"]][,"target"]),
                                          unique(TFsig[[comps[[i]][2]]][["Interaction.pos"]][,"target"]), 
                                          unique(TFsig[[comps[[i]][2]]][["Interaction.neg"]][,"target"]))) }
}
for (i in 1:length(venn)) {
  for (j in 1:length(venn[[i]])) {
    names(venn[[i]][[j]]) = c(paste0(comps[[i]][1]," (+Beta)"),paste0(comps[[i]][1]," (-Beta)"),
                              paste0(comps[[i]][2]," (+Beta)"),paste0(comps[[i]][2]," (-Beta)"))
    venn.diagram(venn[[i]][[j]], paste0("/dcl01/lieber/ajaffe/lab/brain-epigenomics/PWMEnrich/figures/",
                                        names(comps)[i],"_",names(venn[[i]])[j], "Beta_dir_TFoverlap.jpeg"), 
             main=paste0("TF Overlap: ",names(comps)[i], ": ", names(venn[[i]])[j]),
             col = "transparent",
             fill = c("lightpink2","cornflowerblue", "olivedrab2", "khaki1"),
             cat.col = c("palevioletred4", "darkblue", "olivedrab4", "lightgoldenrod4"),
             alpha = 0.50, fontfamily = "Arial", fontface = "bold", cat.fontfamily = "Arial", margin=0.2) 
    }
}


## For UMRs/LMRs
# Correlate TF enrichment between features

df = list("AllLMR:AllUMR" = data.frame(AllLMR = ulTF$AllLMR$transf, AllUMR = ulTF$AllUMR$transf),
          "NeuronsLMR:NeuronsUMR" = data.frame(NeuronsLMR = ulTF$NeuronsLMR$transf, NeuronsUMR = ulTF$NeuronsUMR$transf),
          "GliaLMR:GliaUMR" = data.frame(GliaLMR = ulTF$GliaLMR$transf, GliaUMR = ulTF$GliaUMR$transf),
          "PostnatalLMR:PostnatalUMR" = data.frame(PostnatalLMR = ulTF$PostnatalLMR$transf, PostnatalUMR = ulTF$PostnatalUMR$transf),
          "PrenatalLMR:PrenatalUMR" = data.frame(PrenatalLMR = ulTF$PrenatalLMR$transf, PrenatalUMR = ulTF$PrenatalUMR$transf),
          "Neurons:Glia" = data.frame(Neurons = ulTF$Neurons$transf, Glia = ulTF$Glia$transf),
          "Neurons:Prenatal" = data.frame(Neurons = ulTF$Neurons$transf, Prenatal = ulTF$Prenatal$transf),
          "Neurons:Postnatal" = data.frame(Neurons = ulTF$Neurons$transf, Postnatal = ulTF$Postnatal$transf),
          "Glia:Postnatal" = data.frame(Glia = ulTF$Glia$transf, Postnatal = ulTF$Postnatal$transf),
          "Glia:Prenatal" = data.frame(Glia = ulTF$Glia$transf, Prenatal = ulTF$Prenatal$transf),
          "Postnatal:Prenatal" = data.frame(Postnatal = ulTF$Postnatal$transf, Prenatal = ulTF$PrenatalUMR$transf),
          "All:Neurons" = data.frame(All = ulTF$All$transf, Neurons = ulTF$Neurons$transf),
          "All:Glia" = data.frame(All = ulTF$All$transf, Glia = ulTF$Glia$transf),
          "All:Postnatal" = data.frame(All = ulTF$All$transf, Postnatal = ulTF$Postnatal$transf),
          "All:Prenatal" = data.frame(All = ulTF$All$transf, Prenatal = ulTF$Prenatal$transf),
          "AllUMR:NeuronsUMR" = data.frame(AllUMR = ulTF$AllUMR$transf, NeuronsUMR = ulTF$NeuronsUMR$transf),
          "AllUMR:GliaUMR" = data.frame(AllUMR = ulTF$AllUMR$transf, GliaUMR = ulTF$GliaUMR$transf),
          "AllUMR:PostnatalUMR" = data.frame(AllUMR = ulTF$AllUMR$transf, PostnatalUMR = ulTF$PostnatalUMR$transf),
          "AllUMR:PrenatalUMR" = data.frame(AllUMR = ulTF$AllUMR$transf, PrenatalUMR = ulTF$PrenatalUMR$transf),
          "AllLMR:NeuronsLMR" = data.frame(AllLMR = ulTF$AllLMR$transf, NeuronsLMR = ulTF$NeuronsLMR$transf),
          "AllLMR:GliaLMR" = data.frame(AllLMR = ulTF$AllLMR$transf, GliaLMR = ulTF$GliaLMR$transf),
          "AllLMR:PostnatalLMR" = data.frame(AllLMR = ulTF$AllLMR$transf, PostnatalLMR = ulTF$PostnatalLMR$transf),
          "AllLMR:PrenatalLMR" = data.frame(AllLMR = ulTF$AllLMR$transf, PrenatalLMR = ulTF$PrenatalLMR$transf))
df = lapply(df, function(x) x[is.finite(rowSums(x)),])
cor = lapply(df, function(x) cor.test(x[,1], x[,2]))

x = data.frame(Estimate = unlist(lapply(cor, function(y) y$estimate)),
               Stat = unlist(lapply(cor, function(y) y$statistic)),
               pval = unlist(lapply(cor, function(y) y$p.value)))
x$Comparison = gsub(".cor", "", row.names(x))
write.csv(x, quote=F, row.names=F,file="/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/PWMEnrich/UMR_vs_LMR_vs_shared_vs_specific_TF_correlation.csv")

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/PWMEnrich/figures/UMR_vs_LMR_vs_shared_vs_specific_TF_enrichment.pdf", width = 10, height = 10)
for (i in 1:length(df)) {
  g = ggplot(df[[i]], aes(x = df[[i]][,1], y = df[[i]][,2])) + geom_point() + geom_smooth(method = "lm") + theme_classic() + 
      ylab(colnames(df[[i]])[2]) + xlab(colnames(df[[i]])[1]) +
      ggtitle(paste0("-log10(FDR) TF Enrichment\nBetween Features: ",names(df)[i],"\ncor=", 
                     round(cor[[i]]$estimate,2), ", p=", round(cor[[i]]$p.value,5))) + 
      theme(title = element_text(size = 20)) +
      theme(text = element_text(size = 20), legend.title=element_blank()) + theme(legend.position="bottom")
    print(g)
}
dev.off()

# Venn diagram of significantly detected binding motifs

venn = list()
comps = list(All_vs_Neurons = c("AllLMR","AllUMR","NeuronsLMR","NeuronsUMR"), 
             All_vs_Glia = c("AllLMR","AllUMR","GliaLMR","GliaUMR"), 
             All_vs_Postnatal = c("AllLMR","AllUMR","PostnatalLMR","PostnatalUMR"), ## missing this (accidentally overwrote)
             All_vs_Prenatal = c("AllLMR","AllUMR","PrenatalLMR","PrenatalUMR"), 
             All_groups = c("All","Neurons","Glia","Postnatal","Prenatal"),
             All_groups_UMR = c("AllUMR","NeuronsUMR","GliaUMR","PostnatalUMR","PrenatalUMR"), 
             All_groups_LMR = c("AllLMR","NeuronsLMR","GliaLMR","PostnatalLMR","PrenatalLMR"))
ulTFsig = lapply(ulTF, function(x)  x[which(x$padj<=0.01),])
for (i in 1:length(comps)) {
  if (length(comps[[i]])==4) {
    venn[[i]] = list(unique(ulTFsig[[comps[[i]][1]]][,"target"]), unique(ulTFsig[[comps[[i]][2]]][,"target"]),
                     unique(ulTFsig[[comps[[i]][3]]][,"target"]), unique(ulTFsig[[comps[[i]][4]]][,"target"])) 
    names(venn[[i]]) = comps[[i]]
    venn.diagram(venn[[i]], paste0("/dcl01/lieber/ajaffe/lab/brain-epigenomics/PWMEnrich/figures/",names(comps)[i],"_UMR_LMR_TFoverlap.jpeg"), 
                 main=paste0("TF Overlap: ",names(comps)[i]),
                 col = "transparent",
                 fill = c("lightpink2","cornflowerblue", "olivedrab2", "khaki1"),
                 cat.col = c("palevioletred4", "darkblue", "olivedrab4", "lightgoldenrod4"),
                 alpha = 0.50, fontfamily = "Arial", fontface = "bold", cat.fontfamily = "Arial", margin=0.2)
  } 
  else {
    venn[[i]] = list(unique(ulTFsig[[comps[[i]][1]]][,"target"]), unique(ulTFsig[[comps[[i]][2]]][,"target"]),
                     unique(ulTFsig[[comps[[i]][3]]][,"target"]), unique(ulTFsig[[comps[[i]][4]]][,"target"]), unique(ulTFsig[[comps[[i]][5]]][,"target"]))
    names(venn[[i]]) = comps[[i]]
    venn.diagram(venn[[i]], paste0("/dcl01/lieber/ajaffe/lab/brain-epigenomics/PWMEnrich/figures/",names(comps)[i],"_UMR_LMR_TFoverlap.jpeg"), 
                 main=paste0("TF Overlap: ",names(comps)[i]),
                 col = "transparent",
                 fill = c("lightpink2","cornflowerblue", "olivedrab2", "khaki1", "plum2"),
                 cat.col = c("palevioletred4", "darkblue", "olivedrab4", "lightgoldenrod4", "orchid4"),
                 alpha = 0.50, fontfamily = "Arial", fontface = "bold", cat.fontfamily = "Arial", margin=0.2)
  }
}


## Are motifs for the same target similarly enriched?

splTF = lapply(TF, function(x) lapply(x, function(y) split(y, y$target)))
splTF = lapply(splTF, function(x) lapply(x, function(y) y[elementNROWS(y)>1]))
lapply(splTF$allTF$Age.pos, function(x) x$padj)
# just by eye, some are the same, some are very different


### Identify TFs that are differentially enriched between UMRs and LMRs

## Explore distribution
ulTFdiff = lapply(ulTFdiff, function(x) x$group.bg[which(names(x$group.bg) %in% ulTF$All$id)])

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/PWMEnrich/figures/ulTFdiff_distribution.pdf")
for (i in 1:length(ulTFdiff)) {
  h = hist(ulTFdiff[[i]], main = paste0("Histogram of ",names(ulTFdiff)[i]), xlab = "Difference Statistic")
  print(h)
}
dev.off()

write.csv(do.call(rbind, Map(cbind, lapply(ulTF, function(x) data.frame(max = max(x$raw.score), min = min(x$raw.score), mean = mean(x$raw.score),
          SD = sd(x$raw.score), median = median(x$raw.score))), Group = as.list(names(ulTF)))), quote=F, row.names = F,
          file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/PWMEnrich/ulTFstats.csv")


## Isolate motifs that are in the top and bottom quartile that are at least enriched to +/-2
ulQ = do.call(rbind, lapply(ulTFdiff, function(x) quantile(x, na.rm = T)))

ulDiffQ25 = mapply(function(ul, q) ul[which(ul<=q[2] & abs(ul)>=2)], ulTFdiff, lapply(ulTFdiff, function(x) quantile(x, na.rm = T)))
ulDiffQ75 = mapply(function(ul, q) ul[which(ul>=q[4] & abs(ul)>=2)], ulTFdiff, lapply(ulTFdiff, function(x) quantile(x, na.rm = T)))


## annotate TFs to their entrez IDS

length(unique(TF$allTF$Age.pos$target)) # 673 unique targets
ulDiffQ25 = lapply(ulDiffQ25, function(x) unique(TF$allTF$Age.pos[which(TF$allTF$Age.pos$id %in% names(x)),"target"]))
ulDiffQ75 = lapply(ulDiffQ75, function(x) unique(TF$allTF$Age.pos[which(TF$allTF$Age.pos$id %in% names(x)),"target"]))
ulDiffQ25 = lapply(ulDiffQ25, function(x) targettogeneID[match(x,names(targettogeneID))])
ulDiffQ75 = lapply(ulDiffQ75, function(x) targettogeneID[match(x,names(targettogeneID))])
ulDiffQ25 = lapply(ulDiffQ25, function(x) geneMap[which(geneMap$gencodeID %in% x),])
ulDiffQ75 = lapply(ulDiffQ75, function(x) geneMap[which(geneMap$gencodeID %in% x),])

entrez = mapply(function(x,y) list(Q25 = unique(na.omit(as.character(x$EntrezID))), Q75 = unique(na.omit(as.character(y$EntrezID)))),
                ulDiffQ25, ulDiffQ75, SIMPLIFY = F) 


## Assess enriched terms limiting the gene universe to the terms associated with the master list of TFs

GeneUniverse = unique(na.omit(as.character(geneMap[which(geneMap$gencodeID %in% targettogeneID), "EntrezID"])))
length(GeneUniverse) # 648

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

# Save

save(keggList, goList_MF, goList_BP, goList_CC, goList_DO, compareBP, compareMF, compareCC, compareKegg,
     file="/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/PWMEnrich/UMR_LMR_GO.objects.PWMEnrich.rda")

ulDiffQ25 = mapply(function(ul, q) ul[which(ul<=q[2] & abs(ul)>=2)], ulTFdiff, lapply(ulTFdiff, function(x) quantile(x, na.rm = T)))
ulDiffQ75 = mapply(function(ul, q) ul[which(ul>=q[4] & abs(ul)>=2)], ulTFdiff, lapply(ulTFdiff, function(x) quantile(x, na.rm = T)))

x = list("All" = cbind(ulTF$AllUMR, ulTF$AllLMR[match(ulTF$AllUMR$id, ulTF$AllLMR$id),]), 
         "Prenatal" = cbind(ulTF$PrenatalUMR, ulTF$PrenatalLMR[match(ulTF$PrenatalUMR$id, ulTF$PrenatalLMR$id),]),        
         "Postnatal" = cbind(ulTF$PostnatalUMR, ulTF$PostnatalLMR[match(ulTF$PostnatalUMR$id, ulTF$PostnatalLMR$id),]),
         "Neurons" = cbind(ulTF$NeuronsUMR, ulTF$NeuronsLMR[match(ulTF$NeuronsUMR$id, ulTF$NeuronsLMR$id),]),        
         "Glia" = cbind(ulTF$GliaUMR, ulTF$GliaLMR[match(ulTF$GliaUMR$id, ulTF$GliaLMR$id),]), 
         "pre.post.UMR" = cbind(ulTF$PrenatalUMR, ulTF$PostnatalUMR[match(ulTF$PrenatalUMR$id, ulTF$PostnatalUMR$id),]),
         "pre.neuron.UMR" = cbind(ulTF$PrenatalUMR, ulTF$NeuronsUMR[match(ulTF$PrenatalUMR$id, ulTF$NeuronsUMR$id),]),
         "pre.glia.UMR" = cbind(ulTF$PrenatalUMR, ulTF$GliaUMR[match(ulTF$PrenatalUMR$id, ulTF$GliaUMR$id),]),
         "neuron.glia.UMR" = cbind(ulTF$NeuronsUMR, ulTF$GliaUMR[match(ulTF$NeuronsUMR$id, ulTF$GliaUMR$id),]),
         "pre.post.LMR" = cbind(ulTF$PrenatalLMR, ulTF$PostnatalLMR[match(ulTF$PrenatalLMR$id, ulTF$PostnatalLMR$id),]),
         "pre.neuron.LMR" = cbind(ulTF$PrenatalLMR, ulTF$NeuronsLMR[match(ulTF$PrenatalLMR$id, ulTF$NeuronsLMR$id),]),
         "pre.glia.LMR" = cbind(ulTF$PrenatalLMR, ulTF$GliaLMR[match(ulTF$PrenatalLMR$id, ulTF$GliaLMR$id),]),
         "neuron.glia.LMR" = cbind(ulTF$NeuronsLMR, ulTF$GliaLMR[match(ulTF$NeuronsLMR$id, ulTF$GliaLMR$id),]),
         "all.pre.UMR" = cbind(ulTF$AllUMR, ulTF$PrenatalUMR[match(ulTF$AllUMR$id, ulTF$PrenatalUMR$id),]),
         "all.post.UMR" = cbind(ulTF$AllUMR, ulTF$PostnatalUMR[match(ulTF$AllUMR$id, ulTF$PostnatalUMR$id),]),
         "all.neuron.UMR" = cbind(ulTF$AllUMR, ulTF$NeuronsUMR[match(ulTF$AllUMR$id, ulTF$NeuronsUMR$id),]),
         "all.glia.UMR" = cbind(ulTF$AllUMR, ulTF$GliaUMR[match(ulTF$AllUMR$id, ulTF$GliaUMR$id),]),
         "all.pre.LMR" = cbind(ulTF$AllLMR, ulTF$PrenatalLMR[match(ulTF$AllLMR$id, ulTF$PrenatalLMR$id),]),
         "all.post.LMR" = cbind(ulTF$AllLMR, ulTF$PostnatalLMR[match(ulTF$AllLMR$id, ulTF$PostnatalLMR$id),]),
         "all.neuron.LMR" = cbind(ulTF$AllLMR, ulTF$NeuronsLMR[match(ulTF$AllLMR$id, ulTF$NeuronsLMR$id),]),
         "all.glia.LMR" = cbind(ulTF$AllLMR, ulTF$GliaLMR[match(ulTF$AllLMR$id, ulTF$GliaLMR$id),]))
for (i in 1:length(x)) { 
  x[[i]] = x[[i]][,-c(grep("transf", colnames(x[[i]])),grep("rank", colnames(x[[i]])))]
  colnames(x[[i]]) = list(c(paste0("UMR.",colnames(x[[i]])[1:6]), paste0("LMR.",colnames(x[[i]])[1:6])),
                          c(paste0("UMR.",colnames(x[[i]])[1:6]), paste0("LMR.",colnames(x[[i]])[1:6])),
                          c(paste0("UMR.",colnames(x[[i]])[1:6]), paste0("LMR.",colnames(x[[i]])[1:6])),
                          c(paste0("UMR.",colnames(x[[i]])[1:6]), paste0("LMR.",colnames(x[[i]])[1:6])),
                          c(paste0("UMR.",colnames(x[[i]])[1:6]), paste0("LMR.",colnames(x[[i]])[1:6])),
                          c(paste0("prenatal.UMR.",colnames(x[[i]])[1:6]), paste0("postnatal.UMR.",colnames(x[[i]])[1:6])),
                          c(paste0("prenatal.UMR",colnames(x[[i]])[1:6]), paste0("neuron.UMR.",colnames(x[[i]])[1:6])),
                          c(paste0("prenatal.UMR",colnames(x[[i]])[1:6]), paste0("glia.UMR.",colnames(x[[i]])[1:6])),
                          c(paste0("neuron.UMR",colnames(x[[i]])[1:6]), paste0("glia.UMR.",colnames(x[[i]])[1:6])),
                          c(paste0("prenatal.LMR.",colnames(x[[i]])[1:6]), paste0("postnatal.LMR.",colnames(x[[i]])[1:6])),
                          c(paste0("prenatal.LMR",colnames(x[[i]])[1:6]), paste0("neuron.LMR.",colnames(x[[i]])[1:6])),
                          c(paste0("prenatal.LMR",colnames(x[[i]])[1:6]), paste0("glia.LMR.",colnames(x[[i]])[1:6])),
                          c(paste0("neuron.LMR",colnames(x[[i]])[1:6]), paste0("glia.LMR.",colnames(x[[i]])[1:6])),
                          c(paste0("all.UMR.",colnames(x[[i]])[1:6]), paste0("prenatal.UMR.",colnames(x[[i]])[1:6])),
                          c(paste0("all.UMR.",colnames(x[[i]])[1:6]), paste0("postnatal.UMR.",colnames(x[[i]])[1:6])),
                          c(paste0("all.UMR.",colnames(x[[i]])[1:6]), paste0("neuron.UMR.",colnames(x[[i]])[1:6])),
                          c(paste0("all.UMR.",colnames(x[[i]])[1:6]), paste0("glia.UMR.",colnames(x[[i]])[1:6])),
                          c(paste0("all.LMR.",colnames(x[[i]])[1:6]), paste0("prenatal.LMR.",colnames(x[[i]])[1:6])),
                          c(paste0("all.LMR.",colnames(x[[i]])[1:6]), paste0("postnatal.LMR.",colnames(x[[i]])[1:6])),
                          c(paste0("all.LMR.",colnames(x[[i]])[1:6]), paste0("neuron.LMR.",colnames(x[[i]])[1:6])),
                          c(paste0("all.LMR.",colnames(x[[i]])[1:6]), paste0("glia.LMR.",colnames(x[[i]])[1:6])))[[i]]
  colnames(x[[i]]) = gsub("padj", "FDR", colnames(x[[i]]))
}
                            
ulDiffQ25 = mapply(function(d,x) data.frame(Quantile = "Q25","Difference Statistic" = d, x[match(names(d), x[,2]),]), ulDiffQ25, x, SIMPLIFY = F) 
ulDiffQ75 = mapply(function(d,x) data.frame(Quantile = "Q75","Difference Statistic" = d, x[match(names(d), x[,2]),]), ulDiffQ75, x, SIMPLIFY = F) 

for (i in 1:length(ulDiffQ25)) { 
  write.csv(rbind(ulDiffQ25[[i]], ulDiffQ75[[i]]), quote = F, row.names = F,
            file = paste0("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/PWMEnrich/differentially_enriched_TFs_", names(ulDiffQ75)[i], ".csv"))  }


## plot Gene ontology enrichment

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/PWMEnrich/figures/KEGG_MF_CC_DO_umr_lmr_TFs.pdf", width=12,height=20)
for (i in 1:length(compareKegg)) {
  print(plot(compareKegg[[i]], colorBy= "p.adjust",  showCategory = 400, title= paste0("KEGG Pathway Enrichment: ", names(compareKegg)[i])))
  print(plot(compareMF[[i]], colorBy= "p.adjust",  showCategory = 400, title= paste0("MF Pathway Enrichment: ", names(compareMF)[i])))
  print(plot(compareCC[[i]], colorBy= "p.adjust",  showCategory = 400, title= paste0("CC Pathway Enrichment: ", names(compareCC)[i])))
}
dev.off()
pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/PWMEnrich/figures/BP_umr_lmr_TFs.pdf", width=12,height=80)
for (i in 1:length(compareBP)) {
  print(plot(compareBP[[i]], colorBy= "p.adjust",  showCategory = 700, title= paste0("Biological Process GO: ", names(compareBP)[i])))
}
dev.off()


### Identify TFs that are differentially enriched between DMRs that fall in various features

## Explore distribution of statistics 

TFdiff = lapply(TFdiff, function(x) x$group.bg[which(names(x$group.bg) %in% TF$allTF$CellType.pos$id)])

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/PWMEnrich/figures/TFdiff_distribution.pdf")
for (i in 1:length(TFdiff)) {
  h = hist(TFdiff[[i]], main = paste0("Histogram of ",names(TFdiff)[i]), xlab = "Difference Statistic")
  print(h)
}
dev.off()

write.csv(do.call(rbind, Map(cbind, lapply(TF, function(y)  do.call(rbind, Map(cbind, lapply(y, function(x) data.frame(max = max(x$raw.score), 
                                           min = min(x$raw.score), mean = mean(x$raw.score), SD = sd(x$raw.score), median = median(x$raw.score))), 
                                           Comparison = as.list(names(y))))), Feature = as.list(gsub("TF","",names(TF))))), quote=F, row.names = F,
          file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/PWMEnrich/TFstats.csv")


## Isolate motifs that are in the top and bottom quartile that are at least enriched to +/-2

tfQ = do.call(rbind, lapply(TFdiff, function(x) quantile(x, na.rm = T)))

tfDiffQ25 = mapply(function(tf, q) tf[which(tf<=q[2] & abs(tf)>=2)], TFdiff, lapply(TFdiff, function(x) quantile(x, na.rm = T)))
tfDiffQ75 = mapply(function(tf, q) tf[which(tf>=q[4] & abs(tf)>=2)], TFdiff, lapply(TFdiff, function(x) quantile(x, na.rm = T)))


## annotate TFs to their entrez IDS

length(unique(TF$allTF$Age.pos$target)) # 673 unique targets
tfDiffQ25 = lapply(tfDiffQ25, function(x) unique(TF$allTF$Age.pos[which(TF$allTF$Age.pos$id %in% names(x)),"target"]))
tfDiffQ75 = lapply(tfDiffQ75, function(x) unique(TF$allTF$Age.pos[which(TF$allTF$Age.pos$id %in% names(x)),"target"]))
tfDiffQ25 = lapply(tfDiffQ25, function(x) targettogeneID[match(x,names(targettogeneID))])
tfDiffQ75 = lapply(tfDiffQ75, function(x) targettogeneID[match(x,names(targettogeneID))])
tfDiffQ25 = lapply(tfDiffQ25, function(x) geneMap[which(geneMap$gencodeID %in% x),])
tfDiffQ75 = lapply(tfDiffQ75, function(x) geneMap[which(geneMap$gencodeID %in% x),])

entrez = mapply(function(x,y) list(Q25 = unique(na.omit(as.character(x$EntrezID))), Q75 = unique(na.omit(as.character(y$EntrezID)))),
                tfDiffQ25, tfDiffQ75, SIMPLIFY = F) 


## Assess enriched terms limiting the gene universe to the terms associated with the master list of TFs

GeneUniverse = unique(na.omit(as.character(geneMap[which(geneMap$gencodeID %in% targettogeneID), "EntrezID"])))
length(GeneUniverse) # 648

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

# Save

save(keggList, goList_MF, goList_BP, goList_CC, goList_DO, compareBP, compareMF, compareCC, compareKegg,
     file="/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/PWMEnrich/DMR_GO.objects.PWMEnrich.rda")

tfDiffQ25 = mapply(function(tf, q) tf[which(tf<=q[2] & abs(tf)>=2)], TFdiff, lapply(TFdiff, function(x) quantile(x, na.rm = T)))
tfDiffQ75 = mapply(function(tf, q) tf[which(tf>=q[4] & abs(tf)>=2)], TFdiff, lapply(TFdiff, function(x) quantile(x, na.rm = T)))

x = list("allAge" = cbind(TF$allTF$Age.pos, TF$allTF$Age.neg[match(TF$allTF$Age.pos$id, TF$allTF$Age.neg$id),]),
         "allInt" = cbind(TF$allTF$Interaction.pos, TF$allTF$Interaction.neg[match(TF$allTF$Interaction.pos$id, TF$allTF$Interaction.neg$id),]),                   
         "allCT" = cbind(TF$allTF$CellType.pos, TF$allTF$CellType.neg[match(TF$allTF$CellType.pos$id, TF$allTF$CellType.neg$id),]),
         "promCT" = cbind(TF$promTF$CellType.pos, TF$promTF$CellType.neg[match(TF$promTF$CellType.pos$id, TF$promTF$CellType.neg$id),]),
         "promAge" = cbind(TF$promTF$Age.pos, TF$promTF$Age.neg[match(TF$promTF$Age.pos$id, TF$promTF$Age.neg$id),]),
         "promInt" = cbind(TF$promTF$Interaction.pos, TF$promTF$Interaction.neg[match(TF$promTF$Interaction.pos$id, TF$promTF$Interaction.neg$id),]),
         "intronsCT" = cbind(TF$intronTF$CellType.pos, TF$intronTF$CellType.neg[match(TF$intronTF$CellType.pos$id, TF$intronTF$CellType.neg$id),]),
         "intronsAge" = cbind(TF$intronTF$Age.pos, TF$intronTF$Age.neg[match(TF$intronTF$Age.pos$id, TF$intronTF$Age.neg$id),]),
         "intronsInt" = cbind(TF$intronTF$Interaction.pos, TF$intronTF$Interaction.neg[match(TF$intronTF$Interaction.pos$id, TF$intronTF$Interaction.neg$id),]),
         "interCT" = cbind(TF$intergenicTF$CellType.pos, TF$intergenicTF$CellType.neg[match(TF$intergenicTF$CellType.pos$id, TF$intergenicTF$CellType.neg$id),]),
         "interInt" = cbind(TF$intergenicTF$Interaction.pos, TF$intergenicTF$Interaction.neg[match(TF$intergenicTF$Interaction.pos$id, TF$intergenicTF$Interaction.neg$id),]),
         "prom.intergenic.pos.CT" = cbind(TF$promTF$CellType.pos, TF$intergenicTF$CellType.pos[match(TF$promTF$CellType.pos$id, TF$intergenicTF$CellType.pos$id),]),
         "prom.intron.pos.CT" = cbind(TF$promTF$CellType.pos, TF$intronTF$CellType.pos[match(TF$promTF$CellType.pos$id, TF$intronTF$CellType.pos$id),]),
         "intron.intergenic.pos.CT" = cbind(TF$intronTF$CellType.pos, TF$intergenicTF$CellType.pos[match(TF$intronTF$CellType.pos$id, TF$intergenicTF$CellType.pos$id),]), 
         "prom.intergenic.neg.CT" = cbind(TF$promTF$CellType.neg, TF$intergenicTF$CellType.neg[match(TF$promTF$CellType.neg$id, TF$intergenicTF$CellType.neg$id),]),    
         "prom.intron.neg.CT" = cbind(TF$promTF$CellType.neg, TF$intronTF$CellType.neg[match(TF$promTF$CellType.neg$id, TF$intronTF$CellType.neg$id),]),       
         "intron.intergenic.neg.CT" = cbind(TF$intronTF$CellType.neg, TF$intergenicTF$CellType.neg[match(TF$intronTF$CellType.neg$id, TF$intergenicTF$CellType.neg$id),]),
         "prom.intergenic.pos.Age" = cbind(TF$promTF$Age.pos, TF$intergenicTF$Age.pos[match(TF$promTF$Age.pos$id, TF$intergenicTF$Age.pos$id),]),  
         "prom.intron.pos.Age" = cbind(TF$promTF$Age.pos, TF$intronTF$Age.pos[match(TF$promTF$Age.pos$id, TF$intronTF$Age.pos$id),]),
         "intron.intergenic.pos.Age" = cbind(TF$intronTF$Age.pos, TF$intergenicTF$Age.pos[match(TF$intronTF$Age.pos$id, TF$intergenicTF$Age.pos$id),]),
         "prom.intron.neg.Age" = cbind(TF$promTF$Age.neg, TF$intronTF$Age.neg[match(TF$promTF$Age.neg$id, TF$intronTF$Age.neg$id),]),
         "prom.intergenic.pos.Int" = cbind(TF$promTF$Interaction.pos, TF$intergenicTF$Interaction.pos[match(TF$promTF$Interaction.pos$id, TF$intergenicTF$Interaction.pos$id),]),  
         "prom.intron.pos.Int" = cbind(TF$promTF$Interaction.pos, TF$intronTF$Interaction.pos[match(TF$promTF$Interaction.pos$id, TF$intronTF$Interaction.pos$id),]),
         "intron.intergenic.pos.Int" = cbind(TF$intronTF$Interaction.pos, TF$intergenicTF$Interaction.pos[match(TF$intronTF$Interaction.pos$id, TF$intergenicTF$Interaction.pos$id),]),
         "prom.intergenic.neg.Int" = cbind(TF$promTF$Interaction.neg, TF$intergenicTF$Interaction.neg[match(TF$promTF$Interaction.neg$id, TF$intergenicTF$Interaction.neg$id),]),    
         "prom.intron.neg.Int" = cbind(TF$promTF$Interaction.neg, TF$intronTF$Interaction.neg[match(TF$promTF$Interaction.neg$id, TF$intronTF$Interaction.neg$id),]),      
         "intron.intergenic.neg.Int" = cbind(TF$intronTF$Interaction.neg, TF$intergenicTF$Interaction.neg[match(TF$intronTF$Interaction.neg$id, TF$intergenicTF$Interaction.neg$id),]))

for (i in 1:length(x)) { 
  x[[i]] = x[[i]][,-c(grep("transf", colnames(x[[i]])),grep("rank", colnames(x[[i]])))]
  colnames(x[[i]]) = list(c(paste0("+Beta",colnames(x[[i]])[1:6]), paste0("-Beta",colnames(x[[i]])[1:6])),
                          c(paste0("+Beta",colnames(x[[i]])[1:6]), paste0("-Beta",colnames(x[[i]])[1:6])),
                          c(paste0("+Beta",colnames(x[[i]])[1:6]), paste0("-Beta",colnames(x[[i]])[1:6])),
                          c(paste0("+Beta",colnames(x[[i]])[1:6]), paste0("-Beta",colnames(x[[i]])[1:6])),
                          c(paste0("+Beta",colnames(x[[i]])[1:6]), paste0("-Beta",colnames(x[[i]])[1:6])),
                          c(paste0("+Beta",colnames(x[[i]])[1:6]), paste0("-Beta",colnames(x[[i]])[1:6])),
                          c(paste0("+Beta",colnames(x[[i]])[1:6]), paste0("-Beta",colnames(x[[i]])[1:6])),
                          c(paste0("+Beta",colnames(x[[i]])[1:6]), paste0("-Beta",colnames(x[[i]])[1:6])),
                          c(paste0("+Beta",colnames(x[[i]])[1:6]), paste0("-Beta",colnames(x[[i]])[1:6])),
                          c(paste0("+Beta",colnames(x[[i]])[1:6]), paste0("-Beta",colnames(x[[i]])[1:6])),
                          c(paste0("+Beta",colnames(x[[i]])[1:6]), paste0("-Beta",colnames(x[[i]])[1:6])),
                          c(paste0("Promoter.CellType.+Beta",colnames(x[[i]])[1:6]), paste0("Intergenic.CellType.+Beta",colnames(x[[i]])[1:6])),
                          c(paste0("Promoter.CellType.+Beta",colnames(x[[i]])[1:6]), paste0("Intron.CellType.+Beta",colnames(x[[i]])[1:6])),
                          c(paste0("Intron.CellType.+Beta",colnames(x[[i]])[1:6]), paste0("Intergenic.CellType.+Beta",colnames(x[[i]])[1:6])),
                          c(paste0("Promoter.CellType.-Beta",colnames(x[[i]])[1:6]), paste0("Intergenic.CellType.-Beta",colnames(x[[i]])[1:6])),
                          c(paste0("Promoter.CellType.-Beta",colnames(x[[i]])[1:6]), paste0("Intron.CellType.-Beta",colnames(x[[i]])[1:6])),
                          c(paste0("Intron.CellType.-Beta",colnames(x[[i]])[1:6]), paste0("Intergenic.CellType.-Beta",colnames(x[[i]])[1:6])),
                          c(paste0("Promoter.Age.+Beta",colnames(x[[i]])[1:6]), paste0("Intergenic.Age.+Beta",colnames(x[[i]])[1:6])),
                          c(paste0("Promoter.Age.+Beta",colnames(x[[i]])[1:6]), paste0("Intron.Age.+Beta",colnames(x[[i]])[1:6])),
                          c(paste0("Intron.Age.+Beta",colnames(x[[i]])[1:6]), paste0("Intergenic.Age.+Beta",colnames(x[[i]])[1:6])),
                          c(paste0("Promoter.Age.-Beta",colnames(x[[i]])[1:6]), paste0("Intron.Age.-Beta",colnames(x[[i]])[1:6])),
                          c(paste0("Promoter.Interaction.+Beta",colnames(x[[i]])[1:6]), paste0("Intergenic.Interaction.+Beta",colnames(x[[i]])[1:6])),
                          c(paste0("Promoter.Interaction.+Beta",colnames(x[[i]])[1:6]), paste0("Intron.Interaction.+Beta",colnames(x[[i]])[1:6])),
                          c(paste0("Intron.Interaction.+Beta",colnames(x[[i]])[1:6]), paste0("Intergenic.Interaction.+Beta",colnames(x[[i]])[1:6])),
                          c(paste0("Promoter.Interaction.-Beta",colnames(x[[i]])[1:6]), paste0("Intergenic.Interaction.-Beta",colnames(x[[i]])[1:6])),
                          c(paste0("Promoter.Interaction.-Beta",colnames(x[[i]])[1:6]), paste0("Intron.Interaction.-Beta",colnames(x[[i]])[1:6])),
                          c(paste0("Intron.Interaction.-Beta",colnames(x[[i]])[1:6]), paste0("Intergenic.Interaction.-Beta",colnames(x[[i]])[1:6])))[[i]]
  colnames(x[[i]]) = gsub("padj", "FDR", colnames(x[[i]]))
}

tfDiffQ25 = mapply(function(d,x) data.frame(Quantile = "Q25","Difference Statistic" = d, x[match(names(d), x[,2]),]), tfDiffQ25, x, SIMPLIFY = F) 
tfDiffQ75 = mapply(function(d,x) data.frame(Quantile = "Q75","Difference Statistic" = d, x[match(names(d), x[,2]),]), tfDiffQ75, x, SIMPLIFY = F) 

for (i in 1:length(tfDiffQ25)) { 
  write.csv(rbind(tfDiffQ25[[i]], tfDiffQ75[[i]]), quote = F, row.names = F,
            file = paste0("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/PWMEnrich/differentially_enriched_TFs_", names(tfDiffQ75)[i], ".csv"))  }


## plot Gene ontology enrichment

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/PWMEnrich/figures/KEGG_MF_CC_DO_DMR_TFs.pdf", width=12,height=20)
for (i in 1:length(compareKegg)) {
  print(plot(compareKegg[[i]], colorBy= "p.adjust",  showCategory = 1500, title= paste0("KEGG Pathway Enrichment: ", names(compareKegg)[i])))
  print(plot(compareMF[[i]], colorBy= "p.adjust",  showCategory = 1500, title= paste0("MF Pathway Enrichment: ", names(compareMF)[i])))
  print(plot(compareCC[[i]], colorBy= "p.adjust",  showCategory = 1500, title= paste0("CC Pathway Enrichment: ", names(compareCC)[i])))
}
dev.off()
pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/PWMEnrich/figures/BP_DMR_TFs.pdf", width=12,height=110)
for (i in 1:length(compareBP)) {
  print(plot(compareBP[[i]], colorBy= "p.adjust",  showCategory = 1500, title= paste0("Biological Process GO:\n", names(compareBP)[i])))
}
dev.off()



# excitatory neurons (Mo et al): Egr, Ap-1, Neurod2, Rfx1/3/5, and TBR1
# PV neurons: Mafb/g, Mef2a/c/d
# VIP neurons: developmental role TFs such as Dlx, Pou, Sox family, Arx and Vax2, which are also enriched in fetal and glial hypo-DMRs