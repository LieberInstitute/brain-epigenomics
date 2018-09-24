library(ggplot2)

## load results

birnCG = read.csv("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/DMR/Birnbaum_geneSet_enrichment_DMRs.csv")
birnCH = read.csv("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/non-CpG/Birnbaum_geneSet_enrichment_CpHs.csv")

gwasCG = read.csv("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/DMR/GWAS_fisher_DMR_results.csv")
gwasCGgenes = read.csv("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/DMR/GWASgenes_fisher_DMR_results.csv")

gwasCH = read.csv("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/DMR/GWAS_fisher_CpH_results.csv")
gwasCHgenes = read.csv("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/DMR/GWASgenes_fisher_CpH_results.csv")


## put together gene result table

birnCH = data.frame("X" = birnCH$X, "Model" = birnCH$Model, "GeneSet" = birnCH$GeneSet, P.Value = birnCH$P.value, Odds.Ratio = birnCH$Odds.Ratio, FDR = birnCH$FDR)   

birn = rbind(cbind(C.Context = "CpG-derived DMR", birnCG), cbind(C.Context = "CpH", birnCH)) 
birn = data.frame("C.Context" = birn$C.Context, Disease = NA, "Gene Set" = birn$GeneSet, Group = birn$Model, "Odds Ratio" = birn$Odds.Ratio, "P-value" = birn$P.Value)   
gwasGenes = rbind(cbind(C.Context = "CpG-derived DMR", gwasCGgenes), cbind(C.Context = "CpH", gwasCHgenes))
gwasGenes = data.frame("C.Context" = gwasGenes$C.Context, Disease = gwasGenes$GWAS, "Gene Set" = "Genes within genome-wide association study risk loci", 
                       Group = gwasGenes$Model, "Odds Ratio" = gwasGenes$Odds.Ratio, "P-value" = gwasGenes$P.Value)   
gwasGenes[gwasGenes$Disease=="Alz","Disease"] = "Alzheimers Disease"
gwasGenes[gwasGenes$Disease=="Park","Disease"] = "Parkinsons Disease"
gwasGenes$Disease = gsub("T2D", "Type II Diabetes", gwasGenes$Disease)
gwasGenes$Disease = gsub("SCZ", "Schizophrenia", gwasGenes$Disease)
birn$Gene.Set = as.character(birn$Gene.Set)
birn[grep("ASD", birn$"Gene.Set"),"Disease"] = "Autism Spectrum Disorder"
birn[grep("BPAD",birn$Gene.Set),"Disease"] = "Bipolar Affective Disorder"
birn[grep("ID",birn$Gene.Set),"Disease"] = "Intellectual disability"
birn[grep("NDD",birn$Gene.Set),"Disease"] = "Neurodevelopmental disorder"
birn[grep("Neurodegenerative",birn$Gene.Set),"Disease"] = "Neurodegenerative disorder"
birn[grep("SCZ",birn$Gene.Set),"Disease"] = "Schizophrenia"
birn[grep("CNV", birn$"Gene.Set"),"Gene.Set"] = "Copy number variation studies*"
birn[grep("DATABASE", birn$"Gene.Set"),"Gene.Set"] = "Genes from databases*"
birn[grep("Meta-analysis", birn$"Gene.Set"),"Gene.Set"] = "Meta-analysis*"
birn[grep("ID", birn$"Gene.Set"),"Gene.Set"] = "ID studies*"
birn[grep("NDD", birn$"Gene.Set"),"Gene.Set"] = "NDD studies*"
birn[grep("Neurodegenerative", birn$"Gene.Set"),"Gene.Set"] = "Neurodegenerative disease studies*"
birn[grep("SNV", birn$"Gene.Set"),"Gene.Set"] = "Single nucleotide variation studies*"
birn[grep("GWAS", birn$"Gene.Set"),"Gene.Set"] = "Genes within genome-wide association study risk loci"
birn[birn$Disease=="Bipolar Affective Disorder","Gene.Set"] = "Genes within genome-wide association study risk loci*"

genes = rbind(birn, gwasGenes)
genes$FDR = p.adjust(genes$P.value, method = "fdr")

write.csv(genes, quote = F, file="/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/disease_geneSet_enrichment.csv")

genes[genes$FDR<=0.05,]


## put together GWAS loci table

gwas = rbind(cbind(C.Context = "CpG-derived DMR", gwasCG), cbind(C.Context = "CpH", gwasCH)) 
gwas = data.frame("C.Context" = gwas$C.Context, Disease = gwas$GWAS, Group = gwas$Model, "Odds Ratio" = gwas$OR, "P-value" = gwas$pval)   
gwas$FDR = p.adjust(gwas$P.value, method="fdr")
gwas[gwas$Disease=="Alz","Disease"] = "Alzheimers Disease"
gwas[gwas$Disease=="Park","Disease"] = "Parkinsons Disease"
gwas$Disease = gsub("T2D", "Type II Diabetes", gwas$Disease)
gwas$Disease = gsub("SCZ", "Schizophrenia", gwas$Disease)

write.csv(gwas, quote = F, file="/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/disease_GWAS_loci_enrichment.csv")

gwas[gwas$FDR<=0.05,]

df = gwas[grep("Group", gwas$Group),]
df$Group = gsub(" (", "\n(", df$Group, fixed = T)
df$Disease = gsub(" Disease", "\nDisease", df$Disease)
df$Disease = gsub("Type II Diabetes", "Type II\nDiabetes", df$Disease)

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/DMR/figures/disease_GWAS_loci_enrichment.pdf", height = 8.5, width = 7.5)

ggplot(df, aes(Group, Odds.Ratio, fill = Group)) + 
  geom_col() + scale_fill_brewer(8, palette="Dark2") + facet_grid(Disease ~ ., scales = "free_y") +
  theme_classic() + geom_hline(yintercept=1, linetype="dotted") +
  ylab("Odds Ratio") + 
  xlab("") + 
  ggtitle("Enrichment for GWAS Loci") + 
  theme(title = element_text(size = 20)) + theme(text = element_text(size = 20), legend.title=element_blank()) + theme(legend.position="bottom")

dev.off()





