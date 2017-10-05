library(data.table)
library(ggplot2)


load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/sorted_nuclear_RNA/DE_limma_results_objects.rda")
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/DMR/DMR_objects.rda")

# How many genes are DE by cell type?
dim(nucRNAres[which(nucRNAres$padj.CellTypeNeuron<=0.05),]) # 9994

# How many genes are DE by cell type and greater in neurons?
dim(nucRNAres[which(nucRNAres$padj.CellTypeNeuron<=0.05 & nucRNAres$Coeff.CellTypeNeuron > 0),]) # 6521


### Correlate gene expression and DMR

## by cell type
cell = cbind(DMR$CellType, nucRNAres[match(DMR$CellType$nearestID, nucRNAres$gencodeID),])
dtcell = data.table(cell)

# correlate all DMRs with expression
cor(x = cell[which(cell$fwer<=0.05), "value"], 
    y = cell[which(cell$fwer<=0.05), "Coeff.CellTypeNeuron"], use = "pairwise.complete.obs") # -0.3914629
# correlate promoter DMRs with expression
cor(x = cell[which(cell$fwer<=0.05 & cell$promoter=="Promoter"), "value"], 
    y = cell[which(cell$fwer<=0.05 & cell$promoter=="Promoter"), "Coeff.CellTypeNeuron"], use = "pairwise.complete.obs") # -0.4641855
cor(x = cell[which(cell$fwer<=0.05 & cell$annotation=="Promoter"), "value"], 
    y = cell[which(cell$fwer<=0.05 & cell$annotation=="Promoter"), "Coeff.CellTypeNeuron"], use = "pairwise.complete.obs") # -0.3319433
# correlate promoter or 5'UTR DMRs with expression
cor(x = cell[which(cell$fwer<=0.05 & (cell$promoter=="Promoter" | cell$UTR5=="UTR5")), "value"], 
    y = cell[which(cell$fwer<=0.05 & (cell$promoter=="Promoter" | cell$UTR5=="UTR5")), "Coeff.CellTypeNeuron"], use = "pairwise.complete.obs") #  -0.4665043
cor(x = cell[which(cell$fwer<=0.05 & (cell$annotation=="Promoter" | cell$annotation=="UTR5")), "value"], 
    y = cell[which(cell$fwer<=0.05 & (cell$annotation=="Promoter" | cell$annotation=="UTR5")), "Coeff.CellTypeNeuron"], use = "pairwise.complete.obs") #  -0.4363155
# correlate gene (distance = 0) DMRs with expression
cor(x = cell[which(cell$fwer<=0.05 & cell$distToGene==0), "value"], 
    y = cell[which(cell$fwer<=0.05 & cell$distToGene==0), "Coeff.CellTypeNeuron"], use = "pairwise.complete.obs") #  -0.4449334
# correlate intergenic DMRs with expression
cor(x = cell[which(cell$fwer<=0.05 & cell$annotation=="Other"), "value"], 
    y = cell[which(cell$fwer<=0.05 & cell$annotation=="Other"), "Coeff.CellTypeNeuron"], use = "pairwise.complete.obs") # -0.1755019

length(cell[which(cell$fwer<=0.05 & (cell$promoter=="Promoter" | cell$UTR5=="UTR5")), "value"])
length(na.omit(cell[which(cell$fwer<=0.05 & (cell$promoter=="Promoter" | cell$UTR5=="UTR5")), "value"]))
cell[which(cell$fwer<=0.05 & (cell$promoter=="Promoter" | cell$UTR5=="UTR5")), "Coeff.CellTypeNeuron"]
length(cell[which(cell$fwer<=0.05 & (cell$promoter=="Promoter" | cell$UTR5=="UTR5") & cell$value!="NA" & cell$Coeff.CellTypeNeuron!="NA"), "Coeff.CellTypeNeuron"])

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/sorted_nuclear_RNA/figures/LFC_gene_expression_vs_meanBetadiff_byCellType.pdf", width = 10)
ggplot(cell[which(cell$fwer<=0.05 & (cell$promoter=="Promoter" | cell$UTR5=="UTR5") & cell$value!="NA" & cell$Coeff.CellTypeNeuron!="NA"),], 
       aes(x = value, y = Coeff.CellTypeNeuron)) + geom_point(alpha=1/10) +
  geom_smooth(method=lm) +
  labs(fill="") +
  ylab("Log2(Fold Change)\nGene Expression") + 
  xlab("Mean Difference in Beta") +
  ggtitle("Cell Type DMRs (FWER<0.05) in Promoters or 5'UTR\nvs. Cell Type Gene Expression\nr=-0.47") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
ggplot(cell[which(cell$fwer<=0.05 & cell$annotation=="Other" & cell$value!="NA" & cell$Coeff.CellTypeNeuron!="NA"),], 
       aes(x = value, y = Coeff.CellTypeNeuron)) + geom_point(alpha=1/10) +
  geom_smooth(method=lm) +
  labs(fill="") +
  ylab("Log2(Fold Change)\nGene Expression") + 
  xlab("Mean Difference in Beta") +
  ggtitle("Intergenic Cell Type DMRs (FWER<0.05)\nvs. Nearest Gene Expression\nr=-0.18") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
dev.off()









