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
# correlate Gene Body (ie, CDS | Introns) with expression
cor(x = cell[which(cell$fwer<=0.05 & (cell$cds=="CDS" | cell$intron=="Intron")), "value"], 
    y = cell[which(cell$fwer<=0.05 & (cell$cds=="CDS" | cell$intron=="Intron")), "Coeff.CellTypeNeuron"], use = "pairwise.complete.obs") #  -0.4494576
cor(x = cell[which(cell$fwer<=0.05 & (cell$annotation=="CDS" | cell$annotation=="Intron")), "value"], 
    y = cell[which(cell$fwer<=0.05 & (cell$annotation=="CDS" | cell$annotation=="Intron")), "Coeff.CellTypeNeuron"], use = "pairwise.complete.obs") #  -0.443522
# correlate CDS with expression
cor(x = cell[which(cell$fwer<=0.05 & cell$annotation=="CDS"), "value"], 
    y = cell[which(cell$fwer<=0.05 & cell$annotation=="CDS"), "Coeff.CellTypeNeuron"], use = "pairwise.complete.obs") #  -0.462388


length(cell[which(cell$fwer<=0.05 & (cell$promoter=="Promoter" | cell$UTR5=="UTR5")), "value"])
length(na.omit(cell[which(cell$fwer<=0.05 & (cell$promoter=="Promoter" | cell$UTR5=="UTR5")), "value"]))
cell[which(cell$fwer<=0.05 & (cell$promoter=="Promoter" | cell$UTR5=="UTR5")), "Coeff.CellTypeNeuron"]
length(cell[which(cell$fwer<=0.05 & (cell$promoter=="Promoter" | cell$UTR5=="UTR5") & cell$value!="NA" & cell$Coeff.CellTypeNeuron!="NA"), "Coeff.CellTypeNeuron"])

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/sorted_nuclear_RNA/figures/LFC_gene_expression_vs_meanBetadiff_byCellType.pdf", width = 9, height = 9)
ggplot(cell[which(cell$fwer<=0.05 & (cell$promoter=="Promoter" | cell$UTR5=="UTR5") & cell$value!="NA" & cell$Coeff.CellTypeNeuron!="NA"),], 
       aes(x = value, y = Coeff.CellTypeNeuron)) + geom_point(alpha=1/10) +
  geom_smooth(method=lm) +
  labs(fill="") +
  ylab("Log2(Fold Change)\nGene Expression") + 
  xlab("Mean Difference in Beta") +
  ggtitle("Cell Type DMRs (FWER<0.05) in Promoters or 5'UTR\nvs. Cell Type Gene Expression\nr=-0.47") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
ggplot(cell[which(cell$fwer<=0.05 & cell$cds=="CDS" & cell$value!="NA" & cell$Coeff.CellTypeNeuron!="NA"),], 
       aes(x = value, y = Coeff.CellTypeNeuron)) + geom_point(alpha=1/10) +
  geom_smooth(method=lm) +
  labs(fill="") +
  ylab("Log2(Fold Change)\nGene Expression") + 
  xlab("Mean Difference in Beta") +
  ggtitle("Cell Type DMRs (FWER<0.05) in CDS\nvs. Cell Type Gene Expression\nr=-0.46") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
ggplot(cell[which(cell$fwer<=0.05 & (cell$cds=="CDS" | cell$intron=="Intron") & cell$value!="NA" & cell$Coeff.CellTypeNeuron!="NA"),], 
       aes(x = value, y = Coeff.CellTypeNeuron)) + geom_point(alpha=1/10) +
  geom_smooth(method=lm) +
  labs(fill="") +
  ylab("Log2(Fold Change)\nGene Expression") + 
  xlab("Mean Difference in Beta") +
  ggtitle("Cell Type DMRs (FWER<0.05) in CDS or Introns\nvs. Cell Type Gene Expression\nr=-0.45") +
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

# make contingency tables of quadrants of the above plots
tables = list(Prom.5UTR = data.frame(c(nrow(cell[which(cell$fwer<=0.05 & (cell$promoter=="Promoter" | cell$UTR5=="UTR5") & cell$Coeff.CellTypeNeuron>0 & cell$value<0),]),
                                       nrow(cell[which(cell$fwer<=0.05 & (cell$promoter=="Promoter" | cell$UTR5=="UTR5") & cell$Coeff.CellTypeNeuron<0 & cell$value<0),])),
                                     c(nrow(cell[which(cell$fwer<=0.05 & (cell$promoter=="Promoter" | cell$UTR5=="UTR5") & cell$Coeff.CellTypeNeuron>0 & cell$value>0),]),
                                       nrow(cell[which(cell$fwer<=0.05 & (cell$promoter=="Promoter" | cell$UTR5=="UTR5") & cell$Coeff.CellTypeNeuron<0 & cell$value>0),]))),
              CDSonly = data.frame(c(nrow(cell[which(cell$fwer<=0.05 & cell$cds=="CDS" & cell$Coeff.CellTypeNeuron>0 & cell$value<0),]),
                                     nrow(cell[which(cell$fwer<=0.05 & cell$cds=="CDS" & cell$Coeff.CellTypeNeuron<0 & cell$value<0),])),
                                   c(nrow(cell[which(cell$fwer<=0.05 & cell$cds=="CDS" & cell$Coeff.CellTypeNeuron>0 & cell$value>0),]),
                                     nrow(cell[which(cell$fwer<=0.05 & cell$cds=="CDS" & cell$Coeff.CellTypeNeuron<0 & cell$value>0),]))),
              CDS.Introns = data.frame(c(nrow(cell[which(cell$fwer<=0.05 & (cell$cds=="CDS" | cell$intron=="Intron") & cell$Coeff.CellTypeNeuron>0 & cell$value<0),]),
                                         nrow(cell[which(cell$fwer<=0.05 & (cell$cds=="CDS" | cell$intron=="Intron") & cell$Coeff.CellTypeNeuron<0 & cell$value<0),])),
                                       c(nrow(cell[which(cell$fwer<=0.05 & (cell$cds=="CDS" | cell$intron=="Intron") & cell$Coeff.CellTypeNeuron>0 & cell$value>0),]),
                                         nrow(cell[which(cell$fwer<=0.05 & (cell$cds=="CDS" | cell$intron=="Intron") & cell$Coeff.CellTypeNeuron<0 & cell$value>0),]))),
              Intergenic = data.frame(c(nrow(cell[which(cell$fwer<=0.05 & cell$annotation=="Other" & cell$Coeff.CellTypeNeuron>0 & cell$value<0),]),
                                        nrow(cell[which(cell$fwer<=0.05 & cell$annotation=="Other" & cell$Coeff.CellTypeNeuron<0 & cell$value<0),])),
                                      c(nrow(cell[which(cell$fwer<=0.05 & cell$annotation=="Other" & cell$Coeff.CellTypeNeuron>0 & cell$value>0),]),
                                        nrow(cell[which(cell$fwer<=0.05 & cell$annotation=="Other" & cell$Coeff.CellTypeNeuron<0 & cell$value>0),]))))

fisher = lapply(tables, fisher.test)
p = data.frame(pval.quad = unlist(lapply(fisher, function(x) x$p.value)), 
               odds.quad = unlist(lapply(fisher, function(x) x$estimate)), 
               conf.int1.quad = unlist(lapply(fisher, function(x) x$conf.int[1])),
               conf.int2.quad = unlist(lapply(fisher, function(x) x$conf.int[2])))
#pval.quad odds.quad conf.int1.quad conf.int2.quad
#Prom.5UTR   2.185770e-149  6.500542       5.585876       7.575867
#CDSonly      1.381958e-97  6.083485       5.079201       7.301256
#CDS.Introns  0.000000e+00  6.403200       5.787354       7.088556
#Intergenic   1.297670e-08  2.026202       1.569693       2.631323


# are genes with DMR more likely to be sigDEG?
tables = list(Prom.5UTR = data.frame(c(nrow(cell[which(cell$fwer<=0.05 & (cell$promoter=="Promoter" | cell$UTR5=="UTR5") & cell$fwer<=0.05 & cell$padj.CellTypeNeuron<=0.05),]),
                                       nrow(cell[which(cell$fwer<=0.05 & (cell$promoter=="Promoter" | cell$UTR5=="UTR5") & cell$fwer>0.05 & cell$padj.CellTypeNeuron<=0.05),])),
                                     c(nrow(cell[which(cell$fwer<=0.05 & (cell$promoter=="Promoter" | cell$UTR5=="UTR5") & cell$fwer<=0.05 & cell$padj.CellTypeNeuron>0.05),]),
                                       nrow(cell[which(cell$fwer<=0.05 & (cell$promoter=="Promoter" | cell$UTR5=="UTR5") & cell$fwer>0.05 & cell$padj.CellTypeNeuron>0.05),]))),
              CDSonly = data.frame(c(nrow(cell[which(cell$fwer<=0.05 & cell$cds=="CDS" & cell$fwer<=0.05 & cell$padj.CellTypeNeuron<=0.05),]),
                                     nrow(cell[which(cell$fwer<=0.05 & cell$cds=="CDS" & cell$fwer>0.05 & cell$padj.CellTypeNeuron<=0.05),])),
                                   c(nrow(cell[which(cell$fwer<=0.05 & cell$cds=="CDS" & cell$fwer<=0.05 & cell$padj.CellTypeNeuron>0.05),]),
                                     nrow(cell[which(cell$fwer<=0.05 & cell$cds=="CDS" & cell$fwer>0.05 & cell$padj.CellTypeNeuron>0.05),]))),
              CDS.Introns = data.frame(c(nrow(cell[which(cell$fwer<=0.05 & (cell$cds=="CDS" | cell$intron=="Intron") & cell$fwer<=0.05 & cell$padj.CellTypeNeuron<=0.05),]),
                                         nrow(cell[which(cell$fwer<=0.05 & (cell$cds=="CDS" | cell$intron=="Intron") & cell$fwer>0.05 & cell$padj.CellTypeNeuron<=0.05),])),
                                       c(nrow(cell[which(cell$fwer<=0.05 & (cell$cds=="CDS" | cell$intron=="Intron") & cell$fwer<=0.05 & cell$padj.CellTypeNeuron>0.05),]),
                                         nrow(cell[which(cell$fwer<=0.05 & (cell$cds=="CDS" | cell$intron=="Intron") & cell$fwer>0.05 & cell$padj.CellTypeNeuron>0.05),]))),
              Intergenic = data.frame(c(nrow(cell[which(cell$fwer<=0.05 & cell$annotation=="Other" & cell$fwer<=0.05 & cell$padj.CellTypeNeuron<=0.05),]),
                                        nrow(cell[which(cell$fwer<=0.05 & cell$annotation=="Other" & cell$fwer>0.05 & cell$padj.CellTypeNeuron<=0.05),])),
                                      c(nrow(cell[which(cell$fwer<=0.05 & cell$annotation=="Other" & cell$fwer<=0.05 & cell$padj.CellTypeNeuron>0.05),]),
                                        nrow(cell[which(cell$fwer<=0.05 & cell$annotation=="Other" & cell$fwer>0.05 & cell$padj.CellTypeNeuron>0.05),]))))









