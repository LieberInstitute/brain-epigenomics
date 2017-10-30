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
length(cell[which(cell$fwer<=0.05 & (cell$promoter=="Promoter" | cell$UTR5=="UTR5")), "Coeff.CellTypeNeuron"])
length(cell[which(cell$fwer<=0.05 & (cell$promoter=="Promoter" | cell$UTR5=="UTR5") & cell$value!="NA" & cell$Coeff.CellTypeNeuron!="NA"), "Coeff.CellTypeNeuron"])

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/sorted_nuclear_RNA/figures/LFC_gene_expression_vs_meanBetadiff_byCellType.pdf", width = 9, height = 9)
ggplot(cell[which(cell$fwer<=0.05 & (cell$promoter=="Promoter" | cell$UTR5=="UTR5") & cell$value!="NA" & cell$Coeff.CellTypeNeuron!="NA"),], 
       aes(x = value, y = Coeff.CellTypeNeuron)) + geom_point(alpha=1/10) +
  geom_smooth(method=lm) +
  labs(fill="") +
  ylab("Log2(Fold Change)\nGene Expression") + 
  xlab("Mean Difference in Methylation") +
  ggtitle("Cell Type DMRs (FWER<0.05) in Promoters or 5'UTR\nvs. Cell Type Gene Expression\nr=-0.47") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
ggplot(cell[which(cell$fwer<=0.05 & cell$cds=="CDS" & cell$value!="NA" & cell$Coeff.CellTypeNeuron!="NA"),], 
       aes(x = value, y = Coeff.CellTypeNeuron)) + geom_point(alpha=1/10) +
  geom_smooth(method=lm) +
  labs(fill="") +
  ylab("Log2(Fold Change)\nGene Expression") + 
  xlab("Mean Difference in Methylation") +
  ggtitle("Cell Type DMRs (FWER<0.05) in CDS\nvs. Cell Type Gene Expression\nr=-0.46") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
ggplot(cell[which(cell$fwer<=0.05 & (cell$cds=="CDS" | cell$intron=="Intron") & cell$value!="NA" & cell$Coeff.CellTypeNeuron!="NA"),], 
       aes(x = value, y = Coeff.CellTypeNeuron)) + geom_point(alpha=1/10) +
  geom_smooth(method=lm) +
  labs(fill="") +
  ylab("Log2(Fold Change)\nGene Expression") + 
  xlab("Mean Difference in Methylation") +
  ggtitle("Cell Type DMRs (FWER<0.05) in CDS or Introns\nvs. Cell Type Gene Expression\nr=-0.45") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
ggplot(cell[which(cell$fwer<=0.05 & cell$annotation=="Other" & cell$value!="NA" & cell$Coeff.CellTypeNeuron!="NA"),], 
       aes(x = value, y = Coeff.CellTypeNeuron)) + geom_point(alpha=1/10) +
  geom_smooth(method=lm) +
  labs(fill="") +
  ylab("Log2(Fold Change)\nGene Expression") + 
  xlab("Mean Difference in Methylation") +
  ggtitle("Intergenic Cell Type DMRs (FWER<0.05)\nvs. Nearest Gene Expression\nr=-0.18") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
ggplot(cell[which(cell$fwer<=0.05 & (cell$annotation=="Promoter" | cell$annotation=="UTR5") & cell$value!="NA" & cell$Coeff.CellTypeNeuron!="NA"),], 
       aes(x = value, y = Coeff.CellTypeNeuron)) + geom_point(alpha=1/10) +
  geom_smooth(method=lm) +
  labs(fill="") +
  ylab("Log2(Fold Change)\nGene Expression") + 
  xlab("Mean Difference in Methylation") +
  ggtitle(paste0("Cell Type DMRs (FWER<0.05) in Promoters or 5'UTR\nvs. Cell Type Gene Expression\nr=",round(
    cor(x = cell[which(cell$fwer<=0.05 & (cell$annotation=="Promoter" | cell$annotation=="UTR5")), "value"], 
        y = cell[which(cell$fwer<=0.05 & (cell$annotation=="Promoter" | cell$annotation=="UTR5")), "Coeff.CellTypeNeuron"], use = "pairwise.complete.obs"),3))) +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
ggplot(cell[which(cell$fwer<=0.05 & cell$annotation=="CDS" & cell$value!="NA" & cell$Coeff.CellTypeNeuron!="NA"),], 
       aes(x = value, y = Coeff.CellTypeNeuron)) + geom_point(alpha=1/10) +
  geom_smooth(method=lm) +
  labs(fill="") +
  ylab("Log2(Fold Change)\nGene Expression") + 
  xlab("Mean Difference in Methylation") +
  ggtitle(paste0("Cell Type DMRs (FWER<0.05) in CDS\nvs. Cell Type Gene Expression\nr=",round(
    cor(x = cell[which(cell$fwer<=0.05 & cell$annotation=="CDS"), "value"], 
        y = cell[which(cell$fwer<=0.05 & cell$annotation=="CDS"), "Coeff.CellTypeNeuron"], use = "pairwise.complete.obs"),3))) +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
ggplot(cell[which(cell$fwer<=0.05 & (cell$annotation=="CDS" | cell$annotation=="Intron") & cell$value!="NA" & cell$Coeff.CellTypeNeuron!="NA"),], 
       aes(x = value, y = Coeff.CellTypeNeuron)) + geom_point(alpha=1/10) +
  geom_smooth(method=lm) +
  labs(fill="") +
  ylab("Log2(Fold Change)\nGene Expression") + 
  xlab("Mean Difference in Methylation") +
  ggtitle(paste0("Cell Type DMRs (FWER<0.05) in CDS or Introns\nvs. Cell Type Gene Expression\nr=",round(
    cor(x = cell[which(cell$fwer<=0.05 & (cell$annotation=="CDS" | cell$annotation=="Intron")), "value"], 
        y = cell[which(cell$fwer<=0.05 & (cell$annotation=="CDS" | cell$annotation=="Intron")), "Coeff.CellTypeNeuron"], use = "pairwise.complete.obs"),3))) +
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
#               pval.quad odds.quad conf.int1.quad conf.int2.quad
#Prom.5UTR   2.185770e-149  6.500542       5.585876       7.575867
#CDSonly      1.381958e-97  6.083485       5.079201       7.301256
#CDS.Introns  0.000000e+00  6.403200       5.787354       7.088556
#Intergenic   1.297670e-08  2.026202       1.569693       2.631323


# are genes with DMR more likely to be sigDEG?
tables = list(Prom.5UTR = data.frame(c(length(unique(cell[which(cell$fwer<=0.05 & cell$padj.CellTypeNeuron<=0.05 & (cell$promoter=="Promoter" | cell$UTR5=="UTR5")),"nearestID"])),
                                       length(unique(cell[which(cell$fwer<=0.05 & cell$padj.CellTypeNeuron>0.05 & (cell$promoter=="Promoter" | cell$UTR5=="UTR5")),"nearestID"]))),
                                     c(length(unique(cell[which(cell$fwer>0.05 & cell$padj.CellTypeNeuron<=0.05 & (cell$promoter=="Promoter" | cell$UTR5=="UTR5")),"nearestID"])),
                                       length(unique(cell[which(cell$fwer>0.05 & cell$padj.CellTypeNeuron>0.05 & (cell$promoter=="Promoter" | cell$UTR5=="UTR5")),"nearestID"])))),
              CDSonly = data.frame(c(length(unique(cell[which(cell$fwer<=0.05 & cell$padj.CellTypeNeuron<=0.05 & cell$cds=="CDS"),"nearestID"])),
                                     length(unique(cell[which(cell$fwer<=0.05 & cell$padj.CellTypeNeuron>0.05 & cell$cds=="CDS"),"nearestID"]))),
                                   c(length(unique(cell[which(cell$fwer>0.05 & cell$padj.CellTypeNeuron<=0.05 & cell$cds=="CDS"),"nearestID"])),
                                     length(unique(cell[which(cell$fwer>0.05 & cell$padj.CellTypeNeuron>0.05 & cell$cds=="CDS"),"nearestID"])))),
              CDS.Introns = data.frame(c(length(unique(cell[which(cell$fwer<=0.05 & cell$padj.CellTypeNeuron<=0.05 & (cell$cds=="CDS" | cell$intron=="Intron")),"nearestID"])),
                                         length(unique(cell[which(cell$fwer<=0.05 & cell$padj.CellTypeNeuron>0.05 & (cell$cds=="CDS" | cell$intron=="Intron")),"nearestID"]))),
                                       c(length(unique(cell[which(cell$fwer>0.05 & cell$padj.CellTypeNeuron<=0.05 & (cell$cds=="CDS" | cell$intron=="Intron")),"nearestID"])),
                                         length(unique(cell[which(cell$fwer>0.05 & cell$padj.CellTypeNeuron>0.05 & (cell$cds=="CDS" | cell$intron=="Intron")),"nearestID"])))),
              Intergenic = data.frame(c(length(unique(cell[which(cell$fwer<=0.05 & cell$padj.CellTypeNeuron<=0.05 & cell$annotation=="Other"),"nearestID"])),
                                        length(unique(cell[which(cell$fwer<=0.05 & cell$padj.CellTypeNeuron>0.05 & cell$annotation=="Other"),"nearestID"]))),
                                      c(length(unique(cell[which(cell$fwer>0.05 & cell$padj.CellTypeNeuron<=0.05 & cell$annotation=="Other"),"nearestID"])),
                                        length(unique(cell[which(cell$fwer>0.05 & cell$padj.CellTypeNeuron>0.05 & cell$annotation=="Other"),"nearestID"])))))
fisher = lapply(tables, fisher.test)
p = data.frame(pval.quad = unlist(lapply(fisher, function(x) x$p.value)), 
               odds.quad = unlist(lapply(fisher, function(x) x$estimate)), 
               conf.int1.quad = unlist(lapply(fisher, function(x) x$conf.int[1])),
               conf.int2.quad = unlist(lapply(fisher, function(x) x$conf.int[2])))
#               pval.quad odds.quad conf.int1.quad conf.int2.quad
#Prom.5UTR   7.790004e-78  2.158045       1.990191       2.340004
#CDSonly     9.106769e-31  1.774026       1.607544       1.957919
#CDS.Introns 2.064090e-98  2.032987       1.902478       2.172233
#Intergenic  9.537175e-06  1.327329       1.170376       1.503473

              

### Correlate gene expression and DMR
## by age

nucRNAres$sig.AgeToddler = ifelse(nucRNAres$padj.AgeToddler<=0.05, "sigToddler", "no")
nucRNAres$sig.AgeTeenager = ifelse(nucRNAres$padj.AgeTeenager<=0.05, "sigTeenager", "no")
nucRNAres$LFC.AgeToddler = ifelse(nucRNAres$Coeff.AgeToddler>0, "posToddler", "negToddler")
nucRNAres$LFC.AgeTeenager = ifelse(nucRNAres$Coeff.AgeTeenager>0, "posTeenager", "negTeenager")
nucRNAres$sigAge = paste(nucRNAres$sig.AgeToddler, nucRNAres$sig.AgeTeenager, sep = ":")
nucRNAres$LFCAge = paste(nucRNAres$LFC.AgeToddler, nucRNAres$LFC.AgeTeenager, sep = ":")

x = table(c(nucRNAres$LFCAge=="negToddler:negTeenager" | nucRNAres$LFCAge=="posToddler:posTeenager") == TRUE)
# FALSE  TRUE 
# 11103 35264
round(x/sum(x)*100,3)
# FALSE   TRUE 
# 23.946 76.054 

age = cbind(DMR$Age, nucRNAres[match(DMR$Age$nearestID, nucRNAres$gencodeID),])
dtage = data.table(age)

# correlate all DMRs with expression
cor(x = age[which(age$fwer<=0.05), "value"], 
    y = age[which(age$fwer<=0.05), "Coeff.AgeTeenager"], use = "pairwise.complete.obs") # -0.04910245
# correlate promoter DMRs with expression
cor(x = age[which(age$fwer<=0.05 & age$promoter=="Promoter"), "value"], 
    y = age[which(age$fwer<=0.05 & age$promoter=="Promoter"), "Coeff.AgeTeenager"], use = "pairwise.complete.obs") # -0.06499909
cor(x = age[which(age$fwer<=0.05 & age$annotation=="Promoter"), "value"], 
    y = age[which(age$fwer<=0.05 & age$annotation=="Promoter"), "Coeff.AgeTeenager"], use = "pairwise.complete.obs") # -0.6581951
# correlate promoter or 5'UTR DMRs with expression
cor(x = age[which(age$fwer<=0.05 & (age$promoter=="Promoter" | age$UTR5=="UTR5")), "value"], 
    y = age[which(age$fwer<=0.05 & (age$promoter=="Promoter" | age$UTR5=="UTR5")), "Coeff.AgeTeenager"], use = "pairwise.complete.obs") #  -0.03875779
cor(x = age[which(age$fwer<=0.05 & (age$annotation=="Promoter" | age$annotation=="UTR5")), "value"], 
    y = age[which(age$fwer<=0.05 & (age$annotation=="Promoter" | age$annotation=="UTR5")), "Coeff.AgeTeenager"], use = "pairwise.complete.obs") #  -0.5021648
# correlate gene (distance = 0) DMRs with expression
cor(x = age[which(age$fwer<=0.05 & age$distToGene==0), "value"], 
    y = age[which(age$fwer<=0.05 & age$distToGene==0), "Coeff.AgeTeenager"], use = "pairwise.complete.obs") #  -0.009826338
# correlate intergenic DMRs with expression
cor(x = age[which(age$fwer<=0.05 & age$annotation=="Other"), "value"], 
    y = age[which(age$fwer<=0.05 & age$annotation=="Other"), "Coeff.AgeTeenager"], use = "pairwise.complete.obs") # -0.09765231
# correlate Gene Body (ie, CDS | Introns) with expression
cor(x = age[which(age$fwer<=0.05 & (age$cds=="CDS" | age$intron=="Intron")), "value"], 
    y = age[which(age$fwer<=0.05 & (age$cds=="CDS" | age$intron=="Intron")), "Coeff.AgeTeenager"], use = "pairwise.complete.obs") #  -0.02633665
cor(x = age[which(age$fwer<=0.05 & (age$annotation=="CDS" | age$annotation=="Intron")), "value"], 
    y = age[which(age$fwer<=0.05 & (age$annotation=="CDS" | age$annotation=="Intron")), "Coeff.AgeTeenager"], use = "pairwise.complete.obs") #  -0.0141749
# correlate CDS with expression
cor(x = age[which(age$fwer<=0.05 & age$annotation=="CDS"), "value"], 
    y = age[which(age$fwer<=0.05 & age$annotation=="CDS"), "Coeff.AgeTeenager"], use = "pairwise.complete.obs") #  -0.1903362


length(age[which(age$fwer<=0.05 & (age$promoter=="Promoter" | age$UTR5=="UTR5")), "value"])
length(na.omit(age[which(age$fwer<=0.05 & (age$promoter=="Promoter" | age$UTR5=="UTR5")), "value"]))
length(age[which(age$fwer<=0.05 & (age$promoter=="Promoter" | age$UTR5=="UTR5")), "Coeff.AgeTeenager"])
length(age[which(age$fwer<=0.05 & (age$promoter=="Promoter" | age$UTR5=="UTR5") & age$value!="NA" & age$Coeff.AgeTeenager!="NA"), "Coeff.AgeTeenager"])


pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/sorted_nuclear_RNA/figures/LFC_gene_expression_vs_meanBetadiff_byAge.pdf", width = 9, height = 9)
ggplot(age[which(age$fwer<=0.05 & (age$promoter=="Promoter" | age$UTR5=="UTR5") & age$value!="NA" & age$Coeff.AgeTeenager!="NA"),], 
       aes(x = value, y = Coeff.AgeTeenager)) + geom_point() +
  geom_smooth(method=lm) +
  labs(fill="") +
  ylab("Log2(Fold Change)\nGene Expression") + 
  xlab("Mean Difference in Methylation") +
  ggtitle(paste0("Age DMRs (FWER<0.05) in Promoters or 5'UTR\nvs. Age Gene Expression\nr=",round(
                 cor(x = age[which(age$fwer<=0.05 & (age$promoter=="Promoter" | age$UTR5=="UTR5")), "value"], 
                     y = age[which(age$fwer<=0.05 & (age$promoter=="Promoter" | age$UTR5=="UTR5")), "Coeff.AgeTeenager"], use = "pairwise.complete.obs"),3))) +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
ggplot(age[which(age$fwer<=0.05 & age$cds=="CDS" & age$value!="NA" & age$Coeff.AgeTeenager!="NA"),], 
       aes(x = value, y = Coeff.AgeTeenager)) + geom_point() +
  geom_smooth(method=lm) +
  labs(fill="") +
  ylab("Log2(Fold Change)\nGene Expression") + 
  xlab("Mean Difference in Methylation") +
  ggtitle(paste0("Age DMRs (FWER<0.05) in CDS\nvs. Age Gene Expression\nr=",round(
    cor(x = age[which(age$fwer<=0.05 & age$cds=="CDS"), "value"], 
        y = age[which(age$fwer<=0.05 & age$cds=="CDS"), "Coeff.AgeTeenager"], use = "pairwise.complete.obs"),3))) +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
ggplot(age[which(age$fwer<=0.05 & (age$cds=="CDS" | age$intron=="Intron") & age$value!="NA" & age$Coeff.AgeTeenager!="NA"),], 
       aes(x = value, y = Coeff.AgeTeenager)) + geom_point() +
  geom_smooth(method=lm) +
  labs(fill="") +
  ylab("Log2(Fold Change)\nGene Expression") + 
  xlab("Mean Difference in Methylation") +
  ggtitle(paste0("Age DMRs (FWER<0.05) in CDS or Introns\nvs. Age Gene Expression\nr=",round(
    cor(x = age[which(age$fwer<=0.05 & (age$cds=="CDS" | age$intron=="Intron")), "value"], 
        y = age[which(age$fwer<=0.05 & (age$cds=="CDS" | age$intron=="Intron")), "Coeff.AgeTeenager"], use = "pairwise.complete.obs"),3))) +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
ggplot(age[which(age$fwer<=0.05 & age$annotation=="Other" & age$value!="NA" & age$Coeff.AgeTeenager!="NA"),], 
       aes(x = value, y = Coeff.AgeTeenager)) + geom_point() +
  geom_smooth(method=lm) +
  labs(fill="") +
  ylab("Log2(Fold Change)\nGene Expression") + 
  xlab("Mean Difference in Methylation") +
  ggtitle(paste0("Intergenic Age DMRs (FWER<0.05)\nvs. Nearest Gene Expression\nr=",round(
    cor(x = age[which(age$fwer<=0.05 & age$annotation=="Other"), "value"], 
        y = age[which(age$fwer<=0.05 & age$annotation=="Other"), "Coeff.AgeTeenager"], use = "pairwise.complete.obs"),3))) +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
ggplot(age[which(age$fwer<=0.05 & (age$annotation=="Promoter" | age$annotation=="UTR5") & age$value!="NA" & age$Coeff.AgeTeenager!="NA"),], 
       aes(x = value, y = Coeff.AgeTeenager)) + geom_point() +
  geom_smooth(method=lm) +
  labs(fill="") +
  ylab("Log2(Fold Change)\nGene Expression") + 
  xlab("Mean Difference in Methylation") +
  ggtitle(paste0("Age DMRs (FWER<0.05) in Promoters or 5'UTR\nvs. Age Gene Expression\nr=",round(
    cor(x = age[which(age$fwer<=0.05 & (age$annotation=="Promoter" | age$annotation=="UTR5")), "value"], 
        y = age[which(age$fwer<=0.05 & (age$annotation=="Promoter" | age$annotation=="UTR5")), "Coeff.AgeTeenager"], use = "pairwise.complete.obs"),3))) +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
ggplot(age[which(age$fwer<=0.05 & age$annotation=="CDS" & age$value!="NA" & age$Coeff.AgeTeenager!="NA"),], 
       aes(x = value, y = Coeff.AgeTeenager)) + geom_point() +
  geom_smooth(method=lm) +
  labs(fill="") +
  ylab("Log2(Fold Change)\nGene Expression") + 
  xlab("Mean Difference in Methylation") +
  ggtitle(paste0("Age DMRs (FWER<0.05) in CDS\nvs. Age Gene Expression\nr=",round(
    cor(x = age[which(age$fwer<=0.05 & age$annotation=="CDS"), "value"], 
        y = age[which(age$fwer<=0.05 & age$annotation=="CDS"), "Coeff.AgeTeenager"], use = "pairwise.complete.obs"),3))) +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
ggplot(age[which(age$fwer<=0.05 & (age$annotation=="CDS" | age$annotation=="Intron") & age$value!="NA" & age$Coeff.AgeTeenager!="NA"),], 
       aes(x = value, y = Coeff.AgeTeenager)) + geom_point() +
  geom_smooth(method=lm) +
  labs(fill="") +
  ylab("Log2(Fold Change)\nGene Expression") + 
  xlab("Mean Difference in Methylation") +
  ggtitle(paste0("Age DMRs (FWER<0.05) in CDS or Introns\nvs. Age Gene Expression\nr=",round(
    cor(x = age[which(age$fwer<=0.05 & (age$annotation=="CDS" | age$annotation=="Intron")), "value"], 
        y = age[which(age$fwer<=0.05 & (age$annotation=="CDS" | age$annotation=="Intron")), "Coeff.AgeTeenager"], use = "pairwise.complete.obs"),3))) +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
dev.off()

# make contingency tables of quadrants of the above plots
tables = list(Prom.5UTR = data.frame(c(nrow(age[which(age$fwer<=0.05 & (age$promoter=="Promoter" | age$UTR5=="UTR5") & age$Coeff.AgeTeenager>0 & age$value<0),]),
                                       nrow(age[which(age$fwer<=0.05 & (age$promoter=="Promoter" | age$UTR5=="UTR5") & age$Coeff.AgeTeenager<0 & age$value<0),])),
                                     c(nrow(age[which(age$fwer<=0.05 & (age$promoter=="Promoter" | age$UTR5=="UTR5") & age$Coeff.AgeTeenager>0 & age$value>0),]),
                                       nrow(age[which(age$fwer<=0.05 & (age$promoter=="Promoter" | age$UTR5=="UTR5") & age$Coeff.AgeTeenager<0 & age$value>0),]))),
              CDSonly = data.frame(c(nrow(age[which(age$fwer<=0.05 & age$cds=="CDS" & age$Coeff.AgeTeenager>0 & age$value<0),]),
                                     nrow(age[which(age$fwer<=0.05 & age$cds=="CDS" & age$Coeff.AgeTeenager<0 & age$value<0),])),
                                   c(nrow(age[which(age$fwer<=0.05 & age$cds=="CDS" & age$Coeff.AgeTeenager>0 & age$value>0),]),
                                     nrow(age[which(age$fwer<=0.05 & age$cds=="CDS" & age$Coeff.AgeTeenager<0 & age$value>0),]))),
              CDS.Introns = data.frame(c(nrow(age[which(age$fwer<=0.05 & (age$cds=="CDS" | age$intron=="Intron") & age$Coeff.AgeTeenager>0 & age$value<0),]),
                                         nrow(age[which(age$fwer<=0.05 & (age$cds=="CDS" | age$intron=="Intron") & age$Coeff.AgeTeenager<0 & age$value<0),])),
                                       c(nrow(age[which(age$fwer<=0.05 & (age$cds=="CDS" | age$intron=="Intron") & age$Coeff.AgeTeenager>0 & age$value>0),]),
                                         nrow(age[which(age$fwer<=0.05 & (age$cds=="CDS" | age$intron=="Intron") & age$Coeff.AgeTeenager<0 & age$value>0),]))),
              Intergenic = data.frame(c(nrow(age[which(age$fwer<=0.05 & age$annotation=="Other" & age$Coeff.AgeTeenager>0 & age$value<0),]),
                                        nrow(age[which(age$fwer<=0.05 & age$annotation=="Other" & age$Coeff.AgeTeenager<0 & age$value<0),])),
                                      c(nrow(age[which(age$fwer<=0.05 & age$annotation=="Other" & age$Coeff.AgeTeenager>0 & age$value>0),]),
                                        nrow(age[which(age$fwer<=0.05 & age$annotation=="Other" & age$Coeff.AgeTeenager<0 & age$value>0),]))))

fisher = lapply(tables, fisher.test)
p = data.frame(pval.quad = unlist(lapply(fisher, function(x) x$p.value)), 
               odds.quad = unlist(lapply(fisher, function(x) x$estimate)), 
               conf.int1.quad = unlist(lapply(fisher, function(x) x$conf.int[1])),
               conf.int2.quad = unlist(lapply(fisher, function(x) x$conf.int[2])))
#            pval.quad odds.quad conf.int1.quad conf.int2.quad
#Prom.5UTR   0.2452012 0.5617352      0.1954535       1.573944
#CDSonly     0.1201180 0.4032935      0.1210192       1.272878
#CDS.Introns 0.5499532 0.7361219      0.3112223       1.727350
#Intergenic  1.0000000 0.0000000      0.0000000            Inf


# are genes with DMR more likely to be sigDEG?
tables = list(Prom.5UTR = data.frame(c(length(unique(age[which(age$fwer<=0.05 & age$padj.AgeTeenager<=0.05 & (age$promoter=="Promoter" | age$UTR5=="UTR5")),"nearestID"])),
                                       length(unique(age[which(age$fwer<=0.05 & age$padj.AgeTeenager>0.05 & (age$promoter=="Promoter" | age$UTR5=="UTR5")),"nearestID"]))),
                                     c(length(unique(age[which(age$fwer>0.05 & age$padj.AgeTeenager<=0.05 & (age$promoter=="Promoter" | age$UTR5=="UTR5")),"nearestID"])),
                                       length(unique(age[which(age$fwer>0.05 & age$padj.AgeTeenager>0.05 & (age$promoter=="Promoter" | age$UTR5=="UTR5")),"nearestID"])))),
              CDSonly = data.frame(c(length(unique(age[which(age$fwer<=0.05 & age$padj.AgeTeenager<=0.05 & age$cds=="CDS"),"nearestID"])),
                                     length(unique(age[which(age$fwer<=0.05 & age$padj.AgeTeenager>0.05 & age$cds=="CDS"),"nearestID"]))),
                                   c(length(unique(age[which(age$fwer>0.05 & age$padj.AgeTeenager<=0.05 & age$cds=="CDS"),"nearestID"])),
                                     length(unique(age[which(age$fwer>0.05 & age$padj.AgeTeenager>0.05 & age$cds=="CDS"),"nearestID"])))),
              CDS.Introns = data.frame(c(length(unique(age[which(age$fwer<=0.05 & age$padj.AgeTeenager<=0.05 & (age$cds=="CDS" | age$intron=="Intron")),"nearestID"])),
                                         length(unique(age[which(age$fwer<=0.05 & age$padj.AgeTeenager>0.05 & (age$cds=="CDS" | age$intron=="Intron")),"nearestID"]))),
                                       c(length(unique(age[which(age$fwer>0.05 & age$padj.AgeTeenager<=0.05 & (age$cds=="CDS" | age$intron=="Intron")),"nearestID"])),
                                         length(unique(age[which(age$fwer>0.05 & age$padj.AgeTeenager>0.05 & (age$cds=="CDS" | age$intron=="Intron")),"nearestID"])))),
              Intergenic = data.frame(c(length(unique(age[which(age$fwer<=0.05 & age$padj.AgeTeenager<=0.05 & age$annotation=="Other"),"nearestID"])),
                                        length(unique(age[which(age$fwer<=0.05 & age$padj.AgeTeenager>0.05 & age$annotation=="Other"),"nearestID"]))),
                                      c(length(unique(age[which(age$fwer>0.05 & age$padj.AgeTeenager<=0.05 & age$annotation=="Other"),"nearestID"])),
                                        length(unique(age[which(age$fwer>0.05 & age$padj.AgeTeenager>0.05 & age$annotation=="Other"),"nearestID"])))))
fisher = lapply(tables, fisher.test)
p = data.frame(pval.quad = unlist(lapply(fisher, function(x) x$p.value)), 
               odds.quad = unlist(lapply(fisher, function(x) x$estimate)), 
               conf.int1.quad = unlist(lapply(fisher, function(x) x$conf.int[1])),
               conf.int2.quad = unlist(lapply(fisher, function(x) x$conf.int[2])))
#            pval.quad odds.quad conf.int1.quad conf.int2.quad
#Prom.5UTR   0.3744206  1.344672      0.6534626       2.550465
#CDSonly     0.3745239  1.325453      0.6365305       2.554331
#CDS.Introns 0.1387983  1.500485      0.8535175       2.517915
#Intergenic  0.1049913  3.137352      0.5346431      13.119342










