library(data.table)
library(ggplot2)

load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/sorted_nuclear_RNA/DE_limma_results_objects.rda")
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/bumphunting/bumps_bsseqSmooth_Neuron_cell_250_perm.Rdata")
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/DMR/DMR_objects.rda")
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/misc/rdas/dmr_glia_mean_meth.Rdata")

cell = bumps[[1]]
head(cell)


### Correlate gene expression and DMR in genes split by mean methylation in glial samples

## by cell type
cell = cbind(DMR$CellType, nucRNAres[match(DMR$CellType$nearestID, nucRNAres$gencodeID),], 
             mean_glial_meth = dmr_glia_mean_meth$cell$mean_meth, mean_glial_meth_bin = dmr_glia_mean_meth$cell$mean_meth_bin)
dtcell = data.table(cell)
split = split(cell, cell$mean_glial_meth_bin)

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/sorted_nuclear_RNA/figures/LFC_gene_expression_vs_meanBetadiff_byCellType_glialMethBins.pdf", width = 10, height = 10)
for (i in 1:length(split)){
  tmp = split[[i]]
g = ggplot(tmp[which(tmp$fwer<=0.05 & (tmp$promoter=="Promoter" | tmp$UTR5=="UTR5") & tmp$value!="NA" & tmp$Coeff.CellTypeNeuron!="NA"),], 
       aes(x = value, y = Coeff.CellTypeNeuron)) + geom_point(alpha=1/3) +
  geom_smooth(method=lm) +
  labs(fill="") +
  ylab("Log2(Fold Change)\nGene Expression") + 
  xlab("Mean Difference in Beta") +
  ggtitle(paste0("Cell Type DMRs (FWER<0.05) in Promoters or 5'UTR\nvs. Cell Type Gene Expression: ", names(split)[i], "\nr=", 
      round(cor(x = tmp[which(tmp$fwer<=0.05 & (tmp$promoter=="Promoter" | tmp$UTR5=="UTR5") & tmp$value!="NA" & tmp$Coeff.CellTypeNeuron!="NA"),"value"],
                y = tmp[which(tmp$fwer<=0.05 & (tmp$promoter=="Promoter" | tmp$UTR5=="UTR5") & tmp$value!="NA" & tmp$Coeff.CellTypeNeuron!="NA"),"Coeff.CellTypeNeuron"]),2))) +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
print(g)
h = ggplot(tmp[which(tmp$fwer<=0.05 & tmp$cds=="CDS" & tmp$value!="NA" & tmp$Coeff.CellTypeNeuron!="NA"),], 
           aes(x = value, y = Coeff.CellTypeNeuron)) + geom_point(alpha=1/3) +
  geom_smooth(method=lm) +
  labs(fill="") +
  ylab("Log2(Fold Change)\nGene Expression") + 
  xlab("Mean Difference in Beta") +
  ggtitle(paste0("Cell Type DMRs (FWER<0.05) in CDS\nvs. Cell Type Gene Expression: ", names(split)[i], "\nr=", 
          round(cor(x = tmp[which(tmp$fwer<=0.05 & tmp$cds=="CDS" & tmp$value!="NA" & tmp$Coeff.CellTypeNeuron!="NA"),"value"],
                    y = tmp[which(tmp$fwer<=0.05 & tmp$cds=="CDS" & tmp$value!="NA" & tmp$Coeff.CellTypeNeuron!="NA"),"Coeff.CellTypeNeuron"]),2))) +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
print(h)
a = ggplot(tmp[which(tmp$fwer<=0.05 & (tmp$cds=="CDS" | tmp$intron=="Intron") & tmp$value!="NA" & tmp$Coeff.CellTypeNeuron!="NA"),], 
           aes(x = value, y = Coeff.CellTypeNeuron)) + geom_point(alpha=1/3) +
  geom_smooth(method=lm) +
  labs(fill="") +
  ylab("Log2(Fold Change)\nGene Expression") + 
  xlab("Mean Difference in Beta") +
  ggtitle(paste0("Cell Type DMRs (FWER<0.05) in CDS or Introns\nvs. Cell Type Gene Expression: ", names(split)[i], "\nr=", 
          round(cor(x = tmp[which(tmp$fwer<=0.05 & (tmp$cds=="CDS" | tmp$intron=="Intron") & tmp$value!="NA" & tmp$Coeff.CellTypeNeuron!="NA"),"value"],
                    y = tmp[which(tmp$fwer<=0.05 & (tmp$cds=="CDS" | tmp$intron=="Intron") & tmp$value!="NA" & tmp$Coeff.CellTypeNeuron!="NA"),"Coeff.CellTypeNeuron"]),2))) +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
print(a)
b = ggplot(tmp[which(tmp$fwer<=0.05 & tmp$annotation=="Other" & tmp$value!="NA" & tmp$Coeff.CellTypeNeuron!="NA"),], 
           aes(x = value, y = Coeff.CellTypeNeuron)) + geom_point(alpha=1/3) +
  geom_smooth(method=lm) +
  labs(fill="") +
  ylab("Log2(Fold Change)\nGene Expression") + 
  xlab("Mean Difference in Beta") +
  ggtitle(paste0("Intergenic Cell Type DMRs (FWER<0.05)\nvs. Nearest Gene Expression: ", names(split)[i], "\nr=", 
          round(cor(x = tmp[which(tmp$fwer<=0.05 & tmp$annotation=="Other" & tmp$value!="NA" & tmp$Coeff.CellTypeNeuron!="NA"),"value"],
                    y = tmp[which(tmp$fwer<=0.05 & tmp$annotation=="Other" & tmp$value!="NA" & tmp$Coeff.CellTypeNeuron!="NA"),"Coeff.CellTypeNeuron"]),2))) +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
print(b)
}
dev.off()


# make contingency tables of quadrants of the above plots
tables = list()
for (i in 1:length(split)){
  tmp = split[[i]]
  tables[[i]] = list(Prom.5UTR = data.frame(c(nrow(tmp[which(tmp$fwer<=0.05 & (tmp$promoter=="Promoter" | tmp$UTR5=="UTR5") & tmp$Coeff.CellTypeNeuron>0 & tmp$value<0),]),
                                       nrow(tmp[which(tmp$fwer<=0.05 & (tmp$promoter=="Promoter" | tmp$UTR5=="UTR5") & tmp$Coeff.CellTypeNeuron<0 & tmp$value<0),])),
                                     c(nrow(tmp[which(tmp$fwer<=0.05 & (tmp$promoter=="Promoter" | tmp$UTR5=="UTR5") & tmp$Coeff.CellTypeNeuron>0 & tmp$value>0),]),
                                       nrow(tmp[which(tmp$fwer<=0.05 & (tmp$promoter=="Promoter" | tmp$UTR5=="UTR5") & tmp$Coeff.CellTypeNeuron<0 & tmp$value>0),]))),
              CDSonly = data.frame(c(nrow(tmp[which(tmp$fwer<=0.05 & tmp$cds=="CDS" & tmp$Coeff.CellTypeNeuron>0 & tmp$value<0),]),
                                     nrow(tmp[which(tmp$fwer<=0.05 & tmp$cds=="CDS" & tmp$Coeff.CellTypeNeuron<0 & tmp$value<0),])),
                                   c(nrow(tmp[which(tmp$fwer<=0.05 & tmp$cds=="CDS" & tmp$Coeff.CellTypeNeuron>0 & tmp$value>0),]),
                                     nrow(tmp[which(tmp$fwer<=0.05 & tmp$cds=="CDS" & tmp$Coeff.CellTypeNeuron<0 & tmp$value>0),]))),
              CDS.Introns = data.frame(c(nrow(tmp[which(tmp$fwer<=0.05 & (tmp$cds=="CDS" | tmp$intron=="Intron") & tmp$Coeff.CellTypeNeuron>0 & tmp$value<0),]),
                                         nrow(tmp[which(tmp$fwer<=0.05 & (tmp$cds=="CDS" | tmp$intron=="Intron") & tmp$Coeff.CellTypeNeuron<0 & tmp$value<0),])),
                                       c(nrow(tmp[which(tmp$fwer<=0.05 & (tmp$cds=="CDS" | tmp$intron=="Intron") & tmp$Coeff.CellTypeNeuron>0 & tmp$value>0),]),
                                         nrow(tmp[which(tmp$fwer<=0.05 & (tmp$cds=="CDS" | tmp$intron=="Intron") & tmp$Coeff.CellTypeNeuron<0 & tmp$value>0),]))),
              Intergenic = data.frame(c(nrow(tmp[which(tmp$fwer<=0.05 & tmp$annotation=="Other" & tmp$Coeff.CellTypeNeuron>0 & tmp$value<0),]),
                                        nrow(tmp[which(tmp$fwer<=0.05 & tmp$annotation=="Other" & tmp$Coeff.CellTypeNeuron<0 & tmp$value<0),])),
                                      c(nrow(tmp[which(tmp$fwer<=0.05 & tmp$annotation=="Other" & tmp$Coeff.CellTypeNeuron>0 & tmp$value>0),]),
                                        nrow(tmp[which(tmp$fwer<=0.05 & tmp$annotation=="Other" & tmp$Coeff.CellTypeNeuron<0 & tmp$value>0),]))))
}
