######
library('bsseq')
library('bumphunter')
library('doParallel')
library(limma)

## load BSobj
load('/dcl01/lieber/ajaffe/lab/brain-epigenomics/bumphunting/BSobj_bsseqSmooth_Neuron_minCov_3.Rdata')
grFilter = granges(BSobj)

## full
load("/dcl01/lieber/WGBS/LIBD_Data/bsseqObj/bsseqObj_postNatal_cleaned_CpGonly.rda")
mm = match(grFilter, granges(BSobj))
BSobj = BSobj[mm,]

## extract pheno
pd <- pData(BSobj)

## Get chr coordinates and methylation values
gr <- granges(BSobj)
meth <- getMeth(BSobj, type = 'raw')

#########################
## homogenate age #######
#########################

homIndex = which(pd$Cell.Type == "Homogenate")
fitHom = lmFit(meth[,homIndex], model.matrix(~Age, data=pd[homIndex,]))
ebHom = ebayes(fitHom)

#########################
## cell type specific ###
#########################

cellIndex = which(pd$Cell.Type != "Homogenate")
fitCell = lmFit(meth[,cellIndex], model.matrix(~Age*Cell.Type, data=pd[cellIndex,]))
ebCell = ebayes(fitCell)

stats = data.frame(dm_age_hom = fitHom$coef[,2],
	dm_age_cell = fitCell$coef[,2],
	dm_int_cell = fitCell$coef[,4],
	dm_type_cell = fitCell$coef[,3],
	t_age_hom = ebHom$t[,2],
	t_age_cell = ebCell$t[,2],
	t_int_cell = ebCell$t[,4],
	t_type_cell = ebCell$t[,3],
	pval_age_hom = ebHom$p[,2],
	pval_age_cell = ebCell$p[,2],
	pval_int_cell = ebCell$p[,4],
	pval_type_cell = ebCell$p[,3])
stats$qval_age_hom = p.adjust(stats$pval_age_hom, "fdr")
stats$qval_age_cell = p.adjust(stats$pval_age_cell, "fdr")
stats$qval_int_cell = p.adjust(stats$pval_int_cell, "fdr")
stats$qval_type_cell = p.adjust(stats$pval_type_cell, "fdr")

colSums(stats[,grep("qval", colnames(stats))] < 0.05)	

#########################
## filter ###############
#########################
sigIndex = which(rowSums(stats[,grep("pval", colnames(stats))] < 1e-4) > 0)
methSig = meth[sigIndex,]
grSig = gr[sigIndex]
statsSig = stats[sigIndex,]
rownames(statsSig) = paste0(seqnames(grSig), ":", start(grSig))
save(grSig, pd, methSig, statsSig, file="rdas/signif_stats_1e4_highCov.rda")

## metrics
colSums(statsSig[,grep("pval", colnames(statsSig))] < 1e-4)
colSums(statsSig[,grep("pval", colnames(statsSig))] < 1e-6)
colSums(statsSig[,grep("qval", colnames(statsSig))] < 0.05)

############
## first plot signif in homogenate
load("rdas/signif_stats_1e4_highCov.rda")

statsSig_homSig = statsSig[statsSig$qval_age_hom < 0.05,]

# cell type?
pdf("plots/dm_cellType_vs_homAge_highCov.pdf") # pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/single_CpGs/figures/dm_cellType_vs_homAge_highCov.pdf")
par(mar=c(5,6,2,2),cex.axis=2,cex.lab=2)
plot(dm_type_cell ~ dm_age_hom,data=statsSig,
	subset = statsSig$qval_age_hom < 0.05,pch=21,bg="grey",
	xlim=c(-0.05,0.05),ylim=c(-0.8,0.8),cex=1.4,
	xlab = "Age (Homogenate)",ylab="Cell Type")
abline(h=0,v=0,lty=3)
dev.off()


table(statsSig_homSig$dm_type_cell>0 ,  statsSig_homSig$dm_age_hom<0, dnn = c("Neuron>Glia", "AgeDecr"))
fisher.test(table(statsSig_homSig$dm_type_cell>0 ,  statsSig_homSig$dm_age_hom<0))
fisher.test(table(statsSig_homSig$dm_type_cell>0 ,  statsSig_homSig$dm_age_hom<0))$p.value
# p-value < 2.2e-16
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
 # 6.846840 9.054318
# sample estimates:
# odds ratio
  # 7.865103



# combined
pdf("plots/dm_totalAge_vs_homAge_highCov.pdf") # pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/single_CpGs/figures/dm_totalAge_vs_homAge_highCov.pdf")
par(mar=c(5,6,2,2),cex.axis=2,cex.lab=2)
statsSig$totalAge =(statsSig$dm_age_cell + statsSig$dm_int_cell)
statsSig_homSig$totalAge =(statsSig_homSig$dm_age_cell + statsSig_homSig$dm_int_cell)
plot(totalAge ~ dm_age_hom,data=statsSig,
	subset = statsSig$qval_age_hom < 0.05,
	pch=21,bg="grey",
	xlim=c(-0.05,0.05),ylim=c(-0.05,0.05),
	xlab = "Age (Homogenate)",ylab="Total Cell Type-Specific")
abline(h=0,v=0,lty=3)
dev.off()


table(statsSig_homSig$totalAge<0 ,  statsSig_homSig$dm_age_hom<0, dnn = c("TotalAgeDecr", "HomAgeDecr"))
fisher.test(table(statsSig_homSig$totalAge < 0 ,  statsSig_homSig$dm_age_hom<0))
fisher.test(table(statsSig_homSig$totalAge < 0 ,  statsSig_homSig$dm_age_hom<0))$p.value

# p-value < 2.2e-16
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
 # 30.16703 41.97367
# sample estimates:
# odds ratio
  # 35.50007


cor.test(statsSig_homSig$totalAge,  statsSig_homSig$dm_age_hom)
cor.test(statsSig_homSig$totalAge,  statsSig_homSig$dm_age_hom)$p.value
# t = 82.68, df = 5776, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
 # 0.7241835 0.7478104
# sample estimates:
      # cor
# 0.7362212


## separate
pdf("plots/dm_separateAge_vs_homAge_highCov.pdf")
par(mar=c(5,6,2,2),cex.axis=2,cex.lab=2)
plot(dm_age_cell ~ dm_age_hom,data=statsSig,
	subset = statsSig$qval_age_hom < 0.05,
	pch=21,bg="grey",
	xlim=c(-0.05,0.05),ylim=c(-0.05,0.05),
	xlab = "Age (Homogenate)",ylab="Overall Age (Cell)")
abline(h=0,v=0,lty=3)


table(statsSig_homSig$dm_age_cell<0 ,  statsSig_homSig$dm_age_hom<0, dnn = c("CellAgeDecr", "HomAgeDecr"))
fisher.test(table(statsSig_homSig$dm_age_cell < 0 ,  statsSig_homSig$dm_age_hom<0))
fisher.test(table(statsSig_homSig$dm_age_cell < 0 ,  statsSig_homSig$dm_age_hom<0))$p.value

# p-value < 2.2e-16
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
  # 71.93494 111.73643
# sample estimates:
# odds ratio
   # 89.5082


cor.test(statsSig_homSig$dm_age_cell, statsSig_homSig$dm_age_hom)
# t = 109.88, df = 5776, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
 # 0.8139248 0.8306196
# sample estimates:
      # cor
# 0.8224492


# interaction
plot(dm_int_cell ~ dm_age_hom,data=statsSig,
	subset = statsSig$qval_age_hom < 0.05,pch=21,bg="grey",
	xlim=c(-0.05,0.05),ylim=c(-0.05,0.05),
	xlab = "Age (Homogenate)",ylab="Specific Age (Cell)")
abline(h=0,v=0,lty=3)
dev.off()

fisher.test(data.frame(c(nrow(statsSig[which(statsSig$dm_int_cell>0 & statsSig$dm_age_hom<0 & statsSig$pval_age_hom < 1e-4),]),
                         nrow(statsSig[which(statsSig$dm_int_cell<0 & statsSig$dm_age_hom<0 & statsSig$pval_age_hom < 1e-4),])),
                       c(nrow(statsSig[which(statsSig$dm_int_cell>0 & statsSig$dm_age_hom>0 & statsSig$pval_age_hom < 1e-4),]),
                         nrow(statsSig[which(statsSig$dm_int_cell<0 & statsSig$dm_age_hom>0 & statsSig$pval_age_hom < 1e-4),]))))
#p-value < 2.2e-16
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  2.215964 2.517668
#sample estimates:
#  odds ratio 
#2.361638 

cor.test(x = statsSig[which(statsSig$pval_age_hom < 1e-4),"dm_int_cell"],
         y = statsSig[which(statsSig$pval_age_hom < 1e-4),"dm_age_hom"])
#t = -41.302, df = 22592, p-value < 2.2e-16
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
#  -0.2770483 -0.2528002
#sample estimates:
#  cor 
#-0.2649661 


### just neuronal
neunIndex = which(pd$Cell.Type == "Neuron")
fitNeuron = lmFit(methSig[,neunIndex], 
	model.matrix(~Age, data=pd[neunIndex,]))
ebNeuron = ebayes(fitNeuron)
plot(fitNeuron$coef[,2] ~ dm_age_hom,data=statsSig,
	subset = pval_age_hom < 1e-4,pch=21,bg="grey",
	xlim=c(-0.05,0.05),ylim=c(-0.05,0.05),
	xlab = "Age (Homogenate)",ylab="Age (Neuron)")
abline(h=0,v=0,lty=3)
dev.off()


###
plot(t_age_cell ~ t_age_hom,data=statsSig,
	subset = pval_age_hom < 1e-4,pch=21,bg="grey",
	xlim=c(-10,10),ylim=c(-10,10),
	xlab = "Homogenate",ylab="Cell Type-Specific (Main)")
abline(h=0,v=0,lty=3)
abline(0,1,col="red")
plot(t_int_cell ~ t_age_hom,data=statsSig,
	subset = pval_age_hom < 1e-4,pch=21,bg="grey",
	xlim=c(-10,10),ylim=c(-10,10),
	xlab = "Homogenate",ylab="Cell Type-Specific (Intxn)")
abline(h=0,v=0,lty=3)
abline(0,1,col="red")

#################
### non CpG #####
#################
xx = load("bsseq/bsobj_by_chr/limma_exploration_nonCG_highCov.Rdata")
sigIndex2 = which(rowSums(ebList$interaction$p[,-1] < 1e-4) > 0)
sigStats2 = cbind(ebList$interaction$t[sigIndex2,-1], 
	ebList$interaction$p[sigIndex2,-1])
colnames(sigStats2) = c("t_age", "t_cell", "t_int", 
	"pval_age", "pval_cell", "pval_int")
colSums(sigStats2[,grep("pval",colnames(sigStats2))] < 1e-4)

tsSig$t_age_hom)
abline(0,1,col="red")