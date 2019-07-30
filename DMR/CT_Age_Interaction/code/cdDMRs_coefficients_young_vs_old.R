library(GenomicRanges)
library(GenomicFeatures)
library(data.table)
library(ggplot2)
library(RColorBrewer)
library(bumphunter)
library(limma)


load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/DMR/DMR_objects.rda")
load('/dcl01/lieber/ajaffe/lab/brain-epigenomics/bumphunting/BSobj_bsseqSmooth_Neuron_minCov_3.Rdata')


## Get CpGs and CpHs in each cdDMR

dtinteraction = data.table(DMR$Interaction[DMR$Interaction$sig=="FWER < 0.05",])

gr = granges(BSobj)
meth = bsseq:::getMeth(BSobj, type = 'raw')

oo = findOverlaps(reduce(makeGRangesFromDataFrame(dtinteraction)),gr)

rIndexes = split(subjectHits(oo), queryHits(oo))
splmeth = sapply(rIndexes, function(ii) colMeans(t(t(meth[ii,]))))
meanCG = do.call("rbind", splmeth)

load('/dcl01/lieber/ajaffe/lab/brain-epigenomics/bsseq/bsobj_by_chr/allChrs_postNatal_cleaned_nonCG_noHomogenate_highCov.Rdata')

gr = granges(BSobj)
meth = bsseq:::getMeth(BSobj, type = 'raw')

oo = findOverlaps(reduce(makeGRangesFromDataFrame(dtinteraction)),gr)

rIndexes = split(subjectHits(oo), queryHits(oo))
splmeth = sapply(rIndexes, function(ii) colMeans(t(t(meth[ii,]))))
meanCH = do.call("rbind", splmeth)


# separate young (0-6) from old (6+)

meanCG_young = meanCG[,colnames(meanCG) %in% rownames(pd[pd$Age<6,])]
meanCG_old = meanCG[,colnames(meanCG) %in% rownames(pd[pd$Age>6,])]

meanCH_young = meanCH[,colnames(meanCH) %in% rownames(pd[pd$Age<6,])]
meanCH_old = meanCH[,colnames(meanCH) %in% rownames(pd[pd$Age>6,])]


## Run regression

design <- with(pd[pd$Age<6,], model.matrix(~ Age * Cell.Type))
fit_youngCG <- lmFit(meanCG_young, design)
fit_youngCH <- lmFit(meanCH_young, design)

design <- with(pd[pd$Age>6,], model.matrix(~ Age * Cell.Type))
fit_oldCG <- lmFit(meanCG_old, design)
fit_oldCH <- lmFit(meanCH_old, design)


cg = data.frame(youngCG_Glia = fit_youngCG$coef[,2], oldCG_Glia = fit_oldCG$coef[,2], youngCG_Neuron = fit_youngCG$coef[,4], oldCG_Neuron = fit_oldCG$coef[,4])
ch = data.frame(youngCH_Glia = fit_youngCH$coef[,2], oldCH_Glia = fit_oldCH$coef[,2], youngCH_Neuron = fit_youngCH$coef[,4], oldCH_Neuron = fit_oldCH$coef[,4])

save(cg,ch,fit_youngCG,fit_youngCH,fit_oldCG,fit_oldCH,
     file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/DMR/young_v_old_coefficients.rda")

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/DMR/CT_Age_Interaction/figures/cdDMRs_coefficients_young_vs_old.pdf", 
    height = 5.5, width = 5.5)
ggplot(cg, aes(x = youngCG_Glia, y = oldCG_Glia)) + 
  geom_point() + theme_bw() + geom_hline(aes(yintercept=0)) + 
  geom_vline(aes(xintercept=0)) +
  xlab("Coefficients (0-6 years)") + ylab("Coefficients (6+ years)") + 
  ylim(-0.15,0.15) + xlim(-0.15,0.15) +
  ggtitle("CpG Glia Effect in cdDMRs") + 
  theme(title = element_text(size = 20),
        text = element_text(size = 20))
ggplot(cg, aes(x = youngCG_Neuron, y = oldCG_Neuron)) + 
  geom_point() + theme_bw() + geom_hline(aes(yintercept=0)) + 
  geom_vline(aes(xintercept=0)) +
  xlab("Coefficients (0-6 years)") + ylab("Coefficients (6+ years)") + 
  ylim(-0.15,0.15) + xlim(-0.15,0.15) +
  ggtitle("CpG Neuron Effect in cdDMRs") + 
  theme(title = element_text(size = 20),
        text = element_text(size = 20))
ggplot(ch, aes(x = youngCH_Glia, y = oldCH_Glia)) + 
  geom_point() + theme_bw() + geom_hline(aes(yintercept=0)) + 
  geom_vline(aes(xintercept=0)) +
  xlab("Coefficients (0-6 years)") + ylab("Coefficients (6+ years)") + 
  ylim(-0.06,0.06) + xlim(-0.06,0.06) +
  ggtitle("CpH Glia Effect in cdDMRs") + 
  theme(title = element_text(size = 20),
        text = element_text(size = 20))
ggplot(ch, aes(x = youngCH_Neuron, y = oldCH_Neuron)) + 
  geom_point() + theme_bw() + geom_hline(aes(yintercept=0)) + 
  geom_vline(aes(xintercept=0)) +
  xlab("Coefficients (0-6 years)") + ylab("Coefficients (6+ years)") + 
  ylim(-0.1,0.1) + xlim(-0.1,0.1) +
  ggtitle("CpH Neuron Effect in cdDMRs") + 
  theme(title = element_text(size = 20),
        text = element_text(size = 20))
plot(fit_youngCG$coef[,2], fit_oldCG$coef[,2])
plot(fit_youngCG$coef[,4], fit_oldCG$coef[,4])
plot(fit_youngCH$coef[,2], fit_oldCH$coef[,2])
plot(fit_youngCH$coef[,4], fit_oldCH$coef[,4])
dev.off()


## reviewer requested plots

# separate young (0-3) from old (3+) and run regression

meanCG_young = meanCG[,colnames(meanCG) %in% rownames(pd[pd$Age<3,])]
meanCG_old = meanCG[,colnames(meanCG) %in% rownames(pd[pd$Age>3,])]

meanCH_young = meanCH[,colnames(meanCH) %in% rownames(pd[pd$Age<3,])]
meanCH_old = meanCH[,colnames(meanCH) %in% rownames(pd[pd$Age>3,])]

design <- with(pd[pd$Age<3,], model.matrix(~ Age * Cell.Type))
fit_youngCG <- lmFit(meanCG_young, design)
fit_youngCH <- lmFit(meanCH_young, design)

design <- with(pd[pd$Age>3,], model.matrix(~ Age * Cell.Type))
fit_oldCG <- lmFit(meanCG_old, design)
fit_oldCH <- lmFit(meanCH_old, design)


cg3 = data.frame(youngCG_Glia = fit_youngCG$coef[,2], oldCG_Glia = fit_oldCG$coef[,2], 
                 youngCG_Neuron = fit_youngCG$coef[,4], oldCG_Neuron = fit_oldCG$coef[,4])
ch3 = data.frame(youngCH_Glia = fit_youngCH$coef[,2], oldCH_Glia = fit_oldCH$coef[,2], 
                 youngCH_Neuron = fit_youngCH$coef[,4], oldCH_Neuron = fit_oldCH$coef[,4])

# separate young (0-1) from old (1+) and run regression

meanCG_young = meanCG[,colnames(meanCG) %in% rownames(pd[pd$Age<1,])]
meanCG_old = meanCG[,colnames(meanCG) %in% rownames(pd[pd$Age>1,])]

meanCH_young = meanCH[,colnames(meanCH) %in% rownames(pd[pd$Age<1,])]
meanCH_old = meanCH[,colnames(meanCH) %in% rownames(pd[pd$Age>1,])]

design <- with(pd[pd$Age<1,], model.matrix(~ Age * Cell.Type))
fit_youngCG <- lmFit(meanCG_young, design)
fit_youngCH <- lmFit(meanCH_young, design)

design <- with(pd[pd$Age>1,], model.matrix(~ Age * Cell.Type))
fit_oldCG <- lmFit(meanCG_old, design)
fit_oldCH <- lmFit(meanCH_old, design)


cg1 = data.frame(youngCG_Glia = fit_youngCG$coef[,2], oldCG_Glia = fit_oldCG$coef[,2], 
                 youngCG_Neuron = fit_youngCG$coef[,4], oldCG_Neuron = fit_oldCG$coef[,4])
ch1 = data.frame(youngCH_Glia = fit_youngCH$coef[,2], oldCH_Glia = fit_oldCH$coef[,2], 
                 youngCH_Neuron = fit_youngCH$coef[,4], oldCH_Neuron = fit_oldCH$coef[,4])

path <- "/dcl01/lieber/ajaffe/lab/brain-epigenomics/DMR/CT_Age_Interaction/figures/"
pdf(paste0(path, "cdDMRs_coefficients_young_vs_old_different_ranges.pdf"), 
    height = 5.5, width = 5.5)
ggplot(cg3, aes(x = youngCG_Glia, y = oldCG_Glia)) + 
  geom_point() + theme_bw() + geom_hline(aes(yintercept=0)) + 
  geom_vline(aes(xintercept=0)) +
  xlab("Coefficients (0-3 years)") + ylab("Coefficients (3+ years)") + 
  ylim(-0.15,0.15) + xlim(-0.15,0.15) +
  ggtitle("CpG Glia Effect in cdDMRs") + 
  theme(title = element_text(size = 20),
        text = element_text(size = 20))
ggplot(cg3, aes(x = youngCG_Neuron, y = oldCG_Neuron)) + 
  geom_point() + theme_bw() + geom_hline(aes(yintercept=0)) + 
  geom_vline(aes(xintercept=0)) +
  xlab("Coefficients (0-3 years)") + ylab("Coefficients (3+ years)") + 
  ylim(-0.15,0.15) + xlim(-0.15,0.15) +
  ggtitle("CpG Neuron Effect in cdDMRs") + 
  theme(title = element_text(size = 20),
        text = element_text(size = 20))
ggplot(ch3, aes(x = youngCH_Glia, y = oldCH_Glia)) + 
  geom_point() + theme_bw() + geom_hline(aes(yintercept=0)) + 
  geom_vline(aes(xintercept=0)) +
  xlab("Coefficients (0-3 years)") + ylab("Coefficients (3+ years)") + 
  ylim(-0.06,0.06) + xlim(-0.06,0.06) +
  ggtitle("CpH Glia Effect in cdDMRs") + 
  theme(title = element_text(size = 20),
        text = element_text(size = 20))
ggplot(ch3, aes(x = youngCH_Neuron, y = oldCH_Neuron)) + 
  geom_point() + theme_bw() + geom_hline(aes(yintercept=0)) + 
  geom_vline(aes(xintercept=0)) +
  xlab("Coefficients (0-3 years)") + ylab("Coefficients (3+ years)") + 
  ylim(-0.1,0.1) + xlim(-0.1,0.1) +
  ggtitle("CpH Neuron Effect in cdDMRs") + 
  theme(title = element_text(size = 20),
        text = element_text(size = 20))

ggplot(cg1, aes(x = youngCG_Glia, y = oldCG_Glia)) + 
  geom_point() + theme_bw() + geom_hline(aes(yintercept=0)) + 
  geom_vline(aes(xintercept=0)) +
  xlab("Coefficients (0-1 years)") + ylab("Coefficients (1+ years)") + 
  ylim(-0.15,0.15) + xlim(-0.15,0.15) +
  ggtitle("CpG Glia Effect in cdDMRs") + 
  theme(title = element_text(size = 20),
        text = element_text(size = 20))
ggplot(cg1, aes(x = youngCG_Neuron, y = oldCG_Neuron)) + 
  geom_point() + theme_bw() + geom_hline(aes(yintercept=0)) + 
  geom_vline(aes(xintercept=0)) +
  xlab("Coefficients (0-1 years)") + ylab("Coefficients (1+ years)") + 
  ylim(-0.15,0.15) + xlim(-0.15,0.15) +
  ggtitle("CpG Neuron Effect in cdDMRs") + 
  theme(title = element_text(size = 20),
        text = element_text(size = 20))
ggplot(ch1, aes(x = youngCH_Glia, y = oldCH_Glia)) + 
  geom_point() + theme_bw() + geom_hline(aes(yintercept=0)) + 
  geom_vline(aes(xintercept=0)) +
  xlab("Coefficients (0-1 years)") + ylab("Coefficients (1+ years)") + 
  ylim(-0.06,0.06) + xlim(-0.06,0.06) +
  ggtitle("CpH Glia Effect in cdDMRs") + 
  theme(title = element_text(size = 20),
        text = element_text(size = 20))
ggplot(ch1, aes(x = youngCH_Neuron, y = oldCH_Neuron)) + 
  geom_point() + theme_bw() + geom_hline(aes(yintercept=0)) + 
  geom_vline(aes(xintercept=0)) +
  xlab("Coefficients (0-1 years)") + ylab("Coefficients (1+ years)") + 
  ylim(-0.1,0.1) + xlim(-0.1,0.1) +
  ggtitle("CpH Neuron Effect in cdDMRs") + 
  theme(title = element_text(size = 20),
        text = element_text(size = 20))
dev.off()


