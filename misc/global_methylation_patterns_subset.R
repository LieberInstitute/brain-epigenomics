library('ggplot2')
library(cowplot)

methCH = read.csv("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/non-CpG/proportion_mCH_byCellType_byAge_byTrinucleotideContext.csv")
methCH = read.csv("../../Desktop/BAMS/proportion_mCH_byCellType_byAge_byTrinucleotideContext.csv")

pdf("./non-CpG/figures/nonCG_global_methylation_byTrinucleotideContext_rearranged.pdf", 
    width = 4.5, height = 3)

sub <- methCH[which(methCH$perc=="ten" & methCH$trinucleotide_context %in% c("CAG","CAC")),]

ggplot(sub, aes(x = Age, y = triprop, colour = CellType, shape = trinucleotide_context)) + 
  geom_point() +
  geom_smooth(method = "loess", se = F, aes(linetype = trinucleotide_context)) +  
  theme_classic() + 
  scale_colour_brewer(8, palette="Dark2") +
  scale_shape_manual(values = c(16, 1)) +
  ylab("(mCpH>10%)/CpH") + xlab("Years") +
  theme(title = element_text(size = 18),
        text = element_text(size = 20), 
        legend.title=element_blank(),
        legend.key = element_rect(colour = "transparent")) +
  guides(linetype = guide_legend(override.aes= list(color = "black")))

ggplot(sub, aes(x = Age, y = prop, colour = CellType, shape = trinucleotide_context)) + 
  geom_point() +
  geom_smooth(method = "loess", se = F, aes(linetype = trinucleotide_context)) +  
  theme_classic() + 
  scale_colour_brewer(8, palette="Dark2") +
  scale_shape_manual(values = c(16, 1)) +
  ylab("(mCpH>10%)/CpH") + xlab("Years") +
  theme(title = element_text(size = 18),
        text = element_text(size = 20), 
        legend.title=element_blank(),
        legend.key = element_rect(colour = "transparent")) +
  guides(linetype = guide_legend(override.aes= list(color = "black")))

sub <- methCH[which(methCH$perc=="nonzero" & methCH$trinucleotide_context %in% c("CAG","CAC")),]

ggplot(sub, aes(x = Age, y = triprop, colour = CellType, shape = trinucleotide_context)) + 
  geom_point() +
  geom_smooth(method = "loess", se = F, aes(linetype = trinucleotide_context)) +  
  theme_classic() + 
  scale_colour_brewer(8, palette="Dark2") +
  scale_shape_manual(values = c(16, 1)) +
  ylab("(mCpH>0%)/CpH") + xlab("Years") +
  theme(title = element_text(size = 18),
        text = element_text(size = 20), 
        legend.title=element_blank(),
        legend.key = element_rect(colour = "transparent")) +
  guides(linetype = guide_legend(override.aes= list(color = "black")))

ggplot(sub, aes(x = Age, y = prop, colour = CellType, shape = trinucleotide_context)) + 
  geom_point() +
  geom_smooth(method = "loess", se = F, aes(linetype = trinucleotide_context)) +  
  theme_classic() + 
  scale_colour_brewer(8, palette="Dark2") +
  scale_shape_manual(values = c(16, 1)) +
  ylab("(mCpH>0%)/CpH") + xlab("Years") +
  theme(title = element_text(size = 18),
        text = element_text(size = 20), 
        legend.title=element_blank(),
        legend.key = element_rect(colour = "transparent")) +
  guides(linetype = guide_legend(override.aes= list(color = "black")))

dev.off()
