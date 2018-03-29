library('ggplot2')

methCH = read.csv("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/non-CpG/proportion_mCH_byCellType_byAge_byTrinucleotideContext.csv")

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/non-CpG/figures/nonCG_global_methylation_byTrinucleotideContext_rearranged.pdf", width =7, useDingbats = FALSE)

sub <- methCH[which(methCH$perc=="ten" & methCH$trinucleotide_context %in% c("CAG","CAC")),]

ggplot(sub, 
       aes(x = Age, y = triprop, colour = CellType)) + geom_point() + facet_grid(trinucleotide_context ~ .) +
  geom_smooth(method = "lm") +  
  theme_classic() + scale_colour_brewer(8, palette="Dark2") +
  ylab("mCpH / CpH") + xlab("") +
  ggtitle("Global Proportion of CpHs Methylated >10%\nContext Proportion") + 
  theme(title = element_text(size = 18)) +
  theme(text = element_text(size = 20), legend.title=element_blank()) + theme(legend.position="bottom")
  
  
ggplot(sub, 
       aes(x = Age, y = prop, colour = CellType)) + geom_point() + facet_grid(trinucleotide_context ~ .) +
  geom_smooth(method = "lm") +  
  theme_classic() + scale_colour_brewer(8, palette="Dark2") +
  ylab("mCpH / CpH") + xlab("") +
  ggtitle("Global Proportion of CpHs Methylated >10%\nTotal Proportion") + 
  theme(title = element_text(size = 18)) +
  theme(text = element_text(size = 20), legend.title=element_blank()) + theme(legend.position="bottom")
dev.off()
