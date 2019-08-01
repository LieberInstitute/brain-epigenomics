# Plot LDSC results for 'adjusting for baseline' analyses
# code adapted from Peter Hickey (https://github.com/hansenlab/BrainEpigenomeNN/blob/master/SLDSR/scripts/plot_ldsc.baseline_adjustments.R)
# qrsh -l bluejay,mem_free=100G,h_vmem=100G,h_fsize=500G

### ----------------------------------------------------------------------------
### Packages
###

library(GenomicRanges)
library(readr)
library(dplyr)
library(ggplot2)
library(scales)
library(gplots)
library(cowplot)
library(RColorBrewer)

###-----------------------------------------------------------------------------
### Category sizes
###

# ------------------------------------------------------------------------------
# Brain categories
#

brain_categories <- readRDS("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/ldsc/categories.rds")

brain_categories_df <- data.frame(Category = names(brain_categories),	
                                  Extended = c("Celltype (All)", "Age (All)", "cdDMRs (All)",
                                               "Cell Type (Glia > Neuron)", "Cell Type (Neuron > Glia)", 
                                               "Age (Younger > Older)", "Age (Older > Younger)", 
                                               "Group 1 (Decreasing Glial; Increasing Neuronal)", "Group 2 (Static Glial; Increasing Neuronal)", 
                                               "Group 3 (Static Glial; Decreasing Neuronal)", "Group 4 (Increasing Glial; Static Neuronal)", 
                                               "Group 5 (Increasing Glial; Decreasing Neuronal)", "Group 6 (Decreasing Glial; Static Neuronal)",
                                               "non-DMRs (same N; Cell Type (All))", "non-DMRs (same width; Cell Type (All))",
                                               "non-DMRs (same N; Age (All))", "non-DMRs (same width; Age (All))", 
                                               "non-DMRs (same N; cdDMRs (All))", "non-DMRs (same width; cdDMRs (All))",
                                               "non-DMRs (same N; Cell Type (Glia > Neuron))", "non-DMRs (same width; Cell Type (Glia > Neuron))",
                                               "non-DMRs (same N; Cell Type (Neuron > Glia))", "non-DMRs (same width; Cell Type (Neuron > Glia))", 
                                               "non-DMRs (same N; Age (Younger > Older))", "non-DMRs (same width; Age (Younger > Older))", 
                                               "non-DMRs (same N; Age (Older > Younger))", "non-DMRs (same width; Age (Older > Younger))", 
                                               "non-DMRs (same N; Group 1 cdDMRs)", "non-DMRs (same width; Group 1 cdDMRs)",
                                               "non-DMRs (same N; Group 2 cdDMRs)", "non-DMRs (same width; Group 2 cdDMRs)",
                                               "non-DMRs (same N; Group 3 cdDMRs)", "non-DMRs (same width; Group 3 cdDMRs)",
                                               "non-DMRs (same N; Group 4 cdDMRs)", "non-DMRs (same width; Group 4 cdDMRs)",
                                               "non-DMRs (same N; Group 5 cdDMRs)", "non-DMRs (same width; Group 5 cdDMRs)",
                                               "non-DMRs (same N; Group 6 cdDMRs)", "non-DMRs (same width; Group 6 cdDMRs)", 
                                               "non-DMRs", "All Prenatal LMRs", "All Neuronal LMRs", "All Glial LMRs",
                                               "Prenatal LMRs (Excluding DMRs)", "Neuronal LMRs (Excluding DMRs)","Glial LMRs (Excluding DMRs)", 
                                               "chromHMM (Union)", "chromHMM (Union; Excluding DMRs)", "CNS (LDSC)","CNS (LDSC; Excluding DMRs)"),
                                  Source = c(rep.int("DMRs",13), rep.int("non-DMRs",27), rep.int("Other",8), "Baseline", "Baseline"),
                                  Differential = c(rep.int("TRUE",13), rep.int("FALSE",37)))


# ------------------------------------------------------------------------------
# Record proportion of SNPs and CpGs in each category
#

load('/dcl01/lieber/ajaffe/lab/brain-epigenomics/bumphunting/BSobj_bsseqSmooth_Neuron_minCov_3.Rdata', verbose =T)

gr <- granges(BSobj)
categories <- brain_categories

categories_df <- bind_rows(
  lapply(names(categories), function(x) {
    data_frame(Category = x,
               `Total width (bp)` = sum(width(categories[[x]])),
               `Mean width (bp)` = mean(width(categories[[x]])),
               `Median width (bp)` = median(width(categories[[x]])),
               n = length(categories[[x]]))
  })) %>%
  arrange(`Total width (bp)`) %>%
  inner_join(brain_categories_df)

snp_prop_table <- bind_rows(
  lapply(names(brain_categories), function(bc) {
    x <- read_tsv(paste0("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/ldsc/output/",
                         bc,
                         ".Height.Phase1.results")) %>%
      filter(Category == "L2_0" | Category == "CNS_0")
    data_frame(Category = bc,
               `Prop._SNPs` = unlist(x[, "Prop._SNPs"]))
  }))
colnames(snp_prop_table)[2] <- "Proportion of SNPs"

cpg_prop_table <-
  data_frame(Category = names(categories),
             `Proportion of CpGs` = sapply(categories, function(x) {
               sum(overlapsAny(gr, x)) /
                 length(gr)}))

categories_df <- categories_df %>%
  inner_join(snp_prop_table, c("Category" = "Category")) %>%
  inner_join(cpg_prop_table, c("Category" = "Category")) %>%
  arrange(-Differential, `Total width (bp)`) %>%
  mutate(Category = factor(Category, Category, ordered = TRUE),
         Extended = factor(Extended, levels = Extended,
                                    ordered = TRUE)) %>%
  arrange(Extended)


traits_df <- openxlsx::read.xlsx("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/ldsc/GWAS_summary_statistics_sources.xlsx", 
                                 startRow = 3)
traits_df$`Pretty Trait` = traits_df$Trait
traits_df$Trait = c("ADHD", "Alzheimers_disease","Anorexia_nervosa", "Anxiety_disorder",
                       "Autism_spectrum_disorder", "Bipolar_disorder", "BMI", "Childhood_cognitive_performance",
                       "Cigarettes_per_day", "College_attainment", "Conscientiousness", "Coronary_artery_disease", 
                       "Crohns_disease", "Depressive_symptoms", "Epilepsy", "Ever_smoked", "Extraversion", "Focal_epilepsy", 
                       "Generalized_epilepsy", "Height", "Intracarebral_hemorrhage", "IQ", "Ischemic_stroke", "Major_depressive_disorder",
                       "Neuroticism", "Openness", "PTSD", "Schizophrenia", "Subjective_well-being", "Years_of_education")
traits_df = traits_df[,-which(colnames(traits_df) %in% c("Link.to.GWAS.data", "Filename", "MD5.checksum", "Publication"))]     


## ----------------------------------------------------------------------------
### Load data and construct objects
###

fls <- unlist(lapply(names(brain_categories), function(bc) {
  fls <- list.files("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/ldsc/output",
                    pattern = glob2rx(
                      paste0(bc, ".*Phase1.results")),
                    full.names = TRUE)
  grep("adjusting", fls, invert = TRUE, value = TRUE)
}))

# Read in files, tidy up, and rbind
x <- bind_rows(lapply(fls, function(fl) {
  suppressMessages(read_tsv(fl)) %>%
    filter(Category == "L2_0" | Category == "CNS_0") %>%
    mutate(Category = sapply(strsplit(basename(fl), "\\."), "[[", 1),
           Trait = sapply(strsplit(sub(".Phase1.results", "", basename(fl)),
                                   "\\."),
                          "[[", 2),
           lower = Enrichment - 1.96 * Enrichment_std_error,
           upper = Enrichment + 1.96 * Enrichment_std_error,
           file = fl) %>%
    mutate(Category = factor(Category, Category, ordered = TRUE))
}))

# NOTE: Anttila report these traits "had in sufficient evidence of additive
#       heritability for robust analysis" and excluded them from further
#       analysis
x <- x %>%
  filter(Trait != "Agreeableness",
         Trait != "Cardioembolic_stroke",
         Trait != "Large-vessel_disease",
         Trait != "Small-vessel_disease")
fls <- fls[-c(grep("Agreeableness", fls), grep("Cardioembolic_stroke", fls), 
              grep("Large-vessel_disease", fls), grep("Small-vessel_disease", fls))]

# also use updated summary statistics for 3 psychiatric diseases since Feinberg paper
x <- x %>%
  filter(Trait != "Autism_spectrum_disorder",
         Trait != "Bipolar_disorder",
         Trait != "Major_depressive_disorder")
fls <- fls[-c(grep("Autism_spectrum_disorder.", fls, fixed=T), grep("Bipolar_disorder.", fls, fixed=T), 
              grep("Major_depressive_disorder.", fls, fixed=T))]


# Join munged LDSC output with categories_df and traits_df

x$Trait = gsub("_latest","", x$Trait)
x <- x %>%
  inner_join(categories_df, c("Category" = "Category")) %>%
  inner_join(traits_df, c("Trait" = "Trait")) %>%
  arrange(Extended)

stopifnot(length(fls) == nrow(x))

write.csv(x, file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/ldsc/tables/baseline_ldsc_output_allIncluded.csv", quote=F)


# Add adjusted P-values
x = x[which(x$Category %in% c("Celltype_HypoNeuron", "Celltype_HypoGlia", "Age_Decreasing", "Age_Increasing",
                              "Gr1_cdDMRs", "Gr2_cdDMRs", "Gr3_cdDMRs", "Gr4_cdDMRs", "Gr5_cdDMRs", "Gr6_cdDMRs", 
                              "non_DMRs", "Prenatal_LMRs_all", "Neuronal_LMRs_all", "Glial_LMRs_all", "chromHMM_union", "CNS")),] # not including unsplit DMR groups                                  

x <- x %>%
  group_by(Category, file) %>%
  filter(grepl(Category, file)) %>%
  ungroup() %>%
  mutate(Coefficient_p = pnorm(`Coefficient_z-score`, lower.tail = FALSE)) %>%
  group_by(Trait) %>%
  mutate(Coefficient_holm = p.adjust(Coefficient_p, method = "holm"),
         Enrichment_holm = p.adjust(Enrichment_p, method = "holm"),
         Coefficient_holm_cutoff =
           max(Coefficient_p[Coefficient_holm < 0.05], na.rm = TRUE),
         Enrichment_holm_cutoff =
           max(Enrichment_p[Enrichment_holm < 0.05], na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(sig_coef = Coefficient_holm < 0.05) %>%
  arrange(-Differential, `Total width (bp)`) %>%
  mutate(Extended = factor(Extended,
                                    unique(Extended),
                                    ordered = TRUE),
         `Pretty Trait` = gsub("_", " ", Trait)) %>%
  arrange(Extended)

write.csv(x, file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/ldsc/tables/baseline_ldsc_output.csv", quote=F)

### ----------------------------------------------------------------------------
### Stratify traits by 'Brain-linked (sig)', 'Brain-linked (non-sig)', or
### 'Non-brain-linked' (stratification defined in 'baseline' analysis)
###

x = read.csv("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/ldsc/tables/baseline_ldsc_output.csv")

strata_df <- x %>%
  group_by(Trait) %>%
  summarise(strata = ifelse(any(Coefficient_holm < 0.05),
                            "Brain-linked (sig)",
                            "Brain-linked (non-sig)"),
            Type = unique(Type),
            strata = ifelse(Type == "Additional_phenotype",
                            "Non-brain-linked",
                            strata)) %>%
  dplyr::select(Trait, strata)
x_stratified <- inner_join(x, strata_df)

### ----------------------------------------------------------------------------
### Tables of results
###

x_stratified %>%
  dplyr::select(-Category, -lower, -upper, -n, -`Prop._SNPs`,
                -file, -Source, -Differential,
                -Enrichment_p,
                -Enrichment_holm, -Enrichment_holm_cutoff,
                -`Mean.width..bp.`, -`Median.width..bp.`,
                -`Proportion.of.CpGs`,
                -Coefficient_holm_cutoff, -sig_coef, -Trait,
                -N, -N_cases, -N_controls, -Type) %>%
  dplyr::select(`Pretty.Trait`, strata,
                Extended, `Total.width..bp.`, `Proportion.of.SNPs`,
                starts_with("Coefficient"), everything()) %>%
  rename(Feature = Extended,
         Trait = `Pretty.Trait`,
         `Total width (bp)` = `Total.width..bp.`,
         `Proportion of h2` = `Prop._h2`,
         `Proportion of h2 standard error` = `Prop._h2_std_error`,
         Stratum = `strata`) %>%
  write_csv("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/ldsc/tables/LDSC_results.baseline_adjustments.csv")

x_stratified %>%
  dplyr::select(`Pretty.Trait`, Type, N, N_cases, N_controls) %>%
  distinct() %>%
  rename(Trait = `Pretty.Trait`,
         `Type` = `Type`) %>%
  write_csv("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/ldsc/tables/annotated_traits_df.csv")


### ----------------------------------------------------------------------------
### Plots
###

# ------------------------------------------------------------------------------
# Plot dependent variables by feature, trait
#
x_stratified$Extended = gsub("All ", "", x_stratified$Extended)
x_stratified$Extended = factor(x_stratified$Extended, levels = c("Cell Type (Glia > Neuron)","Cell Type (Neuron > Glia)",
                                                                 "Age (Older > Younger)", "Age (Younger > Older)", 
                                                                 "Group 1 (Decreasing Glial; Increasing Neuronal)",
                                                                 "Group 2 (Static Glial; Increasing Neuronal)",
                                                                 "Group 3 (Static Glial; Decreasing Neuronal)",
                                                                 "Group 4 (Increasing Glial; Static Neuronal)",    
                                                                 "Group 5 (Increasing Glial; Decreasing Neuronal)",
                                                                 "Group 6 (Decreasing Glial; Static Neuronal)",
                                                                 "Prenatal LMRs","Neuronal LMRs","Glial LMRs",
                                                                 "CNS (LDSC)","chromHMM (Union)","non-DMRs"))
x_stratified$Differential = ifelse(x_stratified$Differential==TRUE, "DMR", "Other")
x_stratified$sig_coef = ifelse(x_stratified$sig_coef==TRUE, "Significant","Not significant")
  

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/LDSC/figures/Exploratory_boxplots_byFeature_baseline.pdf",
    height = 4.5, width = 11)

mapply(function(dep, nm) {
  g = ggplot(x_stratified, aes(x = Extended, y = get(dep), fill = Extended)) + 
  geom_jitter() + geom_boxplot(outlier.shape = NA) +
  facet_grid(. ~ Differential, scales ="free") +
  theme_bw() +
  geom_hline(yintercept = 0, lty = 2) +
  labs(fill="") + ylab(nm) + xlab("") +
  theme(title = element_text(size = 20), 
        text = element_text(size = 20),
        legend.title = element_blank(), 
        axis.text.x = element_blank())
  print(g) },
  c("Prop._h2","Enrichment", "Coefficient_z.score","Coefficient_p","Total.width..bp."),
  c("Heritability Explained", "Enrichment", "Coefficient Z Score", 
    "Coefficient P Value", "Total width (bp)"))
ggplot(x_stratified, aes(Coefficient_p, fill = Extended, colour=Extended)) + 
  geom_density(alpha=.2) +
  facet_grid(. ~ Differential, scales ="free") +
  theme_bw() +
  labs(fill="") + ylab("Density") +
  xlab("Coefficient P Value") + xlab("") +
  theme(title = element_text(size = 20), 
        text = element_text(size = 20),
        legend.title = element_blank()) +
  guides(colour=FALSE)
dev.off()


pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/LDSC/figures/Exploratory_boxplots_byTrait_baseline.pdf",
    height = 4.5, width = 16)

mapply(function(dep, nm) {
  g = ggplot(x_stratified, aes(x = Pretty.Trait, y = get(dep))) + 
    geom_jitter() + geom_boxplot(outlier.shape = NA) +
    theme_bw() +
    geom_hline(yintercept = 0, lty = 2) +
    labs(fill="") + ylab(nm) + xlab("") +
    theme(title = element_text(size = 20), 
          text = element_text(size = 20),
          legend.title = element_blank(), 
          axis.text.x = element_text(angle = 45, hjust = .95))
  print(g) },
  c("Prop._h2","Enrichment", "Coefficient_z.score","Coefficient_p"),
  c("Heritability Explained", "Enrichment", "Coefficient Z Score", 
    "Coefficient P Value"))
mapply(function(dep, nm) {
  g = ggplot(x_stratified, aes(x = Pretty.Trait, y = get(dep), fill = Pretty.Trait)) + 
    geom_jitter() + geom_boxplot(outlier.shape = NA) +
    facet_grid(. ~ Differential, scales ="free") +
    theme_bw() +
    geom_hline(yintercept = 0, lty = 2) +
    labs(fill="") + ylab(nm) + xlab("") +
    theme(title = element_text(size = 20), 
          text = element_text(size = 20),
          legend.title = element_blank(), 
          axis.text.x = element_blank())
  print(g) },
  c("Prop._h2","Enrichment", "Coefficient_z.score","Coefficient_p"),
  c("Heritability Explained", "Enrichment", "Coefficient Z Score", 
    "Coefficient P Value"))
ggplot(x_stratified, aes(Coefficient_p, fill = Pretty.Trait, colour=Pretty.Trait)) + 
  geom_density(alpha=.2) +
  facet_grid(. ~ Differential, scales ="free") +
  theme_bw() +
  labs(fill="") + ylab("Density") +
  xlab("Coefficient P Value") + xlab("") +
  theme(title = element_text(size = 20), 
        text = element_text(size = 20),
        legend.title = element_blank()) +
  guides(colour=FALSE)
ggplot(x_stratified, aes(Coefficient_p, fill = Pretty.Trait, colour=Pretty.Trait)) + 
  geom_density(alpha=.2) +
  theme_bw() +
  labs(fill="") + ylab("Density") +
  xlab("Coefficient P Value") + xlab("") +
  theme(title = element_text(size = 20), 
        text = element_text(size = 20),
        legend.title = element_blank()) +
  guides(colour=FALSE)
dev.off()


# ------------------------------------------------------------------------------
# Coefficient Z-score
#

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/LDSC/figures/Coefficient_Z-score.baseline_adjustments.pdf",
    height = 6.5, width = 14)

ggplot(x_stratified, aes(x = Extended, y = Coefficient_z.score, 
                         col = Extended, shape = sig_coef, size = sig_coef)) + 
  geom_point() +
  facet_wrap( ~ `Pretty.Trait`, ncol = 5) +
  theme_bw() +
  xlab("Feature") + ylab("Coefficient Z Score") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.title = element_blank()) +
  geom_hline(yintercept = 0, lty = 2) +
  scale_shape_manual(values = c(1, 16)) +
  scale_size_manual(values = c(2, 3)) +
  scale_colour_manual(values = c("#FC8D62","#66C2A5","#377EB8","#E41A1C",
                                 "#1B9E77","#D95F02","#7570B3","#E7298A","#66A61E","#E6AB02",
                                 "#A65628","#F781BF","#FFFF33","#8DA0CB","#666666","#E5C494"))
dev.off()


pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/LDSC/figures/Coefficient_Z-score.baseline_adjustments.stratified.pdf",
    height = 4, width = 5)

ggplot(data = x_stratified, aes(x = Extended, y = `Coefficient_z.score`,
             col = Extended, shape = sig_coef, size = sig_coef)) +
  geom_jitter(width = 0.3) +
  facet_grid(. ~ strata) +
  theme_bw() + xlab("Feature") + ylab("Coefficient Z Score") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  geom_hline(yintercept = 0, lty = 2) +
  scale_shape_manual(values = c(1, 16)) +
  scale_size_manual(values = c(2, 3)) +
  scale_colour_manual(values = c("#FC8D62","#66C2A5","#377EB8","#E41A1C",
                                 "#1B9E77","#D95F02","#7570B3","#E7298A","#66A61E","#E6AB02",
                                 "#A65628","#F781BF","#FFFF33","#8DA0CB","#666666","#E5C494")) +
  guides(col = FALSE, shape = FALSE, size = FALSE)
dev.off()


# ------------------------------------------------------------------------------
# Enrichment
#

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/LDSC/figures/Enrichment.baseline_adjustments.pdf", height = 8.5, width = 10)

ggplot(data = x_stratified[which(x_stratified$strata=="Brain-linked (sig)"),],
         aes(x = Extended, y = Enrichment, col = Extended,
             shape = sig_coef)) +
  geom_point() +
  geom_pointrange(aes(ymin = Enrichment - 2 * Enrichment_std_error,
                      ymax = Enrichment + 2 * Enrichment_std_error)) +
  facet_wrap( ~ `Pretty.Trait`, ncol = 4, scales="free") +
  theme_bw() + xlab("") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.title = element_blank(), legend.position = c(0.625, 0.115)) +
  scale_shape_manual(values = c(1, 16)) +
  geom_hline(yintercept = 0, lty = 2) +
  scale_colour_manual(values = c("#FC8D62","#66C2A5","#377EB8","#E41A1C",
                                 "#1B9E77","#D95F02","#7570B3","#E7298A","#66A61E","#E6AB02",
                                 "#A65628","#F781BF","#FFFF33","#8DA0CB","#666666","#E5C494")) +
  guides(size=guide_legend(ncol=2),shape=guide_legend(ncol=2),col=guide_legend(ncol=3))

dev.off()


pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/LDSC/figures/Enrichment.baseline_adjustments.sig_stratified.pdf", height = 4, width = 5)

ggplot(data = x_stratified,
       aes(x = Extended, y = Enrichment, col = Extended,
           shape = sig_coef, size = sig_coef)) +
  geom_jitter(width = 0.2) +
  facet_grid(. ~ strata, labeller = labeller(sig = label_both)) +
  theme_bw() + xlab("Feature") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.title = element_blank()) +
  geom_hline(yintercept = 0, lty = 2) +
  scale_shape_manual(values = c(1, 16)) +
  scale_size_manual(values = c(2, 3)) +
  scale_colour_manual(values = c("#FC8D62","#66C2A5","#377EB8","#E41A1C",
                                 "#1B9E77","#D95F02","#7570B3","#E7298A","#66A61E","#E6AB02",
                                 "#A65628","#F781BF","#FFFF33","#8DA0CB","#666666","#E5C494")) +
  guides(col=FALSE, shape = FALSE, size = FALSE)
ggplot(data = x_stratified[which(x_stratified$Enrichment<2000),],
       aes(x = Extended, y = Enrichment, col = Extended,
           shape = sig_coef, size = sig_coef)) +
  geom_jitter(width = 0.2) +
  facet_grid(. ~ strata, labeller = labeller(sig = label_both)) +
  theme_bw() + xlab("Feature") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.title = element_blank()) +
  geom_hline(yintercept = 0, lty = 2) +
  scale_shape_manual(values = c(1, 16)) +
  scale_size_manual(values = c(2, 3)) +
  scale_colour_manual(values = c("#FC8D62","#66C2A5","#377EB8","#E41A1C",
                                 "#1B9E77","#D95F02","#7570B3","#E7298A","#66A61E","#E6AB02",
                                 "#A65628","#F781BF","#FFFF33","#8DA0CB","#666666","#E5C494")) +
  guides(col=FALSE, shape = FALSE, size = FALSE)
ggplot(data = x_stratified,
         aes(x = Extended, y = log(Enrichment), col = Extended,
             shape = sig_coef, size = sig_coef)) +
  geom_jitter(width = 0.2) +
  facet_grid(. ~ strata, labeller = labeller(sig = label_both)) +
  theme_bw() + xlab("Feature") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.title = element_blank()) +
  geom_hline(yintercept = 0, lty = 2) +
  scale_shape_manual(values = c(1, 16)) +
  scale_size_manual(values = c(2, 3)) +
  scale_colour_manual(values = c("#FC8D62","#66C2A5","#377EB8","#E41A1C",
                                 "#1B9E77","#D95F02","#7570B3","#E7298A","#66A61E","#E6AB02",
                                 "#A65628","#F781BF","#FFFF33","#8DA0CB","#666666","#E5C494")) +
  guides(col=FALSE, shape = FALSE, size = FALSE)
dev.off()

### ----------------------------------------------------------------------------
### Create legend used
###

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/LDSC/figures/Feature.color.legend.pdf", height = 6, width = 6)
g = ggplot(data = x_stratified,
         aes(x = Extended, y = -log10(Coefficient_p),
             col = Extended, shape = sig_coef, size = sig_coef)) +
  geom_point() +
  facet_wrap( ~ Trait, ncol = 5) +
  # Holm's cutoff
  geom_hline(aes(yintercept = -log10(Coefficient_holm_cutoff))) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.title = element_blank()) +
  scale_colour_manual(values = c("#FC8D62","#66C2A5","#377EB8","#E41A1C",
                                 "#1B9E77","#D95F02","#7570B3","#E7298A","#66A61E","#E6AB02",
                                 "#A65628","#F781BF","#FFFF33","#8DA0CB","#666666","#E5C494")) +
  scale_shape_manual(values = c(1, 16)) +
  scale_size_manual(values = c(2, 3))
ggdraw(plot_grid(NULL, get_legend(g)))
dev.off()

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/LDSC/figures/Feature.color.legend.longways.pdf", height = 6, width = 22)
x_releveled = x_stratified
x_releveled$Extended = factor(x_releveled$Extended, 
                              levels = c("Prenatal LMRs",                              
                                         "Neuronal LMRs",                                  
                                         "Glial LMRs",                                     
                                         "CNS (LDSC)",                                     
                                         "chromHMM (Union)",                               
                                         "non-DMRs",
                                         "Group 1 (Decreasing Glial; Increasing Neuronal)",
                                         "Group 2 (Static Glial; Increasing Neuronal)",    
                                         "Group 3 (Static Glial; Decreasing Neuronal)",   
                                         "Group 4 (Increasing Glial; Static Neuronal)",    
                                         "Group 5 (Increasing Glial; Decreasing Neuronal)",
                                         "Group 6 (Decreasing Glial; Static Neuronal)",
                                         "Cell Type (Glia > Neuron)",                      
                                         "Cell Type (Neuron > Glia)",                      
                                         "Age (Older > Younger)",                          
                                         "Age (Younger > Older)"))

g = ggplot(data = x_releveled,
           aes(x = Extended, y = -log10(Coefficient_p),
               col = Extended, shape = sig_coef, size = sig_coef)) +
  geom_point() +
  facet_wrap( ~ Trait, ncol = 5) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.title = element_blank(), legend.position = "bottom") +
  scale_colour_manual(values = c("#A65628","#F781BF","#FFFF33","#8DA0CB","#666666","#E5C494",
                                 "#1B9E77","#D95F02","#7570B3","#E7298A","#66A61E","#E6AB02",
                                 "#FC8D62","#66C2A5","#377EB8","#E41A1C")) +
  scale_shape_manual(values = c(1, 16)) +
  scale_size_manual(values = c(2, 3)) +
  guides(size=guide_legend(nrow=2),shape=guide_legend(nrow=2),
         col=guide_legend(ncol=3))
ggdraw(plot_grid(NULL, get_legend(g)))
dev.off()



### ----------------------------------------------------------------------------
### Additional plots
###

x_stratified$Type = gsub("Additional_phenotype", "Other", x_stratified$Type)
x_stratified$Type = factor(x_stratified$Type, levels = c("Psychiatric","Behavioural-cognitive","Neurological","Other"))


pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/LDSC/figures/Coefficient_Z-score.baseline_byType.pdf",
    height = 2.5, width = 6.75)
ggplot(data = x_stratified, aes(x = Extended, y = `Coefficient_z.score`,
                                col = Extended, shape = sig_coef, size = sig_coef)) +
  geom_jitter(width = 0.3) +
  facet_grid(. ~ Type) +
  theme_bw() + xlab("Feature") + ylab("Coefficient Z Score") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  geom_hline(yintercept = 0, lty = 2) +
  scale_shape_manual(values = c(1, 16)) +
  scale_size_manual(values = c(2, 3)) +
  scale_colour_manual(values = c("#FC8D62","#66C2A5","#377EB8","#E41A1C",
                                 "#1B9E77","#D95F02","#7570B3","#E7298A","#66A61E","#E6AB02",
                                 "#A65628","#F781BF","#FFFF33","#8DA0CB","#666666","#E5C494")) +
  guides(col = FALSE, shape = FALSE, size = FALSE)

ggplot(data = x_stratified, aes(x = Extended, y = `Coefficient_z.score`,
                                col = Extended)) +
  geom_boxplot(outlier.shape=NA) + 
  geom_jitter(data = x_stratified[which(x_stratified$sig_coef=="Significant"),],
              width = 0.3, size=3,aes(x = Extended, y = `Coefficient_z.score`)) + 
  facet_grid(. ~ Type) +
  theme_bw() + xlab("Feature") + ylab("Coefficient Z Score") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  geom_hline(yintercept = 0, lty = 2) +
  scale_colour_manual(values = c("#FC8D62","#66C2A5","#377EB8","#E41A1C",
                                 "#1B9E77","#D95F02","#7570B3","#E7298A","#66A61E","#E6AB02",
                                 "#A65628","#F781BF","#FFFF33","#8DA0CB","#666666","#E5C494")) +
  guides(col = FALSE, shape = FALSE, size = FALSE)
dev.off()

x_stratified$three = ifelse(x_stratified$Type=="Psychiatric", "Psychiatric Trait","Other Brain Trait")
x_stratified[x_stratified$Type=="Other","three"] = "Non-Brain Trait"
x_stratified$three = factor(x_stratified$three, levels = c("Psychiatric Trait","Other Brain Trait","Non-Brain Trait"))

pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/LDSC/figures/Coefficient_Z-score.baseline_3panel.pdf",
    height = 2.5, width = 5)
ggplot(data = x_stratified, aes(x = Extended, y = `Coefficient_z.score`,
                                col = Extended, shape = sig_coef, size = sig_coef)) +
  geom_jitter(width = 0.3) +
  facet_grid(. ~ three) +
  theme_bw() + xlab("Feature") + ylab("Coefficient Z Score") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  geom_hline(yintercept = 0, lty = 2) +
  scale_shape_manual(values = c(1, 16)) +
  scale_size_manual(values = c(1, 2)) +
  scale_colour_manual(values = c("#FC8D62","#66C2A5","#377EB8","#E41A1C",
                                 "#1B9E77","#D95F02","#7570B3","#E7298A","#66A61E","#E6AB02",
                                 "#A65628","#F781BF","#FFFF33","#8DA0CB","#666666","#E5C494")) +
  guides(col = FALSE, shape = FALSE, size = FALSE)

ggplot(data = x_stratified[which(x_stratified$sig_coef=="Significant"),],
       aes(x = Extended, y = log(Enrichment), col = Extended)) +
  geom_jitter(width = 0.2, size=2, shape=16) +
  facet_grid(. ~ three, labeller = labeller(sig = label_both)) +
  theme_bw() + xlab("Feature") +
  scale_x_discrete(drop=FALSE) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.title = element_blank()) +
  geom_hline(yintercept = 0, lty = 2) +
  scale_colour_manual(values = c("#FC8D62","#7570B3","#E7298A","#A65628","#F781BF","#8DA0CB","#666666")) +
  guides(col=FALSE, shape = FALSE, size = FALSE)

ggplot(data = x_stratified[which(x_stratified$sig_coef=="Significant"),],
       aes(x = Extended, y = Enrichment, col = Extended)) +
  geom_jitter(width = 0.2, size=2, shape=16) +
  geom_pointrange(aes(ymin = Enrichment - 2 * Enrichment_std_error,
                      ymax = Enrichment + 2 * Enrichment_std_error)) +
  facet_grid(. ~ three, labeller = labeller(sig = label_both)) +
  theme_bw() + xlab("Feature") +
  scale_x_discrete(drop=FALSE) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.title = element_blank()) +
  geom_hline(yintercept = 0, lty = 2) +
  scale_colour_manual(values = c("#FC8D62","#7570B3","#E7298A","#A65628","#F781BF","#8DA0CB","#666666")) +
  guides(col=FALSE, shape = FALSE, size = FALSE)
dev.off()

df = x_stratified[,c("Extended","Coefficient_z.score","sig_coef","Enrichment","three")]
df$Enrichment = log(df$Enrichment)
df = reshape2::melt(df)
df$variable = gsub("Coefficient_z.score", "Coefficient Z Score", df$variable)
df$variable = gsub("Enrichment", "log(Enrichment)", df$variable)
df$variable = factor(df$variable)
df$var2 = paste(df$sig_coef, df$variable, sep=":")
df = df[which(df$var2!="Not significant:log(Enrichment)"),]
head(df)


pdf("/dcl01/lieber/ajaffe/lab/brain-epigenomics/LDSC/figures/Coefficient_Z-score.baseline_3panel_faceted.pdf",
    height = 3.5, width = 5.5)
ggplot(data = df, aes(x = Extended, y = value,
                                col = Extended, shape = sig_coef, size = sig_coef)) +
  geom_jitter(width = 0.3) +
  facet_grid(variable ~ three, scales = "free_y") +
  theme_bw() + xlab("Feature") + ylab("") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  geom_hline(yintercept = 0, lty = 2) +
  scale_shape_manual(values = c(1, 16)) +
  scale_size_manual(values = c(1, 2)) +
  scale_colour_manual(values = c("#FC8D62","#66C2A5","#377EB8","#E41A1C",
                                 "#1B9E77","#D95F02","#7570B3","#E7298A","#66A61E","#E6AB02",
                                 "#A65628","#F781BF","#FFFF33","#8DA0CB","#666666","#E5C494")) +
  guides(col = FALSE, shape = FALSE, size = FALSE)
dev.off()


# Gr 4 neuronal precursor genes include RERE, PROX1, ID2, HHEX, PAX6, ARX
# Gr 4 cell-cell adhesion via plasma-membrane adhesion molecules include "NLGN1"    "PCDHGA1"  "PCDHGA2"  "PCDHGA3"  "PCDHGB1"  "PCDHGA4"  "PCDHGB2"  "PCDHGA5"  "PCDHGB3"  "PCDHGA6"  "PCDHGA7"  "PCDHGB4" 
# "PCDHGA8"  "PCDHGA9"  "PCDHGB6"  "PCDHGA10" "PCDHGB7"  "PCDHGA11" "PCDHGA12" "PCDHGC3"  "PCDHGC4"  "PCDHGC5"  "SDK1"

# Gr 3 modulation of synaptic chemical transmission includes "PRKCZ"  "CLSTN1" "PREPL"  "NISCH"  "CAMK2A" "SYNPO"  "HRH2"   "BAI1"   "GRIN1"  "SYT7"   "IGSF9B" "CLSTN3" "STX1B"  "BAIAP2" "RAB3A" 
# "SNAP25" "SYN1"

library(data.table)

sig = x_stratified[which(x_stratified$sig_coef=="TRUE"),]
sig = sig[,colnames(sig) %in% c("Coefficient_z.score", "Enrichment", "Trait","Category", "Type")]
sig = data.table(sig)
sig$outside = ifelse(sig$Category %in% c("chromHMM_union", "CNS"), "Yes", "No")
sig$typeout = ifelse(sig$Type=="Additional_phenotype", "Yes", "No")
sig

sig[,length(unique(Trait)), by = c("outside","typeout")]
sig[grep("Neuron", sig$Category),length(unique(Trait)),by="typeout"]
sig[grep("Neuron", sig$Category),mean(Enrichment),by=c("typeout", "Category")]
sig[grep("Prenatal", sig$Category),,]
sig[Trait %in% c("Major_depressive_disorder","Neuroticism", "IQ"),,]

sig[Trait=="BMI",,]


x_stratified = data.table(x_stratified)
library(data.table)
sig = data.table(x_stratified[x_stratified$sig_coef=="TRUE",c("Category","Prop._h2","Trait",
                                                      "Enrichment","Type","Coefficient_z.score")])

sig[Trait=="BMI",,]
x_stratified[Trait=="BMI" & Category %in% c("Neuronal_LMRs_all","Prenatal_LMRs_all"),
             c("Category","Prop._h2","Trait","Enrichment","Type",
               "Coefficient_z.score", "Coefficient_holm"),]
sig[grep("Prenatal", Category),,]

unique(x_stratified[,c("Category","Total.width..bp.")])
x_stratified$cattype = ifelse(x_stratified$Category %in% c("Glial_LMRs_all",
                                                           "Neuronal_LMRs_all","Prenatal_LMRs_all"),
                              "LMR", "DMR")
x_stratified[Category=="non_DMRs", "cattype"] = "nonDMR"
x_stratified[Category %in% c("CNS","chromHMM_union"), "cattype"] = "General"
x_stratified[grep("Cell", Category), "cattype"] = "CellType"
x_stratified[grep("Gr", Category), "cattype"] = "cdDMR"

x_stratified[,mean(unique(Total.width..bp.)), by ="cattype"]$V1/499322.7
x_stratified[,mean(unique(Proportion.of.SNPs)), by ="cattype"]$V1/1.692933e-04

sig[grep("Gr", Category),,]
x_stratified[Category=="Gr3_cdDMRs" & Trait=="Schizophrenia",,]

x_stratified[Category=="Gr4_cdDMRs" & Trait=="PTSD",,]


