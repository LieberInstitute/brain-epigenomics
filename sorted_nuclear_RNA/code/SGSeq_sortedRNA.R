library("SGSeq")
library("GenomicFeatures")
library("GenomicRanges")
library("ggplot2")
library("DEXSeq")
library(data.table)
library(reshape2)


# For SGSeq, create sample information dataframe

file_bam = paste0("/dcl01/lieber/ajaffe/Amanda/sorted_nuclear_RNA/merged/", frag_length$sample_name, "_combinedLibraries.bam")
file_bam = paste0("./Desktop/BAMS/", frag_length$sample_name, "_combinedLibraries.bam")

si = data.frame(sample_name = frag_length$sample_name, file_bam = file_bam)
si$sample_name = as.character(si$sample_name)
si$file_bam = as.character(si$file_bam)
si = getBamInfo(si)


# Extract transcript features
gencode = loadDb("/dcl01/lieber/ajaffe/Amanda/sorted_nuclear_RNA/gencode.v25lift37.annotation.sqlite")
gencode = loadDb("./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/SGSeq_out/gencode.v25lift37.annotation.sqlite")

txf = convertToTxFeatures(gencode)
sgf = convertToSGFeatures(txf)
notinheader = c("chr6_ssto_hap7", "chr6_mcf_hap5", "chr6_cox_hap2", "chr6_mann_hap4", "chr6_apd_hap1", "chr6_qbl_hap6", "chr6_dbb_hap3", "chr17_ctg5_hap1", 
                "chr4_ctg9_hap1", "chrUn_gl000225", "chr4_gl000194_random", "chr9_gl000200_random", "chrUn_gl000222", "chrUn_gl000223", "chrUn_gl000224", 
                "chrUn_gl000219", "chrUn_gl000215", "chrUn_gl000216", "chrUn_gl000217", "chrUn_gl000211", "chrUn_gl000213", "chrUn_gl000218", 
                "chr19_gl000209_random", "chrUn_gl000221", "chrUn_gl000214", "chrUn_gl000227", "chr1_gl000191_random", "chr19_gl000208_random", 
                "chr9_gl000198_random", "chrUn_gl000233", "chrUn_gl000237", "chrUn_gl000230", "chrUn_gl000242", "chrUn_gl000243", "chrUn_gl000236", 
                "chrUn_gl000240", "chr17_gl000206_random", "chrUn_gl000232", "chrUn_gl000234", "chrUn_gl000238", "chrUn_gl000244",
                "chrUn_gl000248", "chr8_gl000196_random", "chrUn_gl000249", "chrUn_gl000246", "chr17_gl000203_random", "chr8_gl000197_random", 
                "chrUn_gl000245", "chrUn_gl000247", "chr9_gl000201_random", "chrUn_gl000235", "chrUn_gl000239", "chr21_gl000210_random", "chrUn_gl000231",
                "chrUn_gl000229", "chrUn_gl000226", "chr18_gl000207_random")
seqlevels(txf) = seqlevels(txf)[!seqlevels(txf) %in% notinheader]

# Analyze splice variants in our samples
sgfc = analyzeFeatures(si, features = txf)
save(si, sgfc, file="./Desktop/BAMS/SGSeq_objects_sortedRNA.rda")
sgvc = analyzeVariants(sgfc)
save(si, sgfc, sgvc, file="./Desktop/BAMS/SGSeq_objects_sortedRNA.rda")


## Compare individual splice variant differences by variant type across fraction and age
# Construct the objects needed to test individual splice variant differences
sgv = rowRanges(sgvc)
sgv = getSGVariantCounts(sgv, sample_info = si)
sgv
save(si, sgfc, sgvc, sgv, file="./Desktop/BAMS/SGSeq_objects_sortedRNA.rda")
sgv.counts = counts(sgv)
vid = as.character(variantID(sgv))
eid = as.character(eventID(sgv))

FullModel = ~ sample + exon + Zone:exon + Fetal:exon
ReducedModel = ~ sample + exon + Zone:exon

sampleData = pd[match(colnames(sgv)[1:10], rownames(pd)),]
sampleData = rbind(sampleData, pd[grep("Br5339C1_polyA", rownames(pd)),], pd[grep("Br5340C1_polyA", rownames(pd)),])
rownames(sampleData) = c(rownames(sampleData)[1:10], "Br5339C1_downsamp", "Br5340C1_downsamp")
sampleData$SampleID = rownames(sampleData)

agedxd = DEXSeqDataSet(countData = sgv.counts, featureID = vid, groupID = eid, sampleData = sampleData,
                    design= ageFullModel)
save(agedxd,fracdxd, file="./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/SGSeq_out/DEXSeq_objects.rda")


# Calculate differential splicing
# by Cell Type
dxd = estimateSizeFactors(dxd)
dxd = estimateDispersions(dxd, formula = ageFullModel)
plotDispEsts(agedxd)
agedxd = testForDEU(agedxd, reducedModel = ageReducedModel, fullModel = ageFullModel)
agedxd = estimateExonFoldChanges(agedxd, fitExpToVar="Fetal")
agedxr = DEXSeqResults(agedxd)
save(agedxd,agedxr, file="./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/SGSeq_out/DEXSeq_objects.rda")


# Plot MA

pdf("./Dropbox/sorted_figures/new/github_controlled/intron_retention/figures/SGSeq_out/MA_Plots_Differential_splicing_byFraction_byAge.pdf",height=6,width=6)
plotMA(fracdxr, ylim = c(-15,15), main = "Differential Splicing by Fraction", alpha=0.05)
dev.off()


# Explore results
dexres = dxr[order(dxr$padj),]
mcols = mcols(sgv)
dexres = DataFrame(dexres, vID = as.integer(dexres$featureID))
dexres = cbind(dexres, mcols[match(dexres$vID, mcols$variantID),])
dexres = DataFrame(dexres, more.in.nuc.prenatal = ifelse(dexres[,10]>0, "Yes", "No"))
dexres.sig = dexres[which(dexres$padj<=0.05),]


# Which genes contain splice variants that are differentially included by cell type?

elementNROWS(dexres.sig)
#     Fraction           Age    frac.adult frac.prenatal   age.cytosol   age.nucleus 
#         2158          4608          2512           131          2182          2226 
type = lapply(dexres.sig, function(x) unlist(x$variantType))
type = lapply(type, function(x) data.frame(table(x)[which(names(table(x)) %in% c("SE:S","S2E:S","RI:R","MXE","A5SS:D","A3SS:D","AFE"))]))
type =  Map(cbind, type, comparison = list("By Fraction","By Age","By Fraction\nIn Adult","By Fraction\nIn Prenatal","By Age\nIn Cytosol","By Age\nIn Nucleus"))
type = do.call(rbind, type)
type$comparison = factor(type$comparison, 
                         levels = c("By Age\nIn Nucleus","By Age\nIn Cytosol","By Age","By Fraction\nIn Prenatal","By Fraction\nIn Adult","By Fraction"))


pdf("./Dropbox/sorted_figures/new/github_controlled/intron_retention/figures/SGSeq_out/DSE_counts_byGroup_byFrac_byAge.pdf")
ggplot(type[which(type$comparison!="By Fraction" & type$comparison!="By Age"),], 
       aes(x = comparison, y = Freq, fill = x)) + geom_bar(stat = "identity") +
  coord_flip() +
  labs(fill="") +
  ylab("Count") + xlab("") +
  ggtitle("Differentially Spliced Events\n(FDR < 0.05)") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
ggplot(type[grep("Fraction", type$comparison),], 
       aes(x = comparison, y = Freq, fill = x)) + geom_bar(stat = "identity") +
  coord_flip() +
  labs(fill="") +
  ylab("Count") + xlab("") +
  ggtitle("Differentially Spliced Events\nBy Fraction (FDR < 0.05)") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
ggplot(type[grep("Age", type$comparison),], 
       aes(x = comparison, y = Freq, fill = x)) + geom_bar(stat = "identity") +
  coord_flip() +
  labs(fill="") +
  ylab("Count") + xlab("") +
  ggtitle("Differentially Spliced Events\nBy Age (FDR < 0.05)") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
dev.off()


## Which direction are the log fold changes by variant type?

lapply(dexres,head)
df = Map(cbind, lapply(dexres, function(x) data.frame(x[,c("groupID","featureID","padj")],LFC = x[,10], x[,c("more.in.nuc.prenatal","geneName","variantType","variantName")])),
         comparison = list("By Fraction","By Age","By Fraction\nIn Adult","By Fraction\nIn Prenatal","By Age\nIn Cytosol","By Age\nIn Nucleus"))
df = do.call(rbind, df)
df$geneName = as.character(df$geneName)
df$variantType = as.character(df$variantType)
df$variantName = as.character(df$variantName)
df$threshold = ifelse(df$padj<=0.05, "FDR < 0.05", "FDR > 0.05")
dt = data.table(df)
dt = dt[variantType %in% c("SE:S","S2E:S","RI:R","MXE","A5SS:P","A3SS:P","A5SS:D","A3SS:D","AFE","ALE"),,]
dt = dt[threshold!="NA",,]
dt$variantType = factor(dt$variantType, levels = c("SE:S","S2E:S","RI:R","MXE","A5SS:P","A3SS:P","A5SS:D","A3SS:D","AFE","ALE"))


pdf("./Dropbox/sorted_figures/new/github_controlled/intron_retention/figures/SGSeq_out/volcano_plots_byComparison_byVariantType.pdf",width=30,height=10)
ggplot(dt[grep("Fraction", dt$comparison),,], aes(x=LFC, y=-log10(padj), colour=threshold)) +
  geom_point(alpha=0.4, size=1.75) +
  geom_vline(xintercept=0, linetype="dotted") +
  facet_grid(comparison ~ variantType) +
  xlab("log2 fold change") + ylab("-log10(FDR)") +
  ggtitle("Differential Splicing by Fraction and Variant Type") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 14))
ggplot(dt[grep("Age", dt$comparison),,], aes(x=LFC, y=-log10(padj), colour=threshold)) +
  geom_point(alpha=0.4, size=1.75) + 
  geom_vline(xintercept=0, linetype="dotted") +
  facet_grid(comparison ~ variantType) +
  xlab("log2 fold change") + ylab("-log10(FDR)") +
  ggtitle("Differential Splicing by Age and Variant Type") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 14))
dev.off()


# Calculate the proportion of splice variants that are significant vs not significant

df$rnum = 1:nrow(df) 
df$SE = ifelse(df$rnum %in% grep("SE:S", df$variantType), "SE", "no")
df$S2E = ifelse(df$rnum %in% grep("S2E:S", df$variantType), "S2E", "no")
df$RI = ifelse(df$rnum %in% grep("RI:R", df$variantType), "RI", "no")
df$MXE = ifelse(df$rnum %in% grep("MXE", df$variantType), "MXE", "no")
df$A5SS.P = ifelse(df$rnum %in% grep("A5SS:P", df$variantType), "A5SS:P", "no")
df$A3SS.P = ifelse(df$rnum %in% grep("A3SS:P", df$variantType), "A3SS:P", "no")
df$A5SS.D = ifelse(df$rnum %in% grep("A5SS:D", df$variantType), "A5SS:D", "no")
df$A3SS.D = ifelse(df$rnum %in% grep("A3SS:D", df$variantType), "A3SS:D", "no")
df$AFE = ifelse(df$rnum %in% grep("AFE", df$variantType), "AFE", "no")
df$ALE = ifelse(df$rnum %in% grep("ALE", df$variantType), "ALE", "no")
colnames(df)[12:21]

prop = list(list(),list(),list(),list(),list(),list())
for (i in (1:length(unique(df$comparison)))) {
  for (j in 1:10) {
    prop[[i]][[j]] = data.frame(Anno = c(length(unique(df[which(df$comparison==unique(df$comparison)[i] & df[,colnames(df)[11+j]]!="no" & df$threshold=="FDR < 0.05"),"variantName"])),
                                         length(unique(df[which(df$comparison==unique(df$comparison)[i] & df[,colnames(df)[11+j]]!="no" & df$threshold!="FDR < 0.05"),"variantName"]))),
                                notAnno = c(length(unique(df[which(df$comparison==unique(df$comparison)[i] & df[,colnames(df)[11+j]]=="no" & df$threshold=="FDR < 0.05"),"variantName"])),
                                            length(unique(df[which(df$comparison==unique(df$comparison)[i] & df[,colnames(df)[11+j]]=="no" & df$threshold!="FDR < 0.05"),"variantName"]))),
                                row.names = c("Sig","N.S."))
  }
names(prop[[i]]) = colnames(df)[12:21]
}
names(prop) = c("ByFraction","ByAge","ByFraction:Adult","ByFraction:Prenatal","ByAge:Cytosol","ByAge:Nucleus")
fisher.prop = lapply(prop, function(x) lapply(x, fisher.test))
write.csv(rbind(pvalue = data.frame(lapply(fisher.prop, function(x) unlist(lapply(x, function(y) y$p.value)))),
                data.frame(lapply(fisher.prop, function(x) unlist(lapply(x, function(y) y$estimate))))),quote=F,
          file="./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/SGSeq_out/DSE_byComparison_byVariantType.csv")
names(unlist(lapply(fisher.prop, function(x) unlist(lapply(x, function(y) y$p.value))))
      )[which(unlist(lapply(fisher.prop, function(x) unlist(lapply(x, function(y) y$p.value))))<=0.0008333333)] 
# 34 of 60 tests show a significant relationship between a variant being significantly regulated by fraction and/or age and being a certain type of variant



### is there a relationship between proportion of sig vs non-sig and direction of expression in a variant type?

prop = list(list(),list(),list(),list(),list(),list())
for (i in (1:length(unique(df$comparison)))) {
  for (j in 1:10) {
    prop[[i]][[j]] = data.frame(Pos = c(length(unique(df[which(df$comparison==unique(df$comparison)[i] & df[,colnames(df)[11+j]]!="no" & df$LFC>0 & df$threshold=="FDR < 0.05"),"variantName"])),
                                         length(unique(df[which(df$comparison==unique(df$comparison)[i] & df[,colnames(df)[11+j]]!="no" & df$LFC>0 & df$threshold!="FDR < 0.05"),"variantName"]))),
                                Neg = c(length(unique(df[which(df$comparison==unique(df$comparison)[i] & df[,colnames(df)[11+j]]!="no" & df$LFC<0 & df$threshold=="FDR < 0.05"),"variantName"])),
                                            length(unique(df[which(df$comparison==unique(df$comparison)[i] & df[,colnames(df)[11+j]]!="no" & df$LFC<0 & df$threshold!="FDR < 0.05"),"variantName"]))),
                                row.names = c("Sig","N.S."))
  }
  names(prop[[i]]) = colnames(df)[12:21]
}
names(prop) = c("ByFraction","ByAge","ByFraction:Adult","ByFraction:Prenatal","ByAge:Cytosol","ByAge:Nucleus")
fisher.prop = lapply(prop, function(x) lapply(x, fisher.test))
write.csv(rbind(pvalue = data.frame(lapply(fisher.prop, function(x) unlist(lapply(x, function(y) y$p.value)))),
                data.frame(lapply(fisher.prop, function(x) unlist(lapply(x, function(y) y$estimate))))),quote=F,
          file="./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/SGSeq_out/DSE_byComparison_byVariantType_byLFC_Dir.csv")
names(unlist(lapply(fisher.prop, function(x) unlist(lapply(x, function(y) y$p.value))))
)[which(unlist(lapply(fisher.prop, function(x) unlist(lapply(x, function(y) y$p.value))))<=0.0008333333)]