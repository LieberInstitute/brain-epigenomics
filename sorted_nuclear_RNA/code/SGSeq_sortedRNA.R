library("SGSeq")
library("GenomicFeatures")
library("GenomicRanges")
library("ggplot2")
library("DEXSeq")
library(data.table)
library(clusterProfiler)
library('org.Hs.eg.db')

load("./Desktop/BAMS/DMR_objects.rda")


## For SGSeq, create sample information dataframe

names = c("NeuN_Plus_14","NeuN_Minus_14", "NeuN_Plus_13","NeuN_Minus_13","NeuN_Plus_12","NeuN_Minus_12")
file_bam = paste0("/dcl01/lieber/ajaffe/Amanda/sorted_nuclear_RNA/merged/", frag_length$sample_name, "_combinedLibraries.bam")
file_bam = paste0("./Desktop/BAMS/", frag_length$sample_name, "_combinedLibraries.bam")

si = data.frame(sample_name = frag_length$sample_name, file_bam = file_bam)
si$sample_name = as.character(si$sample_name)
si$file_bam = as.character(si$file_bam)
si = getBamInfo(si)


## Extract transcript features

#gencode = loadDb("/dcl01/lieber/ajaffe/Amanda/sorted_nuclear_RNA/gencode.v25lift37.annotation.sqlite")
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


## Analyze splice variants in our samples

sgfc = analyzeFeatures(si, features = txf)
sgvc = analyzeVariants(sgfc)
sgv = rowRanges(sgvc)
sgv = getSGVariantCounts(sgv, sample_info = si)
save(si, sgfc, sgvc, sgv, file="./Desktop/BAMS/SGSeq_objects_sortedRNA.rda")


### Compare individual splice variant differences by variant type across fraction and age

## Construct the objects needed to test individual splice variant differences

load("./Desktop/BAMS/SGSeq_objects_sortedRNA.rda")
sgv.counts = counts(sgv)
vid = as.character(variantID(sgv))
eid = as.character(eventID(sgv))

FullModel = ~ sample + exon + ageCat:exon + cellType:exon
ReducedModel = ~ sample + exon + ageCat:exon

pd = data.frame(sample = c("NeuN_Plus_12","NeuN_Minus_12","NeuN_Plus_14","NeuN_Minus_14","NeuN_Plus_13","NeuN_Minus_13"),
                cellType = rep.int(c("Neuron","Glia"),3),
                ageNum = c(0.36,0.36,2.71,2.71,14.01,14.01),
                ageCat = c("Infant","Infant","Toddler","Toddler","Teen","Teen"))
sampleData = pd[match(colnames(sgv), pd$sample),]
rownames(sampleData) = sampleData$sample

dxd = DEXSeqDataSet(countData = sgv.counts, featureID = vid, groupID = eid, sampleData = sampleData,
                    design= FullModel)


# Calculate differential splicing by Cell Type
dxd = estimateSizeFactors(dxd)
dxd = estimateDispersions(dxd, formula = FullModel)
plotDispEsts(dxd)
dxd = testForDEU(dxd, reducedModel = ReducedModel, fullModel = FullModel)
dxd = estimateExonFoldChanges(dxd, fitExpToVar="cellType")
dxr = DEXSeqResults(dxd)
save(dxd,dxr, file="./Desktop/BAMS/DEXSeq_objects_sortedRNA.rda")


# Plot MA
pdf("./Desktop/BAMS/MA_Plots_Differential_splicing_byCellType.pdf",height=6,width=6)
plotMA(dxr, ylim = c(-15,15), main = "Differential Splicing by Cell Type\nIn Sorted Nuclear RNA", alpha=0.05)
dev.off()


# Explore results
dexres = dxr[order(dxr$padj),]
mcols = mcols(sgv)
dexres = DataFrame(dexres, vID = as.integer(dexres$featureID))
dexres = cbind(dexres, mcols[match(dexres$vID, mcols$variantID),])
dexres = DataFrame(dexres, more.in.Neuronal = ifelse(dexres$log2fold_Neuron_Glia > 0, "Yes", "No"))
dexres.sig = dexres[which(dexres$padj<=0.05),]
dim(dexres.sig)

# What are the types and numbers of significantly differently regulated splice variants?

table(dexres.sig$more.in.Neuronal)
#No Yes 
#37  33

spl = split(data.frame(dexres.sig), elementNROWS(dexres.sig$variantType)<2)
f = mapply(function(v,m,n) data.frame(variantType = v, more.in.Neuronal = m, variantName = n), 
       dexres.sig$variantType[which(dexres.sig$variantName %in% spl$'FALSE'$variantName)],
       as.list(dexres.sig$more.in.Neuronal[which(dexres.sig$variantName %in% spl$'FALSE'$variantName)]),
       as.list(dexres.sig$variantName[which(dexres.sig$variantName %in% spl$'FALSE'$variantName)]))
fa = list()
for (i in 1:nrow(spl$'FALSE')) { fa[[i]] = f[,i] }
f = do.call(rbind, lapply(fa, function(x) data.frame(variantType = x$variantType, more.in.Neuronal = x$more.in.Neuronal, variantName = x$variantName)))
spl$'TRUE'[grep("OTHER", spl$'TRUE'$variantName),"variantType"] = list("OTHER")
x = data.table(rbind(data.frame(variantType = unlist(spl$'TRUE'$variantType), more.in.Neuronal = spl$'TRUE'$more.in.Neuronal,
                                variantName = spl$'TRUE'$variantName), f))
x = x[,length(unique(variantName)), by = c("variantType","more.in.Neuronal")]
x$more.in.Neuronal = ifelse(x$more.in.Neuronal == "Yes", "Greater in Neurons", "Greater in Glia")

pdf("./Desktop/BAMS/DSEcounts_byVariantType_byCellTypeDirection.pdf", width = 13, height = 5)
ggplot(x, aes(x = variantType, y = V1, fill = more.in.Neuronal)) + geom_bar(stat = "identity") +
  labs(fill="") +
  ylab("Count") + xlab("") +
  ggtitle("Differentially Spliced Events\n(FDR < 0.05)") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
dev.off()


## Which genes contain splice variants that are differentially included by cell type?

# get gene names

spl = split(data.frame(dexres.sig), elementNROWS(dexres.sig$geneName)<2)
f = mapply(function(v,m,n) data.frame(geneName = v, more.in.Neuronal = m, variantName = n), 
           dexres.sig$geneName[which(dexres.sig$variantName %in% spl$'FALSE'$variantName)],
           as.list(dexres.sig$more.in.Neuronal[which(dexres.sig$variantName %in% spl$'FALSE'$variantName)]),
           as.list(dexres.sig$variantName[which(dexres.sig$variantName %in% spl$'FALSE'$variantName)]))
fa = list()
for (i in 1:nrow(spl$'FALSE')) { fa[[i]] = f[,i] }
f = do.call(rbind, lapply(fa, function(x) data.frame(geneName = x$geneName, more.in.Neuronal = x$more.in.Neuronal, variantName = x$variantName)))
spl$'TRUE'[grep("OTHER", spl$'TRUE'$variantName),"variantType"] = list("OTHER")
x = data.table(rbind(data.frame(geneName = unlist(spl$'TRUE'$geneName), more.in.Neuronal = spl$'TRUE'$more.in.Neuronal,
                                variantName = spl$'TRUE'$variantName), f))
x$more.in.Neuronal = ifelse(x$more.in.Neuronal == "Yes", "Greater in Neurons", "Greater in Glia")

sig = split(x$geneName, x$more.in.Neuronal)
sig = lapply(sig, as.character)
entrezID = lapply(sig, function(x) as.character(na.omit(geneMap[which(geneMap$gencodeID %in% x), "EntrezID"])))

# Define universe as all genes expressed in each of the four groups

GeneUniverse = as.character(unique(geneMap[which(geneMap$gencodeID %in% unlist(mcols$geneName)),"EntrezID"]))
GeneUniverse = na.omit(GeneUniverse)

# Find enriched Pathways via KEGG
elementNROWS(entrezID)
keggList = c(lapply(entrezID, function(x) enrichKEGG(as.character(x), organism="human", universe= GeneUniverse, 
                                                   minGSSize=5, pAdjustMethod="BH", qvalueCutoff=1)),
             combined = enrichKEGG(unlist(entrezID), organism="human", universe= GeneUniverse, 
                                   minGSSize=5, pAdjustMethod="BH", qvalueCutoff=1))

# Enriched Molecular Function GOs
goList_MF = c(lapply(entrezID, function(x) enrichGO(as.character(x), ont = "MF", OrgDb = org.Hs.eg.db, 
                                                  universe= GeneUniverse, minGSSize=5, pAdjustMethod="BH",
                                                  qvalueCutoff=1)),
              combined = enrichGO(unlist(entrezID), ont = "MF", OrgDb = org.Hs.eg.db, 
                                  universe= GeneUniverse, minGSSize=5, pAdjustMethod="BH",
                                  qvalueCutoff=1))

# Biological Process GO enrichment
goList_BP = c(lapply(entrezID, function(x) enrichGO(as.character(x), ont = "BP", OrgDb = org.Hs.eg.db, 
                                                  universe= GeneUniverse, minGSSize=5, pAdjustMethod="BH",
                                                  qvalueCutoff=1)),
              combined = enrichGO(unlist(entrezID), ont = "BP", OrgDb = org.Hs.eg.db, 
                                  universe= GeneUniverse, minGSSize=5, pAdjustMethod="BH",
                                  qvalueCutoff=1))

# Cellular Compartment GO enrichment
goList_CC = c(lapply(entrezID, function(x) enrichGO(as.character(x), ont = "CC", OrgDb = org.Hs.eg.db, 
                                                  universe= GeneUniverse, minGSSize=5, pAdjustMethod="BH",
                                                  qvalueCutoff=1)),
              combined = enrichGO(unlist(entrezID), ont = "CC", OrgDb = org.Hs.eg.db, 
                                  universe= GeneUniverse, minGSSize=5, pAdjustMethod="BH",
                                  qvalueCutoff=1))

# Disease Ontology
goList_DO = c(lapply(entrezID, function(x) enrichDO(as.character(x), ont = "DO", universe= GeneUniverse, 
                                                  minGSSize=5, pAdjustMethod="BH", qvalueCutoff=1, readable=TRUE)),
              combined = enrichDO(unlist(entrezID), ont = "DO", universe= GeneUniverse, 
                                  minGSSize=5, pAdjustMethod="BH", qvalueCutoff=1, readable=TRUE))


# Compare the enriched terms between 7 groups
# KEGG
compareKegg = compareCluster(entrezID, fun="enrichKEGG", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
# Biological Process
compareBP = compareCluster(entrezID, fun="enrichGO", ont = "BP", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
# Molecular Function
compareMF = compareCluster(entrezID, fun="enrichGO",  ont = "MF", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
# Cellular Component
compareCC = compareCluster(entrezID, fun="enrichGO",  ont = "CC", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
# Disease Ontology
compareDO = compareCluster(entrezID, fun="enrichDO",  ont = "DO", qvalueCutoff = 0.05, pvalueCutoff = 0.05)

# Save
save(goList_MF, goList_BP, goList_CC, compareBP, compareMF, compareCC, 
     file="./Desktop/BAMS/GO.objects.altSpliced_byCellType.rda")


## plot
pdf("./Desktop/BAMS/BP_altSpliced_byCellType.pdf", width=10,height=12)
plot(compareBP,colorBy="p.adjust",  showCategory = 400, title= "Biological Process GO Enrichment")
dev.off()
pdf("./Desktop/BAMS/CC_altSpliced_byCellType.pdf", width=7,height=5)
plot(compareCC,colorBy="p.adjust",  showCategory = 400, title= "Cellular Compartment GO Enrichment")
dev.off()






