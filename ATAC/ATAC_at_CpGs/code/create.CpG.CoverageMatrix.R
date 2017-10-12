# load matrix of peaks as rows and columns as samples
samps = scan("/dcl01/lieber/ajaffe/Amanda/ATAC/FinalPeaks/scripts/finalizedList-pooled.txt", what = "character")

covMat = list()
for (i in 1:length(samps)){
  tmp = samps[i]
  covMat[[i]] = read.table(paste0("/dcl01/lieber/ajaffe/Amanda/ATAC/Coverage/coverage_at_CpGs/",
                                  samps[i], "_CpGs_coverageBed.txt"),
                           col.names = c("chr","start","end", "Place1", "Place1", "Place3","ReadCount"))
}
names(covMat) = samps
coord = lapply(covMat, function(x) paste0(x$chr, ":", x$start, "-", x$end))
ReadCount = lapply(covMat, function(x) x$ReadCount)
ATACcovCpGs = as.data.frame(ReadCount, row.names = coord[[1]])
ATACcovCpGs = CpGcov[,order(colnames(CpGcov))]
pd = read.table("/dcl01/lieber/ajaffe/Amanda/ATAC/pd.txt", header=TRUE)
rownames(pd) = pd$ID

pd$Race[pd$Race == "CAUC "] <- 'CAUC'
pd$Race[pd$Race == "AA "] <- "AA"
pd$Sex[pd$Sex == " M"] <- 'M'
pd$RIN <- as.numeric(gsub(" ", "", pd$RIN))
pd$Seq.Category = as.character(pd$Seq.Category)
pd$Date.Sorted = as.character(pd$Date.Sorted)
pd$Brain.ID = as.character(pd$Brain.ID)
pd$Age.Bins = as.character(pd$Age.Bins)
pd$Sex = as.character(pd$Sex)
pd$Race = as.character(pd$Race)

SampleID = unique(pd$SampleID)
ATACpd = data.frame(SampleID = SampleID, RNum = 1:length(SampleID), PE_HS3000 = NA, SE_HS3000 = NA, SE_HS2000 = NA, Date.Sorted = NA, BrainID = NA, CellType = NA, Age = NA,
               Age.Bins = NA, PMI = NA, Sex = NA, Race = NA, RIN = NA, pH = NA, Proportion.Neurons = NA, yield.nuclei.mg.tissue = NA, 
               Post.Picard.Markdups.BAM.Path = NA, Filtered.Sorted.BAM.Path = NA, BigWig.Path = NA, TotalReads.Sequenced = NA,          
               TotalAligned.dupsIncl = NA, Alignment.Percent = NA, TotalReads.dupsRemoved = NA, Percent.Duplicates = NA,            
               chrM.final = NA, TotalReads.final = NA, TotalReads.ff = NA, Frag.Size = NA)
ATACpd$CellType = ifelse(ATACpd$RNum %in% grep("Minus", ATACpd$SampleID), "Glia", "Neuron")

for (i in 1:length(SampleID)){
  ATACpd$PE_HS3000[i] = ifelse("PE_HS3000" %in% pd[which(pd$SampleID==ATACpd$SampleID[i]), "Seq.Category"], "PE_HS3000", "no")
  ATACpd$SE_HS3000[i] = ifelse("SE_HS3000" %in% pd[which(pd$SampleID==ATACpd$SampleID[i]), "Seq.Category"], "SE_HS3000", "no")
  ATACpd$SE_HS2000[i] = ifelse("SE_HS2000" %in% pd[which(pd$SampleID==ATACpd$SampleID[i]), "Seq.Category"], "SE_HS2000", "no")
  ATACpd$Date.Sorted[i] = unique(pd[which(pd$SampleID==ATACpd$SampleID[i]), "Date.Sorted"])
  ATACpd$BrainID[i] = unique(pd[which(pd$SampleID==ATACpd$SampleID[i]), "Brain.ID"])
  ATACpd$Age[i] = unique(pd[which(pd$SampleID==ATACpd$SampleID[i]), "Age"])
  ATACpd$Age.Bins[i] = unique(pd[which(pd$SampleID==ATACpd$SampleID[i]), "Age.Bins"])
  ATACpd$PMI[i] = unique(pd[which(pd$SampleID==ATACpd$SampleID[i]), "PMI"])
  ATACpd$Sex[i] = unique(pd[which(pd$SampleID==ATACpd$SampleID[i]), "Sex"])
  ATACpd$Race[i] = unique(pd[which(pd$SampleID==ATACpd$SampleID[i]), "Race"])
  ATACpd$RIN[i] = unique(pd[which(pd$SampleID==ATACpd$SampleID[i]), "RIN"])
  ATACpd$pH[i] = unique(pd[which(pd$SampleID==ATACpd$SampleID[i]), "pH"])
  ATACpd$Proportion.Neurons[i] = unique(pd[which(pd$SampleID==ATACpd$SampleID[i]), "Proportion.Neurons"])
  ATACpd$yield.nuclei.mg.tissue[i] = unique(pd[which(pd$SampleID==ATACpd$SampleID[i]), "yield.nuclei.mg.tissue"])
}


save(ATACcovCpGs, ATACpd, file="/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/ATAC/CpG.single.site.ATAC.coverageMatrix.rda")
