library(RnBeads)

# Made bed files for Chromosome 21 only
# < LIBD_Data/DNAm_Ratios_duplicates_dropped/CpGs/BedFiles/SortedSampleNames.txt parallel -P1 grep -w chr21 
# LIBD_Data/DNAm_Ratios_duplicates_dropped/CpGs/BedFiles/{}_concatenated_duplicatesRemoved.bismark.cov ">" 
# LIBD_Data/DNAm_Ratios_duplicates_dropped/CpGs/BedFiles/{}_chr21.bismark.cov


# Directories for data, annotation file, reports

data.source = "/dcl01/lieber/WGBS/LIBD_Data/DNAm_Ratios_duplicates_dropped/CpGs/BedFiles/"
sample.annotation = "/dcl01/lieber/WGBS/LIBD_Data/DNAm_Ratios_duplicates_dropped/CpGs/BedFiles/sample.anno.chr21.csv"
report.dir = "/dcl01/lieber/ajaffe/Amanda/WGBS/RnBeads_analysis/chr21.reports.tryagain"

# Set Parameters for RnBeads

options(bitmapType='cairo')
rnb.options(analysis.name = "CpGs Test Chr21",
            import.bed.style = "bismarkCov",
            region.aggregation = "coverage.weighted",
            qc.coverage.plots = TRUE,
            qc.coverage.violins = TRUE,
            filtering.greedycut = FALSE,
            filtering.sex.chromosomes.removal = TRUE,
            differential.comparison.columns = c("Cell.Type", "Age", "Sex", "Age.Bin"),
            covariate.adjustment.columns = c("Date", "Data.Yield.PF.Bases.Gb", "Reads.PF.Reads.M", "lessThanOrEqualTo_Q30_percent",
            "Total.ng.DNA", "ng.DNA.per.200K.nuclei.input", "X260.over.280", "Nuclei.Per.mg", "Proportion.Neuron",
            "PMI", "Sex", "Race", "RIN.Ribozero", "pH",
            "avg.Cov", "cov.sDev", "Percent.Duplication", "total_num_trimmed_reads", "total_num_untrimmed_reads", "alignment.efficiency", "Prop.Trimmed.Reads"),
            columns.pairing=c("Cell.Type"="Brain.ID"),
            differential.enrichment = TRUE)

# Run Vanilla Analysis

rnb.run.analysis(dir.reports=report.dir, sample.sheet=sample.annotation, data.source=data.source, data.type="bs.bed.dir")
