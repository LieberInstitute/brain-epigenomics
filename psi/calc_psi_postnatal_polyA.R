## Based on https://github.com/LieberInstitute/BrainRNACompartments/blob/master/intron_retention/SGSeq.R
library('SummarizedExperiment')
library('GenomicFeatures')
library('SGSeq')
library('devtools')

CORES <- 10

## Find BAM files
rse_files <- 
    c('../brainseq_pipeline/polyA_unstranded/rse_gene_polyA_dlpfc_n41.Rdata',
    '../brainseq_pipeline/riboZero_stranded/rse_gene_ribozero_dlpfc_n31.Rdata')
rse_all <- sapply(rse_files, function(x) {
    load(x, verbose = TRUE)
    rowRanges(rse_gene)$meanExprs <- 0
    return(rse_gene)
})
rse <- do.call(cbind, rse_all)

## Keep only the postnatal samples
rse <- rse[, colData(rse)$Age > 0 & colData(rse)$Experiment == 'PolyA']

## Get the BAM info for SGSeq
si <- data.frame(sample_name = colnames(rse), file_bam = rse$bamFile,
    stringsAsFactors = FALSE)
si <- getBamInfo(si, cores = 6)

## More setup for SGSeq
gencode <- makeTxDbFromGFF('ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/GRCh37_mapping/gencode.v25lift37.annotation.gtf.gz', organism='Homo sapiens')
txf <- convertToTxFeatures(gencode)

## Compute psi with SGSeq
sgfc <- analyzeFeatures(si, features = txf, cores = CORES)
colData(sgfc) <- cbind(colData(sgfc), colData(rse))
sgvc10 <- analyzeVariants(sgfc, min_denominator = 10, cores = CORES)

## Save results
dir.create('rda', showWarnings = FALSE)
save(txf, file = 'rda/txf_postnatal_polyA.Rdata')
save(sgfc, file = 'rda/sgfc_postnatal_polyA.Rdata')
save(sgvc10, file = 'rda/sgvc10_postnatal_polyA.Rdata')

## Reproducibility info
proc.time()
message(Sys.time())
options(width = 120)
session_info()
