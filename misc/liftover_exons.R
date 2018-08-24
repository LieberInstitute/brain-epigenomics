####
library(rtracklayer)
library(SummarizedExperiment)

## load exon data
load('/dcl01/lieber/ajaffe/lab/brain-epigenomics/brainseq_pipeline/polyA_unstranded/rse_exon_polyA_dlpfc_n41.Rdata', verbose = TRUE)

## get coordinates
gr = rowRanges(rse_exon)

## liftover
chain = import.chain("hg19ToMm10.over.chain")
lifted = liftOver(gr, chain)

## break back apart
gr$numLiftedPieces = lengths(lifted)
gr$numLiftedPieces_Threshold = gr$numLiftedPieces
gr$numLiftedPieces_Threshold[gr$numLiftedPieces_Threshold > 4] = 5

table(gr$gene_type, gr$numLiftedPieces_Threshold)

## number of multiple chrs
gr$numLiftedChr = lengths(runLength(seqnames(lifted)))
table(gr$gene_type, gr$numLiftedChr > 1)