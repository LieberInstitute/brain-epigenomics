library('bsseq')
library('devtools')
library('bumphunter')
library('GenomicRanges')

## Load the raw data
load('/dcl01/lieber/ajaffe/lab/brain-epigenomics/bumphunting/BSobj_bsseqSmooth_Neuron_minCov_3.Rdata', verbose = TRUE)

models <- c('interaction', 'cell', 'age')

xx <- sapply(models, function(model) {
    # model <- 'interaction'
    
    ## Load the bumphunting results
    load(paste0('/dcl01/lieber/ajaffe/lab/brain-epigenomics/bumphunting/', 
        'bumps_bsseqSmooth_Neuron_', model, '_250_perm.Rdata'), verbose = TRUE)    
    
    ## Build a GRanges object from the table output of bumphunter
    bumps_gr <- GRanges(bumps$table[bumps$table$fwer < 0.05, ])

    ## Next find which bases are part of the DMRs
    idx <- countOverlaps(rowRanges(BSobj), bumps_gr) > 0
    
    ## Subset and then save the subsetted BSobj
    DMR_CpG <- BSobj[idx, ]
    print('Number of DMRs')
    print(length(bumps_gr))
    
    print('Number of CpGs in the DMRs')
    print(nrow(DMR_CpG))
    
    dir.create('rda', showWarnings = FALSE)
    save(DMR_CpG, file = paste0('rda/DMR_CpG_', model, '_CpGbased.Rdata'))
    return(model)
})

## Reproducibility info
proc.time()
message(Sys.time())
options(width = 120)
session_info()
