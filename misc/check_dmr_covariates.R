library('bsseq')
library('GenomicRanges')
library('bumphunter')
library('sessioninfo')
library('SummarizedExperiment')

load('/dcl01/lieber/ajaffe/lab/brain-epigenomics/bumphunting/BSobj_bsseqSmooth_Neuron_minCov_3.Rdata', verbose = TRUE)

pd <- pData(BSobj)

## Add lambda info
## adapted code from 
## https://github.com/LieberInstitute/brain-epigenomics/blob/master/lambda/check_lambda.R
lambda <- read.csv('/dcl01/lieber/WGBS/LIBD_Data/lambda_alignment/MasterList_lambda.csv')
m_lambda <- match(pd$Data.ID, lambda$WGC.ID)
stopifnot(!any(is.na(m_lambda)))
pd$avg_conv_efficiency <- colData(BSobj)$avg_conv_efficiency <- lambda$avg_conv_efficiency[m_lambda]
pd$ratio_trimmed <- colData(BSobj)$ratio_trimmed <- pd$total_num_trimmed_reads/pd$total_num_untrimmed_reads

models <- c('interaction', 'cell', 'age')
# model <- models[2]
rse_mean_meth <- lapply(models, function(model) {
    f <- paste0('/dcl01/lieber/ajaffe/lab/brain-epigenomics/bumphunting/bumps_bsseqSmooth_Neuron_', model, '_250_perm.Rdata')
    message(paste(Sys.time(), 'loading', f))
    load(f, verbose = TRUE)

    ## Build a GRanges object from the table output of bumphunter
    bumps_gr <- GRanges(bumps$table[bumps$table$fwer < 0.05, ])

    ## Next find which bases are part of the DMRs
    idx <- countOverlaps(rowRanges(BSobj), bumps_gr) > 0

    ## Get the methylation only for the bases that are in the DMRs
    ## Use the smoothed values, since that's what we used for
    ## finding the DMRs as in
    ## https://github.com/LieberInstitute/brain-epigenomics/blob/master/bumphunting/find_bumps_bsseqSmooth.R#L124
    message(paste(Sys.time(), 'extracting the methylation values'))
    meth <- getMeth(BSobj[idx, ], type = 'smooth')


    message(paste(Sys.time(), 'computing means across the CpGs in the DMRs'))
    ov <- findOverlaps(bumps_gr, rowRanges(BSobj[idx, ]))
    ov_dmr <- split(subjectHits(ov), queryHits(ov))

    ## Finally compute the mean methylation for the bases in the DMR for the glia samples
    dmr_mean_meth <- do.call(rbind, lapply(ov_dmr, function(dmr_bases) {
        colMeans(meth[dmr_bases, , drop = FALSE])
    }))

    stopifnot(length(bumps_gr) == nrow(dmr_mean_meth))
    rownames(dmr_mean_meth) <- names(bumps_gr)
    final <- SummarizedExperiment(assays = list('mean_meth' = dmr_mean_meth), rowRanges = bumps_gr, colData = pData(BSobj))
    return(final)
})
names(rse_mean_meth) <- models
save(rse_mean_meth, file = 'rdas/rse_mean_meth.Rdata')



covs <- c('Race', 'Sex', 'PMI', 'pH', 'avg.Cov', 'Percent.Duplication', 'alignment.efficiency', 'total_num_trimmed_reads', 'total_num_untrimmed_reads', 'ratio_trimmed', 'avg_conv_efficiency')

lapply(covs, function(cov) {
    # cov <- covs[1]
    lapply(models, function(model) {
        coef <- 2
        if(model == 'cell') {
            design <- model.matrix(~ ., data = pd[, c('Cell.Type', 'Age', cov)])
        } else if (model == 'age') {
            design <- model.matrix(~ ., data = pd[, c('Age', 'Cell.Type', cov)])

        } else if (model == 'interaction') {
            design <- model.matrix(~ ~ Age * Cell.Type + ., data = pd[, c('Age', 'Cell.Type', cov)])
            coef <- grep(':', colnames(design))
        }
        
        message(paste(
            Sys.time(), 'extracting the coefficient for',
            colnames(design)[coef], 'after adjusting for',
            paste(colnames(design)[-c(1, coef)], collapse = ', ')
        ))
        
        rawBeta <- bumphunter:::.getEstimate(mat = assays(rse_mean_meth[[model]])$mean_meth, design = design, coef = coef, full = FALSE)
        
        ## If cov = NULL, with model = 'cell'
        ## then
        # > table(rowRanges(rse_mean_meth[[model]])$value - rawBeta[, 1] < 1e-06)
        #
        #  TRUE
        # 11179
        ## this is based on https://github.com/rafalab/bumphunter/blob/master/R/bumphunter.R#L104-L105
        ## which is the piece of code that gets run when we compute the DMRs with
        ## https://github.com/LieberInstitute/brain-epigenomics/blob/master/bumphunting/find_bumps_bsseqSmooth.R#L137-L141
        
        table(sign(rowRanges(rse_mean_meth[[model]])$value), sign(rawBeta[, 1]))
        
        plot(x = rowRanges(rse_mean_meth[[model]])$value, y = rawBeta[, 1])
        
        plot(x = rank(rowRanges(rse_mean_meth[[model]])$value), y = rank(rawBeta[, 1]))
        
        r1 <- rank(rowRanges(rse_mean_meth[[model]])$value)
        r2 <- rank(rawBeta[, 1])
        hmm <- sapply(seq_len(length(r1)), function(i) { sum(r2[seq_len(i)] <= r1[i] )})
        
        plot(x = seq_len(length(r1)), y = hmm)
    })
})



## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
