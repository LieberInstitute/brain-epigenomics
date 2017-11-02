# qrsh -l bluejay,mem_free=150G,h_vmem=150G,h_fsize=100G
# cd /dcl01/lieber/ajaffe/lab/brain-epigenomics/bumphunting
# module load conda_R/3.4.x
# R
library('limma')
library('GenomicRanges')
library('GGally')
library('devtools')

## Load the raw data
system.time( load('BSobj_bsseqSmooth_Neuron_minCov_3.Rdata') )
cpg <- rowRanges(BSobj)

## Load the DMRs
system.time( load('bumps_bsseqSmooth_Neuron_interaction_250_perm.Rdata') )
dmrs <- GRanges(bumps$table[bumps$table$fwer < 0.05, ])

## Load limma results
system.time( load('limma_Neuron_CpGs_minCov_3.Rdata') )
dim(BSobj)
stopifnot(nrow(BSobj) == nrow(fits[[1]]))
print(object.size(fits), units = 'Mb')

## Check the ordering of the coefficients
sapply(fits[2:3], function(f) colnames(f$coef))

## Specify the age coefficients
coefs <- list('age' = 2, 'age_glia' = 2, 'age_neuron' = c(2, 4), 'age_cell_difference' = 4)

age_coef <- mapply(function(f, coef) {
    if(length(coef) > 1) {
        rowSums(f$coefficients[, coef])
    } else {
        f$coefficients[, coef]
    }
}, fits[c(2, 3, 3, 3)], coefs)
colnames(age_coef) <- c('overall', 'age_glia', 'age_neuron', 'age_cell_difference')
head(age_coef)
summary(age_coef)

## Get age_neuron
f <- contrasts.fit(fits[[3]], c(0, 1, 0, 1))

## Check coefs
## they are basically identical
identical(f$coefficients[, 1], age_coef[, 3])
summary(abs(f$coefficients[, 1] - age_coef[, 3]))

## Get the age t-statistics
coefs <- list('age' = 2, 'age_glia' = 2, 'age_cell_difference' = 4)
age_t <-  mapply(function(f, coef) {
    f$t[, coef]
}, fits[c(2, 3, 3)], coefs)
colnames(age_t) <- c('overall', 'age_glia', 'age_cell_difference')

## add the age_neuron t-stat
age_t <- cbind(age_t, 'age_neuron' = f$t[, 1])

## Get an F for age_neuron
top <- topTable(fits[[3]], coef = c(2, 4), n = nrow(age_coef), sort.by = 'none')
head(top)

age_t <- cbind(age_t, 'age_neuron_F' = top$F)
head(age_t)
summary(age_t)

## Get the age p-values
age_p <-  mapply(function(f, coef) {
    f$p.value[, coef]
}, fits[c(2, 3, 3)], coefs)
colnames(age_p) <- c('overall', 'age_glia', 'age_cell_difference')
age_p <- cbind(age_p, 'age_neuron' = f$p.value[, 1])
age_p <- cbind(age_p, 'age_neuron_F' = top$P.Value)
head(age_p)
summary(age_p)

## F p-value for age and interaction == 0 is not the same as
## t-stat p-value for age_neuron
identical(age_p[, 4], age_p[, 5])
summary(abs(age_p[, 4] - age_p[, 5]))

## Save the results
print(object.size(age_coef), units = 'Mb')
print(object.size(age_t), units = 'Mb')
print(object.size(age_p), units = 'Mb')

limma_age <- list('coef' = age_coef, 't' = age_t, 'pvalue' = age_p)

## Relate DMRs to CpGs
ov <- findOverlaps(cpg, dmrs)

## Save main pieces for later
save(limma_age, dmrs, ov, file = 'limma_Neuron_CpGs_minCov_3_ageInfo.Rdata')

## Check that the CpGs for a given DMR are next to each other
stopifnot(identical(nrun(Rle(subjectHits(ov))), length(unique(subjectHits(ov)))))

age_dmr <- function(age, type) {
    new <- split(data.frame(age[queryHits(ov), ]), subjectHits(ov))
    
    mean <- do.call(rbind, lapply(new, colMeans))[unique(subjectHits(ov)), ]
    sum <- do.call(rbind, lapply(new, colSums))[unique(subjectHits(ov)), ]
    
    colnames(mean) <- paste0(colnames(mean), '_', type, '_mean')
    colnames(sum) <- paste0(colnames(sum), '_', type, '_sum')
    cbind(mean, sum)
}


dmrs$coef <- bumps$coef[which(bumps$table$fwer < 0.05), 1]

mcols(dmrs) <- cbind(mcols(dmrs), age_dmr(age_t[, match(colnames(age_coef), colnames(age_t))], 'tstat'), age_dmr(age_coef, 'coef'))


pdf('age_for_interaction_dmrs.pdf', width = 10, height = 10)
ggpairs(as.data.frame(mcols(dmrs)), columns = grep('tstat_mean', colnames(mcols(dmrs))), upper = list(continuous = 'points'))
ggpairs(as.data.frame(mcols(dmrs)), columns = grep('tstat_sum', colnames(mcols(dmrs))), upper = list(continuous = 'points'))
ggpairs(as.data.frame(mcols(dmrs)), columns = grep('coef_mean', colnames(mcols(dmrs))), upper = list(continuous = 'points'))
ggpairs(as.data.frame(mcols(dmrs)), columns = grep('coef_sum', colnames(mcols(dmrs))), upper = list(continuous = 'points'))

tmp <- as.data.frame(mcols(dmrs))
tmp$age_neuron_tstat_mean <- tmp$age_neuron_tstat_mean * sign(tmp$age_neuron_coef_mean)
tmp$age_neuron_tstat_sum <- tmp$age_neuron_tstat_sum * sign(tmp$age_neuron_coef_mean)
ggpairs(tmp, columns = grep('tstat_mean', colnames(mcols(dmrs))), upper = list(continuous = 'points'))
ggpairs(tmp, columns = grep('tstat_sum', colnames(mcols(dmrs))), upper = list(continuous = 'points'))

dev.off()


pdf('age_for_interaction_dmrs_diff.pdf')

boxplot(dmrs$age_glia_coef_mean, dmrs$age_neuron_coef_mean, names = paste0(c('Glia', 'Neuron'), ' (', round(c(mean(dmrs$age_glia_coef_mean > 0), mean(dmrs$age_neuron_coef_mean > 0)) * 100, 0), '% up)'), ylab = 'Age mean coefficient')
legend('bottomright', legend = paste('p-value <', signif(t.test(dmrs$age_glia_coef_mean, dmrs$age_neuron_coef_mean, paired = TRUE)$p.value, 3)), bty = 'n')


boxplot(abs(dmrs$age_glia_coef_mean), abs(dmrs$age_neuron_coef_mean), names = c('Glia', 'Neuron'), ylab = 'Absolute age mean coefficient')
legend('topright', legend = paste('p-value <', signif(t.test(abs(dmrs$age_glia_coef_mean), abs(dmrs$age_neuron_coef_mean), paired = TRUE)$p.value, 3)), bty = 'n')


d <- abs(dmrs$age_glia_coef_mean) - abs(dmrs$age_neuron_coef_mean) > 0
boxplot(abs(dmrs$age_glia_coef_mean) - abs(dmrs$age_neuron_coef_mean), xlab = paste0('Glia (', sum(d), ') - Neuron (', sum(!d), ')'), ylab = 'Absolute age mean coefficient difference')
legend('topright', legend = paste('p-value <', signif(t.test(abs(dmrs$age_glia_coef_mean) - abs(dmrs$age_neuron_coef_mean))$p.value, 3)), bty = 'n')

dev.off()

save(dmrs, file = 'limma_Neuron_CpGs_minCov_3_ageInfo_dmrs.Rdata')

## Reproducibility info
proc.time()
message(Sys.time())
options(width = 120)
session_info()
