# qrsh -l bluejay,mem_free=150G,h_vmem=150G,h_fsize=100G
# cd /dcl01/lieber/ajaffe/lab/brain-epigenomics/bumphunting
# module load conda_R/3.4.x
# R
library('limma')
library('GenomicRanges')
library('GGally')
library('bsseq')
library('devtools')
library('RColorBrewer')

## Load the raw data
system.time( load('BSobj_bsseqSmooth_Neuron_minCov_3.Rdata') )
cpg <- rowRanges(BSobj)

## Load the DMRs, keep those with FWER < 0.05
system.time( load('bumps_bsseqSmooth_Neuron_interaction_250_perm.Rdata') )
dmrs <- GRanges(bumps$table[bumps$table$fwer < 0.05, ])
print('Number of interaction DMRs')
length(dmrs)

## Load limma results at the CpG level
system.time( load('limma_Neuron_CpGs_minCov_3.Rdata') )
dim(BSobj)
stopifnot(nrow(BSobj) == nrow(fits[[1]]))
print(object.size(fits), units = 'Mb')

## Check the ordering of the coefficients
sapply(fits[2:3], function(f) colnames(f$coef))

## Specify the age coefficients
coefs <- list('age' = 2, 'age_glia' = 2, 'age_neuron' = c(2, 4), 'age_cell_difference' = 4)

## Extract the age coefficients from the different models
age_coef <- mapply(function(f, coef) {
    if(length(coef) > 1) {
        rowSums(f$coefficients[, coef])
    } else {
        f$coefficients[, coef]
    }
}, fits[c(2, 3, 3, 3)], coefs)
colnames(age_coef) <- c('overall', 'age_glia', 'age_neuron', 'age_cell_difference')
print('Summary of the age coefficients at the CpG level. Overall: age while adjusting for cell type. Age glia: age for glia cells. Age neuron: same but for neurons. Age cell difference: the age and cell type interaction.')
head(age_coef)
summary(age_coef)

## Compute t-stats from the Neuron by fitting the a model
## with age adjusted by cell type, with Neurons as the reference group
## instead of Glia (which is what was done before)
pd <- pData(BSobj)
pd$Cell.Type <- relevel(factor(pd$Cell.Type), 'Neuron')
levels(pd$Cell.Type)
meth <- getMeth(BSobj, type = 'raw')
print('Number of CpGs and samples being considered')
dim(meth)

system.time( f_neuron <- eBayes(lmFit(meth, with(pd, model.matrix(~ Age * Cell.Type)))) )

## Check coefficients
## They are the basically the same ones
print('Checking the coefficients from this model vs the original one. Basically, they are the same.')
identical(f_neuron$coefficients[, 'Age'], age_coef[, 'age_neuron'])
summary(abs(f_neuron$coefficients[, 'Age'] - age_coef[, 'age_neuron']))

## Get the age t-statistics
coefs <- list('age' = 2, 'age_glia' = 2, 'age_cell_difference' = 4)
age_t <-  mapply(function(f, coef) {
    f$t[, coef]
}, fits[c(2, 3, 3)], coefs)
colnames(age_t) <- c('overall', 'age_glia', 'age_cell_difference')
age_t <- cbind(age_t, 'age_neuron' = f_neuron$t[, 'Age'])

print('Head and summary of the t-statistics for each of the coefficients. age_neuron_F is an F-statistic for change in age for neuron cells extract from the interaction model.')
head(age_t)
summary(age_t)

## The interaction looks ok though
print('t-statistics basically match for the interaction term in both models (as expected)')
summary(abs(abs(f_neuron$t[, 'Age:Cell.TypeGlia']) - abs(age_t[, 'age_cell_difference'])))


## Get the age p-values
age_p <-  mapply(function(f, coef) {
    f$p.value[, coef]
}, fits[c(2, 3, 3)], coefs)
colnames(age_p) <- c('overall', 'age_glia', 'age_cell_difference')
age_p <- cbind(age_p, 'age_neuron' = f_neuron$p.value[, 'Age'])

print('Head and summary of the p-values for each of the coefficients.')
head(age_p)
summary(age_p)

## Save the results
print('Object sizes for age_coef, age_t and age_p')
print(object.size(age_coef), units = 'Mb')
print(object.size(age_t), units = 'Mb')
print(object.size(age_p), units = 'Mb')

## Combine the coef, t stats and p values into a single list
limma_age <- list('coef' = age_coef, 't' = age_t, 'pvalue' = age_p)

## Relate DMRs to CpGs
ov <- findOverlaps(cpg, dmrs)
print('Number of CpGs that are inside the interaction DMRs')
length(ov)
length(unique(queryHits(ov)))


## Save main pieces for later
save(limma_age, dmrs, ov, file = 'limma_Neuron_CpGs_minCov_3_ageInfo.Rdata')

## Check that the CpGs for a given DMR are next to each other
stopifnot(identical(nrun(Rle(subjectHits(ov))), length(unique(subjectHits(ov)))))


## Compute mean, sum and area of the age values at the CpG level for each DMR
age_dmr <- function(age, type) {
    ## Group CpGs by the DMR they overlap
    new <- split(data.frame(age[queryHits(ov), ]), subjectHits(ov))
    
    ## Compute the mean and sum of the age information
    mean <- do.call(rbind, lapply(new, colMeans))[unique(subjectHits(ov)), ]
    sum <- do.call(rbind, lapply(new, colSums))[unique(subjectHits(ov)), ]
    area <- do.call(rbind, lapply(new, function(x) { colSums(abs(x)) }))[unique(subjectHits(ov)), ]
    
    ## Return the relevant values
    colnames(mean) <- paste0(colnames(mean), '_', type, '_mean')
    colnames(sum) <- paste0(colnames(sum), '_', type, '_sum')
    colnames(area) <- paste0(colnames(area), '_', type, '_area')
    cbind(mean, sum, area)
}

## Extract the coefficient for the DMRs from the bumphunter output
## only for DMRs with FWER < 0.05
dmrs$coef <- bumps$coef[which(bumps$table$fwer < 0.05), 1]

## Add the mean and sum t-stat and coefficients to the dmr object
mcols(dmrs) <- cbind(mcols(dmrs), age_dmr(limma_age$t[, match(colnames(limma_age$coef), colnames(limma_age$t))], 'tstat'), age_dmr(limma_age$coef, 'coef'))
print('head of the dmr info with CpG results appended at the end')
head(dmrs)

print('Compare a bit the signs of some interaction coefficients')
addmargins(table('DMR coef sign' = sign(dmrs$coef), 'interaction mean coef sign' = sign(dmrs$age_cell_difference_coef_mean)))
addmargins(table('DMR value sign' = sign(dmrs$value), 'interaction mean coef sign' = sign(dmrs$age_cell_difference_coef_mean)))

## Group DMRs by glia and neuron mean coefficients
kdf <- data.frame('glia' = dmrs$age_glia_coef_mean,
    'neuron' = dmrs$age_neuron_coef_mean)
k2 <- kmeans(kdf, 2)
k4 <- kmeans(kdf, 4)
k6 <- kmeans(kdf, 6)
k8 <- kmeans(kdf, 8)

kadf <- data.frame('glia' = abs(dmrs$age_glia_coef_mean),
    'neuron' = abs(dmrs$age_neuron_coef_mean))
ka2 <- kmeans(kadf, 2)
ka4 <- kmeans(kadf, 4)
ka6 <- kmeans(kadf, 6)
ka8 <- kmeans(kadf, 8)



grey_lines <- function() {
    abline(h = 0, col = 'grey80')
    abline(v = 0, col = 'grey80')
    abline(a = 0, b = 1, col = 'grey80')
    abline(a = 0, b = -1, col = 'grey80')
}
cols <- brewer.pal(8, 'Dark2')


## Compare glia vs neuron mean coefficient with K-means clusters
pdf('pdf/glia_vs_neuron_mean_coef.pdf')

## Basic plots without clusters
plot(x = dmrs$age_glia_coef_mean, y = dmrs$age_neuron_coef_mean, xlab = 'Glia mean coefficient', ylab = 'Neuron mean coefficient', pch = 20, cex = 0.5)
grey_lines()

plot(x = abs(dmrs$age_glia_coef_mean), y = abs(dmrs$age_neuron_coef_mean), xlab = 'Glia abs mean coefficient', ylab = 'Neuron abs mean coefficient', pch = 20, cex = 0.5)
grey_lines()


## With clusters
plot(x = dmrs$age_glia_coef_mean, y = dmrs$age_neuron_coef_mean, col = cols[k2$cluster], xlab = 'Glia mean coefficient', ylab = 'Neuron mean coefficient', pch = 20, cex = 0.5, main = 'K2')
points(k2$centers, col = 'black', pch = 8, cex = 2)
grey_lines()

plot(x = dmrs$age_glia_coef_mean, y = dmrs$age_neuron_coef_mean, col = cols[k4$cluster], xlab = 'Glia mean coefficient', ylab = 'Neuron mean coefficient', pch = 20, cex = 0.5, main = 'K4')
points(k4$centers, col = 'black', pch = 8, cex = 2)
grey_lines()

plot(x = dmrs$age_glia_coef_mean, y = dmrs$age_neuron_coef_mean, col = cols[k6$cluster], xlab = 'Glia mean coefficient', ylab = 'Neuron mean coefficient', pch = 20, cex = 0.5, main = 'K6')
points(k6$centers, col = 'black', pch = 8, cex = 2)
grey_lines()

plot(x = dmrs$age_glia_coef_mean, y = dmrs$age_neuron_coef_mean, col = cols[k8$cluster], xlab = 'Glia mean coefficient', ylab = 'Neuron mean coefficient', pch = 20, cex = 0.5, main = 'K8')
points(k8$centers, col = 'black', pch = 8, cex = 2)
grey_lines()

## Abs with coefs from non-abs
plot(x = abs(dmrs$age_glia_coef_mean), y = abs(dmrs$age_neuron_coef_mean), xlab = 'Glia abs mean coefficient', ylab = 'Neuron abs mean coefficient', pch = 20, cex = 0.5, col = cols[k2$cluster], main = 'K2 (from non-abs)')
points(abs(k2$centers), col = 'black', pch = 8, cex = 2)
grey_lines()

plot(x = abs(dmrs$age_glia_coef_mean), y = abs(dmrs$age_neuron_coef_mean), xlab = 'Glia abs mean coefficient', ylab = 'Neuron abs mean coefficient', pch = 20, cex = 0.5, col = cols[k4$cluster], main = 'K4 (from non-abs)')
points(abs(k4$centers), col = 'black', pch = 8, cex = 2)
grey_lines()

plot(x = abs(dmrs$age_glia_coef_mean), y = abs(dmrs$age_neuron_coef_mean), xlab = 'Glia abs mean coefficient', ylab = 'Neuron abs mean coefficient', pch = 20, cex = 0.5, col = cols[k6$cluster], main = 'K6 (from non-abs)')
points(abs(k6$centers), col = 'black', pch = 8, cex = 2)
grey_lines()

plot(x = abs(dmrs$age_glia_coef_mean), y = abs(dmrs$age_neuron_coef_mean), xlab = 'Glia abs mean coefficient', ylab = 'Neuron abs mean coefficient', pch = 20, cex = 0.5, col = cols[k8$cluster], main = 'K8 (from non-abs)')
points(abs(k8$centers), col = 'black', pch = 8, cex = 2)
grey_lines()

## Abs with K from abs
plot(x = abs(dmrs$age_glia_coef_mean), y = abs(dmrs$age_neuron_coef_mean), xlab = 'Glia abs mean coefficient', ylab = 'Neuron abs mean coefficient', pch = 20, cex = 0.5, col = cols[ka2$cluster], main = 'K2')
points(abs(ka2$centers), col = 'black', pch = 8, cex = 2)
grey_lines()

plot(x = abs(dmrs$age_glia_coef_mean), y = abs(dmrs$age_neuron_coef_mean), xlab = 'Glia abs mean coefficient', ylab = 'Neuron abs mean coefficient', pch = 20, cex = 0.5, col = cols[ka4$cluster], main = 'K4')
points(abs(ka4$centers), col = 'black', pch = 8, cex = 2)
grey_lines()

plot(x = abs(dmrs$age_glia_coef_mean), y = abs(dmrs$age_neuron_coef_mean), xlab = 'Glia abs mean coefficient', ylab = 'Neuron abs mean coefficient', pch = 20, cex = 0.5, col = cols[ka6$cluster], main = 'K6')
points(abs(ka6$centers), col = 'black', pch = 8, cex = 2)
grey_lines()

plot(x = abs(dmrs$age_glia_coef_mean), y = abs(dmrs$age_neuron_coef_mean), xlab = 'Glia abs mean coefficient', ylab = 'Neuron abs mean coefficient', pch = 20, cex = 0.5, col = cols[ka8$cluster], main = 'K8')
points(abs(ka8$centers), col = 'black', pch = 8, cex = 2)
grey_lines()


## Naive but with clusters from the abs scale
plot(x = dmrs$age_glia_coef_mean, y = dmrs$age_neuron_coef_mean, col = cols[ka2$cluster], xlab = 'Glia mean coefficient', ylab = 'Neuron mean coefficient', pch = 20, cex = 0.5, main = 'K2 (from abs)')
points(ka2$centers, col = 'black', pch = 8, cex = 2)
grey_lines()

plot(x = dmrs$age_glia_coef_mean, y = dmrs$age_neuron_coef_mean, col = cols[ka4$cluster], xlab = 'Glia mean coefficient', ylab = 'Neuron mean coefficient', pch = 20, cex = 0.5, main = 'K4 (from abs)')
points(ka4$centers, col = 'black', pch = 8, cex = 2)
grey_lines()

plot(x = dmrs$age_glia_coef_mean, y = dmrs$age_neuron_coef_mean, col = cols[ka6$cluster], xlab = 'Glia mean coefficient', ylab = 'Neuron mean coefficient', pch = 20, cex = 0.5, main = 'K6 (from abs)')
points(ka6$centers, col = 'black', pch = 8, cex = 2)
grey_lines()

plot(x = dmrs$age_glia_coef_mean, y = dmrs$age_neuron_coef_mean, col = cols[ka8$cluster], xlab = 'Glia mean coefficient', ylab = 'Neuron mean coefficient', pch = 20, cex = 0.5, main = 'K8 (from abs)')
points(ka8$centers, col = 'black', pch = 8, cex = 2)
grey_lines()

dev.off()


## Make a quick plot comparing many of the variables
dir.create('pdf', showWarnings = FALSE)
pdf('pdf/age_for_interaction_dmrs.pdf', width = 10, height = 10)
ggpairs(as.data.frame(mcols(dmrs)), columns = grep('tstat_mean',
    colnames(mcols(dmrs))), upper = list(continuous = 'points'))
ggpairs(as.data.frame(mcols(dmrs)), columns = grep('tstat_sum',
    colnames(mcols(dmrs))), upper = list(continuous = 'points'))
ggpairs(as.data.frame(mcols(dmrs)), columns = grep('tstat_area',
    colnames(mcols(dmrs))), upper = list(continuous = 'points'))
ggpairs(as.data.frame(mcols(dmrs)),
    columns = c(which(colnames(mcols(dmrs)) == 'value'), grep('coef_mean',
    colnames(mcols(dmrs)))), upper = list(continuous = 'points'))
ggpairs(as.data.frame(mcols(dmrs)),
    columns = c(which(colnames(mcols(dmrs)) == 'value'), grep('coef_sum',
    colnames(mcols(dmrs)))), upper = list(continuous = 'points'))
ggpairs(as.data.frame(mcols(dmrs)),
    columns = c(which(colnames(mcols(dmrs)) == 'value'), grep('coef_area',
    colnames(mcols(dmrs)))), upper = list(continuous = 'points'))
dev.off()


if(FALSE) {
    ## Some plots I was playing around with
    plot(x = abs(dmrs$age_glia_coef_mean), y = abs(dmrs$age_neuron_coef_mean))
    abline(a = 0, b = 1, col = 'red')

    plot(dmrs$age_glia_coef_mean / dmrs$age_neuron_coef_mean)
    plot(dmrs$age_glia_coef_mean, dmrs$age_glia_coef_mean / dmrs$age_cell_difference_coef_mean)

    plot(dmrs$age_glia_coef_mean / dmrs$age_neuron_coef_mean, abs(dmrs$age_glia_coef_mean) / abs(dmrs$age_neuron_coef_mean))
    plot(abs(dmrs$age_glia_coef_mean) / abs(dmrs$age_neuron_coef_mean))
}


pdf('pdf/age_for_interaction_dmrs_diff.pdf')

boxplot(dmrs$age_glia_coef_mean, dmrs$age_neuron_coef_mean, names = paste0(c('Glia', 'Neuron'), ' (', round(c(mean(dmrs$age_glia_coef_mean > 0), mean(dmrs$age_neuron_coef_mean > 0)) * 100, 0), '% up)'), ylab = 'Age mean coefficient')
legend('bottomright', legend = paste('p-value <', signif(t.test(dmrs$age_glia_coef_mean, dmrs$age_neuron_coef_mean, paired = TRUE)$p.value, 3)), bty = 'n')

boxplot(abs(dmrs$age_glia_coef_mean), abs(dmrs$age_neuron_coef_mean), names = c('Glia', 'Neuron'), ylab = 'Absolute age mean coefficient')
legend('topright', legend = paste('p-value <', signif(t.test(abs(dmrs$age_glia_coef_mean), abs(dmrs$age_neuron_coef_mean), paired = TRUE)$p.value, 3)), bty = 'n')

d <- abs(dmrs$age_glia_coef_mean) - abs(dmrs$age_neuron_coef_mean) > 0
boxplot(abs(dmrs$age_glia_coef_mean) - abs(dmrs$age_neuron_coef_mean), xlab = paste0('Glia (', sum(d), ') - Neuron (', sum(!d), ')'), ylab = 'Absolute age mean coefficient difference')
legend('topright', legend = paste('p-value <', signif(t.test(abs(dmrs$age_glia_coef_mean) - abs(dmrs$age_neuron_coef_mean))$p.value, 3)), bty = 'n')

boxplot(dmrs$age_glia_coef_area, dmrs$age_neuron_coef_area, names = paste0(c('Glia', 'Neuron'), ' (', round(c(mean(dmrs$age_glia_coef_area > 0), mean(dmrs$age_neuron_coef_area > 0)) * 100, 0), '% up)'), ylab = 'Age coefficient area')
legend('topright', legend = paste('p-value <', signif(t.test(dmrs$age_glia_coef_area, dmrs$age_neuron_coef_area, paired = TRUE)$p.value, 3)), bty = 'n')



d2 <- dmrs$age_glia_coef_area - dmrs$age_neuron_coef_area > 0
boxplot(dmrs$age_glia_coef_area - dmrs$age_neuron_coef_area, xlab = paste0('Glia (', sum(d2), ') - Neuron (', sum(!d2), ')'), ylab = 'Age coefficient area difference')
legend('topright', legend = paste('p-value <', signif(t.test(dmrs$age_glia_coef_area - dmrs$age_neuron_coef_area)$p.value, 3)), bty = 'n')

dev.off()


print("Compare glia coef mean vs neurons")
addmargins(table('Glia > Neuron' = dmrs$age_glia_coef_mean > dmrs$age_neuron_coef_mean, 'abs(glia) > abs(neuron)' = abs(dmrs$age_glia_coef_mean) > abs(dmrs$age_neuron_coef_mean)))

print("Compare the two groups using coef mean and coef area (pages 3 and 5 of pdf/age_for_interaction_dmrs_diff.pdf)")
addmargins(table('abs(glia) > abs(neuron)' = abs(dmrs$age_glia_coef_mean) > abs(dmrs$age_neuron_coef_mean), 'Glia area > neuron area' = dmrs$age_glia_coef_area > dmrs$age_neuron_coef_area))

print('t test comparing absolute glia vs neuron mean coefficients')
t.test(abs(dmrs$age_glia_coef_mean), abs(dmrs$age_neuron_coef_mean), paired = TRUE)
print('Glia mean, then neuron mean, then the ratio (N/G) of the absolute coefficient mean')
mean(abs(dmrs$age_glia_coef_mean))
mean(abs(dmrs$age_neuron_coef_mean))
mean(abs(dmrs$age_neuron_coef_mean)) / mean(abs(dmrs$age_glia_coef_mean))

print('Glia median, then neuron median, then the ratio (N/G) of the absolute coefficient mean')
median(abs(dmrs$age_glia_coef_mean))
median(abs(dmrs$age_neuron_coef_mean))
median(abs(dmrs$age_neuron_coef_mean)) / median(abs(dmrs$age_glia_coef_mean))

print('Glia mean, then neuron mean, then the ratio (N/G) of the coefficient area')
mean(dmrs$age_glia_coef_area)
mean(dmrs$age_neuron_coef_area)
mean(dmrs$age_neuron_coef_area) / mean(dmrs$age_glia_coef_area)

print('Glia median, then neuron median, then the ratio (N/G) of the coefficient area')
median(dmrs$age_glia_coef_area)
median(dmrs$age_neuron_coef_area)
median(dmrs$age_neuron_coef_area) / median(dmrs$age_glia_coef_area)




t_cut <- function(var, cut = 1) {
    rbind(table(abs(mcols(dmrs)[, var]) < cut),
    round(table(abs(mcols(dmrs)[, var]) < cut) / length(dmrs) * 100, 2))
}

tcuts <- lapply(colnames(mcols(dmrs))[grep('tstat_mean', colnames(mcols(dmrs)))], t_cut, cut = quantile(abs(dmrs$age_cell_difference_tstat_mean), 0.025))
names(tcuts) <- colnames(mcols(dmrs))[grep('tstat_mean', colnames(mcols(dmrs)))]
tcuts





## Save results
save(dmrs, file = 'limma_Neuron_CpGs_minCov_3_ageInfo_dmrs.Rdata')


## Reproducibility info
proc.time()
message(Sys.time())
options(width = 120)
session_info()
