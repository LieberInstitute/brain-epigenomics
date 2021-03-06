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
library('ggplot2')
library('ggthemes')

dir.create('rda', showWarnings = FALSE)
dir.create('pdf', showWarnings = FALSE)

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
rm(meth)

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

# ## Relate DMRs to CpGs
# ov <- findOverlaps(cpg, dmrs)
# print('Number of CpGs that are inside the interaction DMRs')
# length(ov)
# length(unique(queryHits(ov)))
#
# ## Check that the CpGs for a given DMR are next to each other
# stopifnot(identical(nrun(Rle(subjectHits(ov))), length(unique(subjectHits(ov)))))


ov2 <- findOverlaps(dmrs, cpg)


## Save main pieces for later
save(limma_age, dmrs, ov2, file = 'rda/limma_Neuron_CpGs_minCov_3_ageInfo.Rdata')




## Compute mean, sum and area of the age values at the CpG level for each DMR
# age_dmr <- function(age, type) {
#     ## Group CpGs by the DMR they overlap
#     new <- split(data.frame(age[queryHits(ov), , drop = FALSE]), subjectHits(ov))
#
#     ## Compute the mean and sum of the age information
#     mean <- do.call(rbind, lapply(new, colMeans))[unique(subjectHits(ov)), , drop = FALSE]
#     sum <- do.call(rbind, lapply(new, colSums))[unique(subjectHits(ov)), , drop = FALSE]
#     area <- do.call(rbind, lapply(new, function(x) { colSums(abs(x)) }))[unique(subjectHits(ov)), , drop = FALSE]
#
#     ## Return the relevant values
#     colnames(mean) <- paste0(colnames(mean), '_', type, '_mean')
#     colnames(sum) <- paste0(colnames(sum), '_', type, '_sum')
#     colnames(area) <- paste0(colnames(area), '_', type, '_area')
#     cbind(mean, sum, area)
# }

age_dmr2 <- function(age, type) {
    ## Group CpGs by the DMR they overlap
    new <- split(data.frame(age[subjectHits(ov2), , drop = FALSE]), queryHits(ov2))
    
    ## Compute the mean and sum of the age information
    mean <- do.call(rbind, lapply(new, colMeans))
    sum <- do.call(rbind, lapply(new, colSums))
    area <- do.call(rbind, lapply(new, function(x) { colSums(abs(x)) }))
    
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
mcols(dmrs) <- cbind(mcols(dmrs), age_dmr2(limma_age$t[, match(colnames(limma_age$coef), colnames(limma_age$t))], 'tstat'), age_dmr2(limma_age$coef, 'coef'))
print('head of the dmr info with CpG results appended at the end')
head(dmrs)

print('Compare a bit the signs of some interaction coefficients')
addmargins(table('DMR coef sign' = sign(dmrs$coef), 'interaction mean coef sign' = sign(dmrs$age_cell_difference_coef_mean)))
addmargins(table('DMR value sign' = sign(dmrs$value), 'interaction mean coef sign' = sign(dmrs$age_cell_difference_coef_mean)))

## Group DMRs by glia and neuron mean coefficients
kdf <- data.frame('glia' = dmrs$age_glia_coef_mean,
    'neuron' = dmrs$age_neuron_coef_mean)
set.seed(20180201)
k2 <- kmeans(kdf, 2, nstart = 100)
k4 <- kmeans(kdf, 4, nstart = 100)
k6 <- kmeans(kdf, 6, nstart = 100)
k8 <- kmeans(kdf, 8, nstart = 100)

kadf <- data.frame('glia' = abs(dmrs$age_glia_coef_mean),
    'neuron' = abs(dmrs$age_neuron_coef_mean))
ka2 <- kmeans(kadf, 2, nstart = 100)
ka4 <- kmeans(kadf, 4, nstart = 100)
ka6 <- kmeans(kadf, 6, nstart = 100)
ka8 <- kmeans(kadf, 8, nstart = 100)



grey_lines <- function() {
    abline(h = 0, col = 'grey80')
    abline(v = 0, col = 'grey80')
    abline(a = 0, b = 1, col = 'grey80')
    abline(a = 0, b = -1, col = 'grey80')
}
cols <- brewer.pal(8, 'Dark2')


## Compare glia vs neuron mean coefficient with K-means clusters
pdf('pdf/glia_vs_neuron_mean_coef.pdf', useDingbats = FALSE)

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


## Main one with letters, just to identify the clusters
plot(x = dmrs$age_glia_coef_mean, y = dmrs$age_neuron_coef_mean, col = cols[k6$cluster], xlab = 'Glia mean coefficient', ylab = 'Neuron mean coefficient', pch = 20, cex = 0.5, main = 'K6 with clusters labeled a, b, ..., f')
points(k6$centers, col = 'black', pch = letters[1:6], cex = 2)
grey_lines()

## Again but with numbers (they don't go counter clockwise like I thought 
## they did)
plot(x = dmrs$age_glia_coef_mean, y = dmrs$age_neuron_coef_mean, col = cols[k6$cluster], xlab = 'Glia mean coefficient', ylab = 'Neuron mean coefficient', pch = 20, cex = 0.5, main = 'K6 with clusters labeled 1, 2, ..., 6')
points(k6$centers, col = 'black', pch = as.character(1:6), cex = 2)
grey_lines()

## Likely final picture
plot(x = dmrs$age_glia_coef_mean, y = dmrs$age_neuron_coef_mean, col = cols[k6$cluster], xlab = 'Glia age mean coefficient', ylab = 'Neuron age mean coefficient', pch = 20, cex = 0.5, main = 'Interaction DMRs grouped by the age coefficients for glia and neuron')
points(k6$centers, col = 'black', pch = 8, cex = 2)
grey_lines()


dev.off()


## Assign the cluster labels
dmrs$k6cluster <- k6$cluster
## These can change order, but setting the seed earlier should fix this
k6_labels <- c('1:G-N+', '2:G0N+', '3:G0N-', '4:G+N0', '5:G+N-', '6:G-N0')
dmrs$k6cluster_label <- factor(k6_labels[k6$cluster], levels = k6_labels)

## Make a quick plot comparing many of the variables
dmrs_df <- as.data.frame(mcols(dmrs))


custom_col <- function(p) {
    ## https://stackoverflow.com/questions/34740210/how-to-change-the-color-palette-for-ggallyggpairs
    for(i in 1:p$nrow) {
      for(j in 1:p$ncol){
        p[i,j] <- p[i,j] + 
            scale_fill_manual(values = cols[1:6]) +
            scale_color_manual(values = cols[1:6]) +
            theme_minimal(base_size = 15)
            # https://cran.r-project.org/web/packages/ggthemes/vignettes/ggthemes.html
      }
    }
    ## Disable the axis label for 1,1 because it's a density label
    ## and doesn't match the rest
    ## https://stackoverflow.com/questions/35090883/remove-all-of-x-axis-labels-in-ggplot
    p[1, 1] <- p[1, 1] + theme(axis.text.y = element_blank())
    return(p)
}

## Plot the t-stats and the coefficient metrics against each other,
## colored by the 6 clusters we labeled already
pdf('pdf/age_for_interaction_dmrs.pdf', width = 10, height = 10, useDingbats = FALSE)
custom_col(ggpairs(dmrs_df, columns = grep('tstat_mean',
    colnames(dmrs_df)), mapping = aes(color = k6cluster_label),
    upper = list(continuous = wrap("cor", size = 4.75, alignPercent = 0.8))))
    
custom_col(ggpairs(dmrs_df, columns = grep('tstat_sum',
    colnames(dmrs_df)), mapping = aes(color = k6cluster_label),
    upper = list(continuous = wrap("cor", size = 4.75, alignPercent = 0.8))))
        
custom_col(ggpairs(dmrs_df, columns = grep('tstat_area',
    colnames(dmrs_df)), mapping = aes(color = k6cluster_label),
    upper = list(continuous = wrap("cor", size = 4.75, alignPercent = 0.8))))
   

custom_col(ggpairs(dmrs_df, columns = grep('coef_mean',
    colnames(dmrs_df)), mapping = aes(color = k6cluster_label),
    upper = list(continuous = wrap("cor", size = 4.75, alignPercent = 0.8))))
    
custom_col(ggpairs(dmrs_df, columns = grep('coef_sum',
    colnames(dmrs_df)), mapping = aes(color = k6cluster_label),
    upper = list(continuous = wrap("cor", size = 4.75, alignPercent = 0.8))))
    
custom_col(ggpairs(dmrs_df, columns = grep('coef_area',
    colnames(dmrs_df)), mapping = aes(color = k6cluster_label),
    upper = list(continuous = wrap("cor", size = 4.75, alignPercent = 0.8))))  
    
custom_col(ggpairs(dmrs_df, columns = c(which(colnames(dmrs_df) == 'value'),
    grep('coef_mean', colnames(dmrs_df))),
    mapping = aes(color = k6cluster_label),
    upper = list(continuous = wrap("cor", size = 4.75, alignPercent = 0.8))))
    
custom_col(ggpairs(dmrs_df, columns = c(which(colnames(dmrs_df) == 'area'),
    grep('coef_area', colnames(dmrs_df))),
    mapping = aes(color = k6cluster_label),
    upper = list(continuous = wrap("cor", size = 4.75, alignPercent = 0.8))))
        
dev.off()  
    

# custom_col(ggpairs(dmrs_df, columns = c(which(colnames(dmrs_df) == 'value'),
#     grep('coef_sum', colnames(dmrs_df))),
#     mapping = aes(color = k6cluster_label),
#     upper = list(continuous = wrap("cor", size = 4.75, alignPercent = 0.8))))
#

    
pdf('pdf/interaction_dmrs_cpg_vs_dmr_metrics.pdf', width = 20, height = 20, useDingbats = FALSE)
dmrs_df2 <- dmrs_df
colnames(dmrs_df2) <- gsub('age_cell_difference_', '', colnames(dmrs_df2))
custom_col(ggpairs(dmrs_df2,
    columns = c(which(colnames(dmrs_df) %in% c('value', 'area', 'coef')),
    grep('age_cell_difference', colnames(dmrs_df))),
    mapping = aes(color = k6cluster_label),
    upper = list(continuous = wrap("cor", size = 4.75, alignPercent = 0.8))))
    
custom_col(ggpairs(dmrs_df2,
    columns = c(which(colnames(dmrs_df) %in% c('value', 'coef')),
    grep('age_cell_difference_coef_[area|mean]', colnames(dmrs_df))),
    mapping = aes(color = k6cluster_label),
    upper = list(continuous = wrap("cor", size = 4.75, alignPercent = 0.8)))) 
        
dev.off()

## Could this be because of the smoothed CpG data?
meth2 <- getMeth(BSobj, type = 'smooth')
print('Number of CpGs and samples being considered')
dim(meth2)
pd2 <- pData(BSobj)
pd2$Cell.Type <- relevel(factor(pd2$Cell.Type), 'Glia')
system.time( f_smooth <- eBayes(lmFit(meth2, with(pd2, model.matrix(~ Age * Cell.Type)))) )
rm(meth2)

print('Compare interaction coefs for all CpGs between using the smoothed and raw meth data')
identical(f_smooth$coef[, 'Age:Cell.TypeNeuron'], limma_age$coef[, 'age_cell_difference'])
summary(abs(f_smooth$coef[, 'Age:Cell.TypeNeuron'] - limma_age$coef[, 'age_cell_difference']))

print('Compare interaction t-stats for all CpGs between using the smoothed and raw meth data')
identical(f_smooth$t[, 'Age:Cell.TypeNeuron'], limma_age$t[, 'age_cell_difference'])
summary(abs(f_smooth$t[, 'Age:Cell.TypeNeuron'] - limma_age$t[, 'age_cell_difference']))

## For plotting
dmrs_int <- cbind(
    age_dmr2(f_smooth$coef[, 'Age:Cell.TypeNeuron', drop = FALSE], 'coef'),
    age_dmr2(f_smooth$t[, 'Age:Cell.TypeNeuron', drop = FALSE], 'tstat'),
    dmrs_df[, which(colnames(dmrs_df) %in% c('value', 'area', 'coef', 'k6cluster_label'))]
)
colnames(dmrs_int) <- gsub('Age.Cell.TypeNeuron', 'intSmooth', colnames(dmrs_int))

## Save for later
save(dmrs_int, file = 'rda/limma_Neuron_CpGs_minCov_3_ageInfo_interaction_smooth.Rdata')

pdf('pdf/interaction_dmrs_cpg_vs_dmr_metrics_with_smooth.pdf', width = 15, height = 15, useDingbats = FALSE)

custom_col(ggpairs(dmrs_int,
    columns = c(which(colnames(dmrs_int) %in% c('value', 'area', 'coef')),
    grep('_mean', colnames(dmrs_int))),
    mapping = aes(color = k6cluster_label),
    upper = list(continuous = wrap("cor", size = 4.75, alignPercent = 0.8))))  
   
custom_col(ggpairs(dmrs_int,
    columns = c(which(colnames(dmrs_int) %in% c('value', 'area', 'coef')),
    grep('_sum', colnames(dmrs_int))),
    mapping = aes(color = k6cluster_label),
    upper = list(continuous = wrap("cor", size = 4.75, alignPercent = 0.8))))
        
custom_col(ggpairs(dmrs_int,
    columns = c(which(colnames(dmrs_int) %in% c('value', 'area', 'coef')),
    grep('_area', colnames(dmrs_int))),
    mapping = aes(color = k6cluster_label),
    upper = list(continuous = wrap("cor", size = 4.75, alignPercent = 0.8))))

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


pdf('pdf/age_for_interaction_dmrs_diff.pdf', useDingbats = FALSE)

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



## Using the absolute t-stat means, compute the number and percent below some
## cutoffs
print('Checking the absolute t-stat means for the overall age changes')
t_cut <- function(var, cut = 1) {
    x <- t(rbind(table(abs(mcols(dmrs)[, var]) < cut),
    round(table(abs(mcols(dmrs)[, var]) < cut) / length(dmrs) * 100, 2)))
    colnames(x) <- c('n', 'percent')
    rownames(x) <- c('Above cut', 'Below cut')
    return(x)
}

cutoff <- quantile(abs(dmrs$age_cell_difference_tstat_mean), 0.025)
print('Cutoff (based on 2.5% quantile of the abs interaction t-stat)')
cutoff
tcuts <- lapply(colnames(mcols(dmrs))[grep('tstat_mean', colnames(mcols(dmrs)))], t_cut, cut = cutoff)
names(tcuts) <- colnames(mcols(dmrs))[grep('tstat_mean', colnames(mcols(dmrs)))]
tcuts

tcuts2 <- lapply(colnames(mcols(dmrs))[grep('tstat_mean', colnames(mcols(dmrs)))], t_cut, cut = 2)
names(tcuts2) <- names(tcuts)
print('Cutoff: 2')
tcuts2


pdf('pdf/age_for_interaction_dmrs_by_overall_mean_tstat.pdf', useDingbats = FALSE)
plot(x = dmrs$overall_tstat_mean, y = dmrs$age_cell_difference_tstat_mean, col = cols[k6$cluster], xlab = 'Overall age mean t-statistic', ylab = 'Interaction mean t-statistic', pch = 20, cex = 0.5, main = 't-statistics for global age vs interaction', sub = 'Lines at interaction quantile 2.5% and 2')
abline(v = c(2, -2), col = 'grey80')
abline(v = c(cutoff, - cutoff), col = 'grey80')


plot(x = dmrs$age_glia_coef_mean, y = dmrs$age_neuron_coef_mean, col = cols[k6$cluster], xlab = 'Glia age mean coefficient', ylab = 'Neuron age mean coefficient', pch = ifelse(abs(dmrs$overall_tstat_mean) < cutoff, 20, 19), cex = 0.5, main = 'Interaction DMRs grouped by the age coefficients for glia and neuron', sub = 'Larger points for age mean t-statistic > interaction quantile 2.5%')
points(k6$centers, col = 'black', pch = 8, cex = 2)
grey_lines()

plot(x = dmrs$age_glia_coef_mean, y = dmrs$age_neuron_coef_mean, col = cols[k6$cluster], xlab = 'Glia age mean coefficient', ylab = 'Neuron age mean coefficient', pch = ifelse(abs(dmrs$overall_tstat_mean) < 2, 20, 19), cex = 0.5, main = 'Interaction DMRs grouped by the age coefficients for glia and neuron', sub = 'Larger points for age mean t-statistic > 2')
points(k6$centers, col = 'black', pch = 8, cex = 2)
grey_lines()
dev.off()

## Save info for later (note reversal of sign: we want the TRUE ones here)
dmrs$global_tcut_cut_quantile <- abs(dmrs$overall_tstat_mean) > cutoff
dmrs$global_tcut_cut_2 <- abs(dmrs$overall_tstat_mean) > 2

## Save results
save(dmrs, file = 'rda/limma_Neuron_CpGs_minCov_3_ageInfo_dmrs.Rdata')


## Reproducibility info
proc.time()
message(Sys.time())
options(width = 120)
session_info()
