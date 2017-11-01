# qrsh -l bluejay,mem_free=150G,h_vmem=150G,h_fsize=100G
# cd /dcl01/lieber/ajaffe/lab/brain-epigenomics/bumphunting
# module load conda_R/3.4.x
# R
library('limma')
library('devtools')

system.time( load('limma_Neuron_CpGs_minCov_3.Rdata') )

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
age_p <- cbind(age_p, 'age_neuron_F' = top$P.value)
head(age_p)
summary(age_p)

identical(age_p[, 4], age_p[, 5])
summary(abs(age_p[, 4] - age_p[, 5]))

## Save the results
print(object.size(age_coef), units = 'Mb')
print(object.size(age_t), units = 'Mb')
print(object.size(age_p), units = 'Mb')

save(age_coef, age_t, age_p, file = 'limma_Neuron_CpGs_minCov_3_ageInfo.Rdata')

## Reproducibility info
proc.time()
message(Sys.time())
options(width = 120)
session_info()
