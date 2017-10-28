library('derfinder')
library('devtools')

bams <- dir('/dcl01/lieber/ajaffe/Amanda/ATAC/BAMs/combined_set/', pattern = '.final.bam$')

samples <- gsub('.final.bam', '', bams)
bams <- file.path('/dcl01/lieber/ajaffe/Amanda/ATAC/BAMs/combined_set/', bams)
names(bams) <- samples

stopifnot(length(unique(samples)) == length(bams))
stopifnot(all(file.exists(bams)))
stopifnot(all(file.exists(paste0(bams, '.bai'))))

atacFullCov <- fullCoverage(files = bams, chrs = paste0('chr', c(1:22, 'X', 'Y', 'M')), mc.cores = 5)

dir.create('rdas', showWarnings = FALSE)
save(atacFullCov, file = 'rdas/atacFullCov.Rdata')

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
