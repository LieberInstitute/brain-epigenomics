## Based on https://github.com/LieberInstitute/brain-epigenomics/blob/master/bsseq/bsobj_by_chr/combine_auto.R

library('devtools')
library('ggplot2')

files <- dir('rda', pattern = '^auto_long', full.names = TRUE)

load_auto <- function(f) {
    message(paste(Sys.time(), 'loading', f))
    load(f)
    return(auto_long)
}

if(!file.exists('auto_long_combined.Rdata')) {
    auto_long <- do.call(rbind, lapply(files, load_auto))
    save(auto_long, file = 'auto_long_combined.Rdata')
} else {
    load('auto_long_combined.Rdata')
}
dim(auto_long)

## This never finished running:
# mod_context_summary <- lapply(unique(auto_long)$model, function(model) {
#     lapply(unique(auto_long)$context, function(context) {
#          sub <- auto_long[auto_long$context == context & auto_long$model == model, ]
#          tapply(sub$acf, sub$lag, summary)
#      })
# })
# mod_context_summary
# save(mod_context_summary, file = 'autocorrelation_summary_by_context_and_model.Rdata')

auto_long$lag <- as.factor(auto_long$lag)

dir.create('pdf', showWarnings = FALSE)
png('pdf/autocorrelation_by_context_and_model.png', width = 480 * 2, height = 480 * 3, type = 'cairo')
ggplot(auto_long, aes(x = lag, y = acf)) + geom_boxplot() + facet_grid(model ~ context) + theme_bw(base_size = 14)
dev.off()

png('pdf/autocorrelation_by_context_and_model_abs.png', width = 480 * 2, height = 480 * 3, type = 'cairo')
ggplot(auto_long, aes(x = lag, y = abs(acf))) + geom_boxplot() + facet_grid(model ~ context) + theme_bw(base_size = 14)
dev.off()

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
