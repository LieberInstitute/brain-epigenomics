library('devtools')
library('ggplot2')

files <- dir('rda', pattern = '^auto_long_chr')

load_auto <- function(f) {
    chr <- gsub('_.*', '', gsub('auto_long_', '', f))
    message(paste(Sys.time(), 'loading', f))
    load(file.path('rda', f))
    auto_long$chr <- chr
    return(auto_long)
}

dir.create('rda', showWarnings = FALSE)
if(!file.exists('rda/auto_long_combined.Rdata')) {
    auto_long <- do.call(rbind, lapply(files, load_auto))
    save(auto_long, file = 'rda/auto_long_combined.Rdata')
} else {
    load('rda/auto_long_combined.Rdata', verbose = TRUE)
}
dim(auto_long)

## This never finished running:
# context_summary <- lapply(unique(auto_long)$context, function(context) {
#     sub <- auto_long[auto_long$context == context, ]
#     tapply(sub$acf, sub$lag, summary)
# })
# context_summary
# save(context_summary, file = 'autocorrelation_summary_by_context.Rdata')

auto_long$lag <- as.factor(auto_long$lag)

png('autocorrelation_by_context.png', width = 480 * 2)
ggplot(auto_long, aes(x = lag, y = acf)) + geom_boxplot() + facet_grid(. ~ context) + theme_bw(base_size = 14)
dev.off()

png('autocorrelation_by_context_abs.png', width = 480 * 2)
ggplot(auto_long, aes(x = lag, y = abs(acf))) + geom_boxplot() + facet_grid(. ~ context) + theme_bw(base_size = 14)
dev.off()

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
