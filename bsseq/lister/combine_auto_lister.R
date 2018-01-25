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
    auto_long$lag <- as.factor(auto_long$lag)
    save(auto_long, file = 'rda/auto_long_combined.Rdata')
} else {
    load('rda/auto_long_combined.Rdata', verbose = TRUE)
}
dim(auto_long)

dir.create('pdf', showWarnings = FALSE)
pdf('pdf/autocorrelation_by_context.pdf', width = 7 * 2, height = 7 * 2)
ggplot(auto_long, aes(x = lag, y = acf_neuron)) + geom_boxplot() + facet_grid(. ~ context) + theme_bw(base_size = 30) + ggtitle('Neuron') + ylab('Auto correlation') + ylim(c(-1, 1))

ggplot(auto_long, aes(x = lag, y = acf_glia)) + geom_boxplot() + facet_grid(. ~ context) + theme_bw(base_size = 30) + ggtitle('Glia') + ylab('Auto correlation')  + ylim(c(-1, 1))
dev.off()

pdf('pdf/autocorrelation_by_context_abs.pdf', width = 7 * 2, height = 7 * 2)
ggplot(auto_long, aes(x = lag, y = abs(acf_neuron))) + geom_boxplot() + facet_grid(. ~ context) + theme_bw(base_size = 30) + ggtitle('Neuron') + ylab('Absolute auto correlation') + ylim(c(0, 1))

ggplot(auto_long, aes(x = lag, y = abs(acf_glia))) + geom_boxplot() + facet_grid(. ~ context) + theme_bw(base_size = 30) + ggtitle('Glia') + ylab('Absolute auto correlation') + ylim(c(0, 1))
dev.off()

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()