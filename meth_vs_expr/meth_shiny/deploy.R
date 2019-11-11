dir.create('www', showWarnings = FALSE)
system('scp e:/dcl01/lieber/ajaffe/lab/brain-epigenomics/meth_vs_expr/rda/meth_data.Rdata .')
system('scp e:/dcl01/lieber/ajaffe/lab/brain-epigenomics/meth_vs_expr/rda/meth_df.Rdata .')
system('scp e:/dcl01/lieber/ajaffe/lab/brain-epigenomics/meth_vs_expr/rda/tf_data.Rdata .')
system('scp e:/dcl01/lieber/ajaffe/lab/brain-epigenomics/meth_vs_expr/rda/meth_summary.Rdata .')

## Add gene_symbol
library('GenomicRanges')
gr <- lapply(meth_data[c('gene', 'exon')], function(x) {
    unlist(GRangesList(lapply(x, function(y) {
        rowRanges(y$expr)
    })))
})
system('scp e:/dcl01/lieber/ajaffe/lab/brain-epigenomics/brainseq_pipeline/polyA_unstranded/rse_gene_polyA_dlpfc_n41.Rdata .')
library('SummarizedExperiment')
load('rse_gene_polyA_dlpfc_n41.Rdata')

load('meth_df.Rdata')
meth_df$symbol <- NA
meth_df$symbol[meth_df$feature != 'exon'] <- rowRanges(rse_gene)$Symbol[match(meth_df$feature_id[meth_df$feature != 'exon'], rowRanges(rse_gene)$gencodeID)]
meth_df$symbol[meth_df$feature == 'exon'] <- gr$exon$Symbol[match(meth_df$feature_id[meth_df$feature == 'exon'], gr$exon$exon_libdID)]
meth_df$symbol <- as.factor(meth_df$symbol)
table(is.na(meth_df$symbol))
save(meth_df, file = 'meth_df_withSymbol.Rdata')

library('rsconnect')
load('.deploy_info.Rdata')
rsconnect::setAccountInfo(name=deploy_info$name, token=deploy_info$token,
    secret=deploy_info$secret)
options(repos = BiocManager::repositories())
rsconnect::deployApp(appFiles = c('ui.R', 'server.R', 'meth_df_withSymbol.Rdata',
    'meth_data.Rdata', 'global.R',
    'tf_data.Rdata', 'google-analytics.js', 'www/LICENSE.txt'),
    appName = 'wgbsExprs', account = 'jhubiostatistics', server = 'shinyapps.io')
Y
