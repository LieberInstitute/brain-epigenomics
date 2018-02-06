library('SGSeq')
library('bsseq')
library('getopt')
library('devtools')
library('limma')
library('gplots')
library('VennDiagram')

## Specify parameters
spec <- matrix(c(
    'feature', 'f', 1, 'character', 'Either: gene, exon, jx or psi',
	'help' , 'h', 0, 'logical', 'Display help'
), byrow=TRUE, ncol=5)
opt <- getopt(spec)

## if help was asked for print a friendly message
## and exit with a non-zero error code
if (!is.null(opt$help)) {
	cat(getopt(spec, usage=TRUE))
	q(status=1)
}

## For testing
if(FALSE) {
    opt <- list('feature' = 'gene')
    opt <- list('feature' = 'exon')
    ## Not yet adapted/designed for the other two
    #opt <- list('feature' = 'jx')
    #opt <- list('feature' = 'psi')
}

stopifnot(opt$feature %in% c('psi', 'gene', 'exon', 'jx'))
cpgs <- c('CpG', 'nonCpG', 'CpGmarg')

## Load annotated meQTL results at FDR 5%
## for the CpG data, I'm using the results near the nonCpG meQTLs
## For the marginal CpG data I'm using the results from 'near'
## just dropping those with infinite statistics
mres <- lapply(cpgs, function(cpg) {
    if(cpg == 'CpG') {
        f <- paste0('rda/me_annotated_FDR5_', cpg, '_', opt$feature,
            '_near_nonCpG_meQTLs.Rdata')
    } else if (cpg == 'nonCpG') {
        f <- paste0('rda/me_annotated_FDR5_', cpg, '_', opt$feature, '.Rdata')
    } else if (cpg == 'CpGmarg') {
        f <- paste0('rda/me_CpG_', opt$feature, '_near_nonCpG_meQTLs.Rdata')
    }
    message(paste(Sys.time(), 'loading the file', f))
    load(f, verbose = TRUE)
    
    ## Just keep those with a finite statistic
    if(cpg == 'CpGmarg') {
        print('Keeping only those with finite statistics')
        print(table(is.finite(me$cis$eqtls$statistic)))
        me_annotated <- list('eqtls' = me$cis$eqtls[is.finite(me$cis$eqtls$statistic), ])
    }
    
    return(me_annotated)
})
names(mres) <- cpgs

load('/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/sorted_nuclear_RNA/DE_limma_results_objects.rda', verbose = TRUE)

## Use the combined results to determine the neuronal genes
if(opt$feature == 'gene') {
    co <- grep('CellType', colnames(fit_gene_combined$p.value))
    print(colnames(fit_gene_combined$p.value)[co])
    top <- topTable(fit_gene_combined, coef = co, sort.by = 'none',
        number = nrow(fit_gene_combined$p.value))
} else if (opt$feature == 'exon') {

    co <- grep('CellType', colnames(fit_exon_combined$p.value))
    print(colnames(fit_exon_combined$p.value)[co])
    top <- topTable(fit_exon_combined, coef = co, sort.by = 'none',
        number = nrow(fit_exon_combined$p.value))
        
    ## Exon names don't match, load the required data
    load("/dcl01/lieber/ajaffe/CellSorting/RNAseq_pipeline/rawCounts_CellSorting_July5_n12.rda", verbose = TRUE)
    load('/dcl01/lieber/ajaffe/lab/brain-epigenomics/brainseq_pipeline/polyA_unstranded/rse_exon_polyA_dlpfc_n41.Rdata', verbose = TRUE)

    ## Note that nrow(top) != nrow(exonMap)
    m1 <- match(rownames(top), rownames(exonMap))
    check_gr <- GRanges(exonMap[m1, ])

    ov <- findOverlaps(check_gr, rowRanges(rse_exon), type = 'equal',
        ignore.strand = FALSE)
    stopifnot(length(unique(queryHits(ov))) == length(unique(subjectHits(ov))))
    ## Names don't match
    table(names(check_gr[queryHits(ov)]) == names(rowRanges(rse_exon)[subjectHits(ov)]))

    ## Drop those not present in the meQTL results
    dim(top)
    top <- top[queryHits(ov), ]
    dim(top)
    ## Finally fix the names
    rownames(top) <- names(rowRanges(rse_exon))[subjectHits(ov)]

    ## Clean up
    rm(check_gr, rse_exon, ov, m1, exonCounts, exonMap, geneCounts, geneMap, jCounts, jMap, metrics, txMap, txNumReads, getRPKM)
        
}
dim(top)
print('Number of DE features at FDR 5%')
table(top$adj.P.Val < 0.05)

## Keep the top 5k DE features (there's enough for genes and exons)
top <- top[top$adj.P.Val < 0.05, ] 

## Keep only those with negative log FC (higher in Neurons)
top5k <- top[sign(top$logFC) == -1, ]
nrow(top5k)
top5k <- top5k[order(top5k$adj.P.Val, decreasing = FALSE), ]
top5k <- head(top5k, 5000)
nrow(top5k)

vinfo <- lapply(mres, function(me) { unique(as.character(me$eqtls$gene)) })
vinfo5k <- lapply(vinfo, function(me) { me[me %in% rownames(top5k)] })

dir.create('pdf', showWarnings = FALSE)
pdf(paste0('pdf/meqtl_venn_', opt$feature, '_using_near.pdf'))
vennres <- venn(vinfo) + title('meQTLs at FDR 5%, CpGs only in proximity to nonCpG')
vennres5k <- venn(vinfo5k) + title(paste0('meQTLs at FDR 5%, CpGs only in proximity to nonCpG\nBased on top ', nrow(top5k), ' ', opt$feature, 's expressed in Neurons'))
dev.off()

## Prettier venn diagrams
pdf(paste0('pdf/meqtl_venn_pretty_', opt$feature, '_using_near.pdf'))
v <- venn.diagram(vinfo, filename = NULL,
    main = 'meQTLs at FDR 5%, CpGs only in proximity to nonCpG',
    col = "transparent", fill = c("lightpink2","cornflowerblue", "olivedrab2"),
    alpha = 0.50, fontface = "bold",
    cat.col = c("palevioletred4", "darkblue", "olivedrab4"),
    margin=0.2)
grid.draw(v)

v5k <- venn.diagram(vinfo5k, filename = NULL,
    main = paste0('meQTLs at FDR 5%, CpGs only in proximity to nonCpG\nBased on top ', nrow(top5k), ' ', opt$feature, 's expressed in Neurons'),
    col = "transparent", fill = c("lightpink2","cornflowerblue", "olivedrab2"),
    alpha = 0.50, ffontface = "bold",
    cat.col = c("palevioletred4", "darkblue", "olivedrab4"),
    margin=0.2)
grid.draw(v5k)
dev.off()


dir.create('rda', showWarnings = FALSE)
save(vennres, vennres5k, top5k, top, v, v5k, file = paste0('rda/meqtl_venn_',
    opt$feature, '_using_near.Rdata'))


message(paste(Sys.time(), 'summarizing the meQTL data by', opt$feature))
m_summary <- do.call(rbind, lapply(1:length(mres), function(i) {
    message(paste(Sys.time(), 'processing', names(mres)[i]))
    me <- mres[[i]]
    gdata <- split(me$eqtls, me$eqtls$gene)
    gdata <- gdata[elementNROWS(gdata) > 0]
    typeres <- do.call(rbind, lapply(gdata, function(g) {
        best <- which.min(g$FDR)
        res <- g[best, , drop = FALSE]
        res$n_meqtls <- nrow(g)
        return(res)
    }))
    typeres$type <- names(mres)[i]
    return(typeres)
}))

m_summary$gene <- as.character(m_summary$gene)
m_summary$snps <- as.character(m_summary$snps)

message(paste(Sys.time(), 'saving summary of the meQTL data by', opt$feature))
save(m_summary, file = paste0('rda/meqtl_summary_', opt$feature,
    '_using_near.Rdata'))


common <- names(table(m_summary$gene)[table(m_summary$gene) == 2])
find_pval <- function(type, t5k = FALSE) {
    res <- -log10(m_summary$FDR[m_summary$gene %in% common & m_summary$type == type])
    if(t5k) {
        res <- res[common %in% rownames(top5k)]
    }
    return(res)
}

pdf(paste0('pdf/scatter_FDR_', opt$feature, '_using_near.pdf'))
plot(x = find_pval('CpG'), y = find_pval('nonCpG'), pch = 20, xlab = 'Best -log10 FDR p-value with CpGs', ylab = 'Best -log10 FDR p-value with non CpGs', main = paste(length(common), 'common', opt$feature, 'meQTLs'))
plot(x = find_pval('CpG', TRUE), y = find_pval('nonCpG', TRUE), pch = 20, xlab = 'Best -log10 FDR p-value with CpGs', ylab = 'Best -log10 FDR p-value with non CpGs', main = paste0(sum(common %in% rownames(top5k)), ' common ', opt$feature, ' meQTLs\nBased on top ', nrow(top5k), ' ', opt$feature, 's expressed in Neurons'))
dev.off()



## Reproducibility info
proc.time()
message(Sys.time())
options(width = 120)
session_info()

