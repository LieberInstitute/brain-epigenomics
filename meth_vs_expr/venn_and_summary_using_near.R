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

## For the CpGmarg
load_dmp <- function(is_cpg) {
    if(is_cpg) {
        load('/dcl01/lieber/ajaffe/lab/brain-epigenomics/bumphunting/BSobj_bsseqSmooth_Neuron_minCov_3.Rdata', verbose = TRUE)
        rowRanges(BSobj)$c_context <- Rle('CG')
        rowRanges(BSobj)$trinucleotide_context <- Rle('CpG')
    } else {
        load('/dcl01/lieber/ajaffe/lab/brain-epigenomics/bsseq/bsobj_by_chr/allChrs_postNatal_cleaned_nonCG_noHomogenate_highCov.Rdata', verbose = TRUE)
    }
    ## Keep Neurons only
    BSobj <- BSobj[, colData(BSobj)$Cell.Type == 'Neuron']
    return(BSobj)
}

## Load data
load_expr <- function(type) {
    if(type == 'psi') {
        load('/dcl01/lieber/ajaffe/lab/brain-epigenomics/psi/rda/sgvc10_postnatal_polyA.Rdata', verbose = TRUE)
        expr <- sgvc10
    } else if (type == 'gene') {
        load('/dcl01/lieber/ajaffe/lab/brain-epigenomics/brainseq_pipeline/polyA_unstranded/rse_gene_polyA_dlpfc_n41.Rdata', verbose = TRUE)
        assays(rse_gene)$norm <- recount::getRPKM(rse_gene, 'Length')
        ## RPKM
        expr <- rse_gene
    } else if (type == 'exon') {
        load('/dcl01/lieber/ajaffe/lab/brain-epigenomics/brainseq_pipeline/polyA_unstranded/rse_exon_polyA_dlpfc_n41.Rdata', verbose = TRUE)
        assays(rse_exon)$norm <- recount::getRPKM(rse_exon, 'Length')
        ## RPKM
        expr <- rse_exon
    } else if (type == 'jx') {
                           load('/dcl01/lieber/ajaffe/lab/brain-epigenomics/brainseq_pipeline/polyA_unstranded/rse_jx_polyA_dlpfc_n41.Rdata', verbose = TRUE)
        rowRanges(rse_jx)$Length <- 100 / 8
        ## RP80m
        assays(rse_jx)$norm <- recount::getRPKM(rse_jx, 'Length')
        expr <- rse_jx
    }
    
    ## Subset to postnatal and PolyA only
    expr <- expr[, colData(expr)$Age > 0 & colData(expr)$Experiment == 'PolyA']
    
    ## Drop low expressed features
    if(type != 'psi') {
        if(!file.exists(paste0('rda/expr_', opt$feature, '_unfiltered.Rdata'))) {
            dir.create('pdf', showWarnings = FALSE)
            pdf(paste0('pdf/suggested_expr_cutoffs_', tolower(opt$feature),
                '.pdf'), width = 12)
            cuts <- expression_cutoff(assays(expr)$norm, seed = 20180119)
            message(paste(cuts, collapse = ' '))
            cut <- max(cuts)
            dev.off()
    
            meanExpr <- rowMeans(assays(expr)$norm)
            rowRanges(expr)$meanExprs <- meanExpr
            rowRanges(expr)$passExprsCut <- meanExpr > cut
            dir.create('rda', showWarnings = FALSE)
            save(expr, file = paste0('rda/expr_', opt$feature, '_unfiltered.Rdata'))
        }
    }
        
    return(expr)
}



## For matching the brain ids
getid <- function(x) {
    as.integer(gsub('Br', '', x))
}

## Subset to around a 1kb window from the nonCpG meQTl results
load_meqtl <- function() {
    load(paste0('rda/me_annotated_FDR5_nonCpG_', opt$feature,
        '.Rdata'), verbose = TRUE)
    return(me_annotated)
}

## Load annotated meQTL results at FDR 5%
## for the CpG data, I'm using the results near the nonCpG meQTLs
## For the marginal CpG data I'm using the results from 'near'
## just dropping those with infinite statistics
## Keep only those with meth_n >= 4
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
        message(paste(Sys.time(), 'Keeping only those with finite statistics'))
        print(table(is.finite(me$cis$eqtls$statistic)))
        me_annotated <- list('eqtls' = me$cis$eqtls[is.finite(me$cis$eqtls$statistic), ], 'meth' = 1, 'expr' = 2)
        
        
        ## For getting the meth_n part
        message(paste(Sys.time(), 'computing the parts for meth_n'))
        expr <- load_expr(opt$feature)
        meqtl <- load_meqtl()
        me_ov <- findOverlaps(resize(rowRanges(meqtl$expr), width(rowRanges(meqtl$expr)) + 2000, fix = 'center'), expr)
        expr <- expr[sort(unique(subjectHits(me_ov))), ]
        rm(me_ov, meqtl)
        BSobj <- load_dmp(TRUE)
        m  <- match(getid(colData(expr)$BrNum), getid(colData(BSobj)$Brain.ID))
        expr <- expr[, which(!is.na(m))]
        colnames(expr) <- paste0('Br', getid(colData(expr)$BrNum))
        BSobj <- BSobj[, m[!is.na(m)]]
        colnames(BSobj) <- paste0('Br', getid(colData(BSobj)$Brain.ID))
        cp_ov <- findOverlaps(resize(rowRanges(expr), width(rowRanges(expr)) + 2000, fix = 'center'), BSobj)
        BSobj <- BSobj[sort(unique(subjectHits(cp_ov))), ]
        
        ## Finally get meth_n
        snp_id <- function(x) { as.integer(gsub('row', '', x)) }
        
        BSobj <- BSobj[snp_id(me_annotated$eqtls$snps), ]
        
        meth_n <- rowSums(getMeth(BSobj, type = 'raw') > 0)
        me_annotated$eqtls <- me_annotated$eqtls[meth_n >= 4, ]
        me_annotated$meth <- BSobj[meth_n >=4, ]
        me_annotated$expr <- expr[match(me_annotated$eqtls$gene, names(rowRanges(expr))), ]
        me_annotated$eqtls$meth_n <- meth_n[meth_n >= 4]
    } else {
        meth_n <- rowSums(getMeth(me_annotated$meth, type = 'raw') > 0)
        me_annotated$eqtls <- me_annotated$eqtls[meth_n >= 4, ]
        me_annotated$meth <- me_annotated$meth[meth_n >= 4, ]
        me_annotated$expr <- me_annotated$expr[meth_n >= 4, ]
        me_annotated$eqtls$meth_n <- meth_n[meth_n >= 4]
    }
    
    message(paste(Sys.time(), 'computing the expression delta'))
    me_annotated$eqtls$expr_delta <- apply(log2(assays(me_annotated$expr)$norm + 1), 1, function(x) diff(range(x)) )
    
    return(me_annotated)
})
names(mres) <- cpgs

dir.create('rda', showWarnings = FALSE)
message(paste(Sys.time(), 'saving the mres object'))
save(mres, file = paste0('rda/meqtl_mres_', opt$feature,
    '_using_near.Rdata'))

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

vennres <- venn(vinfo, show.plot = FALSE)
vennres5k <- venn(vinfo5k, show.plot = FALSE)

## Prettier venn diagrams
pdf(paste0('pdf/meqtl_venn_pretty_', opt$feature, '_using_near.pdf'))
v <- venn.diagram(vinfo, filename = NULL,
    main = 'meQTLs at FDR 5%, CpGs only in proximity to nonCpG',
    col = "transparent", fill = c("lightpink2","cornflowerblue", "olivedrab2"),
    alpha = 0.50, fontface = "bold",
    cat.col = c("palevioletred4", "darkblue", "olivedrab4"),
    margin=0.2)
grid.newpage()
grid.draw(v)

v5k <- venn.diagram(vinfo5k, filename = NULL,
    main = paste0('meQTLs at FDR 5%, CpGs only in proximity to nonCpG\nBased on top ', nrow(top5k), ' ', opt$feature, 's expressed in Neurons'),
    col = "transparent", fill = c("lightpink2","cornflowerblue", "olivedrab2"),
    alpha = 0.50, fontface = "bold",
    cat.col = c("palevioletred4", "darkblue", "olivedrab4"),
    margin=0.2)
grid.newpage()
grid.draw(v5k)
dev.off()


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


common <- names(table(m_summary$gene)[table(m_summary$gene) == 3])
find_pval <- function(type, t5k = FALSE) {
    res <- -log10(m_summary$FDR[m_summary$gene %in% common & m_summary$type == type])
    ## Add the name
    names(res) <- m_summary$gene[m_summary$gene %in% common & m_summary$type == type]
    ## Match the order
    res <- res[match(common, names(res))]
    if(t5k) {
        res <- res[names(res) %in% rownames(top5k)]
    }
    return(res)
}

stopifnot(identical(names(find_pval('CpG')), names(find_pval('nonCpG'))))

## Scatter plot of -log 10 p-values for the CpG and nonCpG data (FDR adjusted p-values)
pdf(paste0('pdf/scatter_FDR_', opt$feature, '_using_near.pdf'))
plot(x = find_pval('CpG'), y = find_pval('nonCpG'), pch = 20, xlab = 'Best -log10 FDR p-value with CpGs', ylab = 'Best -log10 FDR p-value with non CpGs', main = paste(length(common), 'common', opt$feature, 'meQTLs'))
plot(x = find_pval('CpG', TRUE), y = find_pval('nonCpG', TRUE), pch = 20, xlab = 'Best -log10 FDR p-value with CpGs', ylab = 'Best -log10 FDR p-value with non CpGs', main = paste0(sum(common %in% rownames(top5k)), ' common ', opt$feature, ' meQTLs\nBased on top ', nrow(top5k), ' ', opt$feature, 's expressed in Neurons'))
dev.off()

## Make some scatter plots of the methylation vs expr
delta_pval <- data.frame(delta = find_pval('nonCpG') - find_pval('CpG'), gene = names(find_pval('CpG')), nonCpG = find_pval('nonCpG'), CpG = find_pval('CpG'), stringsAsFactors = FALSE)
delta_pval$top5k <- delta_pval$gene %in% rownames(top5k)
delta_pval$snps <- with(subset(m_summary, type == 'nonCpG'), snps[match(delta_pval$gene, gene)])
delta_pval <- delta_pval[order(delta_pval$delta, decreasing = TRUE), ]

delta_pval$i_nonCpG <- mapply(function(g, s) {
    which(mres[['nonCpG']]$eqtls$gene == g & mres[['nonCpG']]$eqtls$snps == s)
}, delta_pval$gene, delta_pval$snps)
# delta_pval$i_CpG <- mapply(function(g, s) {
#     which(mres[['CpG']]$eqtls$gene == g & mres[['CpG']]$eqtls$snps == s)
# }, delta_pval$gene, delta_pval$snps)
# delta_pval$i_CpGmarg <- mapply(function(g, s) {
#     which(mres[['CpGmarg']]$eqtls$gene == g & mres[['CpGmarg']]$eqtls$snps == s)
# }, delta_pval$gene, delta_pval$snps)

save(delta_pval, file = paste0('rda/meqtl_delta_pval_', opt$feature,
    '_using_near.Rdata'))
head(delta_pval)

## Check the meth_n data (note that it's already been filtered to >= 4)
tapply(m_summary$meth_n, m_summary$type, table)
tapply(m_summary$meth_n, m_summary$type, function(x) {
    round(table(x) / length(x) * 100, 2)
})
tapply(m_summary$meth_n, m_summary$type, function(x) {
    cumsum(round(table(x) / length(x) * 100, 2))
})

## Check the expr_delta
tapply(m_summary$expr_delta, m_summary$type, summary)


find_i <- function(n = 10, t5k = FALSE, c_max = 200, nc_max = 200) {
    use <- delta_pval[sign(delta_pval$delta) == 1, ]
    if(t5k) {
        use <- use[use$top5k, ]
    }
    use <- use[use$CpG < c_max, ]
    use <- use[use$nonCpG < nc_max, ]
    head(use$i_nonCpG, n = n)
}


find_i_venn <- function(mtype = 'nonCpG', vennset = 'nonCpG', t5k = FALSE, n = 100) {
    
    use <- subset(m_summary, type == mtype)
    use <- use[order(use$FDR, decreasing = FALSE), ]
    if(t5k) {
        use <- use[use$gene %in% attr(vennres5k, 'intersections')[[vennset]], ]
    } else {
        use <- use[use$gene %in% attr(vennres, 'intersections')[[vennset]], ]
    }
    use <- head(use, n = n)

    ## Get the actual i
    mapply(function(g, s) {
        which(mres[[mtype]]$eqtls$gene == g & mres[[mtype]]$eqtls$snps == s)
    }, use$gene, use$snps)
}


ylab <- ifelse(opt$feature == 'psi', 'PSI', ifelse(opt$feature == 'jx', 'log2 (RP80M + 1)', 'log2 (RPKM + 1)'))

get_y <- function(type, i) {
    if(opt$feature == 'psi') {
        res <- variantFreq(mres[[type]]$expr[i, ])
    } else  {
        res <- log2(assays(mres[[type]]$expr)$norm[i, ] + 1)
    }
    return(res)
}

plotting_code <- function(i, type = 'nonCpG') {
    main <- paste(opt$feature, mres[[type]]$eqtls$gene[i], 'FDR', signif(mres[[type]]$eqtls$FDR[i], 3), '\n',  rowRanges(mres[[type]]$expr[i])$Symbol)
    
    plot(x = jitter(getMeth(mres[[type]]$meth[i, ], type = 'raw'), 0.05), y = jitter(get_y(type, i), 0.05), xlab = 'Methylation', ylab = ylab, main = main, sub = paste(as.vector(seqnames(rowRanges(mres[[type]]$meth)[i])), start(rowRanges(mres[[type]]$meth)[i]), as.vector(strand(rowRanges(mres[[type]]$meth)[i])), as.vector(rowRanges(mres[[type]]$meth)$c_context[i])))
}




## Make meth vs expr scatter plots for different sets
dir.create('pdf', showWarnings = FALSE)


pdf(paste0('pdf/meth_vs_expr_scatter_', 'nonCpG', '_', opt$feature, '_common_high_delta.pdf'))
## Just for checking which am I plotting later
plot(x = delta_pval$CpG, y = delta_pval$nonCpG, pch = 20, xlab = 'Best -log10 FDR p-value with CpGs', ylab = 'Best -log10 FDR p-value with non CpGs', main = paste(length(common), 'common', opt$feature, 'meQTLs'), col = ifelse(delta_pval$gene %in%delta_pval$gene[head(order(delta_pval$delta, decreasing = TRUE), n = 100)], 'red', 'black'))
abline(a = 0, b = 1, col = 'grey80')
for(i in find_i(100)) {
    plotting_code(i)
}
dev.off()

pdf(paste0('pdf/meth_vs_expr_scatter_', 'nonCpG', '_', opt$feature, '_common_high_delta_top5k.pdf'))
## Just for checking which am I plotting later
plot(x = delta_pval$CpG[delta_pval$top5k], y = delta_pval$nonCpG[delta_pval$top5k], pch = 20, xlab = 'Best -log10 FDR p-value with CpGs', ylab = 'Best -log10 FDR p-value with non CpGs', main = paste(sum(delta_pval$top5k), 'common', opt$feature, 'meQTLs\nBased on top ', nrow(top5k), ' ', opt$feature, 's expressed in Neurons'), col = ifelse(1:sum(delta_pval$top5k) %in% 1:100, 'red', 'black'))
abline(a = 0, b = 1, col = 'grey80')
for(i in find_i(100, t5k = TRUE)) {
    plotting_code(i)
}
dev.off()


## Could be adapted to do more venn intersections, check
# names(attr(vennres, 'intersections'))
pdf(paste0('pdf/meth_vs_expr_scatter_venn_', 'nonCpG', '_', opt$feature, '.pdf'))
for(i in find_i_venn()) {
    plotting_code(i)
}
dev.off()

pdf(paste0('pdf/meth_vs_expr_scatter_venn_', 'nonCpG', '_', opt$feature, '_top5k.pdf'))
for(i in find_i_venn(t5k = TRUE)) {
    plotting_code(i)
}
dev.off()


## Reproducibility info
proc.time()
message(Sys.time())
options(width = 120)
session_info()

