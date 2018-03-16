library('SGSeq')
library('bsseq')
library('getopt')
library('devtools')
library('limma')
library('gplots')
library('VennDiagram')
library('RColorBrewer')
library('clusterProfiler')
library('org.Hs.eg.db')
library('ggplot2')

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
    opt <- list('feature' = 'psi')
    opt <- list('feature' = 'jxleft')
    opt <- list('feature' = 'jxright')
    ## Not yet adapted/designed for the other two
    #opt <- list('feature' = 'jx', 'jxside' = 'left')
    #opt <- list('feature' = 'jx', 'jxside' = 'right')
}

stopifnot(opt$feature %in% c('psi', 'gene', 'exon', 'jxleft', 'jxright'))
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
    ## For jxright and jxleft
    type <- gsub('left|right', '', type)
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
        if(!file.exists(paste0('rda/expr_', type, '_unfiltered.Rdata'))) {
            dir.create('pdf', showWarnings = FALSE)
            pdf(paste0('pdf/suggested_expr_cutoffs_', tolower(type),
                '.pdf'), width = 12)
            cuts <- expression_cutoff(assays(expr)$norm, seed = 20180119)
            message(paste(cuts, collapse = ' '))
            cut <- max(cuts)
            dev.off()

            meanExpr <- rowMeans(assays(expr)$norm)
            rowRanges(expr)$meanExprs <- meanExpr
            rowRanges(expr)$passExprsCut <- meanExpr > cut
            dir.create('rda', showWarnings = FALSE)
            save(expr, file = paste0('rda/expr_', type, '_unfiltered.Rdata'))
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
    load(paste0('rda/me_annotated_FDR5_nonCpG_', gsub('left|right', '',
        opt$feature), '.Rdata'), verbose = TRUE)
    return(me_annotated)
}

if(!file.exists(paste0('rda/meqtl_mres_', opt$feature, '_using_near_meth11_proteincoding.Rdata'))) {
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
        f <- paste0('rda/me_annotated_FDR5_', cpg, '_', gsub('left|right', '',
            opt$feature), '.Rdata')
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

        if(grepl('jx', opt$feature)) {
            load('rda/me_annotated_FDR5_nonCpG_jx_only_rowExpr.Rdata', verbose = TRUE)
            ## Subset data to only only 1 side for the junctions
            expr_row <- resize(expr_row, width = 1,
                fix = ifelse(opt$feature == 'jxleft', 'start', 'end'),
                ignore.strand = TRUE)
            rowRanges(expr) <- resize(rowRanges(expr), width = 1,
                fix = ifelse(opt$feature == 'jxleft', 'start', 'end'),
                ignore.strand = TRUE)

            ## Take a window around the end of the jx, then reduce to reduce
            ## mem needed for this step
            expr_row <- reduce(resize(expr_row, width(expr_row) + 1000, fix = 'center'))
            me_ov <- findOverlaps(expr_row, expr)
            expr <- expr[sort(unique(subjectHits(me_ov))), ]
            rm(me_ov, expr_row)
        } else {
            meqtl <- load_meqtl()
            me_ov <- findOverlaps(resize(rowRanges(meqtl$expr), width(rowRanges(meqtl$expr)) + 2000, fix = 'center'), expr)
            expr <- expr[sort(unique(subjectHits(me_ov))), ]
            rm(me_ov, meqtl)
        }

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
        meth_all_n <- rowSums(getMeth(BSobj, type = 'raw') < 1)
        me_annotated$eqtls <- me_annotated$eqtls[meth_n >= 4 & meth_all_n >= 4, ]
        me_annotated$meth <- BSobj[meth_n >= 4 & meth_all_n >= 4, ]
        if(opt$feature == 'psi') {
            names(expr) <- paste0('row', seq_len(nrow(expr)))
            me_annotated$expr <- expr[snp_id(me_annotated$eqtls$gene), ]
        } else {
            me_annotated$expr <- expr[match(me_annotated$eqtls$gene, names(rowRanges(expr))), ]
        }

        me_annotated$eqtls$meth_n <- meth_n[meth_n >= 4 & meth_all_n >= 4]
        me_annotated$eqtls$meth_all_n <- meth_all_n[meth_n >= 4 & meth_all_n >= 4]
    } else {
        meth_n <- rowSums(getMeth(me_annotated$meth, type = 'raw') > 0)
        me_annotated$eqtls <- me_annotated$eqtls[meth_n >= 4 & meth_all_n >= 4, ]
        me_annotated$meth <- me_annotated$meth[meth_n >= 4 & meth_all_n >= 4, ]
        me_annotated$expr <- me_annotated$expr[meth_n >= 4 & meth_all_n >= 4, ]
        me_annotated$eqtls$meth_n <- meth_n[meth_n >= 4 & meth_all_n >= 4]
        me_annotated$eqtls$meth_all_n <- meth_n[meth_n >= 4 & meth_all_n >= 4]
    }

    message(paste(Sys.time(), 'computing the expression delta'))
    if(opt$feature == 'psi') {
        me_annotated$eqtls$expr_delta <- apply(variantFreq(me_annotated$expr), 1, function(x) diff(range(x)) )
    } else {
        me_annotated$eqtls$expr_delta <- apply(log2(assays(me_annotated$expr)$norm + 1), 1, function(x) diff(range(x)) )
    }


    return(me_annotated)
})
names(mres) <- cpgs

dir.create('rda', showWarnings = FALSE)
message(paste(Sys.time(), 'saving the mres object'))
save(mres, file = paste0('rda/meqtl_mres_', opt$feature,
    '_using_near.Rdata'))


## Simplify rowRanges() for the psi info just to make the rest of the
## code work
if(opt$feature == 'psi') {
    gene <- load_expr('gene')

    ## Use the gene-level info to simplify things
    for(cp in names(mres)) {
        newgr <- unlist(range(rowRanges(mres[[cp]]$expr)))
        ## Match using gene id (just first one in cases with more than 1)
        mgr <- match(sapply(mcols(rowRanges(mres[[cp]]$expr))$geneName, '[[', 1), rowRanges(gene)$gencodeID)
        mcols(newgr) <- mcols(gene)[mgr, ]
        names(newgr) <- names(gene)[mgr]
        mres[[cp]]$eqtls$gene <- names(newgr)

        #mres[[cp]]$eqtls$gene_name <- names(gene)[mgr]
        #names(newgr) <- mres[[cp]]$eqtls$gene
        rowRanges(mres[[cp]]$expr) <- newgr
        rm(newgr, mgr)
    }
    rm(cp, gene)
} else if (grepl('jx', opt$feature)) {
    gene <- load_expr('gene')

    ## add some gene-level info
    for(cp in names(mres)) {
        m <- match(rowRanges(mres[[cp]]$expr)$newGeneID,
            rowRanges(gene)$gencodeID)
        rowRanges(mres[[cp]]$expr)$gene_type <- rowRanges(gene)$gene_type[m]

        ## Special case for fusions
        fusion <- grep('-', rowRanges(mres[[cp]]$expr)$newGeneID)
        rowRanges(mres[[cp]]$expr)$gene_type[fusion] <- sapply(
            rowRanges(mres[[cp]]$expr)$newGeneID[fusion], function(x) {
                ifelse('protein_coding' %in%
                    rowRanges(gene)$gene_type[grep(gsub('-', '|', x),
                    rowRanges(gene)$gencodeID)],
                    'protein_coding', 'fusion')
        })
        rm(m, fusion)
    }
    rm(cp, gene)

}

## Filter further to keep only those with meth_n >= 11 & meth_all_n >= 11
## and protein coding
for(cp in names(mres)) {
    pc <- rowRanges(mres[[cp]]$expr)$gene_type == 'protein_coding'
    meth11 <- mres[[cp]]$eqtls$meth_n >= 11
    meth11_all <- mres[[cp]]$eqtls$meth_all_n >= 11
    message(paste(Sys.time(), 'Filtering by meth_n >= 11 and meth_all_n >= 11 and protein_coding genes for set', cp))
    pc_meth_tab <- addmargins(table('protein coding' = pc, 'meth_n >= 11 & meth_all_n >= 11' = meth11 & meth11_all, useNA = 'ifany'))
    print(pc_meth_tab)
    print(round(pc_meth_tab / max(pc_meth_tab) * 100, 2))
    mres[[cp]]$eqtls <- mres[[cp]]$eqtls[which(pc & meth11 & meth11_all), ]
    mres[[cp]]$expr <- mres[[cp]]$expr[which(pc & meth11 & meth11_all), ]
    mres[[cp]]$meth <- mres[[cp]]$meth[which(pc & meth11 & meth11_all), ]
    rm(pc, meth11, meth11_all, pc_meth_tab)
}
rm(cp)

save(mres, file = paste0('rda/meqtl_mres_', opt$feature,
    '_using_near_meth11_proteincoding.Rdata'))
} else {
    load(paste0('rda/meqtl_mres_', opt$feature, '_using_near_meth11_proteincoding.Rdata'), verbose = TRUE)
}

if(!file.exists(paste0('rda/meqtl_venn_', opt$feature, '_using_near.Rdata'))) {
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

} else if (opt$feature %in% c('psi', 'jxleft', 'jxright')) {
    ## Use the gene info here
    load(paste0('rda/meqtl_venn_gene_using_near.Rdata'), verbose = TRUE)
    rm(vennres, vennres5k, vennres5kglia, v, v5k, v5kglia)
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

## Keep only those with positive log FC (higher in Glia)
top5kglia <- top[sign(top$logFC) == 1, ]
nrow(top5kglia)
top5kglia <- top5kglia[order(top5kglia$adj.P.Val, decreasing = FALSE), ]
top5kglia <- head(top5kglia, 5000)

vinfo <- lapply(mres, function(me) { unique(as.character(me$eqtls$gene)) })
if(grepl('jx', opt$feature)) {
    vinfo5k <- lapply(names(vinfo), function(i) {
        me <- vinfo[[i]]
        me[rowRanges(mres[[i]]$expr)[me]$newGeneID %in% rownames(top5k)]
    })
    vinfo5kglia <- lapply(names(vinfo), function(i) {
        me <- vinfo[[i]]
        me[rowRanges(mres[[i]]$expr)[me]$newGeneID %in% rownames(top5kglia)]
    })
    names(vinfo5k) <- names(vinfo5kglia) <- names(vinfo)
} else {
    vinfo5k <- lapply(vinfo, function(me) { me[me %in% rownames(top5k)] })
    vinfo5kglia <- lapply(vinfo, function(me) { me[me %in% rownames(top5kglia)] })
}


dir.create('pdf', showWarnings = FALSE)
pdf(paste0('pdf/meqtl_venn_', opt$feature, '_using_near.pdf'))
vennres <- venn(vinfo) + title('meQTLs at FDR 5%, CpGs only in proximity to nonCpG')
vennres5k <- venn(vinfo5k) + title(paste0('meQTLs at FDR 5%, CpGs only in proximity to nonCpG\nBased on top ', nrow(top5k), ' ', ifelse(opt$feature %in% c('psi', 'jxleft', 'jxright'), 'gene', opt$feature), 's expressed in Neurons'))
vennres5kglia <- venn(vinfo5kglia) + title(paste0('meQTLs at FDR 5%, CpGs only in proximity to nonCpG\nBased on top ', nrow(top5kglia), ' ', ifelse(opt$feature %in% c('psi', 'jxleft', 'jxright'), 'gene', opt$feature), 's expressed in Glia'))
dev.off()

vennres <- venn(vinfo, show.plot = FALSE)
vennres5k <- venn(vinfo5k, show.plot = FALSE)
vennres5kglia <- venn(vinfo5kglia, show.plot = FALSE)

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
    main = paste0('meQTLs at FDR 5%, CpGs only in proximity to nonCpG\nBased on top ', nrow(top5k), ' ', ifelse(opt$feature %in% c('psi', 'jxleft', 'jxright'), 'gene', opt$feature), 's expressed in Neurons'),
    col = "transparent", fill = c("lightpink2","cornflowerblue", "olivedrab2"),
    alpha = 0.50, fontface = "bold",
    cat.col = c("palevioletred4", "darkblue", "olivedrab4"),
    margin=0.2)
grid.newpage()
grid.draw(v5k)


v5kglia <- venn.diagram(vinfo5kglia, filename = NULL,
    main = paste0('meQTLs at FDR 5%, CpGs only in proximity to nonCpG\nBased on top ', nrow(top5kglia), ' ', ifelse(opt$feature %in% c('psi', 'jxleft', 'jxright'), 'gene', opt$feature), 's expressed in Glia'),
    col = "transparent", fill = c("lightpink2","cornflowerblue", "olivedrab2"),
    alpha = 0.50, fontface = "bold",
    cat.col = c("palevioletred4", "darkblue", "olivedrab4"),
    margin=0.2)
grid.newpage()
grid.draw(v5kglia)

dev.off()


save(vennres, vennres5k, vennres5kglia, top, top5k, top5kglia, v, v5k, v5kglia,
    file = paste0('rda/meqtl_venn_', opt$feature, '_using_near.Rdata'))
} else {
    load(paste0('rda/meqtl_venn_', opt$feature, '_using_near.Rdata'), verbose = TRUE)
}

if(!file.exists(paste0('rda/meqtl_summary_', opt$feature, '_using_near.Rdata'))) {
message(paste(Sys.time(), 'summarizing the meQTL data by', opt$feature))
m_summary <- do.call(rbind, lapply(1:length(mres), function(i) {
    message(paste(Sys.time(), 'processing', names(mres)[i]))
    gdata <- split(mres[[i]]$eqtls, mres[[i]]$eqtls$gene)
    gdata <- gdata[elementNROWS(gdata) > 0]

    erow <- elementNROWS(gdata)
    print(table(erow > 1))
    message(paste(Sys.time(), 'creating typeres1'))
    g1 <- sapply(gdata[erow == 1], function(x) { x$gene })
    typeres1 <- mres[[i]]$eqtls[mres[[i]]$eqtls$gene %in% g1, ]
    typeres1$n_meqtls <- 1

    message(paste(Sys.time(), 'creating typeres2'))
    typeres2 <- do.call(rbind, lapply(gdata[erow > 1], function(g) {
        best <- which.min(g$FDR)
        res <- g[best, , drop = FALSE]
        res$n_meqtls <- nrow(g)
        return(res)
    }))

    message(paste(Sys.time(), 'creating typeres'))
    typeres <- rbind(typeres1, typeres2)
    typeres$type <- names(mres)[i]
    return(typeres)
}))

m_summary$gene <- as.character(m_summary$gene)
m_summary$snps <- as.character(m_summary$snps)

message(paste(Sys.time(), 'saving summary of the meQTL data by', opt$feature))
save(m_summary, file = paste0('rda/meqtl_summary_', opt$feature,
    '_using_near.Rdata'))
} else {
    load(paste0('rda/meqtl_summary_', opt$feature, '_using_near.Rdata'), verbose = TRUE)
}


common <- names(table(m_summary$gene)[table(m_summary$gene) == 3])
find_pval <- function(type, t5k = FALSE, t5kg = FALSE) {
    res <- -log10(m_summary$FDR[m_summary$gene %in% common & m_summary$type == type])
    ## Add the name
    names(res) <- m_summary$gene[m_summary$gene %in% common & m_summary$type == type]
    ## Match the order
    res <- res[match(common, names(res))]
    if(t5k) {
        res <- res[names(res) %in% vinfo5k[[type]]]
    } else if (t5kg) {
        res <- res[names(res) %in% vinfo5kglia[[type]]]
    }
    return(res)
}

stopifnot(identical(names(find_pval('CpG')), names(find_pval('nonCpG'))))

## Scatter plot of -log 10 p-values for the CpG and nonCpG data (FDR adjusted p-values)
pdf(paste0('pdf/scatter_FDR_', opt$feature, '_using_near.pdf'))
plot(x = find_pval('CpG'), y = find_pval('nonCpG'), pch = 20, xlab = 'Best -log10 FDR p-value with CpGs', ylab = 'Best -log10 FDR p-value with non CpGs', main = paste(length(common), 'common', opt$feature, 'meQTLs'))
plot(x = find_pval('CpG', TRUE), y = find_pval('nonCpG', TRUE), pch = 20, xlab = 'Best -log10 FDR p-value with CpGs', ylab = 'Best -log10 FDR p-value with non CpGs', main = paste0(sum(common %in% vinfo5k[['nonCpG']]), ' common ', opt$feature, ' meQTLs\nBased on top ', nrow(top5k), ' ', ifelse(opt$feature %in% c('psi', 'jxleft', 'jxright'), 'gene', opt$feature), 's expressed in Neurons'))
plot(x = find_pval('CpG', t5kg = TRUE), y = find_pval('nonCpG', t5kg = TRUE), pch = 20, xlab = 'Best -log10 FDR p-value with CpGs', ylab = 'Best -log10 FDR p-value with non CpGs', main = paste0(sum(common %in% vinfo5kglia[['nonCpG']]), ' common ', opt$feature, ' meQTLs\nBased on top ', nrow(top5kglia), ' ', ifelse(opt$feature %in% c('psi', 'jxleft', 'jxright'), 'gene', opt$feature), 's expressed in Glia'))

## Limiting the x and y axes
xylim <- c(0, round(max(find_pval('nonCpG')) + 0.5, 0))
plot(x = find_pval('CpG'), y = find_pval('nonCpG'), pch = 20, xlab = 'Best -log10 FDR p-value with CpGs', ylab = 'Best -log10 FDR p-value with non CpGs', main = paste(length(common), 'common', opt$feature, 'meQTLs'), xlim = xylim, ylim = xylim)
plot(x = find_pval('CpG', TRUE), y = find_pval('nonCpG', TRUE), pch = 20, xlab = 'Best -log10 FDR p-value with CpGs', ylab = 'Best -log10 FDR p-value with non CpGs', main = paste0(sum(common %in% vinfo5k[['nonCpG']]), ' common ', opt$feature, ' meQTLs\nBased on top ', nrow(top5k), ' ', ifelse(opt$feature %in% c('psi', 'jxleft', 'jxright'), 'gene', opt$feature), 's expressed in Neurons'), xlim = xylim, ylim = xylim)
plot(x = find_pval('CpG', t5kg = TRUE), y = find_pval('nonCpG', t5kg = TRUE), pch = 20, xlab = 'Best -log10 FDR p-value with CpGs', ylab = 'Best -log10 FDR p-value with non CpGs', main = paste0(sum(common %in% vinfo5kglia[['nonCpG']]), ' common ', opt$feature, ' meQTLs\nBased on top ', nrow(top5kglia), ' ', ifelse(opt$feature %in% c('psi', 'jxleft', 'jxright'), 'gene', opt$feature), 's expressed in Glia'), xlim = xylim, ylim = xylim)
dev.off()


## beta instead of FDR
## beta is the coef for change in expr explained by a 1 unit change in meth
find_beta <- function(type, t5k = FALSE, t5kg = FALSE) {
    res <- m_summary$beta[m_summary$gene %in% common & m_summary$type == type]
    ## Add the name
    names(res) <- m_summary$gene[m_summary$gene %in% common & m_summary$type == type]
    ## Match the order
    res <- res[match(common, names(res))]
    if(t5k) {
        res <- res[names(res) %in% vinfo5k[[type]]]
    } else if (t5kg) {
        res <- res[names(res) %in% vinfo5kglia[[type]]]
    }
    return(res)
}

pdf(paste0('pdf/scatter_beta_', opt$feature, '_using_near.pdf'))
plot(x = find_beta('CpG'), y = find_beta('nonCpG'), pch = 20, xlab = 'Beta from best FDR p-value with CpGs', ylab = 'Beta from best FDR p-value with non CpGs', main = paste(length(common), 'common', opt$feature, 'meQTLs'))
abline(v = 0, col = 'grey80')
abline(h = 0, col = 'grey80')
plot(x = find_beta('CpG', TRUE), y = find_beta('nonCpG', TRUE), pch = 20, xlab = 'Beta from best FDR p-value with CpGs', ylab = 'Beta from best FDR p-value with non CpGs', main = paste0(sum(common %in% vinfo5k[['nonCpG']]), ' common ', opt$feature, ' meQTLs\nBased on top ', nrow(top5k), ' ', ifelse(opt$feature %in% c('psi', 'jxleft', 'jxright'), 'gene', opt$feature), 's expressed in Neurons'))
abline(v = 0, col = 'grey80')
abline(h = 0, col = 'grey80')
plot(x = find_beta('CpG', t5kg = TRUE), y = find_beta('nonCpG', t5kg = TRUE), pch = 20, xlab = 'Beta from best FDR p-value with CpGs', ylab = 'Beta from best FDR p-value with non CpGs', main = paste0(sum(common %in% vinfo5kglia[['nonCpG']]), ' common ', opt$feature, ' meQTLs\nBased on top ', nrow(top5kglia), ' ', ifelse(opt$feature %in% c('psi', 'jxleft', 'jxright'), 'gene', opt$feature), 's expressed in Glia'))
abline(v = 0, col = 'grey80')
abline(h = 0, col = 'grey80')

## Limiting axes
xylim <- range(c(find_beta('CpG'), find_beta('nonCpG')))
plot(x = find_beta('CpG'), y = find_beta('nonCpG'), pch = 20, xlab = 'Beta from best FDR p-value with CpGs', ylab = 'Beta from best FDR p-value with non CpGs', main = paste(length(common), 'common', opt$feature, 'meQTLs'), xlim = xylim, ylim = xylim)
abline(v = 0, col = 'grey80')
abline(h = 0, col = 'grey80')
plot(x = find_beta('CpG', TRUE), y = find_beta('nonCpG', TRUE), pch = 20, xlab = 'Beta from best FDR p-value with CpGs', ylab = 'Beta from best FDR p-value with non CpGs', main = paste0(sum(common %in% vinfo5k[['nonCpG']]), ' common ', opt$feature, ' meQTLs\nBased on top ', nrow(top5k), ' ', ifelse(opt$feature %in% c('psi', 'jxleft', 'jxright'), 'gene', opt$feature), 's expressed in Neurons'), xlim = xylim, ylim = xylim)
abline(v = 0, col = 'grey80')
abline(h = 0, col = 'grey80')
plot(x = find_beta('CpG', t5kg = TRUE), y = find_beta('nonCpG', t5kg = TRUE), pch = 20, xlab = 'Beta from best FDR p-value with CpGs', ylab = 'Beta from best FDR p-value with non CpGs', main = paste0(sum(common %in% vinfo5kglia[['nonCpG']]), ' common ', opt$feature, ' meQTLs\nBased on top ', nrow(top5kglia), ' ', ifelse(opt$feature %in% c('psi', 'jxleft', 'jxright'), 'gene', opt$feature), 's expressed in Glia'), xlim = xylim, ylim = xylim)
abline(v = 0, col = 'grey80')
abline(h = 0, col = 'grey80')
dev.off()



## Make some scatter plots of the methylation vs expr
if(!file.exists(paste0('rda/meqtl_delta_pval_', opt$feature, '_using_near.Rdata'))) {
delta_pval <- data.frame(delta = find_pval('nonCpG') - find_pval('CpG'), gene = names(find_pval('CpG')), nonCpG = find_pval('nonCpG'), CpG = find_pval('CpG'), stringsAsFactors = FALSE)
delta_pval$top5k <- delta_pval$gene %in% vinfo5k[['nonCpG']]
delta_pval$top5kglia <- delta_pval$gene %in% vinfo5kglia[['nonCpG']]
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

save(delta_pval, common, file = paste0('rda/meqtl_delta_pval_', opt$feature,
    '_using_near.Rdata'))
} else {
    load(paste0('rda/meqtl_delta_pval_', opt$feature, '_using_near.Rdata'), verbose = TRUE)
}
head(delta_pval)

## Check the meth_n data (note that it's already been filtered to >= 4)
tapply(m_summary$meth_n, m_summary$type, table)
tapply(m_summary$meth_n, m_summary$type, function(x) {
    round(table(x) / length(x) * 100, 2)
})
tapply(m_summary$meth_n, m_summary$type, function(x) {
    cumsum(round(table(x) / length(x) * 100, 2))
})

tapply(m_summary$meth_all_n, m_summary$type, table)
tapply(m_summary$meth_all_n, m_summary$type, function(x) {
    round(table(x) / length(x) * 100, 2)
})
tapply(m_summary$meth_all_n, m_summary$type, function(x) {
    cumsum(round(table(x) / length(x) * 100, 2))
})

## Check the expr_delta
tapply(m_summary$expr_delta, m_summary$type, summary)


find_i <- function(n = 10, t5k = FALSE, c_max = 200, nc_max = 200, t5kg = FALSE) {
    use <- delta_pval[sign(delta_pval$delta) == 1, ]
    if(t5k) {
        use <- use[use$top5k, ]
    } else if (t5kg) {
        use <- use[use$top5kglia, ]
    }
    use <- use[use$CpG < c_max, ]
    use <- use[use$nonCpG < nc_max, ]
    head(use$i_nonCpG, n = n)
}


find_i_venn <- function(mtype = 'nonCpG', vennset = 'nonCpG', t5k = FALSE, n = 100, t5kg = FALSE) {

    use <- subset(m_summary, type == mtype)
    use <- use[order(use$FDR, decreasing = FALSE), ]
    if(t5k) {
        use <- use[use$gene %in% attr(vennres5k, 'intersections')[[vennset]], ]
    } else if (t5kg) {
        use <- use[use$gene %in% attr(vennres5kglia, 'intersections')[[vennset]], ]
    } else {
        use <- use[use$gene %in% attr(vennres, 'intersections')[[vennset]], ]
    }
    use <- head(use, n = n)

    ## Get the actual i
    mapply(function(g, s) {
        which(mres[[mtype]]$eqtls$gene == g & mres[[mtype]]$eqtls$snps == s)
    }, use$gene, use$snps)
}


ylab <- ifelse(opt$feature == 'psi', 'PSI', ifelse(grepl('jx', opt$feature), 'log2 (RP80M + 1)', 'log2 (RPKM + 1)'))

get_y <- function(type, i) {
    if(opt$feature == 'psi') {
        res <- variantFreq(mres[[type]]$expr[i, ])
    } else  {
        res <- log2(assays(mres[[type]]$expr)$norm[i, ] + 1)
    }
    return(res)
}

## Function for getting age colors that match those used in
## bumphunting/plot_bumps_bsseqSmooth.R
get_col <- function(type = 'nonCpG', ag = FALSE) {
    age_group <- factor(ifelse(colData(mres[[type]]$meth)$Age < 0, 'Prenatal',
        ifelse(colData(mres[[type]]$meth)$Age < 1, 'Infant',
        ifelse(colData(mres[[type]]$meth)$Age <= 12, 'Child',
        ifelse(colData(mres[[type]]$meth)$Age <= 17, 'Teen', 'Adult')))),
        levels = c('Infant', 'Child', 'Teen', 'Adult', 'Prenatal'))

    age_group_cell <- factor(paste0(age_group, '_', colData(mres[[type]]$meth)$Cell.Type),
        levels = c(paste0(rep(levels(age_group)[1:4], each = 2),
        '_', c('Glia', 'Neuron')), 'Prenatal_H'))
    if(ag) return(age_group_cell)
    col <- c(brewer.pal(8, "Paired"), 'grey50')[c(5:6, 7:8, 3:4, 1:2, 9)][age_group_cell]
    return(col)
}

if(opt$feature == 'gene') {
    pdf('pdf/meth_vs_expr_scatter_color_labels.pdf', width = 14)
    palette(c(brewer.pal(8, "Paired")[c(5:6, 7:8, 3:4, 1:2)], 'grey50'))
    plot(colData(mres[['nonCpG']]$meth)$Age, type = 'p', pch = 21, ylab = 'Age',
        bg = get_col(ag = TRUE), cex = 3, xlim = c(0, 30))
    legend("topright", levels(get_col(ag = TRUE)), pch = 15, col=1:9, cex=1.4)
    dev.off()
}


set.seed(20180207)
plotting_code <- function(i, type = 'nonCpG') {
    if(length(i) > 1) {
        sapply(i, plotting_code, type = type)
        return(NULL)
    }
    main <- paste(opt$feature, ifelse(opt$feature %in% c('gene', 'psi', 'jxleft', 'jxright'), mres[[type]]$eqtls$gene[i], rowRanges(mres[[type]]$expr[i])$exon_gencodeID), 'FDR', signif(mres[[type]]$eqtls$FDR[i], 3), '\n',  rowRanges(mres[[type]]$expr[i])$Symbol)

    plot(x = jitter(getMeth(mres[[type]]$meth[i, ], type = 'raw'), 0.05), y = jitter(get_y(type, i), 0.05), xlab = 'Methylation', ylab = ylab, main = main, sub = paste(as.vector(seqnames(rowRanges(mres[[type]]$meth)[i])), start(rowRanges(mres[[type]]$meth)[i]), as.vector(strand(rowRanges(mres[[type]]$meth)[i])), as.vector(rowRanges(mres[[type]]$meth)$c_context[i])), col = get_col(type))
}




## Make meth vs expr scatter plots for different sets
dir.create('pdf', showWarnings = FALSE)


pdf(paste0('pdf/meth_vs_expr_scatter_', 'nonCpG', '_', opt$feature, '_common_high_delta.pdf'))
## Just for checking which am I plotting later
plot(x = delta_pval$CpG, y = delta_pval$nonCpG, pch = 20, xlab = 'Best -log10 FDR p-value with CpGs', ylab = 'Best -log10 FDR p-value with non CpGs', main = paste(length(common), 'common', opt$feature, 'meQTLs'), col = ifelse(delta_pval$gene %in% delta_pval$gene[head(order(delta_pval$delta, decreasing = TRUE), n = 100)], 'red', 'black'))
abline(a = 0, b = 1, col = 'grey80')
for(i in find_i(100)) {
    plotting_code(i)
}
dev.off()

pdf(paste0('pdf/meth_vs_expr_scatter_', 'nonCpG', '_', opt$feature, '_common_high_delta_top5k.pdf'))
## Just for checking which am I plotting later
plot(x = delta_pval$CpG[delta_pval$top5k], y = delta_pval$nonCpG[delta_pval$top5k], pch = 20, xlab = 'Best -log10 FDR p-value with CpGs', ylab = 'Best -log10 FDR p-value with non CpGs', main = paste0(sum(delta_pval$top5k), ' common ', opt$feature, ' meQTLs\nBased on top ', nrow(top5k), ' ', opt$feature, 's expressed in Neurons'), col = ifelse(sign(delta_pval$delta[delta_pval$top5k]) == 1, 'red', 'black'))
abline(a = 0, b = 1, col = 'grey80')
for(i in find_i(100, t5k = TRUE)) {
    plotting_code(i)
}
dev.off()


pdf(paste0('pdf/meth_vs_expr_scatter_', 'nonCpG', '_', opt$feature, '_common_high_delta_top5kglia.pdf'))
## Just for checking which am I plotting later
plot(x = delta_pval$CpG[delta_pval$top5kglia], y = delta_pval$nonCpG[delta_pval$top5kglia], pch = 20, xlab = 'Best -log10 FDR p-value with CpGs', ylab = 'Best -log10 FDR p-value with non CpGs', main = paste0(sum(delta_pval$top5kglia), ' common ', opt$feature, ' meQTLs\nBased on top ', nrow(top5kglia), ' ', opt$feature, 's expressed in Glia'), col = ifelse(sign(delta_pval$delta[delta_pval$top5kglia]) == 1, 'red', 'black'))
abline(a = 0, b = 1, col = 'grey80')
for(i in find_i(100, t5kg = TRUE)) {
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

pdf(paste0('pdf/meth_vs_expr_scatter_venn_', 'nonCpG', '_', opt$feature, '_top5kglia.pdf'))
for(i in find_i_venn(t5kg = TRUE)) {
    plotting_code(i)
}
dev.off()



if(!file.exists(paste0('rda/meqtl_venn_go_', opt$feature, '_using_near.Rdata'))) {
## Gene ontology for each of the sets of the main venn diagram
if(opt$feature == 'psi') {
    expr <- load_expr('gene')
} else if (grepl('jx', opt$feature)) {
    expr <- load_expr('gene')
    uni <- unique(rowRanges(expr)$ensemblID[rowRanges(expr)$gene_type == 'protein_coding'])
    expr <- load_expr(opt$feature)
} else {
    expr <- load_expr(opt$feature)
}

if(!grepl('jx', opt$feature)) {
    uni <- unique(rowRanges(expr)$ensemblID[rowRanges(expr)$gene_type == 'protein_coding'])
}

length(uni)

v_symb <- lapply( attr(vennres, 'intersections'), function(vset) {
    if(grepl('jx', opt$feature)) {
        res <- unique(rowRanges(expr)$newGeneID[ names(rowRanges(expr)) %in% vset ])
        ## Make into Ensembl IDs, also deal with fusions
        res <- gsub('\\..*', '', unlist(strsplit(res, '-')))
    } else {
        res <- unique(rowRanges(expr)$ensemblID[ names(rowRanges(expr)) %in% vset ])
    }
    return(res)
})
sapply(v_symb, length)


go_venn_res <- lapply(names(v_symb), function(vname) {
    gores <- enrichGO(v_symb[[vname]], 'org.Hs.eg.db', universe = uni, ont = 'BP', minGSSize = 5, keyType = 'ENSEMBL', readable = TRUE)
    return(gores)
})
names(go_venn_res) <- names(v_symb)

pdf(paste0('pdf/meth_vs_expr_venn_GO_BP_', opt$feature, '.pdf'), width = 14)
lapply(names(v_symb), function(vname) {
    print(dotplot(go_venn_res[[vname]], title = paste('Venn set:', vname, 'with', length(v_symb[[vname]]), 'unique ENSEMBL IDs'), font.size = 18))
    return(NULL)
})
dev.off()


go_cluster_comp <- lapply(c('BP', 'MF', 'CC'), function(bp) {
    compareCluster(v_symb, fun = "enrichGO",
        universe = uni, OrgDb = 'org.Hs.eg.db',
        ont = bp, pAdjustMethod = "BH",
        pvalueCutoff  = 0.1, qvalueCutoff  = 0.05,
        readable= TRUE, keyType = 'ENSEMBL')
})
names(go_cluster_comp) <- c('BP', 'MF', 'CC')

pdf(paste0('pdf/meth_vs_expr_venn_GO_compare_clusters_', opt$feature, '.pdf'), width = 20)
lapply(names(go_cluster_comp), function(bp) {
    print(plot(go_cluster_comp[[bp]], title = paste('ontology:', bp), font.size = 18))
    return(NULL)
})
dev.off()


save(go_venn_res, go_cluster_comp, uni, v_symb, file = paste0('rda/meqtl_venn_go_', opt$feature, '_using_near.Rdata'))
} else {
    load(paste0('rda/meqtl_venn_go_', opt$feature, '_using_near.Rdata'), verbose = TRUE)
}



## Compute coefficient, stats and p-value for methylation vs age
get_age <- function(type) {
    colData(mres[[type]]$meth)$Age
}
get_meth <- function(type, i) {
    as.vector(getMeth(mres[[type]]$meth[i, ], type = 'raw'))
}

if(!file.exists(paste0('rda/meqtl_age_coef_', opt$feature, '_using_near.Rdata'))) {
age_coef <- lapply(names(mres), function(type) {
    message(paste(Sys.time(), 'processing', type))
    age <- get_age(type)
    res <- as.data.frame(t(apply(getMeth(mres[[type]]$meth, type = 'raw'), 1, function(row_i) {
        fit <- lm(age ~ row_i)
        summary(fit)$coef[2, ]
    })))
    res$ageFDR <- p.adjust(res[, 'Pr(>|t|)'], 'fdr')
    return(res)
})
names(age_coef) <- names(mres)
save(age_coef, file = paste0('rda/meqtl_age_coef_', opt$feature, '_using_near.Rdata'))
} else {
    load(paste0('rda/meqtl_age_coef_', opt$feature, '_using_near.Rdata'), verbose = TRUE)
}


beta_common <- do.call(rbind, lapply(seq_len(nrow(delta_pval)), function(i) {
    di <- delta_pval$i_nonCpG[[i]]
    data.frame(
        beta = mres[['nonCpG']]$eqtls$beta[di],
        t = mres[['nonCpG']]$eqtls$statistic[di],
        agebeta = age_coef[['nonCpG']]$Estimate[di],
        aget = age_coef[['nonCpG']][, 't value'][di],
        gtype = ifelse(delta_pval$top5k[i], 'neuron', ifelse(delta_pval$top5kglia[i], 'glia', 'none')),
        stringsAsFactors = FALSE
    )
}))

pdf(paste0('pdf/scatter_beta_vs_agebeta_nonCpGbackground_', opt$feature, '.pdf'), width = 15)
ggplot(beta_common, aes(x = beta, y = agebeta, colour = gtype)) + geom_point() + ylab('Age by methylation beta') + theme_grey(base_size = 18) + xlab('Expr by methylation beta') + scale_colour_discrete(name = 'Feature\ntype') + facet_grid(. ~ gtype) + ggtitle('nonCpGs common with CpGs: best meQTL by FDR per gene')
ggplot(beta_common, aes(x = t, y = aget, colour = gtype)) + geom_point() + ylab('Age by methylation statistic') + theme_grey(base_size = 18) + xlab('Expr by methylation statistic') + scale_colour_discrete(name = 'Feature\ntype') + facet_grid(. ~ gtype) + ggtitle('nonCpGs common with CpGs: best meQTL by FDR per gene')
dev.off()

## Check if methylation still explains the expression
## even after adjusting by age
if(!file.exists(paste0('rda/meqtl_agemeth_coef_', opt$feature, '_using_near.Rdata'))) {
agemeth_coef <- lapply(names(mres), function(type) {
    message(paste(Sys.time(), 'processing', type))
    age <- get_age(type)

    res <- as.data.frame(t(sapply(seq_len(nrow(mres[[type]]$meth)), function(i) {
        fit <- lm(as.vector(get_y(type, i)) ~ get_meth(type, i) + age)
        summary(fit)$coef[2, ]
    })))
    res$ageFDR <- p.adjust(res[, 'Pr(>|t|)'], 'fdr')
    return(res)
})
names(agemeth_coef) <- names(mres)
save(agemeth_coef, file = paste0('rda/meqtl_agemeth_coef_', opt$feature, '_using_near.Rdata'))
} else {
    load(paste0('rda/meqtl_agemeth_coef_', opt$feature, '_using_near.Rdata'), verbose = TRUE)
}

lapply(agemeth_coef, function(x) {
    addmargins(table('FDR < 0.05 after adjusting for age' = x$ageFDR < 0.05, useNA = 'ifany'))
})
lapply(agemeth_coef, function(x) {
    round(addmargins(table('FDR < 0.05 after adjusting for age' = x$ageFDR < 0.05, useNA = 'ifany')) / nrow(x) * 100, 2)
})


if(opt$feature == 'gene') {
    ## Adapted from /dcl01/lieber/ajaffe/lab/brain-epigenomics/sorted_nuclear_RNA/code/check_splicing.R
    expr <- load_expr(opt$feature)

    geneMap = rowRanges(expr)
    genePromoters = GRanges(seqnames(geneMap),
    IRanges(start = ifelse(strand(geneMap) == "+",
    start(geneMap)-2000, end(geneMap)-200),
    end = ifelse(strand(geneMap) == "+",
    start(geneMap)+200, end(geneMap)+2000)),
    strand = strand(geneMap))
    mcols(genePromoters) = mcols(geneMap)
    names(genePromoters) <- names(geneMap)

    geneBody = geneMap
    longIndex = which(width(geneBody) > 250)
    start(geneBody)[longIndex] = ifelse(strand(geneBody) == "+",
    start(geneBody)+200, start(geneBody))[longIndex]
    end(geneBody)[longIndex] = ifelse(strand(geneBody) == "+",
    end(geneBody), end(geneBody)-200)[longIndex]

    geneBodyLeft = geneBodyRight = geneBody
    start(geneBodyLeft) = start(geneBody)-100000
    end(geneBodyLeft) = start(geneBody)
    start(geneBodyRight) = end(geneBody)
    end(geneBodyRight) = end(geneBody) + 100000


    gene_section <- list(promoter = genePromoters, body = geneBody, flanking = c(geneBodyLeft, geneBodyRight))
    save(gene_section, file = 'rda/gene_section.Rdata')
}
load('rda/gene_section.Rdata', verbose = TRUE)

if(!file.exists(paste0('rda/meqtl_c_by_gene_', opt$feature, '_using_near.Rdata'))) {
c_by_gene <- lapply(names(mres), function(type) {
    message(paste(Sys.time(), 'processing', type))
    gr_c <- rowRanges(mres[[type]]$meth)

    typegen <- do.call(cbind, lapply(names(gene_section), function(gsec) {
        message(paste(Sys.time(), 'processing', gsec))
        gr_g <- gene_section[[gsec]]
        ov <- findOverlaps(gr_c, gr_g)
        ov_s <- IntegerList(split(subjectHits(ov), queryHits(ov)))
        gcon <- IntegerList(vector('list', length(gr_c)))
        gcon[as.integer(names(ov_s))] <- ov_s

        ov_gen <- CharacterList(lapply(ov_s, function(x) { gr_g$gencodeID[x] }))
        gencode <- CharacterList(vector('list', length(gr_c)))
        gencode[as.integer(names(ov_s))] <- ov_gen

        res <- DataFrame(geneid = gcon, present = elementNROWS(gcon) > 0, gencodeID = gencode)
        colnames(res) <- paste0(gsec, '_', colnames(res))
        return(res)
    }))
    return(typegen)

})
names(c_by_gene) <- names(mres)
save(c_by_gene, file = paste0('rda/meqtl_c_by_gene_', opt$feature, '_using_near.Rdata'))
} else {
    load(paste0('rda/meqtl_c_by_gene_', opt$feature, '_using_near.Rdata'), verbose = TRUE)
}

## Extract beta and age coef info for the venn groups
if(!file.exists(paste0('rda/meqtl_data_by_venn_', opt$feature, '_using_near.Rdata'))) {
data_by_venn <- do.call(rbind, lapply(c('nonCpG', 'CpGmarg'), function(typeref) {
    message(paste(Sys.time(), 'processing reference', typeref))
    which_v <- grep(typeref, names(attr(vennres, 'intersections')))
    res_ref <- mapply(function(vset, vname) {
        message(paste(Sys.time(), 'processing venn set', vname))
        iset <- which(mres[[typeref]]$eqtls$gene %in% vset)
        ## Put together eQTL info and age coef data
        res <- cbind(DataFrame(mres[[typeref]]$eqtls[iset, ]),
            age_coef[[typeref]][iset, ],
            DataFrame('agemethFDR' = agemeth_coef[[typeref]]$ageFDR[iset],
            'noage' = agemeth_coef[[typeref]]$ageFDR[iset] < 0.05),
            c_by_gene[[typeref]][iset, ]
        )
        ## Fix column names
        colnames(res)[(ncol(mres[[typeref]]$eqtls) + 1):(ncol(mres[[typeref]]$eqtls) + ncol(age_coef[[typeref]]))] <- colnames(age_coef[[typeref]])
        res$vset <- vname
        return(res)
    }, attr(vennres, 'intersections')[which_v], names(attr(vennres, 'intersections')[which_v]))
    res_ref <- do.call(rbind, res_ref)
    res_ref$typeref <- typeref
    return(res_ref)
}))
save(data_by_venn, file = paste0('rda/meqtl_data_by_venn_', opt$feature, '_using_near.Rdata'))
} else {
    load(paste0('rda/meqtl_data_by_venn_', opt$feature, '_using_near.Rdata'), verbose = TRUE)
}


with(data_by_venn, addmargins(table(noage, promoter_present, useNA = 'ifany')))
with(data_by_venn, addmargins(table(noage, body_present, useNA = 'ifany')))
with(data_by_venn, addmargins(table(noage, flanking_present, useNA = 'ifany')))

with(data_by_venn, round(addmargins(table(noage, promoter_present, useNA = 'ifany')) / length(noage) * 100, 2))
with(data_by_venn, round(addmargins(table(noage, body_present, useNA = 'ifany')) / length(noage) * 100, 2))
with(data_by_venn, round(addmargins(table(noage, flanking_present, useNA = 'ifany')) / length(noage) * 100, 2))

## Compare beta, then by venn groups

# From https://stackoverflow.com/questions/3483203/create-a-boxplot-in-r-that-labels-a-box-with-the-sample-size-n
give.n <- function(x){
    return(c(y = mean(x), label = length(x)))
}

venn_5k <- function(sub5k) {
    print(ggplot(sub5k, aes(y = beta, x = vset, fill = vset)) + geom_boxplot() + facet_grid(typeref ~ .) + scale_fill_discrete(name = 'Venn group') + xlab('Venn group') + theme_grey(base_size = 18) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ylab('Beta: expr by methylation') + stat_summary(fun.data = give.n, geom = "text", position = position_nudge(y = min(sub5k$beta))))

    print(ggplot(sub5k, aes(y = Estimate, x = vset, fill = vset)) + geom_boxplot() + facet_grid(typeref ~ .) + scale_fill_discrete(name = 'Venn group') + xlab('Venn group') + theme_grey(base_size = 18) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ylab('Beta: age by methylation') + stat_summary(fun.data = give.n, geom = "text", position = position_nudge(y = min(sub5k$Estimate))))

    print(ggplot(sub5k, aes(y = Estimate, x = beta)) + geom_density_2d() + facet_grid(typeref ~ vset) + xlab('Beta: expr by methylation') + theme_grey(base_size = 18) + ylab('Beta: age by methylation'))

    print(ggplot(sub5k, aes(y = Estimate, x = beta)) + geom_bin2d(bins = 100) + facet_grid(typeref ~ vset) + xlab('Beta: expr by methylation') + theme_grey(base_size = 18) + ylab('Beta: age by methylation'))

    #print(ggplot(sub5k, aes(y = Estimate, x = beta, colour = vset)) + geom_point() + facet_grid(typeref ~ vset) + scale_colour_discrete(name = 'Venn group') + xlab('Beta: expr by methylation') + theme_grey(base_size = 18) + ylab('Beta: age by methylation'))
    return(NULL)
}


pdf(paste0('pdf/meth_vs_expr_venn_beta_', opt$feature, '.pdf'), width = 14, height = 10)
venn_5k(as.data.frame(data_by_venn))
dev.off()

## For top in neurons, then in glia
pdf(paste0('pdf/meth_vs_expr_venn_beta_', opt$feature, '_top5k.pdf'), width = 14, height = 10)
venn_5k(as.data.frame(subset(data_by_venn, gene %in% unlist(vinfo5k))))
dev.off()

pdf(paste0('pdf/meth_vs_expr_venn_beta_', opt$feature, '_top5kglia.pdf'), width = 14, height = 10)
venn_5k(as.data.frame(subset(data_by_venn, gene %in% unlist(vinfo5kglia))))
dev.off()



## Summarize at the gene level, kind of like m_summary
get_summary <- function(var, fun) {
    fun(NumericList(split(var, grp)))
}

summarize_venn <- function(dbv) {
    grp <<- with(dbv, paste0(typeref, '_', vset, '_', gene))
    d_summ <- DataFrame(
        beta_mean = get_summary(dbv$beta, mean),
        beta_median = get_summary(dbv$beta, median),
        beta_prop_pos = get_summary(sign(dbv$beta) == 1, mean),
        beta_sign_mean = get_summary(sign(dbv$beta), mean),
        agebeta_mean = get_summary(dbv$Estimate, mean),
        agebeta_median = get_summary(dbv$Estimate, median),
        agebeta_prop_pos = get_summary(sign(dbv$Estimate) == 1, mean),
        agebeta_sign_mean = get_summary(sign(dbv$Estimate), mean),
        beta_neglog10FDR_mean = get_summary(-log10(dbv$FDR), mean),
        beta_neglog10FDR_median = get_summary(-log10(dbv$FDR), median),
        beta_neglog10pval_mean = get_summary(-log10(dbv$FDR), mean),
        beta_neglog10pval_median = get_summary(-log10(dbv$FDR), median),
        agebeta_neglog10FDR_mean = get_summary(-log10(dbv$ageFDR), mean),
        agebeta_neglog10FDR_median = get_summary(-log10(dbv$ageFDR), median),
        agebeta_neglog10pval_mean = get_summary(-log10(dbv[, 'Pr(>|t|)']), mean),
        agebeta_neglog10pval_median = get_summary(-log10(dbv[, 'Pr(>|t|)']), median),
        meth_n_mean = get_summary(dbv$meth_n, mean),
        meth_n_median = get_summary(dbv$meth_n, median),
        meth_all_n_mean = get_summary(dbv$meth_all_n, mean),
        meth_all_n_median = get_summary(dbv$meth_all_n, median),
        beta_t_mean = get_summary(dbv$statistic, mean),
        beta_t_median = get_summary(dbv$statistic, median),
        agebeta_t_mean = get_summary(dbv[, 't value'], mean),
        agebeta_t_median = get_summary(dbv[, 't value'], median),
        n_meqtls = get_summary(dbv$beta, elementNROWS),
        typeref = sapply(strsplit(unique(grp), '_'), '[[', 1),
        vset = sapply(strsplit(unique(grp), '_'), '[[', 2),
        gene = sapply(strsplit(unique(grp), '_'), function(x) { paste(x[3:length(x)], collapse = '_') })
    )
    d_summ$gtype <- ifelse(d_summ$gene %in% unlist(vinfo5k), 'neuron', ifelse(d_summ$gene %in% unlist(vinfo5kglia), 'glia', 'none'))
    ## https://en.wikipedia.org/wiki/Distance_from_a_point_to_a_line
    d_summ$distance_diag <- with(d_summ, abs( abs(beta_sign_mean) - abs(agebeta_sign_mean)) / sqrt(2))
    return(d_summ)
}

if(!file.exists(paste0('rda/meqtl_data_venn_summ_', opt$feature, '_using_near.Rdata'))) {
data_venn_summ <- summarize_venn(data_by_venn)
save(data_venn_summ, file = paste0('rda/meqtl_data_venn_summ_', opt$feature, '_using_near.Rdata'))
} else {
    load(paste0('rda/meqtl_data_venn_summ_', opt$feature, '_using_near.Rdata'), verbose = TRUE)
}

## Explore briefly
dim(data_venn_summ)
head(data_venn_summ)
summary(as.data.frame(data_venn_summ))
table(data_venn_summ$gtype)


plot_feature_level <- function(dvsumm) {
    d <- as.data.frame(dvsumm)
    print(ggplot(d, aes(y = agebeta_mean, x = beta_mean, colour = gtype)) + geom_point() + facet_grid(typeref ~ vset) + xlab('Expr by methylation mean beta') + theme_grey(base_size = 18) + ylab('Age by methylation mean beta') + scale_colour_discrete(name = 'Feature\ntype'))
    print(ggplot(d, aes(y = agebeta_median, x = beta_median, colour = gtype)) + geom_point() + facet_grid(typeref ~ vset) + xlab('Expr by methylation median beta') + theme_grey(base_size = 18) + ylab('Age by methylation median beta') + scale_colour_discrete(name = 'Feature\ntype'))

    print(ggplot(d, aes(y = agebeta_neglog10FDR_mean, x = beta_neglog10FDR_mean, colour = gtype)) + geom_point() + facet_grid(typeref ~ vset) + xlab('Expr by methylation mean neg log 10 FDR') + theme_grey(base_size = 18) + ylab('Age by methylation mean neg log 10 FDR') + scale_colour_discrete(name = 'Feature\ntype'))
    print(ggplot(d, aes(y = agebeta_neglog10FDR_median, x = beta_neglog10FDR_median, colour = gtype)) + geom_point() + facet_grid(typeref ~ vset) + xlab('Expr by methylation median neg log 10 FDR') + theme_grey(base_size = 18) + ylab('Age by methylation median neg log 10 FDR') + scale_colour_discrete(name = 'Feature\ntype'))
    print(ggplot(d, aes(y = agebeta_neglog10pval_mean, x = beta_neglog10pval_mean, colour = gtype)) + geom_point() + facet_grid(typeref ~ vset) + xlab('Expr by methylation mean neg log 10 pval') + theme_grey(base_size = 18) + ylab('Age by methylation mean neg log 10 pval') + scale_colour_discrete(name = 'Feature\ntype'))
    print(ggplot(d, aes(y = agebeta_t_mean, x = beta_t_mean, colour = gtype)) + geom_point() + facet_grid(typeref ~ vset) + xlab('Expr by methylation mean t') + theme_grey(base_size = 18) + ylab('Age by methylation mean t') + scale_colour_discrete(name = 'Feature\ntype'))
    print(ggplot(d, aes(y = agebeta_t_median, x = beta_t_median, colour = gtype)) + geom_point() + facet_grid(typeref ~ vset) + xlab('Expr by methylation median t') + theme_grey(base_size = 18) + ylab('Age by methylation median t') + scale_colour_discrete(name = 'Feature\ntype'))

    print(ggplot(d, aes(y = agebeta_prop_pos, x = beta_prop_pos, colour = gtype)) + geom_point() + facet_grid(typeref ~ vset) + xlab('Expr by methylation proportion positive') + theme_grey(base_size = 18) + ylab('Age by methylation proportion positive') + scale_colour_discrete(name = 'Feature\ntype'))
    print(ggplot(d, aes(y = agebeta_sign_mean, x = beta_sign_mean, colour = gtype)) + geom_point() + facet_grid(typeref ~ vset) + xlab('Expr by methylation sign mean') + theme_grey(base_size = 18) + ylab('Age by methylation sign mean') + scale_colour_discrete(name = 'Feature\ntype'))
    ## Lots of points on the diagonals in the previous 2 plots (at least at the gene level)
    ## Take absolute on both axes
    print(ggplot(d, aes(y = abs(agebeta_sign_mean), x = abs(beta_sign_mean), colour = gtype)) + geom_point() + facet_grid(typeref ~ vset) + xlab('Expr by methylation sign mean (abs)') + theme_grey(base_size = 18) + ylab('Age by methylation sign mean (abs)') + scale_colour_discrete(name = 'Feature\ntype'))
    print(ggplot(d, aes(y = distance_diag, x = n_meqtls, colour = gtype)) + geom_point() + facet_grid(typeref ~ vset) + ylab('Expr and age by methylation sign distance to diagonal') + theme_grey(base_size = 18) + xlab('Number of meQTLs') + scale_colour_discrete(name = 'Feature\ntype'))
    print(ggplot(d, aes(y = distance_diag, x = beta_mean, colour = gtype)) + geom_point() + facet_grid(typeref ~ vset) + ylab('Expr and age by methylation sign distance to diagonal') + theme_grey(base_size = 18) + xlab('Expr by methylation mean beta') + scale_colour_discrete(name = 'Feature\ntype'))
    print(ggplot(d, aes(y = distance_diag, x = beta_median, colour = gtype)) + geom_point() + facet_grid(typeref ~ vset) + ylab('Expr and age by methylation sign distance to diagonal') + theme_grey(base_size = 18) + xlab('Expr by methylation median beta') + scale_colour_discrete(name = 'Feature\ntype'))
    print(ggplot(d, aes(y = distance_diag, x = agebeta_mean, colour = gtype)) + geom_point() + facet_grid(typeref ~ vset) + ylab('Expr and age by methylation sign distance to diagonal') + theme_grey(base_size = 18) + xlab('Age by methylation mean beta') + scale_colour_discrete(name = 'Feature\ntype'))
    print(ggplot(d, aes(y = distance_diag, x = gtype, colour = gtype)) + geom_boxplot() + facet_grid(typeref ~ vset) + ylab('Expr and age by methylation sign distance to diagonal') + theme_grey(base_size = 18) + xlab('Feature type') + scale_colour_discrete(name = 'Feature\ntype'))

    print(ggplot(d, aes(x = distance_diag < 0.1, y = n_meqtls, colour = gtype)) + geom_boxplot() + facet_grid(typeref ~ vset) + xlab('Expr and age by methylation sign distance to diagonal < 0.1') + theme_grey(base_size = 18) + ylab('Number of meQTLs') + scale_colour_discrete(name = 'Feature\ntype') + stat_summary(fun.data = give.n, geom = "text", position = position_nudge(y = max(d$n_meqtls) * c(0.95, 1, 1.05), x = c(-0.2, 0, 0.2))))
    print(ggplot(d, aes(x = distance_diag < 0.1, y = n_meqtls, colour = gtype)) + geom_boxplot() + facet_grid(typeref ~ vset) + xlab('Expr and age by methylation sign distance to diagonal < 0.1') + theme_grey(base_size = 18) + ylab('Number of meQTLs') + scale_colour_discrete(name = 'Feature\ntype') + stat_summary(fun.data = give.n, geom = "text", position = position_nudge(y = log10(max(d$n_meqtls)) * c(0.65, 0.675, 0.7), x = c(-0.2, 0, 0.2))) + scale_y_log10())

    d$close <- d$distance_diag < 0.1
    print(ggplot(d, aes(y = agebeta_sign_mean, x = beta_sign_mean, colour = gtype)) + geom_point() + facet_grid(typeref + close ~ vset) + xlab('Expr by methylation sign mean') + theme_grey(base_size = 18) + ylab('Age by methylation sign mean') + scale_colour_discrete(name = 'Feature\ntype'))
    return(NULL)
}




## Explore at the feature level
pdf(paste0('pdf/meth_vs_expr_venn_byfeature_', opt$feature, '.pdf'), width = 14, height = 10)
plot_feature_level(data_venn_summ)
dev.off()


## By gene region and whether age affects it or not
pdf(paste0('pdf/meth_vs_expr_venn_byfeature_genebody_hasage_', opt$feature, '.pdf'), width = 14, height = 10)
plot_feature_level(summarize_venn(data_by_venn[intersect(which(data_by_venn$body_present), which(!data_by_venn$noage)), ]))
dev.off()
pdf(paste0('pdf/meth_vs_expr_venn_byfeature_genebody_noage_', opt$feature, '.pdf'), width = 14, height = 10)
plot_feature_level(summarize_venn(data_by_venn[intersect(which(data_by_venn$body_present), which(data_by_venn$noage)), ]))
dev.off()
pdf(paste0('pdf/meth_vs_expr_venn_byfeature_geneflanking_hasage_', opt$feature, '.pdf'), width = 14, height = 10)
plot_feature_level(summarize_venn(data_by_venn[intersect(which(data_by_venn$flanking_present), which(!data_by_venn$noage)), ]))
dev.off()
pdf(paste0('pdf/meth_vs_expr_venn_byfeature_geneflanking_noage_', opt$feature, '.pdf'), width = 14, height = 10)
plot_feature_level(summarize_venn(data_by_venn[intersect(which(data_by_venn$flanking_present), which(data_by_venn$noage)), ]))
dev.off()
pdf(paste0('pdf/meth_vs_expr_venn_byfeature_genepromoter_hasage_', opt$feature, '.pdf'), width = 14, height = 10)
plot_feature_level(summarize_venn(data_by_venn[intersect(which(data_by_venn$promoter_present), which(!data_by_venn$noage)), ]))
dev.off()
pdf(paste0('pdf/meth_vs_expr_venn_byfeature_genepromoter_noage_', opt$feature, '.pdf'), width = 14, height = 10)
plot_feature_level(summarize_venn(data_by_venn[intersect(which(data_by_venn$promoter_present), which(data_by_venn$noage)), ]))
dev.off()


## Reproducibility info
proc.time()
message(Sys.time())
options(width = 120)
session_info()



## For re-loading
saved_files <- c(paste0('rda/meqtl_mres_', opt$feature,
    '_using_near_meth11_proteincoding.Rdata'),
    paste0('rda/meqtl_venn_',opt$feature, '_using_near.Rdata'),
    paste0('rda/meqtl_summary_', opt$feature, '_using_near.Rdata'),
    paste0('rda/meqtl_delta_pval_', opt$feature, '_using_near.Rdata'),
    paste0('rda/meqtl_venn_go_', opt$feature, '_using_near.Rdata'),
    paste0('rda/meqtl_age_coef_', opt$feature, '_using_near.Rdata'),
    paste0('rda/meqtl_data_by_venn_', opt$feature, '_using_near.Rdata'),
    paste0('rda/meqtl_data_venn_summ_', opt$feature, '_using_near.Rdata'),
    'rda/gene_section.Rdata',
    paste0('rda/meqtl_c_by_gene_', opt$feature, '_using_near.Rdata')
)
for(i in saved_files) load(i, verbose = TRUE)
rm(i)

