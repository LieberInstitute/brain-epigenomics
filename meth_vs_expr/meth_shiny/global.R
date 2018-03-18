library('shiny')
library('DT')
library('bsseq')
library('SGSeq')
library('devtools')
library('RColorBrewer')
library('shiny')

load('meth_df_withSymbol.Rdata')
load('meth_data.Rdata')
load('tf_data.Rdata')

get_age <- function(type, feature) {
    colData(meth_data[[feature]][[type]]$meth)$Age
}
get_meth <- function(type, i, feature) {
    as.vector(getMeth(meth_data[[feature]][[type]]$meth[i, ], type = 'raw'))
}

get_ylab <- function(feature) {
    ifelse(feature == 'psi', 'PSI',
           ifelse(grepl('jx', feature), 'log2 (RP80M + 1)','log2 (RPKM + 1)'))
}

get_y <- function(type, i, feature) {
    if(feature == 'psi') {
        res <- variantFreq(meth_data[[feature]][[type]]$expr[i, ])
    } else  {
        res <- log2(assays(meth_data[[feature]][[type]]$expr)$norm[i, ] + 1)
    }
    return(res)
}

## Function for getting age colors that match those used in
## bumphunting/plot_bumps_bsseqSmooth.R
get_col <- function(type = 'nonCpG', feature, ag = FALSE) {
    age <- colData(meth_data[[feature]][[type]]$meth)$Age
    age_group <- factor(ifelse(age < 0, 'Prenatal',
                               ifelse(age < 1, 'Infant',
                                      ifelse(age <= 12, 'Child',
                                             ifelse(age <= 17, 'Teen', 'Adult')))),
                        levels = c('Infant', 'Child', 'Teen', 'Adult', 'Prenatal'))

    age_group_cell <- factor(paste0(age_group, '_',
                                    colData(meth_data[[feature]][[type]]$meth)$Cell.Type),
                             levels = c(paste0(rep(levels(age_group)[1:4], each = 2),
                                               '_', c('Glia', 'Neuron')), 'Prenatal_H'))
    if(ag) return(age_group_cell)
    col <- c(brewer.pal(8, "Paired"), 'grey50')[c(5:6, 7:8, 3:4, 1:2, 9)][age_group_cell]
    return(col)
}

plotting_code <- function(i, type = 'nonCpG', feature) {
    if(length(i) > 1) {
        sapply(i, plotting_code, type = type, feature = feature)
        return(NULL)
    }
    main <- paste(feature,
                  ifelse(feature != 'exon', meth_data[[feature]][[type]]$ids$feature_id[i],
                         rowRanges(meth_data[[feature]][[type]]$expr[i])$exon_gencodeID),
                  'FDR', signif(meth_data[[feature]][[type]]$ids$FDR[i], 3), '\n',
                  ifelse(feature %in% c('gene', 'exon'), rowRanges(meth_data[[feature]][[type]]$expr[i])$Symbol, '')
    )

    plot(
        x = jitter(getMeth(meth_data[[feature]][[type]]$meth[i, ], type = 'raw'), 0.05),
        y = jitter(get_y(type, i, feature), 0.05), xlab = 'Methylation',
        ylab = get_ylab(feature), main = main,
        sub = paste(as.vector(seqnames(rowRanges(meth_data[[feature]][[type]]$meth)[i])),
                    start(rowRanges(meth_data[[feature]][[type]]$meth)[i]),
                    as.vector(strand(rowRanges(meth_data[[feature]][[type]]$meth)[i])),
                    as.vector(rowRanges(meth_data[[feature]][[type]]$meth)$c_context[i])),
        col = get_col(type, feature),
        pch = 21, bg = get_col(type, feature), cex = 2, cex.lab = 1.3, cex.main = 1.5)
}

summarize_meth <- function(d) {
    byfeat <- split(d, d$feature)
    byfeat <- byfeat[elementNROWS(byfeat) > 0]
    final <- do.call(rbind, lapply(names(byfeat), function(feature) {
        bybkg <- split(byfeat[[feature]], byfeat[[feature]]$meth_type)
        bybkg <- bybkg[elementNROWS(bybkg) > 0]
        do.call(rbind, lapply(names(bybkg), function(bkg) {
            bysign <- split(bybkg[[bkg]], sign(bybkg[[bkg]]$meth_statistic))
            names(bysign) <- c('-1' = 'down', '1' = 'up')[names(bysign)]
            bysign <- bysign[elementNROWS(bysign) > 0]
            do.call(rbind, lapply(names(bysign), function(s) {
                bys <- bysign[[s]]
                j <- as.integer(unlist(strsplit(bys$i, ',')))
                res <- data.frame(
                    n = length(j),
                    unique_n_c = length(unique(bys$c_id)),
                    unique_n_feature = length(unique(bys$feature_id)),
                    mean_n_samples_with_meth_not0 = mean(bys$n_samples_with_meth_not0),
                    mean_n_samples_with_meth_not1 = mean(bys$n_samples_with_meth_not1),
                    mean_expr_delta = mean(bys$expr_delta, na.rm = TRUE),
                    prop_age_affected = mean(!bys$meth_adjust_age_FDR_less5percent, na.rm = TRUE),
                    prop_gene_promoter = mean(bys$promoter_present),
                    prop_gene_body = mean(bys$body_present),
                    prop_gene_flanking = mean(bys$flanking_present),
                    feature = feature,
                    meth_type = bkg,
                    direction = s
                )
                context <- factor(rowRanges(meth_data[[feature]][[bkg]]$meth)$c_context[j], levels = c('CG', 'CHG', 'CHH'))

                res_cont <- table(context) / length(context)
                names(res_cont) <- paste0('prop_', names(res_cont))
                res <- cbind(res, as.data.frame(t(as.matrix(res_cont, dimnames = list(names(res_cont), NULL)))))

                gtype <- table(factor(bys$gene_type, levels = c('glia', 'neuron', 'none')))
                res_gtype <-  gtype / sum(gtype)
                names(res_gtype) <- paste0('prop_', names(res_gtype))
                res <- cbind(res, as.data.frame(t(as.matrix(res_gtype, dimnames = list(names(res_gtype), NULL)))))
                return(res)
            }))
        }))
    }))

    ## Make small numbers readable
    for(i in which(sapply(final, class) == 'numeric')) {
        final[, i] <- signif(final[, i], 3)
    }

    ## Reorder
    first_cols <- c('feature', 'meth_type', 'direction')
    final[, c(first_cols, colnames(final)[!colnames(final) %in% first_cols])]
}

make_track_info <- function(gr, name) {
    paste(as.character(seqnames(gr)), start(gr), end(gr), name, 999)
}

#i <- 1
#feature <- 'psi'
#type <- 'CpG'

ucsc_info <- function(i, feature, type) {
    if(is.null(i)) return('')

    expr <- rowRanges(meth_data[[feature]][[type]]$expr[i])
    expr_r <- if(feature == 'psi') range(expr[[1]]) else expr
    cbase <- rowRanges(meth_data[[feature]][[type]]$meth[i])

    textinfo <- c('track name=eqtl description="Methylation-feature pair" visibility=pack',
                  make_track_info(expr_r, ifelse(feature == 'psi', mcols(expr)$variantName, names(expr))),
                  make_track_info(cbase, paste0(type, start(cbase))))

    filename <- paste0(feature, '_', type, '_', i)
    writeLines(textinfo, con = file.path('www', filename))

    res <- paste0('https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=',
           as.character(seqnames(cbase)), ':',
           min(c(start(expr_r), start(cbase))), '-',
           max(c(end(expr_r), end(cbase))),
           '&hub_32407_assembledTx=hide&hubUrl=https://s3.amazonaws.com/LIBD_DLPFC/libdDLPFC/hub.txt&hgt.customText=',
           #file.path(getwd(), 'www', filename)
           filename
    )
    return(res)
}
