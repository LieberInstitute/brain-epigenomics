meth_display <- meth_df[ c('i', 'feature', 'meth_type', 'feature_id', 'symbol', 'gene_type', 'meth_coefficient', 'meth_statistic', 'meth_FDR', 'n_samples_with_meth_not0', 'n_samples_with_meth_not1', 'expr_delta', 'age_coefficient', 'age_statistic', 'age_FDR', 'meth_adjust_age_FDR', 'meth_adjust_age_FDR_less5percent', 'venn_set', 'promoter_present', 'body_present', 'flanking_present', 'promoter_gencodeID',  'body_gencodeID',  'flanking_gencodeID')]

for(i in which(colnames(meth_display) %in% c('feature', 'meth_type', 'gene_type', 'venn_set'))) {
    meth_display[, i] <- as.factor(meth_display[, i])
}

for(i in which(colnames(meth_display) %in% c('meth_coefficient', 'meth_statistic', 'meth_FDR', 'age_coefficient', 'age_statistic', 'age_FDR', 'meth_adjust_age_FDR', 'expr_delta'))) {
    meth_display[, i] <- signif(meth_display[, i], 3)
}


shinyServer(function(input, output, session) {
    output$meth_summary <- DT::renderDataTable(
        meth_display,
        style = 'bootstrap', rownames = FALSE, filter = 'top',
        options = list(
            columnDefs = list(
                list(className = 'dt-center', targets = 1)
            ),
            pageLength = 10,
            lengthMenu = c(5, 10, 25, 50, 100),
            order = list(list(8, 'asc'))
        ),
        selection = 'single'
    )

    output$meth_summarized <- DT::renderDataTable(
        summarize_meth(meth_df[input$meth_summary_rows_all, ]),
        style = 'bootstrap', rownames = FALSE, filter = 'top',
        options = list(
            columnDefs = list(
                list(className = 'dt-center', targets = 1)
            ),
            pageLength = 25,
            lengthMenu = c(25, 50, 100)
        )
    )


    observeEvent(input$meth_type, {
        if(length(input$meth_summary_rows_selected) == 0) output$meth_selector <- renderUI({
            tagList(numericInput('meth_id', label = 'Result (i column in first table)', value = 1,
                min = 1, max = nrow(meth_data[[input$meth_feature]][[input$meth_type]]$ids), step = 1),
            helpText(paste('Max i is', nrow(meth_data[[input$meth_feature]][[input$meth_type]]$ids)))
        )})
    })

    observeEvent(input$meth_feature, {
        if(length(input$meth_summary_rows_selected) == 0) output$meth_selector <- renderUI({
            tagList(numericInput('meth_id', label = 'Result (i column in first table)', value = 1,
                min = 1, max = nrow(meth_data[[input$meth_feature]][[input$meth_type]]$ids), step = 1),
            helpText(paste('Max i is', nrow(meth_data[[input$meth_feature]][[input$meth_type]]$ids)))
        )})
    })

    observeEvent(input$meth_summary_rows_selected, {
        feat <- meth_df$feature[input$meth_summary_rows_selected]
        type <- meth_df$meth_type[input$meth_summary_rows_selected]

        updateSelectInput(session, 'meth_feature', selected = feat)
        updateSelectInput(session, 'meth_type', selected = type)
        output$meth_selector <- renderUI(
            if(length(input$meth_summary_rows_selected) > 0 ) {
                 tagList(numericInput('meth_id', label = 'Result (i column in first table)',
                value = as.integer(strsplit(meth_df$i[input$meth_summary_rows_selected], ',')[[1]][1]),
                min = 1, max = nrow(meth_data[[feat]][[type]]$ids), step = 1),
            helpText(paste('Max i is', nrow(meth_data[[feat]][[type]]$ids))))
            } else {
                tagList(numericInput('meth_id', label = 'Result (i column in first table)', value = 1,
                min = 1, max = nrow(meth_data[[input$meth_feature]][[input$meth_type]]$ids), step = 1),
            helpText(paste('Max i is', nrow(meth_data[[input$meth_feature]][[input$meth_type]]$ids))))
            }
        )
    })

    output$meth_plot <- renderPlot({
        set.seed(20180315)
        if(is.null(input$meth_id)) {
            return(plot(1))
        }
        plotting_code(input$meth_id, input$meth_type, input$meth_feature)
    })

    output$meth_plot_download <-  downloadHandler(
        filename =  function() {
            paste0('LIBD_methylation_explorer_plot_', input$meth_id, '_',
                   input$meth_type, '_', input$meth_feature, '_', Sys.time(), '.pdf')
        },
        content = function(file) {
            pdf(file, useDingbats = FALSE)
            if(is.null(input$meth_id)) {
                return(plot(1))
            }
            plotting_code(input$meth_id, input$meth_type, input$meth_feature)
            dev.off()
        }
    )


    output$expr_info <- renderPrint({
        gr <- rowRanges(meth_data[[input$meth_feature]][[input$meth_type]]$expr[input$meth_id])
        if(input$meth_feature == 'psi') {
            print('----- splicing event information -----')
            print(mcols(gr))
            print('----- ------ ----- ------ ------ -----')
        }
        print(gr)
    }, width = 120)

    output$expr_info_download <-  downloadHandler(
        filename =  function() {
            paste0('LIBD_methylation_expr_info_', input$meth_id, '_',
                   input$meth_type, '_', input$meth_feature, '_', Sys.time(), '.Rdata')
        },
        content = function(file) {
            expr_gr <- rowRanges(meth_data[[input$meth_feature]][[input$meth_type]]$expr[input$meth_id])
            save(expr_gr, file = file)
        }
    )

    output$meth_info <- renderPrint({
        print(rowRanges(meth_data[[input$meth_feature]][[input$meth_type]]$meth[input$meth_id]))
    }, width = 120)

    output$meth_info_download <-  downloadHandler(
        filename =  function() {
            paste0('LIBD_methylation_meth_info_', input$meth_id, '_',
                   input$meth_type, '_', input$meth_feature, '_', Sys.time(), '.Rdata')
        },
        content = function(file) {
            meth_gr <- rowRanges(meth_data[[input$meth_feature]][[input$meth_type]]$meth[input$meth_id])
            save(meth_gr, file = file)
        }
    )

    observeEvent(input$meth_feature, {
        output$tf_info <- renderUI(if(input$meth_feature == 'exon')
            tagList(
                h3('TF information: cytosine level'),
                tableOutput('tf_cbase'),
                h3('TF information: exon level'),
                tableOutput('tf_exon')
            )
        )
    })

    output$tf_cbase <- renderTable({
        if(is.null(input$meth_id)) return('')
        res <- tf_data$cbase[[input$meth_type]]
        j <- which(any(input$meth_id == res$i))
        if(length(j) == 0) return('')
        res <- res[j[1:2]]
        res$i <- paste(res$i, collapse = ',')
        as.data.frame(mcols(res))
    }, digits = -3)

    output$tf_exon <- renderTable({
        if(is.null(input$meth_id)) return(NULL)
        res <- subset(tf_data$exon[[input$meth_type]], i == input$meth_id)
        if(length(res) == 0) return(NULL)
        rownames(res) <- NULL
        res$i <- as.character(res$i)
        return(res[1:2, c('motif', 'raw_score', 'FDR', 'feature_id', 'i')])
    }, digits = -3)



    output$ucsc_link <- renderUI({
      #  print(ucsc_info(input$meth_id, input$meth_feature, input$meth_type))
        tagList(
            HTML(paste0('<a href="',
                    ucsc_info(input$meth_id, input$meth_feature, input$meth_type),
                    '">Export to UCSC Genome Browser</a>'))
    )})

    output$color_guide <- renderPlot({
        palette(c(brewer.pal(8, "Paired")[c(5:6, 7:8, 3:4, 1:2)], 'grey50'))
        plot(colData(meth_data[['gene']][['nonCpG']]$meth)$Age, type = 'p', pch = 21, ylab = 'Age',
             bg = get_col(ag = TRUE, feature = 'gene'), cex = 3, xlim = c(0, 30),
             xlab = 'Sample ID', cex.lab = 1.5)
        legend("topright", levels(get_col(ag = TRUE, feature = 'gene')), pch = 15, col=1:9, cex=1.4)
    })

    output$download_data <- downloadHandler(
        filename = function() { paste0('LIBD_methylation_explorer_selection_', Sys.time(), '.csv') },
        content = function(file) {
            current <- meth_display[input$meth_summary_rows_all, ]
            write.csv(current, file, row.names = FALSE)
        }
    )

    output$download_summary <- downloadHandler(
        filename = function() { paste0('LIBD_methylation_explorer_summary_selection_', Sys.time(), '.csv') },
        content = function(file) {
            current <- summarize_meth(meth_df[input$meth_summary_rows_all, ])
            write.csv(current, file, row.names = FALSE)
        }
    )

    ## Reproducibility info
    output$session_info <- renderPrint(session_info(), width = 120)
})
