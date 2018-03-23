#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

# Define UI for application that draws a histogram
shinyUI(navbarPage(

    title = 'LIBD WGBS Expression explorer',
    tabPanel('home',
        tags$head(
            includeScript("google-analytics.js")
        ),
        p('Welcome to the ', strong('LIBD WGBS Expression explorer!'), 'Here you can explore the methylation and expression associations described by Price et al, 2018. For more details check the documentation tab or simply jump to the next tab.'),
        tags$hr(),
        h3('Main publication'),
        tags$ul(
            tags$li(HTML('<b>Price AJ</b>, Collado-Torres L, Ivanov NA, Xia W, Burke EE, Shin JH, Tao R, Ma L, Jia Y,  Hyde TM, Kleinman JE, Weinberger DR, Jaffe AE. <a href="http://www.sciencemag.org/" onclick="ga(\'send\', \'event\', \'click\', \'link\', \'paper\', 1)">Paper title TBD</a>. <i>paper journal</i>, 2018. doi: todo.'))
        )
    ),
    tabPanel('methylation data',
        p('The following table includes all the CpG and nonCpG associations with nearby features at FDR 5%. The features levels are gene, exon, splicing events that affect percent splicing in (PSI).'),
        helpText('Any filters applied here affect the results of the methylation summary tab.'),
        DT::dataTableOutput('meth_summary'),
        downloadLink('download_data', 'Download methylation data'),
        tags$hr(),
        h1('Association details'),
        p('Here you can explore in further detail a particular methylation and expression association. Choose the feature type, the methylation type (CpG or nonCpG) and select which association id (from column', strong('i'), 'in previous table) you want to explore.'),
        helpText('You can also select a row in the table above and the options will be chosen for you automatically. The i chosen will be the first one in the row.'),
        fillRow(
            selectInput('meth_feature', label = 'Feature', c('gene', 'exon', 'psi'), selected = 'gene'),
            selectInput('meth_type', label = 'Methylation type', c('nonCpG', 'CpG'), selected = 'nonCpG'),
            uiOutput('meth_selector'),
            height = '100px'
        ),
        uiOutput('ucsc_link'),
        tags$hr(),
        h3('Methylation vs expression association plot:'),
        p('Scatter plot comparing methylation vs expression. Samples are colored by age group:'),
        helpText('Red: infant, orange: child, green: teen, blue: adult. Details in documentation tab.'),
        plotOutput('meth_plot', width = '500px'),
        helpText('A small amount of jitter has been added to the data to minimize overplotting.'),
        downloadLink('meth_plot_download', 'Save image to a PDF file'),
        h3('Expression feature information'),
        verbatimTextOutput('expr_info', placeholder = TRUE),
        downloadLink('expr_info_download', 'Download to an Rdata file'),
        h3('Cytosine information'),
        verbatimTextOutput('meth_info', placeholder = TRUE),
        downloadLink('meth_info_download', 'Download to an Rdata file'),
        uiOutput('tf_info'),
        helpText('For more information check the documentation tab.'),
        bookmarkButton(id = 'bookmark1')
    ),
    tabPanel('methylation summary',
        p('This table summarizes the methylation and expression associations that you have filtered in the previous tab.'),
        helpText('Access the previous tab at least once for this table to display. It can take a few seconds for the table to be generated.'),
        DT::dataTableOutput('meth_summarized'),
        downloadLink('download_summary', 'Download methylation summary')
    ),
    tabPanel('documentation',
        h3('Variables in methylation data tab'),
        tags$ul(
            tags$li(p(strong('i'), "id number. Also corresponds to the rank of the association among those from the same feature and methylation type. Could have more than one entry separated by commas as it's typical for PSI events.")),
            tags$li(p(strong('feature'), 'Either gene, exon or event affect the percent spliced in (PSI).')),
            tags$li(p(strong('meth_type'), 'Methylation type. Either CpG or non-CpG')),
            tags$li(p(strong('feature_id'), 'Feature ID. Either GENCODE id or internal exon id.')),
            tags$li(p(strong('symbol'), 'Gene symbol.')),
            tags$li(p(strong('gene_type'), 'Is the feature among those differentially expresed between neurons and glia? Either higher expression in neuron, glia or none of them. For gene and PSI we used gene-level data, for exons we used the list of DE exons. Only up to the top 5,000 DEgenes/exons were considered for each cell type.')),
            tags$li(p(strong('meth_coefficient'), 'Linear regression coefficient for methylation explaining expression changes. Expression is either the log2(normalized expression + 1) or PSI depending on the feature type.')),
            tags$li(p(strong('meth_statistics'), 't-statistic for meth_coefficient.')),
            tags$li(p(strong('meth_FDR'), 'FDR adjusted p-value for meth_coefficient.')),
            tags$li(p(strong('n_samples_with_meth_not0'), 'Number of samples that have non-zero methylation values. Total: 22. Results displayed have been filtered such that the minimum is 11.')),
            tags$li(p(strong('n_samples_with_meth_not1'), 'Number of samples that have non-one methylation values. Results displayed have been filtered such that the minimum is 11.')),
            tags$li(p(strong('expr_delta'), 'Change in expression: max minus min.')),
            tags$li(p(strong('age_coefficient'), 'Linear regression coefficient for age explaining expression changes.')),
            tags$li(p(strong('age_statistics'), 't-statistic for age_coefficient.')),
            tags$li(p(strong('age_FDR'), 'FDR adjusted p-value for age_coefficient.')),
            tags$li(p(strong('meth_adjust_age_FDR'), 'FDR adjuste p-value for meth_coefficient in a multiple linear regression that also adjusts for age.')),
            tags$li(p(strong('meth_adjust_age_FDR_less5percent'), 'Whether meth_adjust_age_FDR is less than 5 percent. If FALSE, then age is confounding the relationship between methylation and expression.')),
            tags$li(p(strong('venn_set'), 'Venn set identifier for the methylation and expression association. The 3 sets are nonCpG FDR <5%, CpG FDR <5% and CpG marg for CpG associations near the nonCpG locations that are marginally significant (not included in the table).')),
            tags$li(p(strong('promoter_present'), 'Is the cytosine present in the promoter of a gene? 2000 to 200 base pairs upstream of a gene.')),
            tags$li(p(strong('body_present'), 'Is the cytosine present in the gene body? For genes > 250 bp, 200 bp upstream to 200 bp downstream of the gene.')),
            tags$li(p(strong('flanking_present'), 'Is the cytosine present in the region flanking a gene? 10kbp upstream and downstream of the gene body as defined above.')),
            tags$li(p(strong('promoter_gencodeID'), 'GENCODE gene ids the cytosine is overlapping. IDs and annotation are from GENCODE v25 on hg19.')),
            tags$li(p(strong('body_gencodeID'), 'GENCODE gene ids the cytosine is overlapping.')),
            tags$li(p(strong('flanking_gencodeID'), 'GENCODE gene ids the cytosine is overlapping.'))
        ),
        tags$hr(),
        h3('Variables in methylation summary tab'),
        p('The table in the methylation summary tab is automatically calculated depending on the filters applied to the methylation data table.'),
        tags$ul(
            tags$li(p(strong('feature'), 'Either gene, exon or event affect the percent spliced in (PSI).')),
            tags$li(p(strong('meth_type'), 'Methylation type. Either CpG or non-CpG')),
            tags$li(p(strong('direction'), 'Does increase in methylation increase expression? Yes: up; no: down.')),
            tags$li(p(strong('n'), 'Number of associations at FDR 5%.')),
            tags$li(p(strong('unique_n_c'), 'Number of unique cytosines in these associations.')),
            tags$li(p(strong('unique_n_feature'), 'Number of unique features in these associations.')),
            tags$li(p(strong('mean_n_samples_with_meth_not0'), 'Mean number of samples that have non-zero methylation values. Total: 22.')),
            tags$li(p(strong('n_samples_with_meth_not1'), 'Mean number of samples that have non-one methylation values.')),
            tags$li(p(strong('mean_expr_delta'), 'Mean change in expression (ignoring NAs).')),
            tags$li(p(strong('prop_age_affected'), 'Proportion of associations that are confounded by age, that is, proportion of associations where meth_adjust_age_FDR_less5percent is FALSE.')),
            tags$li(p(strong('prop_gene_promoter'), 'Proportion of cytosines in these associations that overlap gene promoters. A single cytosone can overlap multiple gene regions for different genes.')),
            tags$li(p(strong('prop_gene_body'), 'Proportion of cytosines in these associations that overlap gene bodies.')),
            tags$li(p(strong('prop_gene_flanking'), 'Proportion of cytosines in these associations that overlap gene flanking regions.')),
            tags$li(p(strong('prop_CG'), 'Proportion of cytosines in these associations that have a CG context.')),
            tags$li(p(strong('prop_CHG'), 'Proportion of cytosines in these associations that have a CHG context.')),
            tags$li(p(strong('prop_CHH'), 'Proportion of cytosines in these associations that have a CHH context.')),
            tags$li(p(strong('prop_glia'), 'Proportion of features in these associations where gene_type is glia (differentially expressed between neurons and glia, with higher expression in glia).')),
            tags$li(p(strong('prop_neuron'), 'Proportion of features in these associations where gene_type is neuron.')),
            tags$li(p(strong('prop_none'), 'Proportion of features in these associations where gene_type is none.'))
        ),
        tags$hr(),
        h3('Information in association details'),
        tags$ul(
            tags$li(p(strong('Expression feature information'), HTML('<a href="http://bioconductor.org/packages/GenomicRanges" onclick="ga(\'send\', \'event\', \'click\', \'link\', \'GenomicRanges\', 1)">GenomicRanges</a> information for the corresponding feature. If the feature type is PSI, it also includes the <a href="http://bioconductor.org/packages/SGSeq" onclick="ga(\'send\', \'event\', \'click\', \'link\', \'GenomicRanges\', 1)">SGSeq</a> data for the PSI event. Uses hg19 coordinates and GENCODE v25 annotation.') )),
            tags$li(p(strong('Cytosine information'), HTML('<a href="http://bioconductor.org/packages/GenomicRanges" onclick="ga(\'send\', \'event\', \'click\', \'link\', \'GenomicRanges\', 1)">GenomicRanges</a> information for the corresponding cytosine. Uses hg19 coordinates.'))),
            tags$li(p(strong('TF information: cytosine level'), '(only for exons) transcription factor motif enrichment in a +- 15 base-pair window around the cytosine. The direction of the methylation association (defined previously), whether the cytosine is inside the exon, the raw score of the TF motif, the FDR adjusted p-value for the TF motif, and the i (defined above) where this cytosine is involved.', HTML('raw-score is as defined by motifEnrichmenth() from <a href="http://bioconductor.org/packages/PWMEnrich" onclick="ga(\'send\', \'event\', \'click\', \'link\', \'PWMEnrich\', 1)">PWMEnrich</a>.'))),
            tags$li(p(strong('TF information: exon level'), '(only for exons) transcription factor motif enrichment inside the exon, similar to the previous section. Exons with a length less than 30 base-pairs were not considered.'))
        ),
        tags$hr(),
        h3('Age colors'),
        p('Below you can see the color palette used for the different age groups. Age cutoffs: prenatal < 0, infant < 1, child <= 12, teen <= 17, adult > 17.'),
        plotOutput('color_guide', width = '800px'),
        tags$hr(),
        h3('Support'),
        p('For issues with this shiny application please get in touch with Leonardo Collado-Torres and Andrew Jaffe at ', HTML('<a href="https://www.libd.org/">libd.org</a>.')),
        p('The following information will be useful to them:'),
        verbatimTextOutput('session_info'),
        p('Also try to include a small reproducible example so they can figure out what went wrong. Thank you!')
    ),
    tags$hr(),
    h3('Data license'),
    p('The data in ', strong('LIBD WGBS Expression explorer'), HTML(' is licensed under CC BY 4.0. The legal text can be found <a href="LICENSE.txt">here</a>.')),
    tags$hr(),
    h3('Acknowledgements'),
    p('This research was supported by NIH R21MH102791-01A1 and the Lieber Institute for Brain Development. We thank the Department of Biostatitics at Johns Hopkins Bloomberg School of Public Health for hosting our application on their shinyapps account.'),
    p(HTML('<a href="https://www.libd.org/" onclick="ga(\'send\', \'event\', \'click\', \'link\', \'LIBD\', 1)"><img src="http://LieberInstitute.github.io/rstatsclub/img/LIBD.jpg" align="left" width = "250"/></a>')),
    p(HTML('<a href="http://www.jhsph.edu/departments/biostatistics/">'), img(src='http://aejaffe.com/media/jhu-bloomberg-logo.jpg', align = 'right', width = '250'), HTML('</a>')),
    tags$br(),
    tags$br(),
    tags$br(),
    tags$br(),
    tags$br(),
    tags$br(),
    tags$br(),
    tags$br()
))
