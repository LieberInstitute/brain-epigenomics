library('clusterProfiler')
library('VennDiagram')
library('GenomicRanges')
library('gplots')


files <- c('/dcl01/lieber/ajaffe/lab/brain-epigenomics/meth_vs_expr/rda/meqtl_venn_go_exon_using_near.Rdata', '/dcl01/lieber/ajaffe/lab/brain-epigenomics/meth_vs_expr/rda/meqtl_venn_go_gene_using_near.Rdata', '/dcl01/lieber/ajaffe/lab/brain-epigenomics/meth_vs_expr/rda/meqtl_venn_go_psi_using_near.Rdata')

go <- lapply(files, function(f) {
    message(paste(Sys.time(), 'loading', f))
    load(f, verbose = TRUE)
    return(go_cluster_comp)
})
names(go) <- c('exon', 'gene', 'psi')


keep <- c('CpH', 'CpG:CpH:CpGmarg', 'CpH:CpGmarg')

pdf('pdf/go_unique.pdf', width = 12, useDingbats = FALSE)
go_all <- lapply(names(go), function(feat) {
    resbp <- lapply(names(go[[1]]), function(bp) {
        gos <- go[[feat]][[bp]]
        comp <- as.data.frame(gos)

        comp = split(comp$Description, comp$Cluster)
        comp = comp[elementNROWS(comp)>0]

        names(comp) <- gsub('nonCpG', 'CpH', names(comp))

        ov <- venn(comp, show.plot = FALSE)
        ov_uniq <- unlist(attr(ov, 'intersections')[keep])

        gos@compareClusterResult$Cluster <- gsub('nonCpG', 'CpH', gos@compareClusterResult$Cluster)

        gos@compareClusterResult = gos@compareClusterResult[which(as.character(gos@compareClusterResult$Cluster) %in% keep),]

        gos@compareClusterResult = gos@compareClusterResult[which(gos@compareClusterResult$Description %in% ov_uniq),]
        if(any(gos@compareClusterResult$Cluster == 'CpG:CpH:CpGmarg')) {
            gos@compareClusterResult$Cluster[gos@compareClusterResult$Cluster == 'CpG:CpH:CpGmarg'] <- 'all'
        }
        if(any(grepl('regulation ', gos@compareClusterResult$Description))) {
            gos@compareClusterResult$Description <- gsub('regulation ', 'reg. ', gos@compareClusterResult$Description)
        }
        if(any(grepl('differentiation', gos@compareClusterResult$Description))) {
            gos@compareClusterResult$Description <- gsub('differentiation', 'diff.', gos@compareClusterResult$Description)
        }
        if(any(grepl('development', gos@compareClusterResult$Description))) {
            gos@compareClusterResult$Description <- gsub('development', 'devel.', gos@compareClusterResult$Description)
        }
        
        print(plot(gos, colorBy="p.adjust", showCategory = 10, title= paste(feat, bp), font.size = 18))
        return(gos)
    })
    names(resbp) <- names(go[[1]])
    return(resbp)
})
names(go_all) <- names(go)

dev.off()