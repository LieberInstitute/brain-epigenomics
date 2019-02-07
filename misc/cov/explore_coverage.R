library('jaffelab')
library('readr')
library('sessioninfo')
library('ggplot2')

## For output files
dir.create('rda', showWarnings = FALSE)

##### Functions:
## Code adapted from https://github.com/LieberInstitute/RNAseq-pipeline/blob/master/sh/create_count_objects-human.R#L209-L217
splitAt <- function(x, pos) unname(split(x, cumsum(seq_along(x) %in% pos)))

compute_fastqc_seqdist <- function(qcData) {
    R <- lapply(qcData, function(x) scan(x, what = "character", sep= "\n", 
    		quiet = TRUE, strip=TRUE) )
    ## Split list into sublists of metric categories
    zz <- lapply(R, function(x) splitAt(x, which(x==">>END_MODULE")+1))

    ## Extract the sequence length distribution
    seqlen <- lapply(zz, function(x) {
        ## access the sequence distribution part of the report
        y <- x[[8]]
        ## Drop header and module separating lines
        y[-c(1:2, length(y))]
    })

    seqdist <- lapply(seqlen, function(dat) {
    
        ## Process the read length section
        len <- ss(dat, "\t", 1)
    
        ## Sometimes there's a range of read lengths
        res <- do.call(rbind, mapply(function(x, hasdash) {
            if(!hasdash) {
                x <- as.numeric(x)
                result <- data.frame(
                    start = x,
                    end = x,
                    stringsAsFactors = FALSE
                )
            } else {
                ## Get the start/end of the range
                result <- data.frame(
                    start = as.numeric(ss(x, '-', 1)),
                    end = as.numeric(ss(x, '-', 2)),
                    stringsAsFactors = FALSE
                )
            }
            ## Compute the median read length based on the range info
            result$median <- median(c(result$start, result$end))
            return(result)
        }, len, grepl('-', len), SIMPLIFY = FALSE))
    
        res$n <- as.numeric(ss(dat, "\t", 2))
    
        return(res)
    })
    return(seqdist)
}
compute_fastqc_auc <- function(qcData) {
    auc <- sapply(compute_fastqc_seqdist(qcData), function(dat) {
        sum(dat$median * dat$n)
    })
    return(auc)
}

## bamcount reader
read_bamcount <- function(f) {
    readr::read_tsv(f, col_names = c('cat', 'n'))$n
}


#### Exploratory code

## Get ids and genome size
ids <- readLines('WGC_IDs_subset.txt')
hg19 <- readr::read_tsv(
    '/dcl01/lieber/ajaffe/lab/brain-epigenomics/misc/cov/test/human.hg19.genome',
    col_names = c('chr', 'length')
)
hg19.size <- sum(hg19$length[hg19$chr %in% paste0('chr', c(1:22, 'X', 'Y', 'M'))])

## Read in fastqc data
get_fastqc_data <- function(id) {
    message(paste(Sys.time(), 'processing sample', id))
    
    ## Build the paths to the FastQC output directories
    id_dirs <- file.path(
        '/dcl02/lieber/WGBS/LIBD_Data/FastQC',
        c('Untrimmed', 'Trimmed'),
    id)
    
    ## Find fastqc data files
    fastqc_files <- dir(id_dirs, pattern = 'data.txt$', recursive = TRUE, full.names = TRUE)
    stopifnot(length(fastqc_files) == 7)
    
    names(fastqc_files) <- gsub('^[[:alnum:]]*_|_fastqc', '', ss(fastqc_files, '/', 9))

    auc_dir <- compute_fastqc_auc(fastqc_files)
    
    res <- data.frame(
        sample = id,
        type = names(auc_dir),
        auc = auc_dir,
        stringsAsFactors = FALSE
    )
    rownames(res) <- NULL
    return(res)
}

cov_fastqc <- do.call(rbind, lapply(ids, get_fastqc_data))
# 2019-02-06 14:34:23 processing sample WGC052316L
# 2019-02-06 14:34:24 processing sample WGC059613L
# 2019-02-06 14:34:32 processing sample WGC059614L
# 2019-02-06 14:34:53 processing sample WGC052317L
# 2019-02-06 14:35:11 processing sample WGC055558L
# 2019-02-06 14:35:33 processing sample WGC055559L
# 2019-02-06 14:35:52 processing sample WGC059596L
# 2019-02-06 14:36:11 processing sample WGC055561L
# 2019-02-06 14:36:33 processing sample WGC052318L
# 2019-02-06 14:36:49 processing sample WGC059588L
# 2019-02-06 14:37:08 processing sample WGC059608L
# 2019-02-06 14:37:30 processing sample WGC059594L
# 2019-02-06 14:37:44 processing sample WGC052309L_reseq
# 2019-02-06 14:38:03 processing sample WGC059592L
# 2019-02-06 14:38:25 processing sample WGC052311L
# 2019-02-06 14:38:42 processing sample WGC059612L
# 2019-02-06 14:38:56 processing sample WGC055562L
# 2019-02-06 14:39:17 processing sample WGC055573L
# 2019-02-06 14:39:30 processing sample WGC052324L
# 2019-02-06 14:39:50 processing sample WGC055575L
# 2019-02-06 14:40:11 processing sample WGC059600L
# 2019-02-06 14:40:33 processing sample WGC052328L
# 2019-02-06 14:40:50 processing sample WGC052314L
# 2019-02-06 14:41:13 processing sample WGC059598L
# 2019-02-06 14:41:31 processing sample WGC059603L
# 2019-02-06 14:41:47 processing sample WGC055577L
# 2019-02-06 14:41:55 processing sample WGC059607L
# 2019-02-06 14:42:01 processing sample WGC052312L
# 2019-02-06 14:42:19 processing sample WGC055570L
# 2019-02-06 14:42:32 processing sample WGC055568L
# 2019-02-06 14:42:50 processing sample WGC059601L
# 2019-02-06 14:43:14 processing sample WGC059604L

save(cov_fastqc, file = 'rda/cov_fastqc.Rdata')

## Read in bamcount output data
get_bamcount_data <- function(id) {
    message(paste(Sys.time(), 'processing sample', id))
    
    ## Find the files
    types <- c('duplicatesRemoved', 'Marked_duplicates')
    bamcount_files <- file.path(
        '/dcl01/lieber/ajaffe/lab/brain-epigenomics/misc/cov/auc_files',
        paste0(id, '_', types, '.auc.tsv')
    )
    stopifnot(all(file.exists(bamcount_files)))
    
    ## Read the data in
    auc_bam <- sapply(bamcount_files, read_bamcount)
    
    res <- data.frame(
        sample = id,
        type = types,
        auc = auc_bam,
        stringsAsFactors = FALSE
    )
    rownames(res) <- NULL
    return(res)
}

cov_bamcount <- do.call(rbind, lapply(ids, get_bamcount_data))
save(cov_bamcount, file = 'rda/cov_bamcount.Rdata')

## Merge both types
cov_merged <- rbind(cov_fastqc, cov_bamcount)
cov_merged$stage <- 'raw'

cov_summarized <- do.call(rbind, lapply(split(cov_merged, cov_merged$sample), function(dat) {
    
    types <- c('initial', 'post-trim', 'mapped', 'no-dups')
    res <- data.frame(
        sample = unique(dat$sample),
        type = factor(types, levels = types),
        auc = c(
            sum(dat$auc[dat$type %in% c('combined_R1', 'combined_R2')]),
            sum(dat$auc[grepl('flash|unpaired', dat$type)]),
            dat$auc[dat$type == 'Marked_duplicates'],
            dat$auc[dat$type == 'duplicatesRemoved']
            ),
        stage = 'summarized',
        stringsAsFactors = FALSE
    )
    return(res)
}))
rownames(cov_summarized) <- NULL

subset(cov_summarized, auc == 0)
#             sample    type auc      stage
# 1 WGC052309L_reseq initial   0 summarized

subset(cov_fastqc, sample == 'WGC052309L_reseq')
#              sample                            type         auc
# 85 WGC052309L_reseq reseq_flashOutput.extendedFrags 39919651703
# 86 WGC052309L_reseq reseq_flashOutput.notCombined_1  5273155607
# 87 WGC052309L_reseq reseq_flashOutput.notCombined_2  4929507357
# 88 WGC052309L_reseq   reseq_output_forward_unpaired 20942023974
# 89 WGC052309L_reseq   reseq_output_reverse_unpaired   303701690
# 90 WGC052309L_reseq               reseq_combined_R1 66313886100
# 91 WGC052309L_reseq               reseq_combined_R2 66313886100
subset(cov_merged, sample == 'WGC052309L_reseq')
#               sample                            type         auc stage
# 85  WGC052309L_reseq reseq_flashOutput.extendedFrags 39919651703   raw
# 86  WGC052309L_reseq reseq_flashOutput.notCombined_1  5273155607   raw
# 87  WGC052309L_reseq reseq_flashOutput.notCombined_2  4929507357   raw
# 88  WGC052309L_reseq   reseq_output_forward_unpaired 20942023974   raw
# 89  WGC052309L_reseq   reseq_output_reverse_unpaired   303701690   raw
# 90  WGC052309L_reseq               reseq_combined_R1 66313886100   raw
# 91  WGC052309L_reseq               reseq_combined_R2 66313886100   raw
# 249 WGC052309L_reseq               duplicatesRemoved 28863207514   raw
# 250 WGC052309L_reseq               Marked_duplicates 59246876569   raw

cov_merged$coverage <- cov_merged$auc / hg19.size
cov_summarized$coverage <- cov_summarized$auc / hg19.size

save(cov_merged, file = 'rda/cov_merged.Rdata')
save(cov_summarized, file = 'rda/cov_summarized.Rdata')

cov_summarized$sample <- as.factor(cov_summarized$sample)
ggplot(cov_summarized, aes(x = type, y = coverage, group = sample)) + geom_line()

boxplot(coverage ~ type, data = cov_summarized)

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
