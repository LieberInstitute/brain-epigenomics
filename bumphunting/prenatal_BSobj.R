## Based on /users/ajaffe/Lieber/Projects/WGBS/Analysis/test_bsseq.R

###
# qrsh -l bluejay,mem_free=20G,h_vmem=20G,h_fsize=200G -pe local 6
# cd /dcl01/lieber/ajaffe/lab/brain-epigenomics/bumphunting
# module load conda_R/3.4.x

library('bsseq')
library('genefilter')
library('devtools')

# read in cleaned phenotype
pd  <- read.csv("/dcl01/lieber/WGBS/LIBD_Data/phenoData/FullMasterList_clean.csv",
    as.is=TRUE)
    
## Keep only prenatal ones
pd <- subset(pd, Age < 0)
dim(pd)

## all report files
pd$reportFiles  <- paste0("/dcl01/lieber/WGBS/LIBD_Data/DNAm_Ratios_duplicates_dropped/Reports/",
	pd$WGC.ID, "/", pd$WGC.ID, ".concatenated.sorted.duplicatesRemoved.CpG_report.txt")
stopifnot(all(file.exists(pd$reportFiles)))


## Read bismark files
bsList <- mclapply(seq_len(nrow(pd)), function(ii) {
    BSobj <- read.bismark(pd$reportFiles[ii], pd$WGC.ID[ii],
        strandCollapse=TRUE, fileType = "cytosineReport")
    cat(".")
	return(BSobj)
}, mc.cores=6, mc.preschedule=FALSE)


BSobj <- combineList(bsList)

## add sample names
sampleNames(BSobj) <- pd$WGC.ID
rownames(pd) <-pd$WGC.ID

## append phenotype data
pData(BSobj) <- DataFrame(pd)


## Load the raw data
foo <- function() {
    load('BSobj_bsseqSmooth_Neuron_minCov_3.Rdata')
    return(BSobj)
}
postnatal <- foo()

## Subset
ov <- findOverlaps(rowRanges(postnatal), rowRanges(BSobj))

save(BSobj, 
	file="")
    
    
## Reproducibility info
proc.time()
message(Sys.time())
options(width = 120)
session_info()
