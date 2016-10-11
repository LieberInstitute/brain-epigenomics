library('bsseq')
library('derfinderHelper')
library('derfinder')
library('limma')
library('devtools')

## Load data
datadir <- '/dcl01/lieber/WGBS/LIBD_Data/DNAm_Ratios_duplicates_dropped/CpGs/BedFiles'
files <- dir(datadir, pattern = 'chr21.*cov', full.names = TRUE)
names(files) <- dir(datadir, pattern = 'chr21.*cov')

## Read phenotype info
pd <- read.csv('/dcl01/lieber/WGBS/LIBD_Data/DNAm_Ratios_duplicates_dropped/CpGs/BedFiles/sample.anno.chr21.csv', stringsAsFactors = FALSE)
map <- match(names(files), pd$Coverage.Files)

## Sort phenotype info by files
pd <- pd[map, ]

## Read in data
message(paste(Sys.time(), 'loading data from bismark cov files'))
bs <- read.bismark(files, sampleNames = pd$Sample.ID, strandCollapse = FALSE, fileType = 'cov')
## Add phenotype info
pData(bs) <- pd


## The average coverage of CpGs on chr21
round(colMeans(getCoverage(bs)), 1)

## Number of CpGs in chr21
length(bs)

## Number of CpGs which are covered by at least 1 read in all samples
sum(rowSums(getCoverage(bs) >= 1) == ncol(bs))

## Number of CpGs with 0 coverage in all samples
sum(rowSums(getCoverage(bs)) == 0)

## 
bs_cov <- getCoverage(bs)


bs_meth <- getMeth(bs, type = 'raw', what = "perBase", confint = FALSE)

logit <- function(p) { p / (1 - p )}
trans <- function(meth, adj = 0.00000000000001) {
    meth[is.nan(meth)] <- adj
    meth[meth == 1] <- 1 - adj
    logit(meth)
}

## Try logit transforming
bs_m <- trans(bs_meth)

## Try simply making the NaNs into a small number
bs_m2 <- bs_meth
bs_m2[is.nan(bs_m2)] <- 0.00000000000001

## Define some models
mod0 <- model.matrix(~ pd$Sex + pd$Race)
mod1 <- model.matrix(~ pd$Age + pd$Sex + pd$Race)


## Calculate the F-stats
f <- fstats.apply(data = bs_m, mod = mod1, mod0 = mod0, method = 'regular')
f2 <- fstats.apply(data = bs_m2, mod = mod1, mod0 = mod0, method = 'regular')

## Fix some NAs manually
f[is.na(f)] <- 0
f2[is.na(f2)] <- 0

## Find a cutoff
cutoff <- derfinder:::.calcFstatCutoff('theoretical', 1e-4, models = list(mod = mod1, mod0 = mod0))

## Get positions
pos <- Rle(seq_len(max(end(bs))) %in% start(bs))

## Find regions
regs <- findRegions(position = pos, fstats = f, cutoff = cutoff, chr = 'chr21')
regs2 <- findRegions(position = pos, fstats = f2, cutoff = cutoff, chr = 'chr21')

## Check the regions
regs
regs2

table(countOverlaps(regs, regs2))
table(countOverlaps(regs2, regs))


## Run analysis with limma

# get the f statistic from 2 lmFit objects
getF <- function(fit, fit0, theData) {
	
	rss1 = rowSums((fitted(fit)-theData)^2)
	df1 = ncol(fit$coef)
	rss0 = rowSums((fitted(fit0)-theData)^2)
	df0 = ncol(fit0$coef)

	fstat = ((rss0-rss1)/(df1-df0))/(rss1/(ncol(theData)-df1))
	f_pval = pf(fstat, df1-df0, ncol(theData)-df1,lower.tail=FALSE)
	fout = cbind(fstat,df1-df0,ncol(theData)-df1,f_pval)
	colnames(fout)[2:3] = c("df1", "df0")
	fout = data.frame(fout, row.names = NULL)
	return(fout)
}

fit <- lmFit(bs_m2, mod1)
eb <- ebayes(fit)
fit0 <- lmFit(bs_m2, mod0)
ff <- getF(fit, fit0, bs_m2)

## Fix NAs
ff$fstat[is.na(ff$fstat)] <- 0
ff$f_pval[is.na(ff$f_pval)] <- 1

regs3 <- findRegions(position = pos, fstats = Rle(ff$f_pval), cutoff = c(1e-04, 1.01), chr = 'chr21')
regs4 <- findRegions(position = pos, fstats = Rle(ff$fstat), cutoff = cutoff, chr = 'chr21')

## Doesn't matter if I select by p-value or not
identical(ranges(regs3, use.names = FALSE), ranges(regs4, use.names = FALSE))

## Check with regs
table(countOverlaps(regs, regs4))
table(countOverlaps(regs4, regs))

## Identical by limma or by derfinderHelper::fstats.apply with this cutoff
identical(ranges(regs2), ranges(regs4))

## Reproducibility info
proc.time()
options(width = 120)
session_info()
