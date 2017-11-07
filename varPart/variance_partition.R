library('variancePartition')
library('bsseq')
library('limma')
library('getop')
library('devtools')
library('doParallel')

## Specify parameters
spec <- matrix(c(
    'set', 's', 1, 'character', 'Either: CpG or Non-CpG',
    'cores', 't', 1, 'integer', 'Number of cores to use',
#    'chromosome', 'c', 1, 'character', 'One of chr1 to chr22, chrX, chrY or chrM'
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
    opt <- list('set' = 'CpG')
}

if(opt$set == 'CpG') {
    load(file.path('/dcl01/lieber/ajaffe/lab/brain-epigenomics',
        'bumphunting', 'BSobj_bsseqSmooth_Neuron_minCov_3.Rdata'))
} else {
    load(file.path('/dcl01/lieber/ajaffe/lab/brain-epigenomics',
        'bsseq/bsobj_by_chr',
        'allChrs_postNatal_cleaned_nonCG_noHomogenate_highCov.Rdata'))
        
    ## Reproduce methylation filter used for the limma exploration:
    ## at least 5 samples have to have methylation greater than 0
    BSobj <- BSobj[rowSums(getMeth(BSobj, type = 'raw') > 0) >= 5, ]
}

## Data used
dim(BSobj)


## extract pheno
pd <- pData(BSobj)

## Fix pheno data
pd$Race[pd$Race == "CAUC "] <- 'CAUC'
pd$Sex[pd$Sex == " M"] <- 'M'
pd$RIN <- as.numeric(gsub(' ', '', pd$RIN))
pd$pg.DNA.nuclei.input <- as.numeric(pd$pg.DNA.nuclei.input)
if(opt$set == 'CpG') pd$Data.Yield. <- as.numeric(pd$Data.Yield.)
pd$Reads <- as.numeric(pd$Reads)
pd$Percent.GreaterThan.Q30 <- as.numeric(pd$Percent.GreaterThan.Q30)

##
cl <- makeCluster(4)
registerDoParallel(cl)

form <- ~ Cell.Type + Age + Age:Cell.Type + (1|Brain.ID) + Sex + Percent.GreaterThan.Q30 + Percent.Duplication

form <- ~ Cell.Type + Age + Age:Cell.Type + Race + Sex + Percent.GreaterThan.Q30 + Percent.Duplication

meth <- getMeth(BSobj[seq_len(1e5), ], type = 'raw')
varPart <- fitExtractVarPartModel(meth[1:10000, ], form, as.data.frame(pd), showWarnings = FALSE)

jo <- fitVarPartModel(meth[1:10000, ], form, as.data.frame(pd), showWarnings = FALSE)
length(jo)
xx <- sapply(jo, colinearityScore)
summary(xx)
colinearityScore(jo[[1]])
plotVarPart(extractVarPart(jo, showWarnings = FALSE))

## Joint modeling
message(paste(Sys.time(), 'performing joint modeling'))
y <- meth[1:10000, ]
system.time( sumSqList <- parallel::mclapply(seq_len(nrow(y)), function(i) {
	if(i %% 1000 == 0) cat(".")
        t(anova(lm(y[i,] ~ Cell.Type + Age + Age:Cell.Type + Race + Sex + Percent.GreaterThan.Q30 + Percent.Duplication, data= as.data.frame(pd)))[2])
}, mc.cores = 1) )

ssOut <- do.call("rbind", sumSqList)
rownames(ssOut) <- NULL
bg <- matrix(rep(rowSums(ssOut), ncol(ssOut)), 
	ncol = ncol(ssOut), nrow = nrow(ssOut))
ssMat <- ssOut / bg

#lab <- c('Age', 'Cell.Type', 'Brain.ID', 'RIN', 'Sex', 'Proportion.of.Neurons', 'yield.nuclei.mg.tissue', 'Race', 'Percent.GreaterThan.Q30', 'Percent.Duplication')
lab <- c('Cell.Type', 'Age', 'Race', 'Sex', 'Percent.GreaterThan.Q30', 'Percent.Duplication', 'Interaction', 'Residuals')
names(lab) <- colnames(ssMat)

par(mar=c(9,5,2,2))
library('RColorBrewer')

palette(brewer.pal(7, "Dark2"))
boxplot(100*ssMat, xaxt="n", ylim = c(0, 100), 
	cex.axis=1.3,cex.lab=1.1, range=2,
	ylab="Percentage variance explained", cex=0.5)
text(seq_len(ncol(ssMat)) + 0.2, y = -8, lab, xpd=TRUE, srt=45, pos=2)


## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
