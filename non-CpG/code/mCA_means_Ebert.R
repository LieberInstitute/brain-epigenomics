library(bsseq)
library(GenomicRanges)


### Collaboration with Dan Ebert and Crystal
### gather the mean mCA levels +/- 3kb around gene body per sample


## load CpH data

load('/dcl01/lieber/ajaffe/lab/brain-epigenomics/bsseq/bsobj_by_chr/allChrs_postNatal_cleaned_nonCG_noHomogenate_highCov.Rdata', verbose=T)
CHgr = granges(BSobj)
meth =getMeth(BSobj, type = 'raw')

# get mCA

filt = meth[which(CHgr$trinucleotide_context %in% c("CAA","CAC","CAG","CAN","CAT")),]
filtgr = CHgr[which(CHgr$trinucleotide_context %in% c("CAA","CAC","CAG","CAN","CAT"))]


## get genes +/- 3kb

load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/DMR/DMR_objects.rda", verbose=T)
geneMapGR = makeGRangesFromDataFrame(geneMap, keep.extra.columns = T)

extend <- function(x, upstream=0, downstream=0)     
{
  if (any(strand(x) == "*"))
    warning("'*' ranges were treated as '+'")
  on_plus <- strand(x) == "+" | strand(x) == "*"
  new_start <- start(x) - ifelse(on_plus, upstream, downstream)
  new_end <- end(x) + ifelse(on_plus, downstream, upstream)
  ranges(x) <- IRanges(new_start, new_end)
  trim(x)
}

extGenes = extend(geneMapGR, upstream = 3000, downstream = 3000)


## Get CAs in each gene

oo = findOverlaps(extGenes, filtgr)
rIndexes = split(subjectHits(oo), queryHits(oo))
names(rIndexes) = extGenes$gencodeID[unique(queryHits(oo))]
splmeth = sapply(rIndexes, function(ii) colMeans(t(t(filt[ii,]))))
meanCA = do.call("rbind", splmeth)

## how many CA within each extended gene region?

numCA = data.frame(ID = names(rIndexes), numCA = elementNROWS(rIndexes))
extGeneMap = as.data.frame(extGenes[unique(queryHits(oo))])
extGeneMap$number.CA = numCA[match(extGeneMap$gencodeID, numCA$ID),"numCA"]

save(meanCA, pd, extGeneMap, 
     file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/non-CpG/meanCA_per3KBExtendedGene_Ebert.rda")



#meanCov = lapply(meanCov, reshape2::melt)
#meanCov = do.call(rbind, Map(cbind, meanCov, Cluster = as.list(names(meanCov))))
#meanCov$Age = factor(pd[match(meanCov$Var2, pd$Data.ID),"Age"])
#meanCov$Cluster = factor(meanCov$Cluster, levels = c("1:G-N+","2:G0N+","3:G0N-","4:G+N0","5:G+N-","6:G-N0"))
#meanCov$CellType = factor(pd[match(meanCov$Var2, pd$Data.ID),"Cell.Type"])














