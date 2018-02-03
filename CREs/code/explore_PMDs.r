library("GenomicFeatures")
library("ggplot2")


load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/PMDs_methylSeekR.rda")
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/PMDs_methylSeekR_prenatal.rda")
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/PMDs_methylSeekR_sortedLister.rda")

# rename Lister samples

listernames = c("GSM1173773" = "ListerNeuron1", "GSM1173774" = "ListerGlia1", "GSM1173776" = "ListerNeuron2", "GSM1173777" = "ListerGlia2")
names(PMDsegments.CG.sortedLister) = listernames[match(names(PMDsegments.CG.sortedLister), names(listernames))]
names(PMDsegments.CG) = paste0("Postnatal-", names(PMDsegments.CG))
names(PMDsegments.CGpren) = paste0("Prenatal-", names(PMDsegments.CGpren))

# What proportion of the genome is covered by each sample?

sLengths = c(lapply(PMDsegments.CG.sortedLister, reduce), lapply(PMDsegments.CG, reduce), 
             lapply(PMDsegments.CGpren, reduce))
pmds_total = c(lapply(PMDsegments.CG.sortedLister, function(x) reduce(x[x$type=="PMD"])),
               lapply(PMDsegments.CG, function(x) reduce(x[x$type=="PMD"])),
               lapply(PMDsegments.CGpren, function(x) reduce(x[x$type=="PMD"])))

propPMDs = mapply(function(s, l) round(sum(as.numeric(width(s)))/sum(as.numeric(width(l))), 2), pmds_total, sLengths)
propPMDs[order(propPMDs)] # 


## What proportion of called PMDs are >100kb?

pmds = c(lapply(PMDsegments.CG.sortedLister, function(x) x[x$type=="PMD"]),lapply(PMDsegments.CG, function(x) x[x$type=="PMD"]),
         lapply(PMDsegments.CGpren, function(x) x[x$type=="PMD"]))

propPMDs = round(elementNROWS(lapply(pmds, function(x) x[which(width(x)>100000)]))/elementNROWS(pmds), 2)
propPMDs[order(propPMDs)] # 3-10% of the reported PMDs are >100kb
# So PMDs aren't really a thing in brain cells...


## is there a difference in prenatal and postnatal?

num = elementNROWS(lapply(pmds, function(x) x[which(width(x)>100000)]))
t.test(num[grep("Pre",names(num))], num[grep("Post",names(num))]) # t = -0.17579, df = 45.961, p-value = 0.8612


## Create new objects reflecting size limit

pmds = c(lapply(PMDsegments.CG.sortedLister, function(x) x[x$type=="PMD"]),lapply(PMDsegments.CG, function(x) x[x$type=="PMD"]),
         lapply(PMDsegments.CGpren, function(x) x[x$type=="PMD"]))
nonpmds = c(lapply(PMDsegments.CG.sortedLister, function(x) x[x$type=="notPMD"]),lapply(PMDsegments.CG, function(x) x[x$type=="notPMD"]),
            lapply(PMDsegments.CGpren, function(x) x[x$type=="notPMD"]))
newnon = lapply(pmds, function(x) x[which(width(x)<100000)])

merge = mapply(function(gr1, gr2) reduce(c(gr1,gr2)), nonpmds, newnon)
oo1 = mapply(function(m, gr1) findOverlaps(m, gr1), merge, nonpmds)
oo2 = mapply(function(m, gr2) findOverlaps(m, gr2), merge, newnon)
grm1 = mapply(function(gr1, o1) split(gr1[subjectHits(o1)], queryHits(o1)), nonpmds, oo1)
grm2 = mapply(function(gr2, o2) split(gr2[subjectHits(o2)], queryHits(o2)), newnon, oo2)

merge2 = rep( list(list()), length(merge)) 
for (j in 1:length(merge)) {
  for (i in 1:length(merge[[j]])) {
    if (as.character(i) %in% names(grm2[[j]])) {
    merge2[[j]][[i]] = c(grm1[[j]][[as.character(i)]], grm2[[j]][[as.character(i)]]) } else {
    merge2[[j]][[i]] = grm1[[j]][[as.character(i)]]
  } } }
nCG = lapply(merge2, function(y) lapply(y, function(x) sum(x$nCG)))
merge = Map(cbind, lapply(merge, as.data.frame), type = "notPMD", nCG = lapply(nCG, unlist))
merge = lapply(merge, function(x) makeGRangesFromDataFrame(x, keep.extra.columns = T))

pmds = lapply(pmds, function(x) x[which(width(x)>=100000)])

total = mapply(function(non, pmd) c(non, pmd), merge, pmds)
total = lapply(total, function(x) x[order(x)])

table(unlist(lapply(c(PMDsegments.CG.sortedLister, PMDsegments.CG, PMDsegments.CGpren), function(x) sum(x$nCG)))==18664892)
table(unlist(lapply(total, function(x) sum(x$nCG)))==18664892)
names(unlist(lapply(total, function(x) sum(x$nCG))))[unlist(lapply(total, function(x) sum(x$nCG)))!=18664892 & 
                                                       unlist(lapply(total, function(x) sum(x$nCG)))!=18308447]
# ??? "Postnatal-WGC052317L" "Postnatal-WGC052318L" "Postnatal-WGC052311L" "Postnatal-WGC059603L"

sLengths = c(lapply(PMDsegments.CG.sortedLister, reduce), lapply(PMDsegments.CG, reduce), 
             lapply(PMDsegments.CGpren, reduce))
red = lapply(total, reduce)
table(elementNROWS(red)==25) # all true
table(unlist(mapply(function(x,y) identical(x,y), red, sLengths))=="TRUE") # all true

save(total, file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/CREs/PMDs_methylSeekR_100kbLimit.rda")


## association with H3K27me3, H3K9me


## association with ATAC peaks








