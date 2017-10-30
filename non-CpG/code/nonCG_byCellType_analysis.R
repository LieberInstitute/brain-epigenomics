load(file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/non-CpG/CH_object.rda")


table(CH$CT.dir=="pos")
# FALSE     TRUE 
# 7630202 33188540 
# 18.7% of CH is greater methylated in glia

table(CH[which(CH$padj.CellType<=0.05),"CT.dir"]=="pos")
#  FALSE    TRUE 
#  79239 7602836
# 99% of all significantly regulated sites by cell type are more methylated in neurons