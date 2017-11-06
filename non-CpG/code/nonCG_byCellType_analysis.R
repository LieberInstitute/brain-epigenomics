load(file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/non-CpG/CH_object.rda")


table(CH$CT.dir=="pos")
# FALSE     TRUE 
# 7630202 33188540 
# 18.7% of CH is greater methylated in glia

table(CH[which(CH$padj.CellType<=0.05),"CT.dir"]=="pos")
#  FALSE    TRUE 
#  79239 7602836
# 99% of all significantly regulated sites by cell type are more methylated in neurons

table(CH$Age.dir=="pos")
#FALSE     TRUE 
#9843075 30975667
#75.9% are increasing over age

table(CH[which(CH$padj.Age<=0.05),"CT.dir"]=="pos")
#FALSE    TRUE 
#24967 3169651 
# 99.2% are increasing over age in significantly dmCH sites.
