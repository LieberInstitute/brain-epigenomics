table1 = read.table("/dcl02/lieber/WGBS/PublicData/Gabel_Greenberg_Nature_2015/FASTQ/ChIPseq/SraRunTable1.txt", header=T,sep="\t")
table2 = read.table("/dcl02/lieber/WGBS/PublicData/Gabel_Greenberg_Nature_2015/FASTQ/ChIPseq/SraRunTable2.txt", header=T,sep="\t")

pairs = list(MECP2.forebrain.WT.8wks = c(ChIP = "SRR1536114", Input = "SRR1536115"),
             MECP2.CB.WT.8wks = c(ChIP = "SRR1536116", Input = "SRR1536117"),
             MECP2.CB.WT.8wks =c(ChIP = "SRR1536118", Input = "SRR1536119"),
             MECP2.FC.WT.6wks = c(ChIP = "SRS884190", Input = "SRS884189"),
             MECP2.FC.WT.6wks = c(ChIP = "SRS884187", Input = "SRS884188"),
             MECP2.Cortex.WT.8wks = c(ChIP = "SRS884186", Input = "SRS884179"))
pairs = data.frame(Name = names(pairs), ChIP = unlist(lapply(pairs, function(x) x[1])), Input = unlist(lapply(pairs, function(x) x[2])))

write.table(pairs[,1], row.names=FALSE, col.names= FALSE, quote=FALSE,
            file = "/dcl02/lieber/WGBS/PublicData/Gabel_Greenberg_Nature_2015/FASTQ/ChIPseq/Test_ids.txt")
write.table(pairs[,2], row.names=FALSE, col.names= FALSE, quote=FALSE,
            file = "/dcl02/lieber/WGBS/PublicData/Gabel_Greenberg_Nature_2015/FASTQ/ChIPseq/ChIP_ids.txt")
write.table(pairs[,3], row.names=FALSE, col.names= FALSE, quote=FALSE,
            file = "/dcl02/lieber/WGBS/PublicData/Gabel_Greenberg_Nature_2015/FASTQ/ChIPseq/Input_ids.txt")
