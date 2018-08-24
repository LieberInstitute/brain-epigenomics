## Find pairs of ChIP and input

tb = read.table("/dcl02/lieber/WGBS/PublicData/Stroud_Greenberg_Cell_2017/FASTQ/ChIPseq/SraRunTable.txt", header=T,sep="\t", fill=T)

pairs = list(polII.phos.Cortex.WT.2wks = c(ChIP = "SRR6129729", Input = "SRR6129696"),
             polII.phos.Cortex.WT.2wks = c(ChIP = "SRR6129730", Input = "SRR6129696"),
             MECP2.Cortex.VIP.8wks = c(ChIP = "SRR6129747", Input = "SRR6129751"),
             MECP2.Cortex.VIP.8wks = c(ChIP = "SRR6129748", Input = "SRR6129751"),
             MECP2.Cortex.PV.8wks = c(ChIP = "SRR6129749", Input = "SRR6129752"),
             MECP2.Cortex.PV.8wks = c(ChIP = "SRR6129750", Input = "SRR6129752"))
pairs = data.frame(Name = names(pairs), ChIP = unlist(lapply(pairs, function(x) x[1])), Input = unlist(lapply(pairs, function(x) x[2])))
                         
write.table(pairs[,1], row.names=FALSE, col.names= FALSE, quote=FALSE,
            file = "/dcl02/lieber/WGBS/PublicData/Stroud_Greenberg_Cell_2017/FASTQ/ChIPseq/Test_ids.txt")
write.table(pairs[,2], row.names=FALSE, col.names= FALSE, quote=FALSE,
            file = "/dcl02/lieber/WGBS/PublicData/Stroud_Greenberg_Cell_2017/FASTQ/ChIPseq/ChIP_ids.txt")
write.table(pairs[,3], row.names=FALSE, col.names= FALSE, quote=FALSE,
            file = "/dcl02/lieber/WGBS/PublicData/Stroud_Greenberg_Cell_2017/FASTQ/ChIPseq/Input_ids.txt")


narrow = list(polII.Cortex.WT.2wks =c(ChIP = "SRR6129731", Input = "SRR6129696"),
              polII.Cortex.WT.2wks = c(ChIP = "SRR6129732", Input = "SRR6129696"))
pairs = data.frame(Name = names(narrow), ChIP = unlist(lapply(narrow, function(x) x[1])), Input = unlist(lapply(narrow, function(x) x[2])))

write.table(pairs[,1], row.names=FALSE, col.names= FALSE, quote=FALSE,
            file = "/dcl02/lieber/WGBS/PublicData/Stroud_Greenberg_Cell_2017/FASTQ/ChIPseq/Test_ids_narrow.txt")
write.table(pairs[,2], row.names=FALSE, col.names= FALSE, quote=FALSE,
            file = "/dcl02/lieber/WGBS/PublicData/Stroud_Greenberg_Cell_2017/FASTQ/ChIPseq/ChIP_ids_narrow.txt")
write.table(pairs[,3], row.names=FALSE, col.names= FALSE, quote=FALSE,
            file = "/dcl02/lieber/WGBS/PublicData/Stroud_Greenberg_Cell_2017/FASTQ/ChIPseq/Input_ids_narrow.txt")
