qrsh

# run this code in /dcl02/lieber/WGBS/PublicData/Stroud_Greenberg_Cell_2017/BAMS/ChIPseq

bash /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/sh/step1-fastqc.sh --experiment "Stroud_Greenberg_ChIPseq" --prefix "Aug212018" --large TRUE
bash /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/sh/step2-trim.sh --experiment "Stroud_Greenberg_ChIPseq" --prefix "Aug212018" --large TRUE --cores 4
