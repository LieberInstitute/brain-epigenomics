
qrsh

# run this code in /dcl02/lieber/WGBS/PublicData/Gabel_Greenberg_Nature_2015/BAMS/WGBS

bash /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/sh/step1-fastqc.sh --experiment "Gabel_Greenberg_WGBS" --prefix "Aug202018" --large TRUE
bash /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/sh/step2-trim.sh --experiment "Gabel_Greenberg_WGBS" --prefix "Aug202018" --large TRUE --cores 4
