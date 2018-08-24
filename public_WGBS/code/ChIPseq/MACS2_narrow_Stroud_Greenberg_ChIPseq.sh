# qsub -cwd -l bluejay,mf=5G,h_vmem=5G,h_fsize=100G,h_stack=256M -t 1 -tc 30 MACS2_narrow_Stroud_Greenberg_ChIPseq.sh

module load macs/2.1.0

cd /dcl02/lieber/WGBS/PublicData/Stroud_Greenberg_Cell_2017/BAMS/ChIPseq

macs2 callpeak -t SRR6129731.bam SRR6129732.bam -c SRR6129696.bam -g mm -n polII.Cortex.WT.2wks --call-summits --outdir /dcl02/lieber/WGBS/PublicData/Stroud_Greenberg_Cell_2017/BAMS/ChIPseq/Peaks
