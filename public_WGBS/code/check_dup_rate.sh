# qsub -cwd -l bluejay,mf=60G,h_vmem=60G,h_fsize=100G,h_stack=256M -t 1 -tc 30 check_dup_rate.sh

module load samtools/1.1
module load java/1.8.0

samtools sort -T /dcl02/lieber/WGBS/PublicData/Stroud_Greenberg_Cell_2017/BAMS/WGBS/SRR6126637.sorted /dcl02/lieber/WGBS/PublicData/Stroud_Greenberg_Cell_2017/BAMS/WGBS/SRR6126637_1_bismark_bt2_pe.bam
java -Xmx50G -jar /dcl01/lieber/WGBS/Software/picard.jar MarkDuplicates \
    I=/dcl02/lieber/WGBS/PublicData/Stroud_Greenberg_Cell_2017/BAMS/WGBS/SRR6126637.sorted.bam \
    M=/dcl02/lieber/WGBS/PublicData/Stroud_Greenberg_Cell_2017/BAMS/WGBS/SRR6126637_duplicateMetrics.txt REMOVE_DUPLICATES=true \
    O=/dcl02/lieber/WGBS/PublicData/Stroud_Greenberg_Cell_2017/BAMS/WGBS/SRR6126637_duplicateRemoved.bam
samtools index /dcl02/lieber/WGBS/PublicData/Stroud_Greenberg_Cell_2017/BAMS/WGBS/SRR6126637_duplicateRemoved.bam

samtools sort -T /dcl02/lieber/WGBS/PublicData/Gabel_Greenberg_Nature_2015/BAMS/WGBS/paired/SRR1930026.sorted /dcl02/lieber/WGBS/PublicData/Gabel_Greenberg_Nature_2015/BAMS/WGBS/paired/SRR1930026_1_bismark_bt2_pe.bam
java -Xmx50G -jar /dcl01/lieber/WGBS/Software/picard.jar MarkDuplicates \
    I=/dcl02/lieber/WGBS/PublicData/Gabel_Greenberg_Nature_2015/BAMS/WGBS/paired/SRR1930026_1_bismark_bt2_pe.bam \
    M=/dcl02/lieber/WGBS/PublicData/Gabel_Greenberg_Nature_2015/BAMS/WGBS/paired/SRR1930026_duplicateMetrics.txt REMOVE_DUPLICATES=true \
    O=/dcl02/lieber/WGBS/PublicData/Gabel_Greenberg_Nature_2015/BAMS/WGBS/paired/SRR1930026_duplicateRemoved.bam
samtools index /dcl02/lieber/WGBS/PublicData/Gabel_Greenberg_Nature_2015/BAMS/WGBS/paired/SRR1930026_duplicateRemoved.bam

samtools sort -T /dcl02/lieber/WGBS/PublicData/Gabel_Greenberg_Nature_2015/BAMS/WGBS/single/SRR1536127.sorted /dcl02/lieber/WGBS/PublicData/Gabel_Greenberg_Nature_2015/BAMS/WGBS/single/SRR1536127_bismark_bt2.bam
java -Xmx50G -jar /dcl01/lieber/WGBS/Software/picard.jar MarkDuplicates \
    I=/dcl02/lieber/WGBS/PublicData/Gabel_Greenberg_Nature_2015/BAMS/WGBS/single/SRR1536127_bismark_bt2.bam \
    M=/dcl02/lieber/WGBS/PublicData/Gabel_Greenberg_Nature_2015/BAMS/WGBS/single/SRR1536127_duplicateMetrics.txt REMOVE_DUPLICATES=true \
    O=/dcl02/lieber/WGBS/PublicData/Gabel_Greenberg_Nature_2015/BAMS/WGBS/single/SRR1536127_duplicateRemoved.bam
samtools index /dcl02/lieber/WGBS/PublicData/Gabel_Greenberg_Nature_2015/BAMS/WGBS/single/SRR1536127_duplicateRemoved.bam
