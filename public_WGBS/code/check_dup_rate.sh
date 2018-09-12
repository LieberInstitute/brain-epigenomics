# qsub -cwd -l bluejay,mf=150G,h_vmem=150G,h_fsize=500G,h_stack=256M -t 1 check_dup_rate.sh

module load samtools/1.1
module load java/1.8.0

#samtools sort /dcl02/lieber/WGBS/PublicData/Stroud_Greenberg_Cell_2017/BAMS/WGBS/SRR6126637_1_bismark_bt2_pe.bam /dcl02/lieber/WGBS/PublicData/Stroud_Greenberg_Cell_2017/BAMS/WGBS/SRR6126637.sorted
java -jar /dcl01/lieber/WGBS/Software/picard.jar ValidateSamFile I=/dcl02/lieber/WGBS/PublicData/Stroud_Greenberg_Cell_2017/BAMS/WGBS/SRR6126637_1_bismark_bt2_pe.bam MODE=SUMMARY IS_BISULFITE_SEQUENCED=TRUE
java -jar /dcl01/lieber/WGBS/Software/picard.jar ValidateSamFile I=/dcl02/lieber/WGBS/PublicData/Stroud_Greenberg_Cell_2017/BAMS/WGBS/SRR6126637.sorted.bam MODE=SUMMARY IS_BISULFITE_SEQUENCED=TRUE

java -Xmx50G -jar /dcl01/lieber/WGBS/Software/picard.jar MarkDuplicates \
    I=/dcl02/lieber/WGBS/PublicData/Stroud_Greenberg_Cell_2017/BAMS/WGBS/SRR6126637.sorted.bam \
    M=/dcl02/lieber/WGBS/PublicData/Stroud_Greenberg_Cell_2017/BAMS/WGBS/SRR6126637_duplicateMetrics.txt REMOVE_DUPLICATES=true \
    O=/dcl02/lieber/WGBS/PublicData/Stroud_Greenberg_Cell_2017/BAMS/WGBS/SRR6126637_duplicateRemoved.bam
samtools index /dcl02/lieber/WGBS/PublicData/Stroud_Greenberg_Cell_2017/BAMS/WGBS/SRR6126637_duplicateRemoved.bam

echo **** finished sample 1 ****

java -jar /dcl01/lieber/WGBS/Software/picard.jar ValidateSamFile I=/dcl02/lieber/WGBS/PublicData/Gabel_Greenberg_Nature_2015/BAMS/WGBS/paired/SRR1930026_1_bismark_bt2_pe.bam MODE=SUMMARY IS_BISULFITE_SEQUENCED=TRUE
samtools sort /dcl02/lieber/WGBS/PublicData/Gabel_Greenberg_Nature_2015/BAMS/WGBS/paired/SRR1930026_1_bismark_bt2_pe.bam /dcl02/lieber/WGBS/PublicData/Gabel_Greenberg_Nature_2015/BAMS/WGBS/paired/SRR1930026.sorted
java -jar /dcl01/lieber/WGBS/Software/picard.jar ValidateSamFile I=/dcl02/lieber/WGBS/PublicData/Gabel_Greenberg_Nature_2015/BAMS/WGBS/paired/SRR1930026.sorted.bam MODE=SUMMARY IS_BISULFITE_SEQUENCED=TRUE

java -Xmx50G -jar /dcl01/lieber/WGBS/Software/picard.jar MarkDuplicates \
    I=/dcl02/lieber/WGBS/PublicData/Gabel_Greenberg_Nature_2015/BAMS/WGBS/paired/SRR1930026.sorted.bam \
    M=/dcl02/lieber/WGBS/PublicData/Gabel_Greenberg_Nature_2015/BAMS/WGBS/paired/SRR1930026_duplicateMetrics.txt REMOVE_DUPLICATES=true \
    O=/dcl02/lieber/WGBS/PublicData/Gabel_Greenberg_Nature_2015/BAMS/WGBS/paired/SRR1930026_duplicateRemoved.bam
samtools index /dcl02/lieber/WGBS/PublicData/Gabel_Greenberg_Nature_2015/BAMS/WGBS/paired/SRR1930026_duplicateRemoved.bam

echo **** finished sample 2 ****

#samtools sort /dcl02/lieber/WGBS/PublicData/Gabel_Greenberg_Nature_2015/BAMS/WGBS/single/SRR1536127_bismark_bt2.bam /dcl02/lieber/WGBS/PublicData/Gabel_Greenberg_Nature_2015/BAMS/WGBS/single/SRR1536127.sorted
java -jar /dcl01/lieber/WGBS/Software/picard.jar ValidateSamFile I=/dcl02/lieber/WGBS/PublicData/Gabel_Greenberg_Nature_2015/BAMS/WGBS/single/SRR1536127_bismark_bt2.bam MODE=SUMMARY IS_BISULFITE_SEQUENCED=TRUE
java -jar /dcl01/lieber/WGBS/Software/picard.jar ValidateSamFile I=/dcl02/lieber/WGBS/PublicData/Gabel_Greenberg_Nature_2015/BAMS/WGBS/single/SRR1536127.sorted.bam MODE=SUMMARY IS_BISULFITE_SEQUENCED=TRUE

java -Xmx50G -jar /dcl01/lieber/WGBS/Software/picard.jar MarkDuplicates \
    I=/dcl02/lieber/WGBS/PublicData/Gabel_Greenberg_Nature_2015/BAMS/WGBS/single/SRR1536127.sorted.bam \
    M=/dcl02/lieber/WGBS/PublicData/Gabel_Greenberg_Nature_2015/BAMS/WGBS/single/SRR1536127_duplicateMetrics.txt REMOVE_DUPLICATES=true \
    O=/dcl02/lieber/WGBS/PublicData/Gabel_Greenberg_Nature_2015/BAMS/WGBS/single/SRR1536127_duplicateRemoved.bam
samtools index /dcl02/lieber/WGBS/PublicData/Gabel_Greenberg_Nature_2015/BAMS/WGBS/single/SRR1536127_duplicateRemoved.bam

echo **** finished sample 3 ****
