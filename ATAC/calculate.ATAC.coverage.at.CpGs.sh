# qsub -tc 245 -l mf=10G,h_vmem=10G,h_fsize=15G,h_stack=256M -m e -M amanda.joy.price@gmail.com -t 1-34 calculate.ATAC.coverage.at.CpGs.sh

FILELIST=/dcl01/lieber/ajaffe/Amanda/ATAC/FinalPeaks/scripts/finalizedList-pooled.txt
ID=$(awk "NR==$SGE_TASK_ID" $FILELIST)

PEhs3=/dcl01/lieber/ajaffe/Amanda/ATAC/FinalPeaks/individual_sample_peaks/${ID}.final.10bp-added2.sorted.bed
SEhs3=/dcl01/lieber/ajaffe/Amanda/ATAC/BAMs/unmerged_remaining_SE/HS3000/${ID}_SEhs3.final.filtered.sorted.bed
SEhs2=/dcl01/lieber/ajaffe/Amanda/ATAC/BAMs/unmerged_remaining_SE/HS2000/${ID}_SEhs2.final.filtered.sorted.bed

cd /dcl01/lieber/ajaffe/Amanda/ATAC/Coverage/coverage_at_CpGs
module load bedtools/2.24.0

coverageBed -sorted -g ../../genome_files/hg19.genome.reordered -counts -a /dcl01/lieber/ajaffe/Amanda/ATAC/Coverage/coverage_at_CpGs/{{{CPGS_FILE}}}.bed -b $PEhs3 > ${ID}_PEhs3_CpGs_coverageBed.txt
coverageBed -sorted -g ../../genome_files/hg19.genome.reordered -counts -a /dcl01/lieber/ajaffe/Amanda/ATAC/Coverage/coverage_at_CpGs/{{{CPGS_FILE}}}.bed -b $SEhs3 > ${ID}_SEhs3_CpGs_coverageBed.txt
coverageBed -sorted -g ../../genome_files/hg19.genome.reordered -counts -a /dcl01/lieber/ajaffe/Amanda/ATAC/Coverage/coverage_at_CpGs/{{{CPGS_FILE}}}.bed -b $SEhs2 > ${ID}_SEhs2_CpGs_coverageBed.txt
