# qsub -tc 245 -l mf=5G,h_vmem=5G,h_stack=256M -m e -M amanda.joy.price@gmail.com -t 1 sort.CpG.bed.sh

cd /dcl01/lieber/ajaffe/lab/brain-epigenomics/bumphunting/ 

sort -k1,1 -k2,2n BSobj_bsseqSmooth_Neuron_minCov_3.bed > /dcl01/lieber/ajaffe/Amanda/ATAC/Coverage/coverage_at_CpGs/BSobj_bsseqSmooth_Neuron_minCov_3_sorted.bed
