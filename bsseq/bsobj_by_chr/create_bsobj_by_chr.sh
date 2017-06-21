#!/bin/bash

## Usage:
# sh create_bsobj_by_chr.sh

mkdir -p logs

for chrnum in {1..22} X Y M
do 

chr="chr${chrnum}"
SHORT="bsseq_${chr}"

# Construct shell file
echo "Creating script for chromosome ${chr}"
cat > .${SHORT}.sh <<EOF
#!/bin/bash
#$ -cwd
#$ -l mem_free=20G,h_vmem=25G
#$ -N ${SHORT}
#$ -pe local 8
#$ -o ./logs/${SHORT}.txt
#$ -e ./logs/${SHORT}.txt
#$ -m e

echo "**** Job starts ****"
date

Rscript create_bsobj_by_chr.R -c ${chr} -t 8

echo "**** Job ends ****"
date
EOF
    call="qsub .${SHORT}.sh"
    echo $call
    $call
done
