#!/bin/bash

## Usage:
# sh explore_acf.sh

mkdir -p logs

#for chrnum in {1..22} X Y M
#for chrnum in M Y 21 22
for chrnum in {2..20} X
#for chrnum in {2..22} X Y M
do 

    chr="chr${chrnum}"
    
    for context in all nonCG CG CHG CHH
    do

    
    SHORT="acf_${chr}_context_${context}"

    # Construct shell file
    echo "Creating script for chromosome ${chr}"
    cat > .${SHORT}.sh <<EOF
#!/bin/bash
#$ -cwd
#$ -l mem_free=35G,h_vmem=40G,h_fsize=100G
#$ -N ${SHORT}
#$ -o ./logs/${SHORT}.txt
#$ -e ./logs/${SHORT}.txt
#$ -m e

echo "**** Job starts ****"
date

Rscript explore_acf.R -c ${chr} -t 1 -x ${context}

echo "**** Job ends ****"
date
EOF
        call="qsub .${SHORT}.sh"
        echo $call
        $call
    done
done
