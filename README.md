# brain-epigenomics
[![DOI](https://zenodo.org/badge/68875556.svg)](https://zenodo.org/badge/latestdoi/68875556)

Located in JHPCE at `/dcl01/lieber/ajaffe/lab/brain-epigenomics`

This repository contains the code used in the analyses in the following citation:

Price AJ, Collado-Torres L, Ivanov NA, Xia W, Burke EE, Shin JH, Tao R, Ma L, Jia Y,  Hyde TM, Kleinman JE, Weinberger DR, Jaffe AE. "Divergent neuronal DNA methylation patterns across human cortical development: Critical periods and a unique role of CpH methylation." *bioRxiv*, DOI 10.1101/428391. 2018.

Pre-print URL: https://www.biorxiv.org/content/early/2018/09/29/428391.

If you use anything from this repository please cite the above.

License Attribution-NonCommercial: CC BY-NC
This license lets others remix, tweak, and build upon our work non-commercially as long as they acknowledge our work.

Annotation of the included files:

| Directory | Description |
|:------------- |:-------------|
| brainseq_pipeline | processing the homogenate RNA-seq data. Uses the LIBD RNA-seq pipeline developed by EE Burke, L Collado-Torres, and AE Jaffe. | 
| BSobj_subsets | subsetting the BSSeq objects based on differential methylation at individual cytosines | 
| bsseq | creating the BSSeq objects used in later analyses | 
| bumphunting | identifying and visualizing differentially methylated regions of CpGs | 
| CREs | identifying, exploring and visualizing DNA methylation features like UMRs, LMRs, etc. | 
| DMR | exploring and visualizing the DMRs by cell type, age, and the interaction between cell type and age | 
| DMR_acf | assessing the autocorrelation of methylation between neighboring cytosines within the DMRs | 
| homogenate_RNA | explore RNA-seq from homogenate DLPFC | 
| lambda | check the lambda spike-in of the whole genome bisulfite sequencing samples | 
| meth_vs_expr | analysis of the relationship between cytosine methylation levels, splicing and gene expression | 
| misc | miscellanious analysis including global methylation patterns | 
| non-CpG | exploration and visualization of non-CpG methylation patterns | 
| psi | calculate the "percent spliced in" of splicing variants in homogenate RNA-seq data | 
| PWMEnrich | explore enrichment for canonical splice donor and acceptor sequences and transcription factors | 
| single_CpGs | analyzing single CpGs, including comparing to Lister *et al* (*Science*, 2013) data and comparing detection of age effects in homogenate and cell type-specific samples and deconvolution of neuronal subtypes | 
| sorted_nuclear_RNA | analyzing the sorted neuronal (NeuN+) and glial (NeuN-) nuclear RNAseq samples |  
| README.md | README file for the github repository | 
