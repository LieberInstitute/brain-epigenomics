# brain-epigenomics
[![DOI](https://zenodo.org/badge/68875556.svg)](https://zenodo.org/badge/latestdoi/68875556)

Located in JHPCE at `/dcl01/lieber/ajaffe/lab/brain-epigenomics`

This repository contains the code used in the analyses in the following citation:

Price AJ &ast;, Collado-Torres L &ast;, Ivanov NA, Xia W, Burke EE, Shin JH, Tao R, Ma L, Jia Y,  Hyde TM, Kleinman JE, Weinberger DR, Jaffe AE &dagger;. Divergent neuronal DNA methylation patterns across human cortical development reveal critical periods and a unique role of CpH methylation. Genome Biology. 2019. DOI: [10.1186/s13059-019-1805-1](https://doi.org/10.1186/s13059-019-1805-1). Pre-print: bioRxiv, 2018. DOI: [10.1101/428391](https://doi.org/10.1101/428391).

If you use anything from this repository please cite the above. Here's a bibtex entry:

```
@article{price_divergent_2019,
	title = {Divergent neuronal {DNA} methylation patterns across human cortical development reveal critical periods and a unique role of {CpH} methylation},
	volume = {20},
	issn = {1474-760X},
	url = {https://doi.org/10.1186/s13059-019-1805-1},
	doi = {10.1186/s13059-019-1805-1},
	abstract = {DNA methylation (DNAm) is a critical regulator of both development and cellular identity and shows unique patterns in neurons. To better characterize maturational changes in DNAm patterns in these cells, we profile the DNAm landscape at single-base resolution across the first two decades of human neocortical development in NeuN+ neurons using whole-genome bisulfite sequencing and compare them to non-neurons (primarily glia) and prenatal homogenate cortex.},
	number = {1},
	urldate = {2019-09-26},
	journal = {Genome Biology},
	author = {Price, Amanda J. and Collado-Torres, Leonardo and Ivanov, Nikolay A. and Xia, Wei and Burke, Emily E. and Shin, Joo Heon and Tao, Ran and Ma, Liang and Jia, Yankai and Hyde, Thomas M. and Kleinman, Joel E. and Weinberger, Daniel R. and Jaffe, Andrew E.},
	month = sep,
	year = {2019},
	pages = {196}
}
```

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
