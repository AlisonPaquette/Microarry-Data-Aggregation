# Objective

The **objective** of the code presented in this repository is to (A) *preprocess* microarray data to reduce disperate probes to genes, then (B) *aggregate* this uniformly microarray data to create a large dataset composed of indivdiual microarray datasets, while adjusting for batch effects using comBat

###Step 1: Preprocess Microarray Data across Platforms
In my analysis, I analyzed data across illumina and affymetrix array platforms.  Each microarray is designed differently, with different probes which capture unique genomic features.I performed RMA normalization on affymetrix data, and Quantile Normalization on Illumina data, as recommended by manufacturers. To aggregate data, probes must be reduced to a gene level.  The subtleys of this approach and challenges are discussed in detail in [Key Issues in Conducting a Meta-Analysis of Gene Expression Microarray Datasets, Ramasamy et al, 2008](http://journals.plos.org/plosmedicine/article?id=10.1371/journal.pmed.0050184)

The R Markdown file entitled "Part1_ProcessMicroarrayData.rmd" contains functions to reduce probes to genes and remove self hybridizing probes.

**Dependencies:** Oligo (Bioconductor), Lumi (Bioconductor)

####Step 2: Aggregate Data Analysis
For my analysis, I created an aggregate dataset from merged datasets. I corrected for batch effects using [comBat](https://www.bu.edu/jlab/wp-assets/ComBat/Abstract.html), with extensive quality control checks (included), to evaluate how variability in data was assocaited with study of origin. Next,I created a series of linear regression models to identify assocations with gene expression and outcome of interest using [LIMMA](http://bioinf.wehi.edu.au/limma/) (Linear Regression for Microarray Data). I then prepared the data for downstream analysis using [DIRAC](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1000792) (Differential Rank Conservation)

**Dependencies**: LIMMA (Bioconductor)




