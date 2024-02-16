# `funomics` Package

The `funomics` R package is a collection of functions for aggregating omics data into higher-level functional representations such as pathways, protein complexes, and cellular locations. The package provides a tool for aggregating omics data from high-throughput experiments (e.g. transcriptomics, metabolomics, proteomics) into higher-level functional activity scores that can then be used for further analysis and modeling. 

The package provides different pooling operators, such as aggregation statistics (mean, median, standard deviation, min, max), dimension-reduction derived scores (pca, nmf, mds, pathifier), or test statistics (t-test, Wilcoxon test, Kolmogorov–Smirnov test) with options for adjusting parameters and settings to suit specific research questions and data types. The package is also well-documented, with detailed descriptions of each function and an example of usage.

## Installation

The package will soon be submitted to the BioConductor repository. Until its release, you can use the following command to install the `funomics` package:

```R
devtools::install_github("elisagdelope/funomics") 
```


## Usage

To use the Funomics R package, you can load the package using the following command:

```R
library(funomics)
```

You can then access the main function provided by the package, _summarize_pathway_level_, with the type of pooling operator desired to be applied for each molecular set. The available aggregation operators and other parameters options are described in detail in the package documentation.

This function has several options for adjusting parameters and settings to suit specific research questions and data types. Here is an example usage, with a simulated gene expression matrix `X` of dimensions `p*n`, and a list of 100 gene sets `pathways`. Pathway activity is summarized using the mean pooling aggregation for those sets containing at least 12 genes. Note that you can adjust the minsize and type of aggregation as desired. 

```R
# Example usage:
p <- 10000
n <- 20
X <- matrix(rnorm(p * n), nrow = p, dimnames = list(paste0("g", 1:p), paste0("s", 1:n)))
pathways <- as.list(sample(10:100, size = 100, replace = TRUE))
pathways <- lapply(pathways, function(n, p) paste0("g", sample(1:p, size = n, replace = FALSE)), p)
names(pathways) <- paste0("pathways", 1:length(pathways))
pathway_activity <- summarize_pathway_level(X, pathways, type = "mean", minsize = 12)
```
This example mimics gene expression data and pathway gene sets, but `funomics` can be used to aggregate other types of omics data and molecular sets. For example, it can be similarly applied to gene expression data and gene sets of GO terms or protein complexes of the CORUM database. It can also be applied to a metabolomics matrix `X` and KEGG metabolic pathways.

If you have any questions or issues with the funomics R package, please contact <elisa.gomezdelope@uni.lu>. I welcome feedback and suggestions for improving the package.


## Disclaimer

The R package `funomics` implements functions for aggregating omics data into higher-level functional representations.

Copyright &copy; 2024 Elisa Gómez de Lope, University of Luxembourg, Luxembourg Centre for Systems Biomedicine (LCSB), Biomedical Data Science (BDS)

This program is free software: you can redistribute it and/or modify it under the terms of the MIT License. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; including but not limited to the warranties of MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, and NONINFRINGEMENT. See the terms of the MIT License for more details.

