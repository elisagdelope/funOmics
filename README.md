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


## Contact

If you have any questions or issues with the funomics R package, please contact <elisa.gomezdelope@uni.lu>. I welcome feedback and suggestions for improving the package.


## Disclaimer

The R package `funomics` implements functions for aggregating omics data into higher-level functional representations.

Copyright &copy; 2024 Elisa Gómez de Lope, University of Luxembourg, Luxembourg Centre for Systems Biomedicine (LCSB), Biomedical Data Science (BDS)

This program is free software: you can redistribute it and/or modify it under the terms of the MIT License. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; including but not limited to the warranties of MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, and NONINFRINGEMENT. See the terms of the MIT License for more details.

