
![STExplorer](vignettes/images/STExplorer_logo_hex.png)

# STExplorer v1

<!-- badges: start -->
![GitHub R package version](https://img.shields.io/github/r-package/v/LefterisZ/STExplorer)
![GitHub branch status](https://img.shields.io/github/checks-status/LefterisZ/STExplorer/main)
![Libraries.io dependency status for GitHub repo](https://img.shields.io/librariesio/github/LefterisZ/STExplorer)
![CRAN Version](https://www.r-pkg.org/badges/version/STExplorer)
<!-- badges: end -->

The goal of STExplorer is to utilise Geography's spatial data analysis techniques to perform spatially-aware analysis of spatial transcriptomics data. 

## Installation

### GitHub Development version
You can install the development version of STExplorer from [GitHub](https://github.com/) with:

``` r
# To install the package run the below:
# install.packages("devtools")
devtools::install_github("LefterisZ/STExplorer")
```

### Bioconductor stable version
This option is still not available. Please use the GitHub installation from above. The package is Bioconductor-compatible and we are working towards submitting it to Bioconductor as soon as possible.
You will be able to install the stable version of STExplorer from [Bioconductor](https://bioconductor.org/) with:

``` r
# To install the package run the below:
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("STExplorer")
```
## Loading

To load the package in R (or RStudio) follow the usual route:

``` r
# Load package to use it
library(STExplorer)
```

## Example

The package includes a comprehensive set of vignettes with a variety of examples and use cases. Please refer to the vignettes.

# Release Notes
## Version 1.0.0
This is the release version of STExplorer (*as of /10/2024*)

In this version STExplorer includes the below functionalities: 

1. Analysis methods
    - Geographically Weighted PCA (GWPCA)
    - Fuzzy Geographically Weighted Clustering (FGWC)
    - Geographically Weighted Regression (GWR)
    - Spatial Autocorrelation (SA) statistics [Global and Local]
      - Moran's I
      - Getis & Ord's G
      - Geary's C
2. Infrastructure
    - Utilises [SpatialFeatureExperiment](https://bioconductor.org/packages/devel/bioc/html/SpatialFeatureExperiment.html) (SFE) package for the data structure
    - Introduces the MetaSpatialFeatureExperiment S4 class objects which are essentially lists of SFE objects expected to include one sample each. This is done to allow smooth operation of STExplorer's geospatial methods.
3. Dependencies on geospatial packages
    - [spdep](https://cran.r-project.org/web/packages/spdep/index.html)
    - [sf](https://cran.r-project.org/web/packages/sf/index.html)
    - [GWmodel](https://cran.r-project.org/web/packages/GWmodel/index.html)
4. Dependencies on other single cell/spatial omics packages
    - [scran](https://bioconductor.org/packages/release/bioc/html/scran.html) used in modelling gene variance and selecting high variable genes
    - [scater](https://bioconductor.org/packages/release/bioc/html/scater.html) used in adding QC metrics in col- and rowData
5. Interoperability
    - Vignette explaining how STExplorer can work around [Seurat](https://satijalab.org/seurat/).
6. Examples of use
    - Vignette 1: STExplorer structure
    - Vignette 2: STExplorer pipeline with comments on the process
    - Vignette 3: STExplorer pipeline without comments (for those who are busy)
    - Vignette 4: GWPCA in detail
    - Vignette 5: FGWC in detail
    - Vignette 6: GWR in detail
    - Vignette 7: SA in detail
5. Future/upcoming work
    - Submission to Bioconductor (however the package is compatible with other Bioconductor packages)
    - Introduce more vignettes analysing other data types (VisiumHD, Xenium, CosMx, etc.)


