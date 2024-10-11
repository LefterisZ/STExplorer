# STExplorer 1.0.0 (2024-10-)
This is the release version of STExplorer. Below is an overview of STExplorer's main components.

## Analysis methods
- Geographically Weighted PCA (GWPCA)
- Fuzzy Geographically Weighted Clustering (FGWC)
- Geographically Weighted Regression (GWR)
- Spatial Autocorrelation (SA) statistics [Global and Local]
  - Moran's I
  - Getis & Ord's G
  - Geary's C

## Infrastructure
- Utilises [SpatialFeatureExperiment](https://bioconductor.org/packages/devel/bioc/html/SpatialFeatureExperiment.html) (SFE) package for the data structure
- Introduces the MetaSpatialFeatureExperiment S4 class objects which are essentially lists of SFE objects expected to include one sample each. This is done to allow smooth operation of STExplorer's geospatial methods.

## Dependencies on geospatial packages
- [spdep](https://cran.r-project.org/web/packages/spdep/index.html)
- [sf](https://cran.r-project.org/web/packages/sf/index.html)
- [GWmodel](https://cran.r-project.org/web/packages/GWmodel/index.html)

## Dependencies on other single cell/spatial omics packages
- [scran](https://bioconductor.org/packages/release/bioc/html/scran.html) used in modelling gene variance and selecting high variable genes
- [scater](https://bioconductor.org/packages/release/bioc/html/scater.html) used in adding QC metrics in col- and rowData

## Interoperability
- Vignette explaining how STExplorer can work around [Seurat](https://satijalab.org/seurat/).

## Examples of use
- Vignette 1: STExplorer structure
- Vignette 2: STExplorer pipeline with comments on the process
- Vignette 3: STExplorer pipeline without comments (for those who are busy)
- Vignette 4: GWPCA in detail
- Vignette 5: FGWC in detail
- Vignette 6: GWR in detail
- Vignette 7: SA in detail

## Future/upcoming work
- Submission to Bioconductor (however the package is compatible with other Bioconductor packages)
- Introduce more vignettes analysing other data types (VisiumHD, Xenium, CosMx, etc.)
