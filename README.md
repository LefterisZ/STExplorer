
<!-- README.md is generated from README.Rmd. Please edit that file -->

# STExplorer

<!-- badges: start -->

[![GitHub
issues](https://img.shields.io/github/issues/LefterisZ/STExplorer)](https://github.com/LefterisZ/STExplorer/issues)
[![GitHub
pulls](https://img.shields.io/github/issues-pr/LefterisZ/STExplorer)](https://github.com/LefterisZ/STExplorer/pulls)
<!-- badges: end -->

The goal of `STExplorer` is to utilise Geography’s spatial data analysis
techniques to perform spatially-aware analysis of spatial
transcriptomics data.

## Installation instructions

Get the latest stable `R` release from
[CRAN](http://cran.r-project.org/). Then install `STExplorer` from
[Bioconductor](http://bioconductor.org/) using the following code:

This option is still **NOT** available. Please use the GitHub
installation from below. The package is Bioconductor-compatible and we
are working towards submitting it to Bioconductor as soon as possible.

``` r
# if (!requireNamespace("BiocManager", quietly = TRUE)) {
#     install.packages("BiocManager")
# }
# 
# BiocManager::install("STExplorer")
```

And the development version from
[GitHub](https://github.com/LefterisZ/STExplorer) with:

``` r
BiocManager::install("LefterisZ/STExplorer")
```

## Example

The package includes a comprehensive set of vignettes with a variety of
examples and use cases. Please refer to the vignettes.

## Citation

Below is the citation output from using `citation('STExplorer')` in R.
Please run this yourself to check for any updates on how to cite
**STExplorer**.

``` r
print(citation('STExplorer'), bibtex = TRUE)
#> <0-length citation>
```

Please note that the `STExplorer` was only made possible thanks to many
other R and bioinformatics software authors, which are cited either in
the vignettes and/or the paper(s) describing this package.

## Code of Conduct

Please note that the `STExplorer` project is released with a
[Contributor Code of
Conduct](http://bioconductor.org/about/code-of-conduct/). By
contributing to this project, you agree to abide by its terms.
