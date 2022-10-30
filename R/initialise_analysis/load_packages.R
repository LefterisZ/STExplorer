#----------------------------------------------------#
## LOAD/ INSTALL PACKAGES ----
#----------------------------------------------------#
## 1 Bioconductor ----
pkgBio <- c("Spaniel", "scater", "biomaRt",
            "batchelor", "scran", "DESeq2")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

## Check if packages are installed and load them or install&load them if not.
pkg.check <- lapply(
  pkgBio, 
  FUN <- function(x) {
    if (!require(x, character.only = TRUE)) {
      BiocManager::install(x, update = FALSE)
      library(x, character.only = TRUE)
    } 
  }
)


## 2 Cran ----
pkgCRAN <- c("Seurat", "cowplot", 
             "RColorBrewer",
             "harmony", "dplyr", 
             "spdep", "sf", "jsonlite",
             "tidyverse", "GWmodel", 
             "gridExtra", "ggbeeswarm",
             "egg", "ggpubr", "scales", 
             "pheatmap", "rlist")

## Check if packages are installed and load them or install&load them if not.
pkg.check <- lapply(
  pkgCRAN, 
  FUN <- function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)

## 3 GitHub ----
pkgGit <- c("RachelQueen1/SCFunctionsV3",
            "eddelbuettel/rbenchmark",
            "mtennekes/cols4all")

if (!require("devtools", quietly = TRUE))
  install.packages("devtools")

## Check if packages are installed and load them or install&load them if not.
pkg.check <- lapply(
  pkgGit, 
  FUN <- function(x) {
    pkg.name <- sub(".*/", "", x)
    if (!require(pkg.name, character.only = TRUE)) {
      devtools::install_git(paste0("https://github.com/", x),
                            force = TRUE)
      library(pkg.name, character.only = TRUE)
    }
  }
)

## 4 source scripts ----
# List all scripts in R folder and remove the ones located in the core folder.
file.sources <- list.files(path = "./R/core",
                           pattern = "*.R$",
                           recursive = TRUE,
                           full.names = TRUE,
                           ignore.case = TRUE)

sapply(file.sources, source, .GlobalEnv)

