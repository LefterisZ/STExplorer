% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sa_moranI.R
\name{.int_moranITest}
\alias{.int_moranITest}
\title{Internal Function: Calculate Moran's I with statistical testing}
\usage{
.int_moranITest(
  gene,
  sfe,
  assay,
  listw,
  randomisation,
  alternative,
  rank,
  na.action,
  spChk,
  adjust.n,
  drop.EI2,
  zero.policy
)
}
\description{
This internal function calculates the I statistic for a vector of gene
expression alongside a statistical significance testing. For the parameter
arguments check the \code{link[STExplorer]{moranGlobalITest}} function.
}
\author{
Eleftherios (Lefteris) Zormpas
}
