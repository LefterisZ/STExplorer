% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sa_gearyC.R
\name{.int_gearyCTest}
\alias{.int_gearyCTest}
\title{Internal Function: Calculate Geary's C with statistical testing}
\usage{
.int_gearyCTest(
  gene,
  sfe,
  assay,
  listw,
  randomisation,
  alternative,
  spChk,
  adjust.n,
  zero.policy
)
}
\description{
This internal function calculates the C statistic for a vector of gene
expression alongside a statistical significance testing. For the parameter
arguments check the \code{link[STExplorer]{gearyGlobalCTest}} function.
}
\author{
Eleftherios (Lefteris) Zormpas
}
