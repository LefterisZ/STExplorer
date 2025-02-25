% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fgwc_plot.R
\name{plotFGWC_multiHeatmap}
\alias{plotFGWC_multiHeatmap}
\title{Plot FGWC Membership percentages Heatmap}
\usage{
plotFGWC_multiHeatmap(
  fgwc,
  loc_annot = c("annotation", "cluster", "both", "none"),
  order_rows = c("annotation", "cluster", "none"),
  subset_row = NULL,
  colour = NULL,
  ...
)
}
\arguments{
\item{fgwc}{An object of class \code{fgwc} containing the results of FGWC
analysis.}

\item{loc_annot}{A character string specifying the type of annotations to
include in the heatmap. Possible values are:
\itemize{
\item \code{"annotation"}: Include only the annotation column.
\item \code{"cluster"}: Include only the cluster column.
\item \code{"both"}: Include both annotation and cluster columns.
\item \code{"none"}: Include no annotations.
}
Defaults to \code{"both"}.}

\item{order_rows}{order rows based on annotation or cluster. Possible values
are:
\itemize{
\item \code{"annotation"}: Order by the annotation column (alphabetical).
\item \code{"cluster"}: Order by the cluster column (alphabetical).
\item \code{"none"}: No ordering annotations (rows are clustered then).
}
Defaults to \code{"none"}.}

\item{subset_row}{subset rows (locations) using a boolean vector of length
equal to the number of locations.}

\item{colour}{The colour parameter passed to the \code{colour} argument in
\code{pheatmap}.}

\item{...}{Arguments passed to \code{pheatmap}.}
}
\value{
A heatmap plot of cluster membership \% with optional annotations.
}
\description{
This function generates a clustered heatmap of membership \% obtained from
FGWC (Fuzzy Geographically Weighted Clustering) analysis. It allows for
optional annotations based on clusters or annotations.
}
\examples{
\dontrun{
# Assuming `fgwc_object` is a valid FGWC object
plotFGWC_multiHeatmap(fgwc_object, loc_annot = "both")
plotFGWC_multiHeatmap(fgwc_object, loc_annot = "annotation")
plotFGWC_multiHeatmap(fgwc_object, loc_annot = "none")
}

}
\author{
Eleftherios (Lefteris) Zormpas=
}
